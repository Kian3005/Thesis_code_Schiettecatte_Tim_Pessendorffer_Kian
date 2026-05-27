# ============================================================
# Standard RSF prediction framework for first-goal timing
# - Separate home and away marginal RSF models
# - Expanding-window development validation
# - Final held-out evaluation on season 2024-2025
# - Marginal C-index and IPCW Brier scores
# - Derived football-event probabilities
# - Out-of-fold isotonic recalibration diagnostic
# - Final variable importance and Brier plot
# ============================================================

# ----------------------------
# 0) Packages and seed
# ----------------------------
needed_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "readr",
  "survival", "randomForestSRC", "tibble", "purrr"
)

for (p in needed_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
invisible(lapply(needed_pkgs, library, character.only = TRUE))

set.seed(123)

# ----------------------------
# 1) Settings
# ----------------------------
data_path <- "SurvDataShortFormat.csv"
output_dir <- "rsf_prediction_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

date_var <- "Date"
league_var <- "Round_grouped"
season_var <- "season"
row_id_var <- "..row_id"
league_feature_var <- "league_feature"

final_test_season <- 2425

home_time_var <- "home_first_goal_minute"
home_status_var <- "home_scored"
away_time_var <- "away_first_goal_minute"
away_status_var <- "away_scored"

base_cat_vars <- c(
  "derby_type",
  "Round",
  "home_formation_grouped",
  "away_formation_grouped",
  league_feature_var
)

roll_pattern <- "_roll"
corr_threshold <- 0.95

# Development expanding-window design
n_folds <- 8
anchor_quantile_start <- 0.40
anchor_quantile_end <- 0.88
anchor_jitter_days <- 28
test_horizon_range_days <- c(21, 49)
min_gap_days <- 7
min_train_n <- 2500
min_test_n <- 150
min_test_leagues <- 4

# RSF tuning grid used for the reported standard RSF analysis
ntree_search <- 800
nodesize_grid <- c(10, 15, 20, 25)
mtry_mode_grid <- c("sqrt", "p4", "p3", "p2")
nodedepth_grid <- list(NULL, 8, 12)

early_stop_tolerance <- 0.0015
early_stop_min_trees <- 100
early_stop_step <- 25

surv_eval_times <- c(15, 30, 45, 60, 75, 90)
minute_grid <- 0:90

calibration_targets <- c(
  "p_home_goal_by_45",
  "p_away_goal_by_45",
  "p_any_goal_by_45",
  "p_home_first",
  "p_away_first",
  "p_no_goal_by_90",
  "p_any_goal_second_half_cond_no_goal_45"
)

# ----------------------------
# 2) General helper functions
# ----------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

load_data_any <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file format. Use .csv or .rds.")
  }
}

safe_as_date <- function(x) {
  out <- as.Date(x)
  if (all(is.na(out))) out <- as.Date(readr::parse_date(as.character(x)))
  out
}

drop_nuisance_cols <- function(df) {
  nuisance <- c(
    "Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.2",
    "Unnamed.0", "Unnamed.0.1", "Unnamed.0.2",
    "index", "id", "match_id"
  )
  nuisance <- nuisance[nuisance %in% names(df)]
  if (length(nuisance) > 0) df <- df %>% select(-all_of(nuisance))
  df
}

get_roll_vars <- function(df, pattern = "_roll") {
  names(df)[grepl(pattern, names(df))]
}

drop_zero_var_cols <- function(df) {
  if (ncol(df) == 0) return(list(data = df, removed = character(0)))
  sds <- vapply(df, function(x) sd(x, na.rm = TRUE), numeric(1))
  keep <- !(is.na(sds) | sds == 0)
  list(data = df[, keep, drop = FALSE], removed = names(df)[!keep])
}

find_high_corr_to_drop <- function(df_num, threshold = 0.95) {
  if (ncol(df_num) <= 1) return(character(0))
  C <- cor(as.matrix(df_num), use = "pairwise.complete.obs")
  C[is.na(C)] <- 0
  drop_vars <- character(0)

  repeat {
    C_work <- C
    diag(C_work) <- 0
    max_corr <- max(abs(C_work))
    if (!is.finite(max_corr) || max_corr < threshold) break

    ij <- which(abs(C_work) == max_corr, arr.ind = TRUE)[1, ]
    i <- ij[1]
    j <- ij[2]

    mean_i <- mean(abs(C[i, -i]), na.rm = TRUE)
    mean_j <- mean(abs(C[j, -j]), na.rm = TRUE)
    drop_idx <- if (mean_i >= mean_j) i else j
    drop_vars <- c(drop_vars, colnames(C)[drop_idx])

    keep <- setdiff(colnames(C), drop_vars)
    if (length(keep) <= 1) break
    C <- C[keep, keep, drop = FALSE]
  }

  unique(drop_vars)
}

align_factor_levels_safe <- function(train_df, test_df, factor_vars, new_level = "__NEW__") {
  if (length(factor_vars) == 0) return(list(train = train_df, test = test_df))

  for (v in factor_vars) {
    train_chr <- as.character(train_df[[v]])
    test_chr <- as.character(test_df[[v]])

    train_levels <- sort(unique(train_chr[!is.na(train_chr)]))
    all_levels <- c(train_levels, new_level)
    test_chr[!is.na(test_chr) & !(test_chr %in% train_levels)] <- new_level

    train_df[[v]] <- factor(train_chr, levels = all_levels)
    test_df[[v]] <- factor(test_chr, levels = all_levels)
  }

  list(train = train_df, test = test_df)
}

# ----------------------------
# 3) Expanding-window folds
# ----------------------------
build_random_expanding_folds <- function(data,
                                         date_var,
                                         league_var,
                                         n_folds,
                                         anchor_quantile_start,
                                         anchor_quantile_end,
                                         anchor_jitter_days,
                                         test_horizon_range_days,
                                         min_gap_days,
                                         min_train_n,
                                         min_test_n,
                                         min_test_leagues,
                                         seed = 123) {
  set.seed(seed)
  dd <- data %>% arrange(.data[[date_var]])
  unique_dates <- sort(unique(dd[[date_var]]))
  unique_dates <- unique_dates[!is.na(unique_dates)]

  probs <- seq(anchor_quantile_start, anchor_quantile_end, length.out = n_folds)
  target_dates <- as.Date(
    as.numeric(stats::quantile(unique_dates, probs = probs, type = 1)),
    origin = "1970-01-01"
  )

  folds <- vector("list", n_folds)
  last_test_end <- min(unique_dates) - 1

  for (i in seq_len(n_folds)) {
    cand <- unique_dates[
      unique_dates >= target_dates[i] - anchor_jitter_days &
        unique_dates <= target_dates[i] + anchor_jitter_days &
        unique_dates > last_test_end + min_gap_days
    ]
    if (length(cand) == 0) cand <- unique_dates[unique_dates > last_test_end + min_gap_days]
    if (length(cand) == 0) next

    valid_fold_found <- FALSE
    attempt <- 1

    while (!valid_fold_found && attempt <= 200) {
      anchor <- sample(cand, size = 1)
      horizon_days <- sample(seq(test_horizon_range_days[1], test_horizon_range_days[2]), size = 1)
      test_end <- anchor + horizon_days

      train_idx <- dd[[date_var]] <= anchor
      test_idx <- dd[[date_var]] > anchor & dd[[date_var]] <= test_end

      n_train <- sum(train_idx)
      n_test <- sum(test_idx)
      test_leagues <- dd %>% filter(test_idx) %>% distinct(.data[[league_var]]) %>% nrow()

      if (n_train >= min_train_n && n_test >= min_test_n && test_leagues >= min_test_leagues) {
        folds[[i]] <- list(
          fold_id = i,
          anchor_date = anchor,
          test_end_date = test_end,
          horizon_days = horizon_days,
          train_idx = which(train_idx),
          test_idx = which(test_idx),
          n_train = n_train,
          n_test = n_test,
          n_test_leagues = test_leagues
        )
        last_test_end <- test_end
        valid_fold_found <- TRUE
      }

      attempt <- attempt + 1
    }
  }

  folds <- folds[!vapply(folds, is.null, logical(1))]
  if (length(folds) == 0) stop("No valid expanding-window folds could be created.")
  folds
}

# ----------------------------
# 4) Design preparation
# ----------------------------
prepare_fold_design <- function(train_df,
                                test_df,
                                roll_vars,
                                cat_vars,
                                outcome_time,
                                outcome_status,
                                date_var,
                                league_var,
                                row_id_var,
                                corr_threshold = 0.95) {
  core_keep <- unique(c(
    row_id_var, date_var, league_var,
    outcome_time, outcome_status,
    cat_vars, roll_vars
  ))

  train_df <- train_df %>% select(all_of(core_keep))
  test_df <- test_df %>% select(all_of(core_keep))

  zv <- drop_zero_var_cols(train_df[, roll_vars, drop = FALSE])
  roll_keep <- colnames(zv$data)

  high_corr_drop <- find_high_corr_to_drop(train_df[, roll_keep, drop = FALSE], threshold = corr_threshold)
  roll_keep <- setdiff(roll_keep, high_corr_drop)

  aligned <- align_factor_levels_safe(
    train_df = train_df[, cat_vars, drop = FALSE],
    test_df = test_df[, cat_vars, drop = FALSE],
    factor_vars = cat_vars
  )

  meta_cols <- c(row_id_var, date_var, league_var, outcome_time, outcome_status)

  train_out <- bind_cols(
    train_df[, meta_cols, drop = FALSE],
    aligned$train,
    train_df[, roll_keep, drop = FALSE]
  )

  test_out <- bind_cols(
    test_df[, meta_cols, drop = FALSE],
    aligned$test,
    test_df[, roll_keep, drop = FALSE]
  )

  list(
    train = train_out,
    test = test_out,
    kept_roll_vars = roll_keep,
    dropped_zero_var = zv$removed,
    dropped_high_corr = high_corr_drop
  )
}

get_model_data <- function(df, row_id_var, date_var, league_var) {
  df %>% select(-all_of(c(row_id_var, date_var, league_var)))
}

# ----------------------------
# 5) RSF fitting, tuning and prediction
# ----------------------------
get_mtry_value <- function(p, mode) {
  if (mode == "sqrt") return(max(2, floor(sqrt(p))))
  if (mode == "p4") return(max(2, floor(p / 4)))
  if (mode == "p3") return(max(2, floor(p / 3)))
  if (mode == "p2") return(max(2, floor(p / 2)))
  stop("Unknown mtry mode: ", mode)
}

choose_ntree_from_oob <- function(err_rate, tolerance, min_trees, step) {
  err <- as.numeric(err_rate)
  candidate_idx <- seq(min_trees, length(err), by = step)
  if (length(candidate_idx) == 0) candidate_idx <- seq_along(err)

  err_sub <- err[candidate_idx]
  best_err <- min(err_sub, na.rm = TRUE)
  acceptable <- candidate_idx[err_sub <= best_err + tolerance]

  if (length(acceptable) == 0) candidate_idx[which.min(err_sub)] else min(acceptable)
}

fit_rsf_model <- function(train_df,
                          outcome_time,
                          outcome_status,
                          ntree,
                          nodesize,
                          mtry,
                          nodedepth = NULL,
                          importance = FALSE,
                          seed = 123) {
  mod_df <- get_model_data(train_df, row_id_var, date_var, league_var)

  randomForestSRC::rfsrc(
    formula = as.formula(paste0("Surv(", outcome_time, ", ", outcome_status, ") ~ .")),
    data = mod_df,
    ntree = ntree,
    nodesize = nodesize,
    mtry = mtry,
    nodedepth = nodedepth,
    importance = importance,
    block.size = 1,
    seed = seed,
    forest = TRUE,
    save.memory = TRUE
  )
}

tune_rsf_margin <- function(train_df,
                            outcome_time,
                            outcome_status,
                            ntree_search,
                            nodesize_grid,
                            mtry_mode_grid,
                            nodedepth_grid,
                            seed = 123) {
  p <- ncol(train_df) - 5
  candidates <- expand.grid(
    nodesize = nodesize_grid,
    mtry_mode = mtry_mode_grid,
    depth_id = seq_along(nodedepth_grid),
    stringsAsFactors = FALSE
  )

  best <- NULL
  best_err <- Inf

  for (k in seq_len(nrow(candidates))) {
    nodesize <- candidates$nodesize[k]
    mtry_mode <- candidates$mtry_mode[k]
    mtry <- get_mtry_value(p, mtry_mode)
    nodedepth <- nodedepth_grid[[candidates$depth_id[k]]]

    fit <- fit_rsf_model(
      train_df = train_df,
      outcome_time = outcome_time,
      outcome_status = outcome_status,
      ntree = ntree_search,
      nodesize = nodesize,
      mtry = mtry,
      nodedepth = nodedepth,
      importance = FALSE,
      seed = seed + k
    )

    ntree_opt <- choose_ntree_from_oob(
      fit$err.rate,
      tolerance = early_stop_tolerance,
      min_trees = early_stop_min_trees,
      step = early_stop_step
    )

    oob_best <- min(fit$err.rate, na.rm = TRUE)

    if (oob_best < best_err) {
      best_err <- oob_best
      best <- list(
        nodesize = nodesize,
        mtry_mode = mtry_mode,
        mtry = mtry,
        nodedepth = nodedepth,
        ntree = ntree_opt,
        oob_best = oob_best
      )
    }
  }

  best
}

predict_rsf_survival <- function(fit, newdata_df) {
  newx <- get_model_data(newdata_df, row_id_var, date_var, league_var)
  predict(fit, newdata = newx)
}

risk_from_rfsrc_prediction <- function(pred_obj) {
  if (!is.null(pred_obj$mortality)) return(as.numeric(pred_obj$mortality))
  if (!is.null(pred_obj$predicted)) return(as.numeric(pred_obj$predicted))
  if (!is.null(pred_obj$chf)) return(rowMeans(pred_obj$chf, na.rm = TRUE))
  stop("Could not extract a risk score from the RSF prediction object.")
}

calc_cindex <- function(time, status, risk) {
  survival::concordance(
    survival::Surv(time, status) ~ risk,
    reverse = TRUE
  )$concordance
}

# ----------------------------
# 6) Survival probabilities and football-event probabilities
# ----------------------------
surv_matrix_at_times <- function(pred_obj, eval_times) {
  surv_mat <- pred_obj$survival
  pred_times <- pred_obj$time.interest
  if (is.null(surv_mat) || is.null(pred_times)) stop("Prediction object missing survival curves.")

  idx <- vapply(eval_times, function(tt) {
    j <- suppressWarnings(max(which(pred_times <= tt)))
    if (!is.finite(j)) NA_integer_ else j
  }, integer(1))

  out <- matrix(1, nrow = nrow(surv_mat), ncol = length(eval_times))
  colnames(out) <- paste0("t", eval_times)

  for (k in seq_along(eval_times)) {
    if (!is.na(idx[k])) out[, k] <- surv_mat[, idx[k]]
  }

  out
}

derive_match_probabilities <- function(home_pred, away_pred, minute_grid = 0:90) {
  Sh <- surv_matrix_at_times(home_pred, minute_grid)
  Sa <- surv_matrix_at_times(away_pred, minute_grid)

  idx45 <- which(minute_grid == 45)
  idx90 <- which(minute_grid == 90)

  Sh45 <- Sh[, idx45]
  Sa45 <- Sa[, idx45]
  Sh90 <- Sh[, idx90]
  Sa90 <- Sa[, idx90]

  # Discrete minute-wise event probabilities on 1:90.
  pH <- Sh[, -ncol(Sh), drop = FALSE] - Sh[, -1, drop = FALSE]
  pA <- Sa[, -ncol(Sa), drop = FALSE] - Sa[, -1, drop = FALSE]

  # Strict first-scorer convention: same-minute first goals are not counted as either team first.
  p_home_first <- rowSums(pH * Sa[, -1, drop = FALSE], na.rm = TRUE)
  p_away_first <- rowSums(pA * Sh[, -1, drop = FALSE], na.rm = TRUE)

  p_home_goal_by_45 <- 1 - Sh45
  p_away_goal_by_45 <- 1 - Sa45
  p_any_goal_by_45 <- 1 - Sh45 * Sa45
  p_no_goal_by_90 <- Sh90 * Sa90

  p_any_goal_second_half_cond_no_goal_45 <-
    1 - (Sh90 / pmax(Sh45, 1e-8)) * (Sa90 / pmax(Sa45, 1e-8))

  tibble(
    p_home_goal_by_45 = p_home_goal_by_45,
    p_away_goal_by_45 = p_away_goal_by_45,
    p_any_goal_by_45 = p_any_goal_by_45,
    p_home_first = p_home_first,
    p_away_first = p_away_first,
    p_no_goal_by_90 = p_no_goal_by_90,
    p_any_goal_second_half_cond_no_goal_45 = p_any_goal_second_half_cond_no_goal_45
  )
}

observed_match_outcomes <- function(df) {
  h_time <- df[[home_time_var]]
  a_time <- df[[away_time_var]]
  h_sc <- df[[home_status_var]]
  a_sc <- df[[away_status_var]]

  h_goal_45 <- as.integer(h_sc == 1 & h_time <= 45)
  a_goal_45 <- as.integer(a_sc == 1 & a_time <= 45)
  any_goal_45 <- as.integer(h_goal_45 == 1 | a_goal_45 == 1)

  home_first <- as.integer(h_sc == 1 & (a_sc == 0 | h_time < a_time))
  away_first <- as.integer(a_sc == 1 & (h_sc == 0 | a_time < h_time))

  no_goal_90 <- as.integer(h_sc == 0 & a_sc == 0)
  no_goal_45 <- 1L - any_goal_45

  home_goal_second_half <- as.integer(h_sc == 1 & h_time > 45 & h_time <= 90)
  away_goal_second_half <- as.integer(a_sc == 1 & a_time > 45 & a_time <= 90)
  any_goal_second_half <- as.integer(home_goal_second_half == 1 | away_goal_second_half == 1)

  tibble(
    !!row_id_var := df[[row_id_var]],
    obs_home_goal_by_45 = h_goal_45,
    obs_away_goal_by_45 = a_goal_45,
    obs_any_goal_by_45 = any_goal_45,
    obs_home_first = home_first,
    obs_away_first = away_first,
    obs_no_goal_by_90 = no_goal_90,
    obs_any_goal_second_half_cond_no_goal_45 = ifelse(no_goal_45 == 1, any_goal_second_half, NA_integer_)
  )
}

# ----------------------------
# 7) Performance metrics
# ----------------------------
km_censor_fit <- function(time, status) {
  survival::survfit(survival::Surv(time, 1 - status) ~ 1)
}

km_surv_at <- function(km_fit, t) {
  s <- summary(km_fit, times = t, extend = TRUE)$surv
  if (length(s) == 0) return(1)
  s[is.na(s)] <- 1
  pmax(s, 1e-6)
}

brier_ipcw_at_time <- function(time, status, surv_prob_t, eval_time) {
  km_cens <- km_censor_fit(time, status)
  G_t <- km_surv_at(km_cens, eval_time)

  obs <- as.integer(time > eval_time)
  weights <- rep(0, length(time))

  event_by_t <- time <= eval_time & status == 1
  survive_past_t <- time > eval_time

  if (any(event_by_t)) {
    G_event <- vapply(time[event_by_t], function(tt) km_surv_at(km_cens, tt - 1e-8), numeric(1))
    weights[event_by_t] <- 1 / G_event
  }

  if (any(survive_past_t)) {
    weights[survive_past_t] <- 1 / G_t
  }

  mean(weights * (obs - surv_prob_t)^2, na.rm = TRUE)
}

margin_brier_summary <- function(test_df, pred_obj, outcome_time, outcome_status, eval_times) {
  surv_at_t <- surv_matrix_at_times(pred_obj, eval_times)
  tibble(
    time = eval_times,
    brier = sapply(seq_along(eval_times), function(k) {
      brier_ipcw_at_time(
        time = test_df[[outcome_time]],
        status = test_df[[outcome_status]],
        surv_prob_t = surv_at_t[, k],
        eval_time = eval_times[k]
      )
    })
  )
}

clip_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

binary_brier <- function(y, p) {
  idx <- !is.na(y) & !is.na(p)
  if (!any(idx)) return(NA_real_)
  mean((y[idx] - p[idx])^2)
}

binary_logloss <- function(y, p) {
  idx <- !is.na(y) & !is.na(p)
  if (!any(idx)) return(NA_real_)
  p <- clip_prob(p[idx])
  y <- y[idx]
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

score_event_target <- function(df, p_name) {
  obs_name <- sub("^p_", "obs_", p_name)
  y <- df[[obs_name]]
  p <- df[[p_name]]
  idx <- !is.na(y) & !is.na(p)

  if (!any(idx)) {
    return(tibble(
      target = p_name,
      event_rate = NA_real_,
      baseline_brier = NA_real_,
      model_brier = NA_real_,
      model_logloss = NA_real_,
      n = 0L
    ))
  }

  y <- y[idx]
  p <- p[idx]
  p0 <- mean(y)

  tibble(
    target = p_name,
    event_rate = p0,
    baseline_brier = binary_brier(y, rep(p0, length(y))),
    model_brier = binary_brier(y, p),
    model_logloss = binary_logloss(y, p),
    n = length(y)
  )
}

# ----------------------------
# 8) Isotonic calibration diagnostic
# ----------------------------
fit_isotonic_calibrator <- function(p, y) {
  df <- tibble(p = p, y = y) %>% filter(!is.na(p), !is.na(y)) %>% arrange(p)
  if (nrow(df) < 20) return(NULL)

  iso <- isoreg(df$p, df$y)
  list(fun = approxfun(x = iso$x, y = iso$yf, method = "linear", ties = "ordered", rule = 2))
}

predict_isotonic_calibrator <- function(cal_obj, p) {
  if (is.null(cal_obj)) return(p)
  clip_prob(cal_obj$fun(clip_prob(p)))
}

fit_calibrators_from_oof <- function(oof_df, targets) {
  calibrators <- list()
  for (nm in targets) {
    obs_nm <- sub("^p_", "obs_", nm)
    if (nm %in% names(oof_df) && obs_nm %in% names(oof_df)) {
      calibrators[[nm]] <- fit_isotonic_calibrator(oof_df[[nm]], oof_df[[obs_nm]])
    }
  }
  calibrators
}

apply_calibrators <- function(df, calibrators) {
  out <- df
  for (nm in names(calibrators)) {
    if (nm %in% names(out)) {
      out[[paste0(nm, "_cal")]] <- predict_isotonic_calibrator(calibrators[[nm]], out[[nm]])
    }
  }
  out
}

compare_raw_calibrated <- function(df, p_name) {
  obs_name <- sub("^p_", "obs_", p_name)
  p_cal_name <- paste0(p_name, "_cal")

  tibble(
    target = p_name,
    raw_brier = binary_brier(df[[obs_name]], df[[p_name]]),
    raw_logloss = binary_logloss(df[[obs_name]], df[[p_name]]),
    cal_brier = if (p_cal_name %in% names(df)) binary_brier(df[[obs_name]], df[[p_cal_name]]) else NA_real_,
    cal_logloss = if (p_cal_name %in% names(df)) binary_logloss(df[[obs_name]], df[[p_cal_name]]) else NA_real_
  )
}

# ----------------------------
# 9) Fit/predict one margin
# ----------------------------
fit_predict_one_margin <- function(train_df,
                                   test_df,
                                   outcome_time,
                                   outcome_status,
                                   seed,
                                   importance = FALSE,
                                   tuned = NULL) {
  if (is.null(tuned)) {
    tuned <- tune_rsf_margin(
      train_df = train_df,
      outcome_time = outcome_time,
      outcome_status = outcome_status,
      ntree_search = ntree_search,
      nodesize_grid = nodesize_grid,
      mtry_mode_grid = mtry_mode_grid,
      nodedepth_grid = nodedepth_grid,
      seed = seed
    )
  }

  fit <- fit_rsf_model(
    train_df = train_df,
    outcome_time = outcome_time,
    outcome_status = outcome_status,
    ntree = tuned$ntree,
    nodesize = tuned$nodesize,
    mtry = tuned$mtry,
    nodedepth = tuned$nodedepth,
    importance = importance,
    seed = seed + 1000
  )

  pred <- predict_rsf_survival(fit, test_df)
  risk <- risk_from_rfsrc_prediction(pred)

  cindex <- calc_cindex(
    time = test_df[[outcome_time]],
    status = test_df[[outcome_status]],
    risk = risk
  )

  brier_tbl <- margin_brier_summary(
    test_df = test_df,
    pred_obj = pred,
    outcome_time = outcome_time,
    outcome_status = outcome_status,
    eval_times = surv_eval_times
  )

  list(
    fit = fit,
    pred = pred,
    risk = risk,
    cindex = cindex,
    tuned = tuned,
    brier_tbl = brier_tbl,
    test_df = test_df
  )
}

# ----------------------------
# 10) One expanding-window fold
# ----------------------------
run_one_fold <- function(fold_i, dev_df) {
  train_fold <- dev_df[fold_i$train_idx, , drop = FALSE] %>% arrange(.data[[date_var]])
  test_fold <- dev_df[fold_i$test_idx, , drop = FALSE] %>% arrange(.data[[date_var]])

  home_design <- prepare_fold_design(
    train_fold, test_fold, roll_vars, cat_vars,
    home_time_var, home_status_var, date_var, league_var, row_id_var, corr_threshold
  )

  away_design <- prepare_fold_design(
    train_fold, test_fold, roll_vars, cat_vars,
    away_time_var, away_status_var, date_var, league_var, row_id_var, corr_threshold
  )

  if (!identical(home_design$test[[row_id_var]], away_design$test[[row_id_var]])) {
    stop("Home and away fold test sets are not aligned.")
  }

  home_fit <- fit_predict_one_margin(
    train_df = home_design$train,
    test_df = home_design$test,
    outcome_time = home_time_var,
    outcome_status = home_status_var,
    seed = 100 + fold_i$fold_id,
    importance = FALSE
  )

  away_fit <- fit_predict_one_margin(
    train_df = away_design$train,
    test_df = away_design$test,
    outcome_time = away_time_var,
    outcome_status = away_status_var,
    seed = 200 + fold_i$fold_id,
    importance = FALSE
  )

  observed <- observed_match_outcomes(test_fold)
  probs <- derive_match_probabilities(home_fit$pred, away_fit$pred, minute_grid)

  pred_tbl <- bind_cols(
    home_design$test %>% select(all_of(c(row_id_var, date_var, league_var, home_time_var, home_status_var))),
    away_design$test %>% select(all_of(c(away_time_var, away_status_var))),
    probs,
    observed %>% select(-all_of(row_id_var))
  )

  derived_scores <- bind_rows(lapply(calibration_targets, function(x) score_event_target(pred_tbl, x)))

  list(
    fold_summary = tibble(
      fold_id = fold_i$fold_id,
      anchor_date = fold_i$anchor_date,
      test_end_date = fold_i$test_end_date,
      horizon_days = fold_i$horizon_days,
      home_cindex = home_fit$cindex,
      away_cindex = away_fit$cindex,
      home_nodesize = home_fit$tuned$nodesize,
      away_nodesize = away_fit$tuned$nodesize,
      home_mtry = home_fit$tuned$mtry,
      away_mtry = away_fit$tuned$mtry,
      home_ntree = home_fit$tuned$ntree,
      away_ntree = away_fit$tuned$ntree,
      home_nodedepth = ifelse(is.null(home_fit$tuned$nodedepth), NA, home_fit$tuned$nodedepth),
      away_nodedepth = ifelse(is.null(away_fit$tuned$nodedepth), NA, away_fit$tuned$nodedepth)
    ),
    home_brier = home_fit$brier_tbl %>% mutate(fold_id = fold_i$fold_id),
    away_brier = away_fit$brier_tbl %>% mutate(fold_id = fold_i$fold_id),
    derived_scores = derived_scores %>% mutate(fold_id = fold_i$fold_id),
    oof_predictions = pred_tbl
  )
}

# ----------------------------
# 11) Load data and create analysis dataset
# ----------------------------
log_msg("Loading and preparing data")

df <- load_data_any(data_path) %>% drop_nuisance_cols()

required_cols <- c(
  date_var, league_var, season_var,
  home_time_var, home_status_var,
  away_time_var, away_status_var
)

missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) stop("Missing required columns: ", paste(missing_required, collapse = ", "))

df <- df %>%
  mutate(
    !!row_id_var := seq_len(n()),
    !!date_var := safe_as_date(.data[[date_var]]),
    !!home_time_var := as.numeric(.data[[home_time_var]]),
    !!away_time_var := as.numeric(.data[[away_time_var]]),
    !!home_status_var := as.integer(.data[[home_status_var]]),
    !!away_status_var := as.integer(.data[[away_status_var]]),
    !!league_feature_var := .data[[league_var]]
  )

cat_vars <- base_cat_vars[base_cat_vars %in% names(df)]
cat_vars <- setdiff(cat_vars, league_var)
roll_vars <- get_roll_vars(df, roll_pattern)
if (length(roll_vars) == 0) stop("No rolling variables found.")

analysis_df <- df %>%
  filter(
    !is.na(.data[[date_var]]),
    !is.na(.data[[season_var]]),
    is.finite(.data[[home_time_var]]),
    is.finite(.data[[away_time_var]]),
    .data[[home_time_var]] >= 0,
    .data[[away_time_var]] >= 0,
    .data[[home_status_var]] %in% c(0, 1),
    .data[[away_status_var]] %in% c(0, 1)
  ) %>%
  select(all_of(c(
    row_id_var, date_var, league_var, season_var,
    home_time_var, home_status_var,
    away_time_var, away_status_var,
    cat_vars, roll_vars
  )))

dev_df <- analysis_df %>% filter(.data[[season_var]] != final_test_season)
final_test_df <- analysis_df %>% filter(.data[[season_var]] == final_test_season)

log_msg("Development rows: ", nrow(dev_df))
log_msg("Final test rows: ", nrow(final_test_df))

# ----------------------------
# 12) Development expanding-window validation
# ----------------------------
log_msg("Building expanding-window folds")

folds <- build_random_expanding_folds(
  data = dev_df,
  date_var = date_var,
  league_var = league_var,
  n_folds = n_folds,
  anchor_quantile_start = anchor_quantile_start,
  anchor_quantile_end = anchor_quantile_end,
  anchor_jitter_days = anchor_jitter_days,
  test_horizon_range_days = test_horizon_range_days,
  min_gap_days = min_gap_days,
  min_train_n = min_train_n,
  min_test_n = min_test_n,
  min_test_leagues = min_test_leagues,
  seed = 123
)

fold_overview <- bind_rows(lapply(folds, function(x) {
  tibble(
    fold_id = x$fold_id,
    anchor_date = x$anchor_date,
    test_end_date = x$test_end_date,
    horizon_days = x$horizon_days,
    n_train = x$n_train,
    n_test = x$n_test,
    n_test_leagues = x$n_test_leagues
  )
}))

write.csv(fold_overview, file.path(output_dir, "fold_overview_dev.csv"), row.names = FALSE)

log_msg("Running development folds")
fold_results <- vector("list", length(folds))
for (i in seq_along(folds)) {
  log_msg("Fold ", i, " / ", length(folds))
  fold_results[[i]] <- run_one_fold(folds[[i]], dev_df)
}

fold_summary_df <- bind_rows(lapply(fold_results, `[[`, "fold_summary"))
home_brier_dev <- bind_rows(lapply(fold_results, `[[`, "home_brier"))
away_brier_dev <- bind_rows(lapply(fold_results, `[[`, "away_brier"))
derived_dev_scores <- bind_rows(lapply(fold_results, `[[`, "derived_scores"))
oof_df <- bind_rows(lapply(fold_results, `[[`, "oof_predictions"))

write.csv(fold_summary_df, file.path(output_dir, "dev_fold_summary.csv"), row.names = FALSE)
write.csv(home_brier_dev, file.path(output_dir, "home_dev_brier.csv"), row.names = FALSE)
write.csv(away_brier_dev, file.path(output_dir, "away_dev_brier.csv"), row.names = FALSE)
write.csv(derived_dev_scores, file.path(output_dir, "derived_dev_scores.csv"), row.names = FALSE)
write.csv(oof_df, file.path(output_dir, "oof_predictions_raw.csv"), row.names = FALSE)

# ----------------------------
# 13) OOF isotonic calibration diagnostic
# ----------------------------
log_msg("Fitting isotonic calibrators from OOF predictions")

calibrators <- fit_calibrators_from_oof(oof_df, calibration_targets)
oof_df_cal <- apply_calibrators(oof_df, calibrators)
oof_calibration_comparison <- bind_rows(lapply(calibration_targets, function(x) compare_raw_calibrated(oof_df_cal, x)))

write.csv(oof_df_cal, file.path(output_dir, "oof_predictions_calibrated.csv"), row.names = FALSE)
write.csv(oof_calibration_comparison, file.path(output_dir, "oof_raw_vs_calibrated.csv"), row.names = FALSE)
saveRDS(calibrators, file.path(output_dir, "isotonic_calibrators.rds"))

# ----------------------------
# 14) Lock final hyperparameter settings from development folds
# ----------------------------
lock_margin_settings <- function(fold_summary_df, margin = c("home", "away")) {
  margin <- match.arg(margin)

  if (margin == "home") {
    fold_summary_df %>%
      group_by(home_nodesize, home_mtry, home_nodedepth) %>%
      summarise(
        n_folds = n(),
        mean_cindex = mean(home_cindex, na.rm = TRUE),
        median_ntree = median(home_ntree, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_folds), desc(mean_cindex)) %>%
      slice(1) %>%
      transmute(
        nodesize = home_nodesize,
        mtry = home_mtry,
        nodedepth = home_nodedepth,
        ntree = as.integer(round(median_ntree)),
        n_folds = n_folds,
        mean_cindex = mean_cindex
      )
  } else {
    fold_summary_df %>%
      group_by(away_nodesize, away_mtry, away_nodedepth) %>%
      summarise(
        n_folds = n(),
        mean_cindex = mean(away_cindex, na.rm = TRUE),
        median_ntree = median(away_ntree, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_folds), desc(mean_cindex)) %>%
      slice(1) %>%
      transmute(
        nodesize = away_nodesize,
        mtry = away_mtry,
        nodedepth = away_nodedepth,
        ntree = as.integer(round(median_ntree)),
        n_folds = n_folds,
        mean_cindex = mean_cindex
      )
  }
}

home_locked <- lock_margin_settings(fold_summary_df, "home")
away_locked <- lock_margin_settings(fold_summary_df, "away")

write.csv(home_locked, file.path(output_dir, "locked_home_settings.csv"), row.names = FALSE)
write.csv(away_locked, file.path(output_dir, "locked_away_settings.csv"), row.names = FALSE)

# ----------------------------
# 15) Final refit and held-out evaluation
# ----------------------------
log_msg("Final held-out test: season ", final_test_season)

final_home_design <- prepare_fold_design(
  dev_df, final_test_df, roll_vars, cat_vars,
  home_time_var, home_status_var, date_var, league_var, row_id_var, corr_threshold
)

final_away_design <- prepare_fold_design(
  dev_df, final_test_df, roll_vars, cat_vars,
  away_time_var, away_status_var, date_var, league_var, row_id_var, corr_threshold
)

if (!identical(final_home_design$test[[row_id_var]], final_away_design$test[[row_id_var]])) {
  stop("Final home and away test sets are not aligned.")
}

home_tuned <- list(
  nodesize = home_locked$nodesize,
  mtry = home_locked$mtry,
  nodedepth = if (is.na(home_locked$nodedepth)) NULL else home_locked$nodedepth,
  ntree = home_locked$ntree
)

away_tuned <- list(
  nodesize = away_locked$nodesize,
  mtry = away_locked$mtry,
  nodedepth = if (is.na(away_locked$nodedepth)) NULL else away_locked$nodedepth,
  ntree = away_locked$ntree
)

final_home_fit <- fit_predict_one_margin(
  train_df = final_home_design$train,
  test_df = final_home_design$test,
  outcome_time = home_time_var,
  outcome_status = home_status_var,
  seed = 999,
  importance = TRUE,
  tuned = home_tuned
)

final_away_fit <- fit_predict_one_margin(
  train_df = final_away_design$train,
  test_df = final_away_design$test,
  outcome_time = away_time_var,
  outcome_status = away_status_var,
  seed = 1999,
  importance = TRUE,
  tuned = away_tuned
)

final_observed <- observed_match_outcomes(final_test_df)
final_probs <- derive_match_probabilities(final_home_fit$pred, final_away_fit$pred, minute_grid)

final_predictions_raw <- bind_cols(
  final_home_design$test %>% select(all_of(c(row_id_var, date_var, league_var, home_time_var, home_status_var))),
  final_away_design$test %>% select(all_of(c(away_time_var, away_status_var))),
  final_probs,
  final_observed %>% select(-all_of(row_id_var))
)

final_predictions_cal <- apply_calibrators(final_predictions_raw, calibrators)

final_test_summary <- tibble(
  season = final_test_season,
  home_cindex = final_home_fit$cindex,
  away_cindex = final_away_fit$cindex,
  home_ntree = final_home_fit$tuned$ntree,
  away_ntree = final_away_fit$tuned$ntree,
  home_nodesize = final_home_fit$tuned$nodesize,
  away_nodesize = final_away_fit$tuned$nodesize,
  home_mtry = final_home_fit$tuned$mtry,
  away_mtry = final_away_fit$tuned$mtry,
  home_nodedepth = ifelse(is.null(final_home_fit$tuned$nodedepth), NA, final_home_fit$tuned$nodedepth),
  away_nodedepth = ifelse(is.null(final_away_fit$tuned$nodedepth), NA, final_away_fit$tuned$nodedepth)
)

final_home_brier <- final_home_fit$brier_tbl
final_away_brier <- final_away_fit$brier_tbl

final_event_scores_raw <- bind_rows(lapply(calibration_targets, function(x) score_event_target(final_predictions_raw, x)))
final_event_scores_cal <- bind_rows(lapply(calibration_targets, function(x) compare_raw_calibrated(final_predictions_cal, x)))

final_home_vimp <- tibble(
  variable = names(final_home_fit$fit$importance),
  importance = as.numeric(final_home_fit$fit$importance)
) %>% arrange(desc(importance))

final_away_vimp <- tibble(
  variable = names(final_away_fit$fit$importance),
  importance = as.numeric(final_away_fit$fit$importance)
) %>% arrange(desc(importance))

write.csv(final_test_summary, file.path(output_dir, "final_test_summary.csv"), row.names = FALSE)
write.csv(final_home_brier, file.path(output_dir, "final_home_brier.csv"), row.names = FALSE)
write.csv(final_away_brier, file.path(output_dir, "final_away_brier.csv"), row.names = FALSE)
write.csv(final_predictions_raw, file.path(output_dir, "final_predictions_raw.csv"), row.names = FALSE)
write.csv(final_predictions_cal, file.path(output_dir, "final_predictions_calibrated.csv"), row.names = FALSE)
write.csv(final_event_scores_raw, file.path(output_dir, "final_event_scores_raw.csv"), row.names = FALSE)
write.csv(final_event_scores_cal, file.path(output_dir, "final_raw_vs_calibrated.csv"), row.names = FALSE)
write.csv(final_home_vimp, file.path(output_dir, "final_home_vimp.csv"), row.names = FALSE)
write.csv(final_away_vimp, file.path(output_dir, "final_away_vimp.csv"), row.names = FALSE)

saveRDS(final_home_fit$fit, file.path(output_dir, "final_home_model.rds"))
saveRDS(final_away_fit$fit, file.path(output_dir, "final_away_model.rds"))

# ----------------------------
# 16) Brier plot reported in the thesis
# ----------------------------
final_brier_long <- bind_rows(
  final_home_brier %>% mutate(margin = "Home"),
  final_away_brier %>% mutate(margin = "Away")
)

p_final_brier <- ggplot(final_brier_long, aes(x = time, y = brier, linetype = margin, group = margin)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Held-out IPCW Brier scores by margin",
    x = "Time (minutes)",
    y = "IPCW Brier score",
    linetype = "Margin"
  )

print(p_final_brier)
ggsave(file.path(output_dir, "brier_plot_RSF.png"), p_final_brier, width = 7, height = 5, dpi = 300)

# ----------------------------
# 17) Print thesis-relevant outputs
# ----------------------------
cat("\n===== Final model settings and C-index =====\n")
print(final_test_summary)

cat("\n===== Final marginal IPCW Brier scores: home =====\n")
print(final_home_brier)

cat("\n===== Final marginal IPCW Brier scores: away =====\n")
print(final_away_brier)

cat("\n===== Final event-probability scores: raw RSF =====\n")
print(final_event_scores_raw)

cat("\n===== Final raw versus calibrated scores =====\n")
print(final_event_scores_cal)

cat("\n===== Top 20 final home VIMP =====\n")
print(head(final_home_vimp, 20))

cat("\n===== Top 20 final away VIMP =====\n")
print(head(final_away_vimp, 20))

log_msg("Done. Outputs saved to: ", output_dir)
