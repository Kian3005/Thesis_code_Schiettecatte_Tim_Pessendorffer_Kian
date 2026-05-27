# ============================================================
# Two-stage Cure Random Survival Forest prediction framework
# - Stage 1: binary Random Forest for scoring incidence
# - Stage 2: conditional event-time Random Survival Forest among scoring teams
# - Held-out evaluation on season 2024-2025
# - Development-fold tuning and OOF isotonic calibration diagnostic
# ============================================================

# ----------------------------
# 0) Packages and seed
# ----------------------------
needed_pkgs <- c(
  "dplyr", "tidyr", "readr", "tibble", "purrr",
  "survival", "ranger", "pROC"
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
output_dir <- "cure_rsf_output_clean"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

final_test_season <- 2425
cure_time_threshold <- 120

row_id_var <- "..row_id"
date_var <- "Date"
league_var <- "Round_grouped"
season_var <- "season"
league_feature_var <- "league_feature"

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

# Expanding-window development validation
n_folds <- 8
anchor_quantile_start <- 0.40
anchor_quantile_end <- 0.88
anchor_jitter_days <- 28
test_horizon_range_days <- c(21, 49)
min_gap_days <- 7
min_train_n <- 2500
min_test_n <- 150
min_test_leagues <- 4

# Evaluation grid
surv_eval_times <- c(15, 30, 45, 60, 75, 90)
minute_grid <- 0:90

# Hyperparameter grids used for the reported Cure RSF extension
incidence_ntree <- 600
incidence_nodesize_grid <- c(10, 15, 20)
incidence_mtry_mode_grid <- c("sqrt", "p3")

latency_ntree <- 400
latency_nodesize_grid <- c(20, 25, 30)
latency_mtry_mode_grid <- c("sqrt", "p4", "p3")
latency_nodedepth_grid <- list(8, 10, 12)

importance_mode_final <- TRUE

calibration_targets <- c(
  "p_home_scores",
  "p_away_scores",
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
  cat(sprintf("[%s] %s
", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

load_data_any <- function(path) {
  if (grepl("\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else if (grepl("\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file type. Use .csv or .rds.")
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
    "index", "id", "match_id",
    "home_formation", "away_formation", "Score", "HT_score", "notes"
  )
  nuisance <- nuisance[nuisance %in% names(df)]
  if (length(nuisance) > 0) df <- df %>% select(-all_of(nuisance))
  df
}

get_roll_vars <- function(df, pattern = "_roll") {
  names(df)[grepl(pattern, names(df))]
}

drop_zero_var_cols <- function(df) {
  if (ncol(df) == 0) return(list(data = df, dropped = character(0)))
  sds <- vapply(df, function(x) sd(x, na.rm = TRUE), numeric(1))
  keep <- !(is.na(sds) | sds == 0)
  list(data = df[, keep, drop = FALSE], dropped = names(df)[!keep])
}

find_high_corr_to_drop <- function(df_num, threshold = 0.95) {
  if (ncol(df_num) <= 1) return(character(0))

  C <- cor(as.matrix(df_num), use = "pairwise.complete.obs")
  C[is.na(C)] <- 0
  drop_vars <- character(0)

  repeat {
    C_work <- abs(C)
    diag(C_work) <- 0

    max_corr <- max(C_work)
    if (!is.finite(max_corr) || max_corr < threshold) break

    ij <- which(C_work == max_corr, arr.ind = TRUE)[1, ]
    var1 <- rownames(C)[ij[1]]
    var2 <- colnames(C)[ij[2]]

    mean_corr_1 <- mean(abs(C[var1, setdiff(colnames(C), var1)]), na.rm = TRUE)
    mean_corr_2 <- mean(abs(C[var2, setdiff(colnames(C), var2)]), na.rm = TRUE)
    drop_var <- if (mean_corr_1 >= mean_corr_2) var1 else var2

    drop_vars <- c(drop_vars, drop_var)
    keep <- setdiff(colnames(C), drop_vars)
    if (length(keep) <= 1) break
    C <- C[keep, keep, drop = FALSE]
  }

  unique(drop_vars)
}

align_factor_levels_safe <- function(train_df, test_df, factor_vars, new_level = "__NEW__") {
  for (v in factor_vars) {
    if (!(v %in% names(train_df)) || !(v %in% names(test_df))) next

    train_chr <- as.character(train_df[[v]])
    test_chr <- as.character(test_df[[v]])

    train_levels <- sort(unique(train_chr[!is.na(train_chr)]))
    test_chr[!is.na(test_chr) & !(test_chr %in% train_levels)] <- new_level

    all_levels <- c(train_levels, new_level)
    train_df[[v]] <- factor(train_chr, levels = all_levels)
    test_df[[v]] <- factor(test_chr, levels = all_levels)
  }

  list(train = train_df, test = test_df)
}

get_mtry_value <- function(p, mode = "sqrt") {
  if (mode == "sqrt") return(max(2, floor(sqrt(p))))
  if (mode == "p4") return(max(2, floor(p / 4)))
  if (mode == "p3") return(max(2, floor(p / 3)))
  if (mode == "p2") return(max(2, floor(p / 2)))
  stop("Unknown mtry mode: ", mode)
}

get_modal_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

# ----------------------------
# 3) Expanding-window folds
# ----------------------------
build_random_expanding_folds <- function(data,
                                         date_var,
                                         league_var,
                                         n_folds = 8,
                                         anchor_quantile_start = 0.40,
                                         anchor_quantile_end = 0.88,
                                         anchor_jitter_days = 28,
                                         test_horizon_range_days = c(21, 49),
                                         min_gap_days = 7,
                                         min_train_n = 2500,
                                         min_test_n = 150,
                                         min_test_leagues = 4,
                                         seed = 123) {
  set.seed(seed)

  dd <- data %>% arrange(.data[[date_var]])
  unique_dates <- sort(unique(dd[[date_var]]))
  unique_dates <- unique_dates[!is.na(unique_dates)]

  if (length(unique_dates) < n_folds + 20) {
    stop("Not enough unique dates to construct folds.")
  }

  probs <- seq(anchor_quantile_start, anchor_quantile_end, length.out = n_folds)
  target_dates <- as.Date(
    as.numeric(stats::quantile(unique_dates, probs = probs, type = 1)),
    origin = "1970-01-01"
  )

  folds <- vector("list", n_folds)
  last_test_end <- min(unique_dates) - 1

  for (i in seq_len(n_folds)) {
    valid_fold_found <- FALSE
    attempt <- 1

    while (!valid_fold_found && attempt <= 200) {
      anchor <- target_dates[i] + sample(-anchor_jitter_days:anchor_jitter_days, 1)
      horizon_days <- sample(seq(test_horizon_range_days[1], test_horizon_range_days[2]), 1)

      test_start <- max(anchor + min_gap_days, last_test_end + min_gap_days)
      test_end <- test_start + horizon_days

      train_idx <- which(dd[[date_var]] < test_start)
      test_idx <- which(dd[[date_var]] >= test_start & dd[[date_var]] <= test_end)
      n_test_leagues <- length(unique(dd[[league_var]][test_idx]))

      if (
        length(train_idx) >= min_train_n &&
          length(test_idx) >= min_test_n &&
          n_test_leagues >= min_test_leagues
      ) {
        folds[[i]] <- list(
          fold_id = i,
          train_idx = train_idx,
          test_idx = test_idx,
          test_start = test_start,
          test_end = test_end,
          n_train = length(train_idx),
          n_test = length(test_idx),
          n_test_leagues = n_test_leagues
        )
        last_test_end <- test_end
        valid_fold_found <- TRUE
      }

      attempt <- attempt + 1
    }
  }

  folds <- folds[!vapply(folds, is.null, logical(1))]
  if (length(folds) == 0) stop("No valid development folds could be created.")
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
                                corr_threshold = 0.95) {
  keep_cols <- unique(c(
    row_id_var, date_var, league_var,
    outcome_time, outcome_status,
    cat_vars, roll_vars
  ))

  train_df <- train_df %>% select(all_of(keep_cols))
  test_df <- test_df %>% select(all_of(keep_cols))

  zv <- drop_zero_var_cols(train_df[, roll_vars, drop = FALSE])
  roll_keep <- colnames(zv$data)

  high_corr_drop <- find_high_corr_to_drop(train_df[, roll_keep, drop = FALSE], corr_threshold)
  roll_keep <- setdiff(roll_keep, high_corr_drop)

  aligned_cat <- align_factor_levels_safe(
    train_df = train_df[, cat_vars, drop = FALSE],
    test_df = test_df[, cat_vars, drop = FALSE],
    factor_vars = cat_vars
  )

  meta_cols <- c(row_id_var, date_var, league_var, outcome_time, outcome_status)

  train_out <- bind_cols(
    train_df[, meta_cols, drop = FALSE],
    aligned_cat$train,
    train_df[, roll_keep, drop = FALSE]
  )

  test_out <- bind_cols(
    test_df[, meta_cols, drop = FALSE],
    aligned_cat$test,
    test_df[, roll_keep, drop = FALSE]
  )

  list(
    train = train_out,
    test = test_out,
    kept_roll_vars = roll_keep,
    dropped_zero_var = zv$dropped,
    dropped_high_corr = high_corr_drop
  )
}

# ----------------------------
# 5) Evaluation functions
# ----------------------------
surv_matrix_at_times <- function(pred_obj, eval_times) {
  surv_mat <- pred_obj$survival
  pred_times <- pred_obj$time.interest

  if (is.null(surv_mat) || is.null(pred_times)) {
    stop("Prediction object must contain survival and time.interest.")
  }

  out <- matrix(1, nrow = nrow(surv_mat), ncol = length(eval_times))
  colnames(out) <- paste0("t", eval_times)

  for (k in seq_along(eval_times)) {
    tt <- eval_times[k]
    idx <- suppressWarnings(max(which(pred_times <= tt)))
    if (is.finite(idx)) out[, k] <- surv_mat[, idx]
  }

  out
}

km_censor_fit <- function(time, status) {
  survival::survfit(survival::Surv(time, 1 - status) ~ 1)
}

km_surv_at <- function(km_fit, t) {
  s <- summary(km_fit, times = t, extend = TRUE)$surv
  if (length(s) == 0 || is.na(s)) return(1)
  pmax(s, 1e-6)
}

brier_ipcw_at_time <- function(time, status, surv_prob_t, eval_time) {
  km_cens <- km_censor_fit(time, status)
  G_t <- km_surv_at(km_cens, eval_time)

  obs_survival <- as.integer(time > eval_time)
  weights <- rep(0, length(time))

  event_by_t <- time <= eval_time & status == 1
  known_survival_past_t <- time > eval_time

  if (any(event_by_t)) {
    G_event <- vapply(
      time[event_by_t],
      function(tt) km_surv_at(km_cens, tt - 1e-8),
      numeric(1)
    )
    weights[event_by_t] <- 1 / G_event
  }

  if (any(known_survival_past_t)) {
    weights[known_survival_past_t] <- 1 / G_t
  }

  mean(weights * (obs_survival - surv_prob_t)^2, na.rm = TRUE)
}

margin_brier_summary <- function(test_df, pred_obj, outcome_time, outcome_status, eval_times) {
  surv_at_t <- surv_matrix_at_times(pred_obj, eval_times)

  tibble(
    time = eval_times,
    brier = vapply(seq_along(eval_times), function(k) {
      brier_ipcw_at_time(
        time = test_df[[outcome_time]],
        status = test_df[[outcome_status]],
        surv_prob_t = surv_at_t[, k],
        eval_time = eval_times[k]
      )
    }, numeric(1))
  )
}

calc_cindex <- function(time, status, risk) {
  survival::concordance(
    survival::Surv(time, status) ~ risk,
    reverse = TRUE
  )$concordance
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

score_binary_targets <- function(df, targets, model_col_suffix = "") {
  bind_rows(lapply(targets, function(target) {
    pred_col <- paste0(target, model_col_suffix)
    obs_col <- sub("^p_", "obs_", target)

    if (!(pred_col %in% names(df)) || !(obs_col %in% names(df))) {
      return(tibble(
        target = target,
        event_rate = NA_real_,
        baseline_brier = NA_real_,
        model_brier = NA_real_,
        model_logloss = NA_real_,
        n = 0L
      ))
    }

    y <- df[[obs_col]]
    p <- df[[pred_col]]
    idx <- !is.na(y) & !is.na(p)

    if (!any(idx)) {
      return(tibble(
        target = target,
        event_rate = NA_real_,
        baseline_brier = NA_real_,
        model_brier = NA_real_,
        model_logloss = NA_real_,
        n = 0L
      ))
    }

    y_use <- y[idx]
    p_use <- p[idx]
    p0 <- mean(y_use)

    tibble(
      target = target,
      event_rate = p0,
      baseline_brier = binary_brier(y_use, rep(p0, length(y_use))),
      model_brier = binary_brier(y_use, p_use),
      model_logloss = binary_logloss(y_use, p_use),
      n = length(y_use)
    )
  }))
}

# ----------------------------
# 6) Isotonic calibration diagnostic
# ----------------------------
fit_isotonic_calibrator <- function(p, y) {
  df <- tibble(p = p, y = y) %>%
    filter(!is.na(p), !is.na(y)) %>%
    arrange(p)

  if (nrow(df) < 20) return(NULL)

  iso <- isoreg(df$p, df$y)

  list(
    x = iso$x,
    yhat = iso$yf,
    fun = approxfun(iso$x, iso$yf, method = "linear", rule = 2, ties = mean)
  )
}

predict_isotonic_calibrator <- function(cal_obj, p) {
  if (is.null(cal_obj)) return(p)
  clip_prob(cal_obj$fun(clip_prob(p)))
}

fit_calibrators_from_oof <- function(oof_df, targets) {
  calibrators <- list()

  for (target in targets) {
    obs_col <- sub("^p_", "obs_", target)
    if (target %in% names(oof_df) && obs_col %in% names(oof_df)) {
      calibrators[[target]] <- fit_isotonic_calibrator(oof_df[[target]], oof_df[[obs_col]])
    }
  }

  calibrators
}

apply_calibrators <- function(df, calibrators) {
  out <- df

  for (target in names(calibrators)) {
    if (target %in% names(out)) {
      out[[paste0(target, "_cal")]] <- predict_isotonic_calibrator(calibrators[[target]], out[[target]])
    }
  }

  out
}

# ----------------------------
# 7) Stage 1: scoring incidence RF
# ----------------------------
fit_incidence_rf <- function(train_df,
                             outcome_status,
                             outcome_time,
                             ntree,
                             mtry_val,
                             nodesize_val,
                             importance_mode = FALSE,
                             seed = 123) {
  set.seed(seed)

  mod_df <- train_df %>% select(-all_of(c(row_id_var, date_var, league_var, outcome_time)))
  y <- factor(mod_df[[outcome_status]], levels = c(0, 1))
  X <- mod_df %>% select(-all_of(outcome_status))

  ranger::ranger(
    x = X,
    y = y,
    num.trees = ntree,
    mtry = mtry_val,
    min.node.size = nodesize_val,
    probability = TRUE,
    classification = TRUE,
    importance = if (importance_mode) "impurity" else "none",
    seed = seed
  )
}

tune_incidence_rf <- function(train_df,
                              outcome_status,
                              outcome_time,
                              ntree,
                              nodesize_grid,
                              mtry_mode_grid,
                              seed = 123) {
  p <- ncol(train_df) - 5

  candidates <- expand.grid(
    nodesize = nodesize_grid,
    mtry_mode = mtry_mode_grid,
    stringsAsFactors = FALSE
  )

  best <- NULL
  best_error <- Inf

  for (k in seq_len(nrow(candidates))) {
    ns <- candidates$nodesize[k]
    mtry <- get_mtry_value(p, candidates$mtry_mode[k])

    fit <- tryCatch(
      fit_incidence_rf(
        train_df = train_df,
        outcome_status = outcome_status,
        outcome_time = outcome_time,
        ntree = ntree,
        mtry_val = mtry,
        nodesize_val = ns,
        importance_mode = FALSE,
        seed = seed + k
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) next
    err <- fit$prediction.error
    if (!is.finite(err)) next

    if (err < best_error) {
      best_error <- err
      best <- list(nodesize = ns, mtry = mtry, mtry_mode = candidates$mtry_mode[k], oob_error = err)
    }
  }

  if (is.null(best)) stop("Incidence tuning failed.")
  best
}

predict_incidence_rf <- function(fit, newdata_df, outcome_time, outcome_status) {
  newx <- newdata_df %>% select(-all_of(c(row_id_var, date_var, league_var, outcome_time, outcome_status)))
  pred <- predict(fit, data = newx)$predictions

  if (!("1" %in% colnames(pred))) {
    stop("Incidence model prediction matrix has no class '1' column.")
  }

  as.numeric(pred[, "1"])
}

# ----------------------------
# 8) Stage 2: latency RSF among scoring teams
# ----------------------------
fit_latency_rsf <- function(train_df,
                            outcome_time,
                            outcome_status,
                            ntree,
                            nodesize_val,
                            mtry_val,
                            nodedepth_val,
                            importance_mode = FALSE,
                            seed = 123) {
  set.seed(seed)

  train_scorers <- train_df %>% filter(.data[[outcome_status]] == 1)
  if (nrow(train_scorers) < 50) stop("Too few scoring observations to fit latency model.")

  mod_df <- train_scorers %>% select(-all_of(c(row_id_var, date_var, league_var, outcome_status)))
  y_time <- mod_df[[outcome_time]]
  y_status <- rep(1, nrow(mod_df))
  X <- mod_df %>% select(-all_of(outcome_time))

  ranger::ranger(
    x = X,
    y = survival::Surv(y_time, y_status),
    num.trees = ntree,
    mtry = mtry_val,
    min.node.size = nodesize_val,
    max.depth = nodedepth_val,
    importance = if (importance_mode) "impurity" else "none",
    seed = seed
  )
}

tune_latency_rsf <- function(train_df,
                             outcome_time,
                             outcome_status,
                             ntree,
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
  best_error <- Inf

  for (k in seq_len(nrow(candidates))) {
    ns <- candidates$nodesize[k]
    mtry <- get_mtry_value(p, candidates$mtry_mode[k])
    depth <- nodedepth_grid[[candidates$depth_id[k]]]

    fit <- tryCatch(
      fit_latency_rsf(
        train_df = train_df,
        outcome_time = outcome_time,
        outcome_status = outcome_status,
        ntree = ntree,
        nodesize_val = ns,
        mtry_val = mtry,
        nodedepth_val = depth,
        importance_mode = FALSE,
        seed = seed + k
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) next
    err <- fit$prediction.error
    if (!is.finite(err)) next

    if (err < best_error) {
      best_error <- err
      best <- list(
        nodesize = ns,
        mtry = mtry,
        mtry_mode = candidates$mtry_mode[k],
        nodedepth = depth,
        ntree = ntree,
        prediction_error = err
      )
    }
  }

  if (is.null(best)) stop("Latency tuning failed.")
  best
}

predict_latency_rsf <- function(fit, newdata_df, outcome_time, outcome_status, time_grid = 0:90) {
  newx <- newdata_df %>% select(-all_of(c(row_id_var, date_var, league_var, outcome_time, outcome_status)))
  pred <- predict(fit, data = newx)

  surv_raw <- pred$survival
  raw_times <- pred$unique.death.times

  surv_interp <- matrix(1, nrow = nrow(surv_raw), ncol = length(time_grid))

  for (i in seq_len(nrow(surv_raw))) {
    surv_interp[i, ] <- approx(
      x = c(0, raw_times),
      y = c(1, surv_raw[i, ]),
      xout = time_grid,
      method = "constant",
      f = 0,
      rule = 2
    )$y

    surv_interp[i, ] <- pmin(pmax(surv_interp[i, ], 0), 1)
    surv_interp[i, ] <- cummin(surv_interp[i, ])
  }

  list(
    survival = surv_interp,
    chf = -log(pmax(surv_interp, 1e-10)),
    time.interest = time_grid
  )
}

# ----------------------------
# 9) Cure survival prediction
# ----------------------------
predict_cure_survival <- function(incidence_fit,
                                  latency_fit,
                                  newdata_df,
                                  outcome_time,
                                  outcome_status,
                                  time_grid = 0:90,
                                  latency_pred = NULL) {
  p_score <- predict_incidence_rf(
    fit = incidence_fit,
    newdata_df = newdata_df,
    outcome_time = outcome_time,
    outcome_status = outcome_status
  )

  if (is.null(latency_pred)) {
    latency_pred <- predict_latency_rsf(
      fit = latency_fit,
      newdata_df = newdata_df,
      outcome_time = outcome_time,
      outcome_status = outcome_status,
      time_grid = time_grid
    )
  }

  pi_cure <- 1 - p_score
  surv_cure <- sweep(latency_pred$survival, 1, 1 - pi_cure, "*")
  surv_cure <- sweep(surv_cure, 1, pi_cure, "+")
  surv_cure <- pmin(pmax(surv_cure, 0), 1)

  list(
    survival = surv_cure,
    chf = -log(pmax(surv_cure, 1e-10)),
    time.interest = latency_pred$time.interest,
    p_score = p_score,
    pi_cure = pi_cure,
    latency_survival = latency_pred$survival
  )
}

latency_diagnostics <- function(test_df, latency_pred, outcome_time, outcome_status) {
  scorer_idx <- test_df[[outcome_status]] == 1
  if (!any(scorer_idx)) {
    return(tibble(mae = NA_real_, rmse = NA_real_, n_scorers = 0L))
  }

  actual <- test_df[[outcome_time]][scorer_idx]
  surv_mat <- latency_pred$survival[scorer_idx, , drop = FALSE]
  times <- latency_pred$time.interest

  pred_median <- vapply(seq_len(nrow(surv_mat)), function(i) {
    below <- which(surv_mat[i, ] < 0.5)
    if (length(below) == 0) return(max(times))
    times[min(below)]
  }, numeric(1))

  tibble(
    mae = mean(abs(actual - pred_median), na.rm = TRUE),
    rmse = sqrt(mean((actual - pred_median)^2, na.rm = TRUE)),
    n_scorers = sum(scorer_idx)
  )
}

# ----------------------------
# 10) Match-level probabilities and observed outcomes
# ----------------------------
derive_match_probabilities_cure <- function(home_pred, away_pred, minute_grid = 0:90) {
  Sh <- surv_matrix_at_times(home_pred, minute_grid)
  Sa <- surv_matrix_at_times(away_pred, minute_grid)

  idx45 <- which(minute_grid == 45)
  idx90 <- which(minute_grid == 90)

  Sh45 <- Sh[, idx45]
  Sa45 <- Sa[, idx45]
  Sh90 <- Sh[, idx90]
  Sa90 <- Sa[, idx90]

  pH <- Sh[, -ncol(Sh), drop = FALSE] - Sh[, -1, drop = FALSE]
  pA <- Sa[, -ncol(Sa), drop = FALSE] - Sa[, -1, drop = FALSE]

  p_home_first <- rowSums(pH * Sa[, -1, drop = FALSE], na.rm = TRUE)
  p_away_first <- rowSums(pA * Sh[, -1, drop = FALSE], na.rm = TRUE)

  p_home_second_half <- pmax(0, Sh45 - Sh90)
  p_away_second_half <- pmax(0, Sa45 - Sa90)

  tibble(
    p_home_scores = home_pred$p_score,
    p_away_scores = away_pred$p_score,
    p_home_first = p_home_first,
    p_away_first = p_away_first,
    p_home_goal_by_45 = 1 - Sh45,
    p_away_goal_by_45 = 1 - Sa45,
    p_no_goal_by_45 = Sh45 * Sa45,
    p_any_goal_by_45 = 1 - Sh45 * Sa45,
    p_no_goal_by_90 = Sh90 * Sa90,
    p_any_goal_by_90 = 1 - Sh90 * Sa90,
    p_home_goal_second_half = p_home_second_half,
    p_away_goal_second_half = p_away_second_half,
    p_any_goal_second_half = 1 - (1 - p_home_second_half) * (1 - p_away_second_half),
    p_any_goal_second_half_cond_no_goal_45 = 1 - (Sh90 / pmax(Sh45, 1e-8)) * (Sa90 / pmax(Sa45, 1e-8))
  )
}

observed_match_outcomes_cure <- function(df) {
  h_time <- df[[home_time_var]]
  a_time <- df[[away_time_var]]
  h_sc <- df[[home_status_var]]
  a_sc <- df[[away_status_var]]

  h_goal_45 <- as.integer(h_sc == 1 & h_time <= 45)
  a_goal_45 <- as.integer(a_sc == 1 & a_time <= 45)
  any_goal_45 <- as.integer(h_goal_45 == 1 | a_goal_45 == 1)
  no_goal_45 <- 1L - any_goal_45

  home_first <- as.integer(h_sc == 1 & (a_sc == 0 | h_time < a_time))
  away_first <- as.integer(a_sc == 1 & (h_sc == 0 | a_time < h_time))

  home_second_half <- as.integer(h_sc == 1 & h_time > 45 & h_time <= 90)
  away_second_half <- as.integer(a_sc == 1 & a_time > 45 & a_time <= 90)
  any_second_half <- as.integer(home_second_half == 1 | away_second_half == 1)

  tibble(
    !!row_id_var := df[[row_id_var]],
    obs_home_scores = h_sc,
    obs_away_scores = a_sc,
    obs_home_first = home_first,
    obs_away_first = away_first,
    obs_home_goal_by_45 = h_goal_45,
    obs_away_goal_by_45 = a_goal_45,
    obs_no_goal_by_45 = no_goal_45,
    obs_any_goal_by_45 = any_goal_45,
    obs_no_goal_by_90 = as.integer(h_sc == 0 & a_sc == 0),
    obs_any_goal_by_90 = as.integer(h_sc == 1 | a_sc == 1),
    obs_home_goal_second_half = home_second_half,
    obs_away_goal_second_half = away_second_half,
    obs_any_goal_second_half = any_second_half,
    obs_any_goal_second_half_cond_no_goal_45 = ifelse(no_goal_45 == 1, any_second_half, NA_integer_)
  )
}

# ----------------------------
# 11) Fit/predict one Cure RSF margin
# ----------------------------
fit_predict_cure_margin <- function(train_df,
                                    test_df,
                                    outcome_time,
                                    outcome_status,
                                    incidence_settings = NULL,
                                    latency_settings = NULL,
                                    tune = TRUE,
                                    importance_mode = FALSE,
                                    seed = 123) {
  if (tune) {
    incidence_settings <- tune_incidence_rf(
      train_df = train_df,
      outcome_status = outcome_status,
      outcome_time = outcome_time,
      ntree = incidence_ntree,
      nodesize_grid = incidence_nodesize_grid,
      mtry_mode_grid = incidence_mtry_mode_grid,
      seed = seed
    )

    latency_settings <- tune_latency_rsf(
      train_df = train_df,
      outcome_time = outcome_time,
      outcome_status = outcome_status,
      ntree = latency_ntree,
      nodesize_grid = latency_nodesize_grid,
      mtry_mode_grid = latency_mtry_mode_grid,
      nodedepth_grid = latency_nodedepth_grid,
      seed = seed
    )
  }

  incidence_fit <- fit_incidence_rf(
    train_df = train_df,
    outcome_status = outcome_status,
    outcome_time = outcome_time,
    ntree = incidence_ntree,
    mtry_val = incidence_settings$mtry,
    nodesize_val = incidence_settings$nodesize,
    importance_mode = importance_mode,
    seed = seed
  )

  latency_fit <- fit_latency_rsf(
    train_df = train_df,
    outcome_time = outcome_time,
    outcome_status = outcome_status,
    ntree = latency_settings$ntree,
    nodesize_val = latency_settings$nodesize,
    mtry_val = latency_settings$mtry,
    nodedepth_val = latency_settings$nodedepth,
    importance_mode = importance_mode,
    seed = seed
  )

  latency_pred <- predict_latency_rsf(
    fit = latency_fit,
    newdata_df = test_df,
    outcome_time = outcome_time,
    outcome_status = outcome_status,
    time_grid = minute_grid
  )

  cure_pred <- predict_cure_survival(
    incidence_fit = incidence_fit,
    latency_fit = latency_fit,
    newdata_df = test_df,
    outcome_time = outcome_time,
    outcome_status = outcome_status,
    time_grid = minute_grid,
    latency_pred = latency_pred
  )

  risk <- cure_pred$chf[, ncol(cure_pred$chf)]

  latency_diag <- latency_diagnostics(
    test_df = test_df,
    latency_pred = latency_pred,
    outcome_time = outcome_time,
    outcome_status = outcome_status
  )

  list(
    incidence_fit = incidence_fit,
    latency_fit = latency_fit,
    incidence_settings = incidence_settings,
    latency_settings = latency_settings,
    pred = cure_pred,
    test_df = test_df,
    cindex = calc_cindex(test_df[[outcome_time]], test_df[[outcome_status]], risk),
    brier_tbl = margin_brier_summary(test_df, cure_pred, outcome_time, outcome_status, surv_eval_times),
    incidence_auc = as.numeric(pROC::auc(test_df[[outcome_status]], cure_pred$p_score, quiet = TRUE)),
    pi_cure_mean = mean(cure_pred$pi_cure),
    latency_mae = latency_diag$mae,
    latency_rmse = latency_diag$rmse,
    latency_n_scorers = latency_diag$n_scorers
  )
}

# ----------------------------
# 12) Development fold runner
# ----------------------------
run_one_fold_cure <- function(fold, dev_df, roll_vars, cat_vars) {
  log_msg("Running development fold ", fold$fold_id)

  train_fold <- dev_df[fold$train_idx, , drop = FALSE] %>% arrange(.data[[date_var]])
  test_fold <- dev_df[fold$test_idx, , drop = FALSE] %>% arrange(.data[[date_var]])

  home_design <- prepare_fold_design(
    train_df = train_fold,
    test_df = test_fold,
    roll_vars = roll_vars,
    cat_vars = cat_vars,
    outcome_time = home_time_var,
    outcome_status = home_status_var,
    corr_threshold = corr_threshold
  )

  away_design <- prepare_fold_design(
    train_df = train_fold,
    test_df = test_fold,
    roll_vars = roll_vars,
    cat_vars = cat_vars,
    outcome_time = away_time_var,
    outcome_status = away_status_var,
    corr_threshold = corr_threshold
  )

  if (!identical(home_design$test[[row_id_var]], away_design$test[[row_id_var]])) {
    stop("Home and away fold test sets are not aligned.")
  }

  home_fit <- fit_predict_cure_margin(
    train_df = home_design$train,
    test_df = home_design$test,
    outcome_time = home_time_var,
    outcome_status = home_status_var,
    tune = TRUE,
    importance_mode = FALSE,
    seed = 100 + fold$fold_id
  )

  away_fit <- fit_predict_cure_margin(
    train_df = away_design$train,
    test_df = away_design$test,
    outcome_time = away_time_var,
    outcome_status = away_status_var,
    tune = TRUE,
    importance_mode = FALSE,
    seed = 200 + fold$fold_id
  )

  observed <- observed_match_outcomes_cure(test_fold)
  derived_probs <- derive_match_probabilities_cure(home_fit$pred, away_fit$pred, minute_grid)

  pred_tbl <- bind_cols(
    test_fold %>% select(all_of(c(row_id_var, date_var, league_var))),
    observed %>% select(-all_of(row_id_var)),
    derived_probs
  )

  derived_scores <- score_binary_targets(pred_tbl, calibration_targets)

  fold_summary <- list(
    fold_id = fold$fold_id,
    n_train = fold$n_train,
    n_test = fold$n_test,
    home_cindex = home_fit$cindex,
    away_cindex = away_fit$cindex,
    home_incidence_auc = home_fit$incidence_auc,
    away_incidence_auc = away_fit$incidence_auc,
    home_pi_cure = home_fit$pi_cure_mean,
    away_pi_cure = away_fit$pi_cure_mean,
    home_latency_mae = home_fit$latency_mae,
    home_latency_rmse = home_fit$latency_rmse,
    away_latency_mae = away_fit$latency_mae,
    away_latency_rmse = away_fit$latency_rmse,
    home_incidence_settings = home_fit$incidence_settings,
    away_incidence_settings = away_fit$incidence_settings,
    home_latency_settings = home_fit$latency_settings,
    away_latency_settings = away_fit$latency_settings,
    home_brier_tbl = home_fit$brier_tbl,
    away_brier_tbl = away_fit$brier_tbl,
    derived_scores = derived_scores
  )

  fold_dir <- file.path(output_dir, paste0("dev_fold_", fold$fold_id))
  dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(fold_summary, file.path(fold_dir, "fold_summary_light.rds"))
  write.csv(pred_tbl, file.path(fold_dir, "predictions_oof.csv"), row.names = FALSE)
  write.csv(home_fit$brier_tbl, file.path(fold_dir, "home_margin_brier.csv"), row.names = FALSE)
  write.csv(away_fit$brier_tbl, file.path(fold_dir, "away_margin_brier.csv"), row.names = FALSE)
  write.csv(derived_scores, file.path(fold_dir, "derived_scores.csv"), row.names = FALSE)

  list(fold_summary = fold_summary, predictions = pred_tbl)
}

# ----------------------------
# 13) Hyperparameter locking
# ----------------------------
lock_incidence_settings <- function(fold_results, margin = c("home", "away")) {
  margin <- match.arg(margin)
  param_name <- paste0(margin, "_incidence_settings")
  params <- lapply(fold_results, function(x) x$fold_summary[[param_name]])

  list(
    nodesize = as.numeric(get_modal_value(sapply(params, `[[`, "nodesize"))),
    mtry = as.numeric(get_modal_value(sapply(params, `[[`, "mtry")))
  )
}

lock_latency_settings <- function(fold_results, margin = c("home", "away")) {
  margin <- match.arg(margin)
  param_name <- paste0(margin, "_latency_settings")
  params <- lapply(fold_results, function(x) x$fold_summary[[param_name]])

  list(
    nodesize = as.numeric(get_modal_value(sapply(params, `[[`, "nodesize"))),
    mtry = as.numeric(get_modal_value(sapply(params, `[[`, "mtry"))),
    nodedepth = as.numeric(get_modal_value(sapply(params, function(x) x$nodedepth))),
    ntree = as.numeric(get_modal_value(sapply(params, `[[`, "ntree")))
  )
}

# ----------------------------
# 14) Load and prepare data
# ----------------------------
log_msg("Loading and preparing data")

df <- load_data_any(data_path) %>%
  drop_nuisance_cols() %>%
  mutate(
    !!date_var := safe_as_date(.data[[date_var]]),
    !!season_var := as.integer(.data[[season_var]]),
    !!league_feature_var := .data[[league_var]],
    !!home_time_var := as.numeric(.data[[home_time_var]]),
    !!away_time_var := as.numeric(.data[[away_time_var]]),
    !!home_status_var := as.integer(.data[[home_status_var]]),
    !!away_status_var := as.integer(.data[[away_status_var]]),
    !!row_id_var := row_number()
  )

required_cols <- c(
  date_var, league_var, season_var,
  home_time_var, home_status_var,
  away_time_var, away_status_var
)
missing_required <- setdiff(required_cols, names(df))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

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
    .data[[home_status_var]] %in% c(0, 1),
    .data[[away_status_var]] %in% c(0, 1)
  ) %>%
  filter(if_all(all_of(c(cat_vars, roll_vars)), ~ !is.na(.))) %>%
  mutate(
    !!home_time_var := ifelse(.data[[home_status_var]] == 0, cure_time_threshold, .data[[home_time_var]]),
    !!away_time_var := ifelse(.data[[away_status_var]] == 0, cure_time_threshold, .data[[away_time_var]])
  )

log_msg("Rows after cleaning: ", nrow(analysis_df))

if (!(final_test_season %in% unique(analysis_df[[season_var]]))) {
  stop("Final test season not found in data.")
}

dev_df <- analysis_df %>% filter(.data[[season_var]] != final_test_season)
final_test_df <- analysis_df %>% filter(.data[[season_var]] == final_test_season)

log_msg("Development rows: ", nrow(dev_df))
log_msg("Held-out rows: ", nrow(final_test_df))

# ----------------------------
# 15) Development folds
# ----------------------------
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
    test_start = x$test_start,
    test_end = x$test_end,
    n_train = x$n_train,
    n_test = x$n_test,
    n_test_leagues = x$n_test_leagues
  )
}))
write.csv(fold_overview, file.path(output_dir, "fold_overview.csv"), row.names = FALSE)

fold_results <- vector("list", length(folds))
all_oof_predictions <- vector("list", length(folds))

for (i in seq_along(folds)) {
  fold_results[[i]] <- run_one_fold_cure(folds[[i]], dev_df, roll_vars, cat_vars)
  all_oof_predictions[[i]] <- fold_results[[i]]$predictions
  invisible(gc())
}

oof_df <- bind_rows(all_oof_predictions)
write.csv(oof_df, file.path(output_dir, "oof_predictions_raw.csv"), row.names = FALSE)

fold_summary_df <- bind_rows(lapply(fold_results, function(x) {
  fs <- x$fold_summary
  tibble(
    fold_id = fs$fold_id,
    n_train = fs$n_train,
    n_test = fs$n_test,
    home_cindex = fs$home_cindex,
    away_cindex = fs$away_cindex,
    home_incidence_auc = fs$home_incidence_auc,
    away_incidence_auc = fs$away_incidence_auc,
    home_pi_cure = fs$home_pi_cure,
    away_pi_cure = fs$away_pi_cure,
    home_latency_mae = fs$home_latency_mae,
    home_latency_rmse = fs$home_latency_rmse,
    away_latency_mae = fs$away_latency_mae,
    away_latency_rmse = fs$away_latency_rmse
  )
}))
write.csv(fold_summary_df, file.path(output_dir, "fold_summary.csv"), row.names = FALSE)

home_brier_dev <- bind_rows(lapply(fold_results, function(x) {
  x$fold_summary$home_brier_tbl %>% mutate(fold_id = x$fold_summary$fold_id)
}))
away_brier_dev <- bind_rows(lapply(fold_results, function(x) {
  x$fold_summary$away_brier_tbl %>% mutate(fold_id = x$fold_summary$fold_id)
}))
write.csv(home_brier_dev, file.path(output_dir, "home_dev_brier.csv"), row.names = FALSE)
write.csv(away_brier_dev, file.path(output_dir, "away_dev_brier.csv"), row.names = FALSE)

calibrators <- fit_calibrators_from_oof(oof_df, calibration_targets)
oof_df_cal <- apply_calibrators(oof_df, calibrators)
write.csv(oof_df_cal, file.path(output_dir, "oof_predictions_calibrated.csv"), row.names = FALSE)
saveRDS(calibrators, file.path(output_dir, "isotonic_calibrators.rds"))

# ----------------------------
# 16) Final held-out fit with locked CV hyperparameters
# ----------------------------
log_msg("Fitting final Cure RSF models on all development seasons")

locked_home_incidence <- lock_incidence_settings(fold_results, "home")
locked_away_incidence <- lock_incidence_settings(fold_results, "away")
locked_home_latency <- lock_latency_settings(fold_results, "home")
locked_away_latency <- lock_latency_settings(fold_results, "away")

final_home_design <- prepare_fold_design(
  train_df = dev_df,
  test_df = final_test_df,
  roll_vars = roll_vars,
  cat_vars = cat_vars,
  outcome_time = home_time_var,
  outcome_status = home_status_var,
  corr_threshold = corr_threshold
)

final_away_design <- prepare_fold_design(
  train_df = dev_df,
  test_df = final_test_df,
  roll_vars = roll_vars,
  cat_vars = cat_vars,
  outcome_time = away_time_var,
  outcome_status = away_status_var,
  corr_threshold = corr_threshold
)

if (!identical(final_home_design$test[[row_id_var]], final_away_design$test[[row_id_var]])) {
  stop("Final home and away test sets are not aligned.")
}

final_home_fit <- fit_predict_cure_margin(
  train_df = final_home_design$train,
  test_df = final_home_design$test,
  outcome_time = home_time_var,
  outcome_status = home_status_var,
  incidence_settings = locked_home_incidence,
  latency_settings = locked_home_latency,
  tune = FALSE,
  importance_mode = importance_mode_final,
  seed = 999
)

final_away_fit <- fit_predict_cure_margin(
  train_df = final_away_design$train,
  test_df = final_away_design$test,
  outcome_time = away_time_var,
  outcome_status = away_status_var,
  incidence_settings = locked_away_incidence,
  latency_settings = locked_away_latency,
  tune = FALSE,
  importance_mode = importance_mode_final,
  seed = 1999
)

final_test_aligned <- final_test_df %>%
  filter(.data[[row_id_var]] %in% final_home_design$test[[row_id_var]]) %>%
  arrange(match(.data[[row_id_var]], final_home_design$test[[row_id_var]]))

final_observed <- observed_match_outcomes_cure(final_test_aligned)
final_probs <- derive_match_probabilities_cure(final_home_fit$pred, final_away_fit$pred, minute_grid)

final_predictions <- bind_cols(
  final_test_aligned %>% select(all_of(c(row_id_var, date_var, league_var))),
  final_observed %>% select(-all_of(row_id_var)),
  final_probs
)

final_predictions_cal <- apply_calibrators(final_predictions, calibrators)

# ----------------------------
# 17) Thesis-relevant output tables
# ----------------------------
final_marginal_summary <- tibble(
  metric = c(
    "Held-out C-index",
    "Incidence AUC",
    "Estimated cure fraction",
    "Observed no-goal proportion"
  ),
  home_margin = c(
    final_home_fit$cindex,
    final_home_fit$incidence_auc,
    final_home_fit$pi_cure_mean,
    mean(final_test_aligned[[home_status_var]] == 0)
  ),
  away_margin = c(
    final_away_fit$cindex,
    final_away_fit$incidence_auc,
    final_away_fit$pi_cure_mean,
    mean(final_test_aligned[[away_status_var]] == 0)
  )
)

final_event_scores_raw <- score_binary_targets(final_predictions, calibration_targets) %>%
  mutate(
    target_label = recode(
      target,
      p_home_scores = "Home scores",
      p_away_scores = "Away scores",
      p_home_first = "Home scores first",
      p_away_first = "Away scores first",
      p_home_goal_by_45 = "Home goal by 45",
      p_away_goal_by_45 = "Away goal by 45",
      p_any_goal_by_45 = "Any goal by 45",
      p_no_goal_by_90 = "No goal by 90",
      p_any_goal_second_half_cond_no_goal_45 = "Any goal in second half, conditional on 0--0 at halftime"
    )
  ) %>%
  select(target_label, event_rate, baseline_brier, cure_rsf_brier = model_brier, cure_rsf_logloss = model_logloss, n)

final_event_scores_calibrated <- score_binary_targets(
  final_predictions_cal,
  calibration_targets,
  model_col_suffix = "_cal"
) %>%
  mutate(
    target_label = recode(
      target,
      p_home_scores = "Home scores",
      p_away_scores = "Away scores",
      p_home_first = "Home scores first",
      p_away_first = "Away scores first",
      p_home_goal_by_45 = "Home goal by 45",
      p_away_goal_by_45 = "Away goal by 45",
      p_any_goal_by_45 = "Any goal by 45",
      p_no_goal_by_90 = "No goal by 90",
      p_any_goal_second_half_cond_no_goal_45 = "Any goal in second half, conditional on 0--0 at halftime"
    )
  ) %>%
  select(target_label, event_rate, baseline_brier, cure_rsf_brier_cal = model_brier, cure_rsf_logloss_cal = model_logloss, n)

final_home_brier <- final_home_fit$brier_tbl
final_away_brier <- final_away_fit$brier_tbl

final_home_incidence_vimp <- tibble(
  variable = names(final_home_fit$incidence_fit$variable.importance),
  importance = as.numeric(final_home_fit$incidence_fit$variable.importance)
) %>% arrange(desc(importance))

final_away_incidence_vimp <- tibble(
  variable = names(final_away_fit$incidence_fit$variable.importance),
  importance = as.numeric(final_away_fit$incidence_fit$variable.importance)
) %>% arrange(desc(importance))

final_home_latency_vimp <- tibble(
  variable = names(final_home_fit$latency_fit$variable.importance),
  importance = as.numeric(final_home_fit$latency_fit$variable.importance)
) %>% arrange(desc(importance))

final_away_latency_vimp <- tibble(
  variable = names(final_away_fit$latency_fit$variable.importance),
  importance = as.numeric(final_away_fit$latency_fit$variable.importance)
) %>% arrange(desc(importance))

# ----------------------------
# 18) Save outputs
# ----------------------------
write.csv(final_predictions, file.path(output_dir, "final_predictions_raw.csv"), row.names = FALSE)
write.csv(final_predictions_cal, file.path(output_dir, "final_predictions_calibrated.csv"), row.names = FALSE)
write.csv(final_marginal_summary, file.path(output_dir, "final_marginal_summary.csv"), row.names = FALSE)
write.csv(final_event_scores_raw, file.path(output_dir, "final_event_scores_raw.csv"), row.names = FALSE)
write.csv(final_event_scores_calibrated, file.path(output_dir, "final_event_scores_calibrated.csv"), row.names = FALSE)
write.csv(final_home_brier, file.path(output_dir, "final_home_brier.csv"), row.names = FALSE)
write.csv(final_away_brier, file.path(output_dir, "final_away_brier.csv"), row.names = FALSE)
write.csv(final_home_incidence_vimp, file.path(output_dir, "final_home_incidence_vimp.csv"), row.names = FALSE)
write.csv(final_away_incidence_vimp, file.path(output_dir, "final_away_incidence_vimp.csv"), row.names = FALSE)
write.csv(final_home_latency_vimp, file.path(output_dir, "final_home_latency_vimp.csv"), row.names = FALSE)
write.csv(final_away_latency_vimp, file.path(output_dir, "final_away_latency_vimp.csv"), row.names = FALSE)

saveRDS(
  list(
    final_home_fit = final_home_fit,
    final_away_fit = final_away_fit,
    locked_home_incidence = locked_home_incidence,
    locked_away_incidence = locked_away_incidence,
    locked_home_latency = locked_home_latency,
    locked_away_latency = locked_away_latency,
    final_home_design = final_home_design,
    final_away_design = final_away_design,
    calibrators = calibrators
  ),
  file.path(output_dir, "cure_rsf_objects.rds")
)

cat("
===== Final marginal Cure RSF summary =====
")
print(final_marginal_summary)

cat("
===== Final Cure RSF event scores, raw probabilities =====
")
print(final_event_scores_raw, n = Inf)

cat("
===== Final Cure RSF event scores, calibrated diagnostic =====
")
print(final_event_scores_calibrated, n = Inf)

log_msg("Done. Outputs saved to: ", output_dir)
