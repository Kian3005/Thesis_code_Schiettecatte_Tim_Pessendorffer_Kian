# ============================================================
# Weibull AFT benchmark for first-goal prediction
# - Home and away marginal Weibull AFT models
# - Train on all seasons before 2024-2025
# - Evaluate on held-out season 2024-2025
# - Derive the same football-event probabilities as the RSF framework
# - Compute binary Brier score and log-loss
# ============================================================

# ----------------------------
# 0) Packages and seed
# ----------------------------
needed_pkgs <- c("dplyr", "tidyr", "readr", "tibble", "survival")

for (p in needed_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
invisible(lapply(needed_pkgs, library, character.only = TRUE))

set.seed(123)

# ----------------------------
# 1) Settings
# ----------------------------
data_path <- "SurvDataShortFormat.csv"
output_dir <- "weibull_aft_benchmark_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

final_test_season <- 2425

row_id_var <- "..row_id"
date_var <- "Date"
league_var <- "Round_grouped"
season_var <- "season"

home_time_var <- "home_first_goal_minute"
home_status_var <- "home_scored"
away_time_var <- "away_first_goal_minute"
away_status_var <- "away_scored"

minute_grid <- 0:90

# Predictors retained in the inferential Weibull AFT models.
home_weibull_vars <- c(
  "Home_Number_passes_for_roll5",
  "Away_Number_passes_for_roll5",
  "Home_xG_for_roll5",
  "Away_xG_against_roll5",
  "Away_goals_for_roll5",
  "Home_goals_for_roll5",
  "Home_shots_for_roll5",
  "Home_PSxG_15_for_roll5",
  "Away_xG_for_roll5",
  "Away_Shots_on_Target_against_roll5",
  "Home_PSxG_15_against_roll5",
  "Home_xG_15_against_roll5",
  "Home_Passing_accuracy_for_roll5",
  "home_formation_grouped",
  "Round_grouped"
)

away_weibull_vars <- c(
  "Away_Number_passes_for_roll5",
  "Away_xG_for_roll5",
  "Away_goals_for_roll5",
  "Away_Shots_on_Target_for_roll5",
  "Away_Corners_for_roll5",
  "Away_PSxG_90_for_roll5",
  "Away_Passing_accuracy_against_roll5",
  "Away_goals_against_roll5",
  "Home_xG_against_roll5",
  "Home_shots_against_roll5",
  "Home_Number_passes_for_roll5",
  "Home_xG_for_roll5",
  "Home_possesion_against_roll5",
  "Home_xG_15_for_roll5",
  "Home_xG_90_for_roll5",
  "Home_PSxG_15_against_roll5",
  "Home_PSxG_90_against_roll5",
  "Round_grouped"
)

evaluation_targets <- c(
  "p_home_first",
  "p_away_first",
  "p_home_goal_by_45",
  "p_away_goal_by_45",
  "p_any_goal_by_45",
  "p_no_goal_by_90",
  "p_any_goal_second_half_cond_no_goal_45"
)

# ----------------------------
# 2) Helper functions
# ----------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

safe_as_date <- function(x) {
  out <- as.Date(x)
  if (all(is.na(out))) out <- as.Date(readr::parse_date(as.character(x)))
  out
}

get_numeric_and_factor_vars <- function(vars) {
  list(
    numeric = vars[!grepl("grouped$", vars)],
    factor = vars[grepl("grouped$", vars)]
  )
}

make_standardization_params <- function(train_df, numeric_vars) {
  tibble(
    variable = numeric_vars,
    mean = sapply(numeric_vars, function(v) mean(train_df[[v]], na.rm = TRUE)),
    sd = sapply(numeric_vars, function(v) sd(train_df[[v]], na.rm = TRUE))
  ) %>%
    mutate(sd = ifelse(is.na(sd) | sd == 0, 1, sd))
}

apply_standardization <- function(df, standardization_params) {
  out <- df
  for (i in seq_len(nrow(standardization_params))) {
    v <- standardization_params$variable[i]
    out[[v]] <- (out[[v]] - standardization_params$mean[i]) / standardization_params$sd[i]
  }
  out
}

prepare_weibull_data <- function(train_df,
                                 test_df,
                                 outcome_time,
                                 outcome_status,
                                 predictor_vars) {
  var_types <- get_numeric_and_factor_vars(predictor_vars)

  keep_cols <- unique(c(row_id_var, outcome_time, outcome_status, predictor_vars))

  train_out <- train_df %>% select(all_of(keep_cols))
  test_out <- test_df %>% select(all_of(keep_cols))

  # Standardise numeric predictors using training-set means and SDs only.
  std_params <- make_standardization_params(train_out, var_types$numeric)
  train_out <- apply_standardization(train_out, std_params)
  test_out <- apply_standardization(test_out, std_params)

  # Align factor levels using the training set as reference.
  for (v in var_types$factor) {
    train_levels <- sort(unique(as.character(train_out[[v]])))
    train_out[[v]] <- factor(as.character(train_out[[v]]), levels = train_levels)
    test_out[[v]] <- factor(as.character(test_out[[v]]), levels = train_levels)
  }

  # Remove rows with missing outcome or predictor values.
  train_out <- train_out %>% filter(if_all(all_of(c(outcome_time, outcome_status, predictor_vars)), ~ !is.na(.)))
  test_out <- test_out %>% filter(if_all(all_of(c(outcome_time, outcome_status, predictor_vars)), ~ !is.na(.)))

  list(train = train_out, test = test_out, standardization = std_params)
}

predict_weibull_survival <- function(model, newdata, time_grid = 0:90) {
  lp <- predict(model, newdata = newdata, type = "lp")
  sigma <- model$scale

  surv_matrix <- sapply(time_grid, function(tt) {
    exp(-((tt / exp(lp))^(1 / sigma)))
  })

  surv_matrix <- as.matrix(surv_matrix)
  surv_matrix <- pmin(pmax(surv_matrix, 0), 1)

  list(
    survival = surv_matrix,
    time.interest = time_grid
  )
}

derive_match_probabilities <- function(home_pred, away_pred, minute_grid = 0:90) {
  Sh <- home_pred$survival
  Sa <- away_pred$survival

  idx45 <- which(minute_grid == 45)
  idx90 <- which(minute_grid == 90)

  Sh45 <- Sh[, idx45]
  Sa45 <- Sa[, idx45]
  Sh90 <- Sh[, idx90]
  Sa90 <- Sa[, idx90]

  # Minute-wise event probabilities on 1:90.
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
    p_home_first = p_home_first,
    p_away_first = p_away_first,
    p_home_goal_by_45 = p_home_goal_by_45,
    p_away_goal_by_45 = p_away_goal_by_45,
    p_any_goal_by_45 = p_any_goal_by_45,
    p_no_goal_by_90 = p_no_goal_by_90,
    p_any_goal_second_half_cond_no_goal_45 = p_any_goal_second_half_cond_no_goal_45
  )
}

observed_match_outcomes <- function(df) {
  h_time <- df[[home_time_var]]
  a_time <- df[[away_time_var]]
  h_sc <- df[[home_status_var]]
  a_sc <- df[[away_status_var]]

  home_goal_45 <- as.integer(h_sc == 1 & h_time <= 45)
  away_goal_45 <- as.integer(a_sc == 1 & a_time <= 45)
  any_goal_45 <- as.integer(home_goal_45 == 1 | away_goal_45 == 1)
  no_goal_45 <- 1L - any_goal_45

  home_first <- as.integer(h_sc == 1 & (a_sc == 0 | h_time < a_time))
  away_first <- as.integer(a_sc == 1 & (h_sc == 0 | a_time < h_time))

  no_goal_90 <- as.integer(h_sc == 0 & a_sc == 0)

  home_goal_second_half <- as.integer(h_sc == 1 & h_time > 45 & h_time <= 90)
  away_goal_second_half <- as.integer(a_sc == 1 & a_time > 45 & a_time <= 90)
  any_goal_second_half <- as.integer(home_goal_second_half == 1 | away_goal_second_half == 1)

  tibble(
    !!row_id_var := df[[row_id_var]],
    obs_home_first = home_first,
    obs_away_first = away_first,
    obs_home_goal_by_45 = home_goal_45,
    obs_away_goal_by_45 = away_goal_45,
    obs_any_goal_by_45 = any_goal_45,
    obs_no_goal_by_90 = no_goal_90,
    obs_any_goal_second_half_cond_no_goal_45 = ifelse(no_goal_45 == 1, any_goal_second_half, NA_integer_)
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

score_target <- function(df, p_name) {
  obs_name <- sub("^p_", "obs_", p_name)

  if (!(p_name %in% names(df)) || !(obs_name %in% names(df))) {
    return(tibble(
      target = p_name,
      event_rate = NA_real_,
      baseline_brier = NA_real_,
      weibull_brier = NA_real_,
      weibull_logloss = NA_real_,
      n = 0L
    ))
  }

  y <- df[[obs_name]]
  p <- df[[p_name]]
  idx <- !is.na(y) & !is.na(p)

  if (!any(idx)) {
    return(tibble(
      target = p_name,
      event_rate = NA_real_,
      baseline_brier = NA_real_,
      weibull_brier = NA_real_,
      weibull_logloss = NA_real_,
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
    weibull_brier = binary_brier(y, p),
    weibull_logloss = binary_logloss(y, p),
    n = length(y)
  )
}

# ----------------------------
# 3) Load and prepare data
# ----------------------------
log_msg("Loading data")

df <- read.csv(data_path, stringsAsFactors = FALSE) %>%
  mutate(
    !!date_var := safe_as_date(.data[[date_var]]),
    !!season_var := as.integer(.data[[season_var]]),
    !!home_time_var := as.numeric(.data[[home_time_var]]),
    !!away_time_var := as.numeric(.data[[away_time_var]]),
    !!home_status_var := as.integer(.data[[home_status_var]]),
    !!away_status_var := as.integer(.data[[away_status_var]]),
    !!row_id_var := row_number()
  ) %>%
  filter(
    !is.na(.data[[date_var]]),
    !is.na(.data[[season_var]]),
    is.finite(.data[[home_time_var]]),
    is.finite(.data[[away_time_var]]),
    .data[[home_status_var]] %in% c(0, 1),
    .data[[away_status_var]] %in% c(0, 1)
  )

required_predictors <- unique(c(home_weibull_vars, away_weibull_vars))
missing_predictors <- setdiff(required_predictors, names(df))
if (length(missing_predictors) > 0) {
  stop("Missing required predictors: ", paste(missing_predictors, collapse = ", "))
}

# Keep rows that are complete for all predictors needed for the benchmark.
df <- df %>% filter(if_all(all_of(required_predictors), ~ !is.na(.)))

train_df <- df %>% filter(.data[[season_var]] != final_test_season)
test_df <- df %>% filter(.data[[season_var]] == final_test_season)

log_msg("Training rows: ", nrow(train_df))
log_msg("Test rows: ", nrow(test_df))

# ----------------------------
# 4) Prepare home and away design matrices
# ----------------------------
log_msg("Preparing Weibull AFT design matrices")

home_data <- prepare_weibull_data(
  train_df = train_df,
  test_df = test_df,
  outcome_time = home_time_var,
  outcome_status = home_status_var,
  predictor_vars = home_weibull_vars
)

away_data <- prepare_weibull_data(
  train_df = train_df,
  test_df = test_df,
  outcome_time = away_time_var,
  outcome_status = away_status_var,
  predictor_vars = away_weibull_vars
)

if (!identical(home_data$test[[row_id_var]], away_data$test[[row_id_var]])) {
  stop("Home and away test sets are not aligned by row_id.")
}

# ----------------------------
# 5) Fit Weibull AFT models
# ----------------------------
log_msg("Fitting home and away Weibull AFT models")

home_formula <- as.formula(paste(
  "Surv(", home_time_var, ",", home_status_var, ") ~",
  paste(home_weibull_vars, collapse = " + ")
))

away_formula <- as.formula(paste(
  "Surv(", away_time_var, ",", away_status_var, ") ~",
  paste(away_weibull_vars, collapse = " + ")
))

weibull_home <- survival::survreg(home_formula, data = home_data$train, dist = "weibull")
weibull_away <- survival::survreg(away_formula, data = away_data$train, dist = "weibull")

# ----------------------------
# 6) Predict survival curves and derive football probabilities
# ----------------------------
log_msg("Predicting held-out survival curves and football-event probabilities")

weibull_home_pred <- predict_weibull_survival(weibull_home, home_data$test, minute_grid)
weibull_away_pred <- predict_weibull_survival(weibull_away, away_data$test, minute_grid)

weibull_probs <- derive_match_probabilities(weibull_home_pred, weibull_away_pred, minute_grid)

observed <- observed_match_outcomes(
  test_df %>%
    filter(.data[[row_id_var]] %in% home_data$test[[row_id_var]]) %>%
    arrange(match(.data[[row_id_var]], home_data$test[[row_id_var]]))
)

weibull_predictions <- bind_cols(
  home_data$test %>% select(all_of(c(row_id_var, home_time_var, home_status_var))),
  away_data$test %>% select(all_of(c(away_time_var, away_status_var))),
  observed %>% select(-all_of(row_id_var)),
  weibull_probs
)

# ----------------------------
# 7) Evaluate held-out event probabilities
# ----------------------------
log_msg("Computing held-out Brier scores and log-loss values")

weibull_scores <- bind_rows(lapply(evaluation_targets, function(x) score_target(weibull_predictions, x)))

# Thesis table order and labels.
weibull_scores_thesis <- weibull_scores %>%
  mutate(
    target_label = recode(
      target,
      p_home_first = "Home scores first",
      p_away_first = "Away scores first",
      p_home_goal_by_45 = "Home goal by 45",
      p_away_goal_by_45 = "Away goal by 45",
      p_any_goal_by_45 = "Any goal by 45",
      p_no_goal_by_90 = "No goal by 90",
      p_any_goal_second_half_cond_no_goal_45 = "Any goal in second half, conditional on 0--0 at halftime"
    )
  ) %>%
  select(target_label, event_rate, baseline_brier, weibull_brier, weibull_logloss, n)

# ----------------------------
# 8) Save outputs
# ----------------------------
write.csv(weibull_predictions, file.path(output_dir, "weibull_predictions.csv"), row.names = FALSE)
write.csv(weibull_scores, file.path(output_dir, "weibull_scores.csv"), row.names = FALSE)
write.csv(weibull_scores_thesis, file.path(output_dir, "weibull_scores_thesis_table.csv"), row.names = FALSE)
write.csv(home_data$standardization, file.path(output_dir, "home_standardization_params.csv"), row.names = FALSE)
write.csv(away_data$standardization, file.path(output_dir, "away_standardization_params.csv"), row.names = FALSE)

saveRDS(
  list(
    weibull_home = weibull_home,
    weibull_away = weibull_away,
    home_data = home_data,
    away_data = away_data
  ),
  file.path(output_dir, "weibull_aft_benchmark_objects.rds")
)

# ----------------------------
# 9) Print thesis-relevant outputs
# ----------------------------
cat("\n===== Weibull AFT benchmark scores =====\n")
print(weibull_scores_thesis, n = Inf)

log_msg("Done. Outputs saved to: ", output_dir)
