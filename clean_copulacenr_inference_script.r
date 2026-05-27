# ============================================================
# CopulaCenR one-stage copula survival model
# - Pooled focal-team representation with two rows per match
# - Home-team row: focal team = home team
# - Away-team row: focal team = away team
# - Outcome: focal team's own first-goal time, censored at 90 if no goal
# - Model grid over copula and marginal distributions
# - Final coefficient, dependence, and AIC summary tables
# ============================================================

# ----------------------------
# 0) Packages
# ----------------------------
needed_pkgs <- c("readr", "dplyr", "tibble", "CopulaCenR")

for (p in needed_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
invisible(lapply(needed_pkgs, library, character.only = TRUE))

# ----------------------------
# 1) Settings
# ----------------------------
input_path <- "df_final_kian_matchup_long.csv"
output_prefix <- "copulacenr_pooled_team_opp_home_indicator"

horizon_min <- 90

fit_method <- "BFGS"
fit_iter <- 200
fit_stepsize <- 1e-5

boundary_tol <- 1e-3
pd_tol <- 1e-10

copulas_to_fit <- c("Frank", "Clayton", "Gumbel", "Joe")
margins_to_fit <- c("Weibull", "Gompertz", "Loglogistic")

continuous_covariates <- c(
  "team_Number_passes_for_roll5",
  "opp_Number_passes_for_roll5",
  "team_xG_for_roll5",
  "opp_xG_for_roll5",
  "team_xG_against_roll5",
  "opp_xG_against_roll5",
  "opp_goals_for_roll5",
  "team_goals_for_roll5",
  "team_Passing_accuracy_for_roll5",
  "opp_Passing_accuracy_for_roll5",
  "opp_PSxG_90_for_roll5",
  "opp_PSxG_90_against_roll5",
  "team_PSxG_90_for_roll5",
  "team_PSxG_90_against_roll5",
  "is_home"
)

categorical_covariates <- c("Round_grouped")

# ----------------------------
# 2) General helpers
# ----------------------------
to_binary01 <- function(x, var_name = "variable") {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))

  x_chr <- trimws(tolower(as.character(x)))
  out <- dplyr::case_when(
    x_chr %in% c("1", "true", "t", "yes", "y", "home") ~ 1L,
    x_chr %in% c("0", "false", "f", "no", "n", "away") ~ 0L,
    is.na(x_chr) | x_chr == "" ~ NA_integer_,
    TRUE ~ NA_integer_
  )

  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    stop(
      "Could not convert ", var_name, " to 0/1. Problem values: ",
      paste(unique(as.character(x[bad])), collapse = ", ")
    )
  }

  out
}

same_or_all_na <- function(x) {
  ux <- unique(x)
  ux <- ux[!is.na(ux)]
  length(ux) <= 1
}

stars_from_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

clean_param_label <- function(x) {
  x %>%
    gsub("^Round_grouped", "", .) %>%
    gsub("_", " ", .) %>%
    trimws()
}

# ----------------------------
# 3) Copula/dependence helpers
# ----------------------------
debye1 <- function(theta) {
  if (!is.finite(theta) || abs(theta) < 1e-8) return(1)

  integrate(
    f = function(t) ifelse(abs(t) < 1e-8, 1, t / (exp(t) - 1)),
    lower = 0,
    upper = theta,
    subdivisions = 1000,
    rel.tol = 1e-8
  )$value / theta
}

effective_eta <- function(copula, eta, tol = boundary_tol) {
  if (!is.finite(eta)) return(NA_real_)

  if (copula %in% c("Gumbel", "Joe")) {
    if (eta < 1 && eta >= 1 - tol) return(1)
    return(eta)
  }

  if (copula == "Clayton") {
    if (eta < 0 && eta >= -tol) return(0)
    return(eta)
  }

  eta
}

is_admissible_eta <- function(copula, eta, tol = boundary_tol) {
  if (!is.finite(eta)) return(FALSE)

  if (copula == "Frank") return(TRUE)
  if (copula == "Clayton") return(eta >= -tol)
  if (copula %in% c("Gumbel", "Joe")) return(eta >= 1 - tol)

  FALSE
}

kendall_tau <- function(copula, eta) {
  if (!is.finite(eta)) return(NA_real_)
  eta_eff <- effective_eta(copula, eta)

  if (copula == "Frank") {
    if (abs(eta_eff) < 1e-8) return(0)
    return(1 + 4 * (debye1(eta_eff) - 1) / eta_eff)
  }

  if (copula == "Clayton") {
    if (eta_eff <= 0) return(0)
    return(eta_eff / (eta_eff + 2))
  }

  if (copula == "Gumbel") {
    if (eta_eff <= 1) return(0)
    return(1 - 1 / eta_eff)
  }

  if (copula == "Joe") {
    if (eta_eff <= 1) return(0)

    out <- tryCatch({
      integrand <- function(t) {
        A <- 1 - (1 - t)^eta_eff
        phi <- -log(A)
        dphi <- -eta_eff * (1 - t)^(eta_eff - 1) / A
        phi / dphi
      }
      1 + 4 * integrate(
        integrand,
        lower = 1e-8,
        upper = 1 - 1e-8,
        subdivisions = 1000,
        rel.tol = 1e-8
      )$value
    }, error = function(e) NA_real_)

    return(out)
  }

  NA_real_
}

independence_test <- function(copula, eta, se) {
  if (!is.finite(eta) || !is.finite(se) || se <= 0) {
    return(c(null_eta = NA_real_, z = NA_real_, p = NA_real_))
  }

  null_eta <- if (copula %in% c("Frank", "Clayton")) 0 else 1
  z <- (eta - null_eta) / se
  p <- 2 * (1 - pnorm(abs(z)))

  c(null_eta = null_eta, z = z, p = p)
}

tau_ci_from_eta_ci <- function(copula, eta_ci) {
  tau_vals <- vapply(eta_ci, function(x) kendall_tau(copula, x), numeric(1))
  sort(tau_vals)
}

# ----------------------------
# 4) Coefficient extraction helpers
# ----------------------------
make_coefmat_auto <- function(fit) {
  v <- fit$summary
  if (is.null(v) || !is.numeric(v) || length(v) %% 4 != 0) return(NULL)

  mats <- list(
    byrow = matrix(v, ncol = 4, byrow = TRUE),
    bycol = matrix(v, ncol = 4, byrow = FALSE)
  )

  score_mat <- function(M) {
    sum(vapply(seq_len(4), function(j) {
      x <- M[, j]
      all(is.finite(x)) && all(x >= 0 & x <= 1)
    }, logical(1)))
  }

  scores <- vapply(mats, score_mat, numeric(1))
  M <- mats[[names(scores)[which.max(scores)]]]

  p_col <- which.max(vapply(seq_len(4), function(j) {
    mean(M[, j] >= 0 & M[, j] <= 1, na.rm = TRUE)
  }, numeric(1)))

  remaining <- setdiff(seq_len(4), p_col)
  nonnegative_cols <- remaining[vapply(remaining, function(j) all(M[, j] >= 0, na.rm = TRUE), logical(1))]

  if (length(nonnegative_cols) > 0) {
    se_col <- nonnegative_cols[which.min(vapply(nonnegative_cols, function(j) median(M[, j], na.rm = TRUE), numeric(1)))]
  } else {
    se_col <- remaining[which.min(vapply(remaining, function(j) median(abs(M[, j]), na.rm = TRUE), numeric(1)))]
  }

  remaining <- setdiff(remaining, se_col)
  statistic_col <- remaining[which.max(vapply(remaining, function(j) median(abs(M[, j]), na.rm = TRUE), numeric(1)))]
  estimate_col <- setdiff(remaining, statistic_col)

  out <- cbind(
    estimate = M[, estimate_col],
    SE = M[, se_col],
    statistic = M[, statistic_col],
    pvalue = M[, p_col]
  )

  if (!is.null(fit$estimates) && length(names(fit$estimates)) == nrow(out)) {
    rownames(out) <- names(fit$estimates)
  } else {
    rownames(out) <- paste0("par_", seq_len(nrow(out)))
  }

  out
}

coef_table_from_fit <- function(fit) {
  if (is.null(fit)) return(NULL)

  if (!is.null(fit$summary) && (is.matrix(fit$summary) || is.data.frame(fit$summary))) {
    out <- as.data.frame(fit$summary)
    out$parameter <- rownames(out)
    out <- out %>% relocate(parameter)

    cleaned_names <- tolower(gsub("[^a-z]", "", names(out)))
    names(out)[cleaned_names %in% c("estimate", "estimates", "est")] <- "estimate"
    names(out)[cleaned_names %in% c("se", "stderr", "stdse", "standarderror")] <- "SE"
    names(out)[cleaned_names %in% c("z", "zvalue", "stat", "statistic", "wald")] <- "statistic"
    names(out)[cleaned_names %in% c("p", "pvalue", "prz", "pr") | grepl("pvalue", cleaned_names)] <- "pvalue"

    if (!"estimate" %in% names(out) && !is.null(fit$estimates)) {
      out$estimate <- as.numeric(fit$estimates[out$parameter])
    }

    return(out)
  }

  cm <- make_coefmat_auto(fit)
  if (!is.null(cm)) {
    return(as.data.frame(cm) %>% rownames_to_column("parameter"))
  }

  if (!is.null(fit$estimates)) {
    return(tibble(parameter = names(fit$estimates), estimate = as.numeric(fit$estimates)))
  }

  NULL
}

extract_eta <- function(fit) {
  ct <- coef_table_from_fit(fit)
  if (is.null(ct) || !"parameter" %in% names(ct)) {
    return(c(eta = NA_real_, se = NA_real_))
  }

  eta_row <- ct %>% filter(parameter == "eta")
  if (nrow(eta_row) == 0) {
    return(c(eta = NA_real_, se = NA_real_))
  }

  eta <- if ("estimate" %in% names(eta_row)) as.numeric(eta_row$estimate[1]) else NA_real_
  se <- if ("SE" %in% names(eta_row)) as.numeric(eta_row$SE[1]) else NA_real_

  c(eta = eta, se = se)
}

is_stable_fit <- function(fit, pd_tol = pd_tol) {
  if (is.null(fit)) return(FALSE)
  if (is.null(fit$estimates) || any(!is.finite(fit$estimates))) return(FALSE)

  ct <- coef_table_from_fit(fit)
  if (is.null(ct) || !"SE" %in% names(ct)) return(FALSE)

  se <- as.numeric(ct$SE)
  if (any(!is.finite(se)) || any(se <= 0)) return(FALSE)

  if (is.null(fit$inv_info)) return(FALSE)

  eigs <- tryCatch(
    eigen(fit$inv_info, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  if (any(!is.finite(eigs)) || any(eigs <= pd_tol)) return(FALSE)

  TRUE
}

# ----------------------------
# 5) Load data and construct focal-team outcome
# ----------------------------
df_raw <- readr::read_csv(input_path, show_col_types = FALSE)

required_cols <- c(
  "match_id",
  "is_home",
  "home_first_goal_minute",
  "away_first_goal_minute",
  "home_scored",
  "away_scored",
  continuous_covariates,
  categorical_covariates
)

missing_required <- setdiff(required_cols, names(df_raw))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

df <- df_raw %>%
  mutate(
    match_id = as.character(match_id),
    is_home = to_binary01(is_home, "is_home"),
    home_scored = to_binary01(home_scored, "home_scored"),
    away_scored = to_binary01(away_scored, "away_scored"),
    home_first_goal_minute = as.numeric(home_first_goal_minute),
    away_first_goal_minute = as.numeric(away_first_goal_minute),
    ind = if_else(is_home == 1L, 1L, 2L),
    focal_scored = if_else(is_home == 1L, home_scored, away_scored),
    focal_first_goal_minute = if_else(is_home == 1L, home_first_goal_minute, away_first_goal_minute),
    obs_time = if_else(focal_scored == 1L, pmin(focal_first_goal_minute, horizon_min), horizon_min),
    status = focal_scored
  )

if (any(is.na(df$is_home))) stop("is_home contains missing values after conversion.")
if (any(is.na(df$home_scored))) stop("home_scored contains missing values after conversion.")
if (any(is.na(df$away_scored))) stop("away_scored contains missing values after conversion.")

if (any(df$home_scored == 1L & is.na(df$home_first_goal_minute))) {
  stop("Some rows have home_scored = 1 but missing home_first_goal_minute.")
}

if (any(df$away_scored == 1L & is.na(df$away_first_goal_minute))) {
  stop("Some rows have away_scored = 1 but missing away_first_goal_minute.")
}

if (any(df$obs_time <= 0, na.rm = TRUE)) {
  stop("Some obs_time values are <= 0. Check first-goal-minute coding.")
}

# ----------------------------
# 6) Validate paired structure and actual outcome construction
# ----------------------------
match_check <- df %>%
  group_by(match_id) %>%
  summarise(
    n_rows = n(),
    n_home = sum(is_home == 1L, na.rm = TRUE),
    n_away = sum(is_home == 0L, na.rm = TRUE),
    ind_values = paste(sort(unique(ind)), collapse = ","),
    home_scored_consistent = same_or_all_na(home_scored),
    away_scored_consistent = same_or_all_na(away_scored),
    home_time_consistent = same_or_all_na(home_first_goal_minute),
    away_time_consistent = same_or_all_na(away_first_goal_minute),
    competition_consistent = same_or_all_na(Round_grouped),
    .groups = "drop"
  )

bad_matches <- match_check %>%
  filter(
    n_rows != 2L |
      n_home != 1L |
      n_away != 1L |
      ind_values != "1,2" |
      !home_scored_consistent |
      !away_scored_consistent |
      !home_time_consistent |
      !away_time_consistent |
      !competition_consistent
  )

if (nrow(bad_matches) > 0) {
  print(head(bad_matches, 20))
  stop("Some matches are not valid two-row home/away pairs.")
}

outcome_check <- df %>%
  transmute(
    match_id,
    is_home,
    ind,
    status,
    expected_status = if_else(is_home == 1L, home_scored, away_scored),
    obs_time,
    expected_obs_time = if_else(
      expected_status == 1L,
      pmin(if_else(is_home == 1L, home_first_goal_minute, away_first_goal_minute), horizon_min),
      horizon_min
    )
  ) %>%
  summarise(
    n_rows = n(),
    status_mismatches = sum(status != expected_status, na.rm = TRUE),
    time_mismatches = sum(abs(obs_time - expected_obs_time) > 1e-8, na.rm = TRUE),
    missing_status = sum(is.na(status)),
    missing_time = sum(is.na(obs_time))
  )

print(match_check %>% count(n_rows, n_home, n_away, ind_values))
print(outcome_check)

if (
  outcome_check$status_mismatches > 0 ||
    outcome_check$time_mismatches > 0 ||
    outcome_check$missing_status > 0 ||
    outcome_check$missing_time > 0
) {
  stop("Reconstructed status/obs_time does not match actual focal-team scoring variables.")
}

event_pattern <- df %>%
  group_by(match_id) %>%
  summarise(
    n_events = sum(status, na.rm = TRUE),
    n_distinct_times = n_distinct(obs_time),
    .groups = "drop"
  )

print(event_pattern %>% count(n_events, n_distinct_times))

if (!any(event_pattern$n_events == 2L)) {
  stop("No matches with two scoring teams. This does not look like actual home/away first-goal timing data.")
}

# ----------------------------
# 7) Prepare CopulaCenR model input
# ----------------------------
df_cc <- df %>%
  mutate(
    id_original = match_id,
    id = as.integer(factor(match_id)),
    ind = as.integer(ind),
    obs_time = as.numeric(obs_time),
    status = as.integer(status),
    is_home = as.numeric(is_home),
    Round_grouped = as.character(Round_grouped),
    Round_grouped = if_else(is.na(Round_grouped), "Missing", Round_grouped),
    Round_grouped = factor(Round_grouped)
  ) %>%
  arrange(id, ind)

if ("Bundesliga" %in% levels(df_cc$Round_grouped)) {
  df_cc <- df_cc %>% mutate(Round_grouped = relevel(Round_grouped, ref = "Bundesliga"))
}

league_dummies <- model.matrix(~ Round_grouped, data = df_cc)
league_dummies <- league_dummies[, colnames(league_dummies) != "(Intercept)", drop = FALSE]
colnames(league_dummies) <- make.names(colnames(league_dummies), unique = TRUE)
league_dummies <- as.data.frame(league_dummies)

model_data <- bind_cols(
  df_cc %>% select(id, id_original, ind, obs_time, status, all_of(continuous_covariates)),
  league_dummies
)

var_list <- c(continuous_covariates, colnames(league_dummies))

# Keep only complete match pairs.
model_cols <- c("obs_time", "status", var_list)
model_data$.row_complete <- complete.cases(model_data[, model_cols])

ids_keep <- model_data %>%
  group_by(id) %>%
  summarise(
    n_rows = n(),
    has_home_and_away = setequal(ind, c(1L, 2L)),
    complete_pair = all(.row_complete),
    .groups = "drop"
  ) %>%
  filter(n_rows == 2L, has_home_and_away, complete_pair) %>%
  select(id)

n_ids_before <- n_distinct(model_data$id)
model_data <- model_data %>%
  semi_join(ids_keep, by = "id") %>%
  select(-.row_complete) %>%
  arrange(id, ind)
n_ids_after <- n_distinct(model_data$id)

cat("\nComplete-pair filtering\n")
cat("IDs before:", n_ids_before, "\n")
cat("IDs after :", n_ids_after, "\n")
cat("Dropped IDs:", n_ids_before - n_ids_after, "\n")

stopifnot(nrow(model_data) == 2L * n_distinct(model_data$id))
stopifnot(all(table(model_data$id) == 2L))

# Drop zero-variance covariates.
zero_var <- var_list[
  vapply(var_list, function(v) {
    x <- model_data[[v]]
    all(is.na(x)) || isTRUE(sd(x, na.rm = TRUE) == 0)
  }, logical(1))
]

if (length(zero_var) > 0) {
  warning("Dropping zero-variance covariates: ", paste(zero_var, collapse = ", "))
  var_list <- setdiff(var_list, zero_var)
  continuous_covariates <- setdiff(continuous_covariates, zero_var)
  model_data <- model_data %>% select(-all_of(zero_var))
}

# Standardise continuous football covariates, but not is_home or league dummies.
scale_cols <- setdiff(continuous_covariates, "is_home")
model_data <- model_data %>% mutate(across(all_of(scale_cols), ~ as.numeric(scale(.x))))

stopifnot(all(vapply(model_data[var_list], is.numeric, logical(1))))
stopifnot(all(table(model_data$id) == 2L))

cat("\nFinal CopulaCenR input\n")
print(model_data %>% summarise(n_rows = n(), n_ids = n_distinct(id), rows_per_id = n_rows / n_ids))
print(model_data %>% group_by(id) %>% summarise(n_events = sum(status), n_distinct_times = n_distinct(obs_time), .groups = "drop") %>% count(n_events, n_distinct_times))

write.csv(model_data, paste0(output_prefix, "_model_input.csv"), row.names = FALSE)

# ----------------------------
# 8) Fit grid of copula-margin models
# ----------------------------
fit_copula_model <- function(data, var_list, copula, margin) {
  d <- as.data.frame(data)
  d$id <- as.integer(factor(d$id))
  d$ind <- as.integer(d$ind)
  d <- d[order(d$id, d$ind), ]

  warning_messages <- character(0)

  fit <- tryCatch(
    withCallingHandlers(
      {
        rc_par_copula(
          data = d,
          var_list = var_list,
          copula = copula,
          m.dist = margin,
          method = fit_method,
          iter = fit_iter,
          stepsize = fit_stepsize
        )
      },
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      ok = FALSE,
      copula = copula,
      margin = margin,
      logLik = NA_real_,
      AIC = NA_real_,
      npar = NA_integer_,
      n_warnings = length(warning_messages),
      warning_examples = paste(unique(head(warning_messages, 3)), collapse = " | "),
      error = conditionMessage(fit),
      fit = NULL
    ))
  }

  list(
    ok = TRUE,
    copula = copula,
    margin = margin,
    logLik = suppressWarnings(as.numeric(fit$llk)),
    AIC = suppressWarnings(as.numeric(fit$AIC)),
    npar = length(fit$estimates),
    n_warnings = length(warning_messages),
    warning_examples = paste(unique(head(warning_messages, 3)), collapse = " | "),
    error = "",
    fit = fit
  )
}

fit_results <- list()
k <- 1L

for (margin in margins_to_fit) {
  for (copula in copulas_to_fit) {
    cat("\nFitting ", copula, " copula with ", margin, " margins\n", sep = "")
    fit_results[[k]] <- fit_copula_model(model_data, var_list, copula, margin)
    k <- k + 1L
  }
}

results_table <- bind_rows(lapply(seq_along(fit_results), function(i) {
  r <- fit_results[[i]]
  eta_info <- if (isTRUE(r$ok)) extract_eta(r$fit) else c(eta = NA_real_, se = NA_real_)
  indep <- independence_test(r$copula, eta_info["eta"], eta_info["se"])

  tibble(
    idx = i,
    copula = r$copula,
    margin = r$margin,
    ok = r$ok,
    stable = if (isTRUE(r$ok)) is_stable_fit(r$fit) else FALSE,
    logLik = r$logLik,
    AIC = r$AIC,
    npar = r$npar,
    n_warnings = r$n_warnings,
    warning_examples = r$warning_examples,
    error = r$error,
    eta = eta_info["eta"],
    eta_se = eta_info["se"],
    eta_used = effective_eta(r$copula, eta_info["eta"]),
    admissible = is_admissible_eta(r$copula, eta_info["eta"]),
    kendall_tau = kendall_tau(r$copula, eta_info["eta"]),
    independence_null_eta = indep["null_eta"],
    independence_z = indep["z"],
    independence_p = indep["p"]
  )
}))

results_table_sorted <- results_table %>%
  arrange(desc(ok), desc(stable), desc(admissible), AIC)

write.csv(results_table_sorted, paste0(output_prefix, "_aic_results_all.csv"), row.names = FALSE)
print(results_table_sorted)

# ----------------------------
# 9) Select final reportable model
# ----------------------------
successful_rows <- which(results_table$ok & is.finite(results_table$AIC) & is.finite(results_table$logLik))
if (length(successful_rows) == 0) stop("No successful finite-AIC fits were obtained.")

best_overall_idx <- successful_rows[which.min(results_table$AIC[successful_rows])]
best_overall <- results_table %>% filter(idx == best_overall_idx)
best_overall_fit <- fit_results[[best_overall_idx]]$fit

reportable_rows <- which(
  results_table$ok &
    results_table$stable &
    results_table$admissible &
    is.finite(results_table$AIC) &
    is.finite(results_table$logLik)
)

if (length(reportable_rows) == 0) {
  stop("No reportable models were obtained. A reportable model must be successful, stable, admissible, and have finite AIC/logLik.")
}

best_reportable_idx <- reportable_rows[which.min(results_table$AIC[reportable_rows])]
best_reportable <- results_table %>% filter(idx == best_reportable_idx)
best_reportable_fit <- fit_results[[best_reportable_idx]]$fit

reportable_results <- results_table %>%
  filter(ok, stable, admissible, is.finite(AIC), is.finite(logLik)) %>%
  arrange(AIC)

write.csv(reportable_results, paste0(output_prefix, "_aic_results_reportable.csv"), row.names = FALSE)

cat("\nBest overall model, possibly non-reportable\n")
print(best_overall)

cat("\nBest reportable model\n")
print(best_reportable)

# ----------------------------
# 10) Coefficient table for final reportable model
# ----------------------------
best_coef_table <- coef_table_from_fit(best_reportable_fit)

if (is.null(best_coef_table) || !"estimate" %in% names(best_coef_table)) {
  stop("Could not extract a usable coefficient table from the final reportable model.")
}

if (!"SE" %in% names(best_coef_table)) best_coef_table$SE <- NA_real_
if (!"pvalue" %in% names(best_coef_table)) best_coef_table$pvalue <- NA_real_
if (!"statistic" %in% names(best_coef_table)) best_coef_table$statistic <- NA_real_

best_coef_table <- best_coef_table %>%
  mutate(
    estimate = as.numeric(estimate),
    SE = as.numeric(SE),
    statistic = as.numeric(statistic),
    pvalue = as.numeric(pvalue),
    stars = stars_from_p(pvalue)
  )

alpha_hat <- best_coef_table %>%
  filter(parameter %in% c("alpha", "k", "shape")) %>%
  slice(1) %>%
  pull(estimate)

if (length(alpha_hat) == 0 || !is.finite(alpha_hat[1])) alpha_hat <- NA_real_ else alpha_hat <- alpha_hat[1]

selected_margin <- best_reportable$margin[1]
beta_rows <- best_coef_table$parameter %in% var_list
weibull_selected <- selected_margin == "Weibull"

best_coef_table <- best_coef_table %>%
  mutate(
    AF = if_else(beta_rows & weibull_selected & is.finite(alpha_hat), exp(estimate / alpha_hat), NA_real_),
    section = case_when(
      parameter %in% c("lambda", "alpha", "k", "mu", "sigma", "shape", "scale") ~ "Model parameters",
      parameter == "eta" ~ "Copula dependence",
      parameter %in% colnames(league_dummies) ~ "Competition",
      parameter %in% continuous_covariates ~ "Covariates",
      TRUE ~ "Other"
    ),
    label = case_when(
      parameter == "is_home" ~ "Home-team indicator",
      parameter %in% colnames(league_dummies) ~ clean_param_label(parameter),
      TRUE ~ clean_param_label(parameter)
    )
  ) %>%
  select(section, parameter, label, estimate, AF, SE, pvalue, stars, statistic, everything())

write.csv(best_coef_table, paste0(output_prefix, "_best_reportable_coef_table.csv"), row.names = FALSE)
print(best_coef_table)

# ----------------------------
# 11) Dependence summary for final reportable model
# ----------------------------
best_eta <- best_reportable$eta[1]
best_eta_se <- best_reportable$eta_se[1]
best_copula <- best_reportable$copula[1]
best_margin <- best_reportable$margin[1]

eta_ci_raw <- c(NA_real_, NA_real_)
eta_ci_used <- c(NA_real_, NA_real_)
tau_ci <- c(NA_real_, NA_real_)

if (is.finite(best_eta) && is.finite(best_eta_se) && best_eta_se > 0) {
  eta_ci_raw <- best_eta + c(-1, 1) * qnorm(0.975) * best_eta_se

  if (best_copula %in% c("Gumbel", "Joe")) {
    eta_ci_used <- c(max(1, eta_ci_raw[1]), eta_ci_raw[2])
  } else if (best_copula == "Clayton") {
    eta_ci_used <- c(max(0, eta_ci_raw[1]), eta_ci_raw[2])
  } else {
    eta_ci_used <- eta_ci_raw
  }

  tau_ci <- tau_ci_from_eta_ci(best_copula, eta_ci_used)
}

dependence_summary <- tibble(
  copula = best_copula,
  margin = best_margin,
  eta = best_eta,
  eta_used_for_interpretation = best_reportable$eta_used[1],
  eta_se = best_eta_se,
  kendall_tau = best_reportable$kendall_tau[1],
  independence_null_eta = best_reportable$independence_null_eta[1],
  independence_z = best_reportable$independence_z[1],
  independence_p = best_reportable$independence_p[1],
  eta_ci_raw_lower = eta_ci_raw[1],
  eta_ci_raw_upper = eta_ci_raw[2],
  eta_ci_used_lower = eta_ci_used[1],
  eta_ci_used_upper = eta_ci_used[2],
  tau_ci_lower = tau_ci[1],
  tau_ci_upper = tau_ci[2]
)

write.csv(dependence_summary, paste0(output_prefix, "_best_reportable_dependence_summary.csv"), row.names = FALSE)
print(dependence_summary)

# ----------------------------
# 12) Save compact run summary
# ----------------------------
capture.output(
  {
    cat("CopulaCenR pooled focal-team copula model\n")
    cat("==========================================\n\n")

    cat("Outcome definition\n")
    cat("- Two rows per match: home-team focal row and away-team focal row.\n")
    cat("- status = 1 if the focal team scored.\n")
    cat("- obs_time = focal team's first-goal minute if scored, otherwise 90.\n")
    cat("- Stoppage-time first goals are collapsed to minute 90.\n\n")

    cat("Fitting settings\n")
    cat("fit_method =", fit_method, "\n")
    cat("fit_iter =", fit_iter, "\n")
    cat("fit_stepsize =", fit_stepsize, "\n")
    cat("boundary_tol =", boundary_tol, "\n\n")

    cat("Pairing and outcome checks\n")
    print(match_check %>% count(n_rows, n_home, n_away, ind_values))
    print(outcome_check)
    print(event_pattern %>% count(n_events, n_distinct_times))

    cat("\nFinal model input\n")
    print(model_data %>% summarise(n_rows = n(), n_ids = n_distinct(id), rows_per_id = n_rows / n_ids))

    cat("\nCovariates passed to CopulaCenR\n")
    print(var_list)

    cat("\nAll model grid results\n")
    print(results_table_sorted)

    cat("\nBest overall model\n")
    print(best_overall)

    cat("\nBest reportable model\n")
    print(best_reportable)

    cat("\nBest reportable coefficient table\n")
    print(best_coef_table)

    cat("\nBest reportable dependence summary\n")
    print(dependence_summary)
  },
  file = paste0(output_prefix, "_best_reportable_summary.txt")
)

cat("\nDone. Saved files:\n")
cat("- ", paste0(output_prefix, "_model_input.csv"), "\n", sep = "")
cat("- ", paste0(output_prefix, "_aic_results_all.csv"), "\n", sep = "")
cat("- ", paste0(output_prefix, "_aic_results_reportable.csv"), "\n", sep = "")
cat("- ", paste0(output_prefix, "_best_reportable_coef_table.csv"), "\n", sep = "")
cat("- ", paste0(output_prefix, "_best_reportable_dependence_summary.csv"), "\n", sep = "")
cat("- ", paste0(output_prefix, "_best_reportable_summary.txt"), "\n", sep = "")
