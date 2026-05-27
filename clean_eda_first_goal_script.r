# ============================================================
# Exploratory Data Analysis for time to first goal
# - Home and away right-censored first-goal margins
# - Figures and tables used in the EDA chapter
# - Censoring summaries, KM curves, league heterogeneity,
#   Weibull diagnostics, dependence EDA, multicollinearity,
#   goal-time correlations, and formation frequencies
# ============================================================

# ----------------------------
# 0) Packages and settings
# ----------------------------
needed_pkgs <- c(
  "dplyr", "tidyr", "readr", "ggplot2", "tibble",
  "survival", "survminer", "ggcorrplot", "muhaz"
)

for (p in needed_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
invisible(lapply(needed_pkgs, library, character.only = TRUE))

# Input data used for the EDA chapter.
data_path <- "df_final_with_rolling_averages.rds"

# Output directory for all EDA tables and figures.
output_dir <- "eda_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Main variables.
league_var <- "Round_grouped"
home_time_var <- "home_first_goal_minute"
home_status_var <- "home_scored"
away_time_var <- "away_first_goal_minute"
away_status_var <- "away_scored"

# Goal-count variables used in the correlation table.
# Change these names here if your dataset uses different names.
home_goals_var <- "Home_goals"
away_goals_var <- "Away_goals"

time_horizon <- 90
corr_threshold <- 0.70

# Le Coz formation frequencies used for the comparison table in the thesis.
lecoz_formations <- tibble::tribble(
  ~formation, ~lecoz_percent,
  "4-2-3-1", 28.9,
  "4-3-3",   23.8,
  "4-4-2",   22.7,
  "3-5-2",    8.7,
  "3-4-3",    6.3,
  "4-5-1",    9.2,
  "5-4-1",    0.3,
  "5-3-2",    0.1
)

# ----------------------------
# 1) Helper functions
# ----------------------------
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

save_plot <- function(plot, filename, width = 8, height = 5, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

save_survplot <- function(survplot, filename, width = 8, height = 6, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = print(survplot),
    width = width,
    height = height,
    dpi = dpi
  )
}

tidy_percent <- function(x, digits = 3) {
  round(as.numeric(x), digits)
}

rename_for_table <- function(x) {
  x %>%
    gsub("_roll5$", "", .) %>%
    gsub("_", " ", .)
}

# Safely convert common binary encodings to 0/1.
to_binary01 <- function(x, var_name = "variable") {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))

  x_chr <- trimws(tolower(as.character(x)))
  out <- dplyr::case_when(
    x_chr %in% c("1", "true", "t", "yes", "y") ~ 1L,
    x_chr %in% c("0", "false", "f", "no", "n") ~ 0L,
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

extract_km_medians <- function(sf_league) {
  med_tbl <- as.data.frame(summary(sf_league)$table)
  med_tbl$league_raw <- rownames(med_tbl)
  rownames(med_tbl) <- NULL

  med_tbl %>%
    transmute(
      league = gsub(paste0("^", league_var, "="), "", league_raw),
      median = median,
      lower_95 = `0.95LCL`,
      upper_95 = `0.95UCL`
    ) %>%
    arrange(median, league)
}

weibull_diag_df <- function(fit) {
  if (is.null(fit$strata)) {
    stop("Use a grouped survfit object, for example survfit(Surv(t, d) ~ side, data = df_long).")
  }

  n_per_stratum <- as.integer(fit$strata)
  idx_end <- cumsum(n_per_stratum)
  idx_start <- c(1, head(idx_end, -1) + 1)
  strata_names <- names(fit$strata)

  out <- vector("list", length(n_per_stratum))

  for (j in seq_along(n_per_stratum)) {
    ii <- idx_start[j]:idx_end[j]
    out[[j]] <- tibble(
      time = fit$time[ii],
      survival = fit$surv[ii],
      strata = strata_names[j]
    )
  }

  bind_rows(out) %>%
    filter(is.finite(time), is.finite(survival), time > 0, survival > 0, survival < 1) %>%
    mutate(
      log_time = log(time),
      log_minus_log_survival = log(-log(survival)),
      strata = gsub("side=", "", strata)
    )
}

cor_test_row <- function(data, x, y, label_x, label_y) {
  dd <- data %>%
    select(all_of(c(x, y))) %>%
    filter(if_all(everything(), ~ is.finite(.)))

  test <- stats::cor.test(dd[[x]], dd[[y]], method = "pearson")

  tibble(
    time = label_x,
    goals = label_y,
    correlation = unname(test$estimate),
    ci_lower = test$conf.int[1],
    ci_upper = test$conf.int[2],
    p_value = test$p.value,
    n = nrow(dd)
  )
}

# ----------------------------
# 2) Load and clean data
# ----------------------------
log_msg("Loading data")

df_raw <- readRDS(data_path)

required_cols <- c(
  league_var,
  home_time_var, home_status_var,
  away_time_var, away_status_var,
  "home_formation_grouped", "away_formation_grouped"
)

missing_required <- setdiff(required_cols, names(df_raw))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

df0 <- df_raw %>%
  mutate(
    !!league_var := as.factor(.data[[league_var]]),
    !!home_status_var := to_binary01(.data[[home_status_var]], home_status_var),
    !!away_status_var := to_binary01(.data[[away_status_var]], away_status_var),
    !!home_time_var := as.numeric(.data[[home_time_var]]),
    !!away_time_var := as.numeric(.data[[away_time_var]])
  ) %>%
  filter(
    !is.na(.data[[league_var]]),
    is.finite(.data[[home_time_var]]),
    is.finite(.data[[away_time_var]]),
    .data[[home_time_var]] >= 0,
    .data[[away_time_var]] >= 0,
    .data[[home_status_var]] %in% c(0, 1),
    .data[[away_status_var]] %in% c(0, 1)
  )

log_msg("Rows after EDA cleaning: ", nrow(df0))

# Long format for survival plots.
df_long <- df0 %>%
  transmute(
    match_row = row_number(),
    Round_grouped = .data[[league_var]],
    t_home = .data[[home_time_var]],
    d_home = .data[[home_status_var]],
    t_away = .data[[away_time_var]],
    d_away = .data[[away_status_var]]
  ) %>%
  pivot_longer(
    cols = c(t_home, d_home, t_away, d_away),
    names_to = c(".value", "side"),
    names_pattern = "([td])_(home|away)"
  ) %>%
  mutate(
    side = factor(side, levels = c("home", "away"), labels = c("Home", "Away")),
    interval_15 = cut(
      t,
      breaks = c(0, 15, 30, 45, 60, 75, Inf),
      right = TRUE,
      include.lowest = TRUE,
      labels = c("0--15", "16--30", "31--45", "46--60", "61--75", "76--90")
    )
  )

# ----------------------------
# 3) Censoring summaries
# ----------------------------
eda_censor_overall <- df0 %>%
  summarise(
    home_censor_rate = mean(.data[[home_status_var]] == 0),
    away_censor_rate = mean(.data[[away_status_var]] == 0),
    both_censor_rate = mean(.data[[home_status_var]] == 0 & .data[[away_status_var]] == 0),
    n_matches = n()
  )

eda_league_counts_censor <- df0 %>%
  group_by(Round_grouped = .data[[league_var]]) %>%
  summarise(
    n_matches = n(),
    home_censor_rate = mean(.data[[home_status_var]] == 0),
    away_censor_rate = mean(.data[[away_status_var]] == 0),
    both_censor_rate = mean(.data[[home_status_var]] == 0 & .data[[away_status_var]] == 0),
    .groups = "drop"
  ) %>%
  arrange(desc(n_matches))

write.csv(eda_censor_overall, file.path(output_dir, "eda_censor_overall.csv"), row.names = FALSE)
write.csv(eda_league_counts_censor, file.path(output_dir, "eda_league_counts_censor.csv"), row.names = FALSE)

cat("\n===== Overall censoring summary =====\n")
print(eda_censor_overall)

cat("\n===== League censoring summary =====\n")
print(eda_league_counts_censor)

# ----------------------------
# 4) Distribution of observed first-goal times
# ----------------------------
p_hist_side <- ggplot(df_long, aes(x = t)) +
  geom_histogram(
    breaks = seq(0, time_horizon, by = 15),
    boundary = 0,
    closed = "right"
  ) +
  facet_wrap(~side, ncol = 1) +
  labs(
    x = "Observed first-goal time (minutes)",
    y = "Count",
    title = "Distribution of observed first-goal times by margin"
  ) +
  theme_minimal()

print(p_hist_side)
save_plot(p_hist_side, "eda_hist_side.png", width = 8, height = 6)

# Optional 15-minute interval count table.
eda_hist_15_table <- df_long %>%
  count(side, interval_15, name = "count") %>%
  group_by(side) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()
write.csv(eda_hist_15_table, file.path(output_dir, "eda_hist_15_table.csv"), row.names = FALSE)

# ----------------------------
# 5) Kaplan-Meier curves and log-rank test
# ----------------------------
sf_overall <- survival::survfit(survival::Surv(t, d) ~ 1, data = df_long)
sf_side <- survival::survfit(survival::Surv(t, d) ~ side, data = df_long)
sf_league <- survival::survfit(survival::Surv(t, d) ~ Round_grouped, data = df_long)

p_km_side <- survminer::ggsurvplot(
  sf_side,
  data = df_long,
  conf.int = FALSE,
  risk.table = TRUE,
  xlab = "Minutes",
  ylab = "Survival probability",
  title = "Kaplan--Meier curves for home and away first-goal times",
  legend.title = "Margin",
  legend.labs = c("Home", "Away")
)
print(p_km_side)
save_survplot(p_km_side, "eda_km_side.png", width = 9, height = 7)

logrank_side <- survival::survdiff(survival::Surv(t, d) ~ side, data = df_long)
logrank_summary <- tibble(
  chisq = unname(logrank_side$chisq),
  df = length(logrank_side$n) - 1,
  p_value = stats::pchisq(logrank_side$chisq, df = length(logrank_side$n) - 1, lower.tail = FALSE)
)
write.csv(logrank_summary, file.path(output_dir, "eda_logrank_home_away.csv"), row.names = FALSE)

cat("\n===== Log-rank test, home vs away =====\n")
print(logrank_summary)

p_km_league <- survminer::ggsurvplot(
  sf_league,
  data = df_long,
  conf.int = FALSE,
  risk.table = FALSE,
  xlab = "Minutes",
  ylab = "Survival probability",
  title = "Kaplan--Meier curves by league"
)
print(p_km_league)
save_survplot(p_km_league, "eda_km_by_league.png", width = 9, height = 6)

eda_median_by_league <- extract_km_medians(sf_league)
write.csv(eda_median_by_league, file.path(output_dir, "eda_median_by_league.csv"), row.names = FALSE)

cat("\n===== KM median time by league =====\n")
print(eda_median_by_league)

# ----------------------------
# 6) Home-away observed-time dependence plot
# ----------------------------
dep_df <- df0 %>%
  transmute(
    home_first_goal_minute = .data[[home_time_var]],
    away_first_goal_minute = .data[[away_time_var]],
    home_scored = .data[[home_status_var]],
    away_scored = .data[[away_status_var]]
  ) %>%
  filter(
    is.finite(home_first_goal_minute),
    is.finite(away_first_goal_minute)
  )

tau_raw <- cor(
  dep_df$home_first_goal_minute,
  dep_df$away_first_goal_minute,
  method = "kendall",
  use = "complete.obs"
)

eda_raw_kendall <- tibble(raw_kendall_tau_observed_times = tau_raw, n = nrow(dep_df))
write.csv(eda_raw_kendall, file.path(output_dir, "eda_raw_kendall_home_away_observed_times.csv"), row.names = FALSE)

p_scatter_times <- ggplot(dep_df, aes(x = home_first_goal_minute, y = away_first_goal_minute)) +
  geom_count(alpha = 0.65) +
  labs(
    x = "Observed home first-goal time (minutes)",
    y = "Observed away first-goal time (minutes)",
    size = "Number of matches",
    title = "Observed home and away first-goal times"
  ) +
  theme_minimal()

print(p_scatter_times)
save_plot(p_scatter_times, "eda_scatter_home_away_times.png", width = 8, height = 6)

cat("\n===== Raw Kendall tau on observed home/away times =====\n")
print(eda_raw_kendall)

# ----------------------------
# 7) Weibull plausibility diagnostics
# ----------------------------
# Smoothed hazard estimates by margin, shown together in one figure.
make_hazard_df <- function(time, status, side_label) {
  dd <- tibble(time = time, status = status) %>%
    filter(is.finite(time), is.finite(status), time > 0)

  hz <- muhaz::muhaz(
    times = dd$time,
    delta = dd$status,
    min.time = min(dd$time, na.rm = TRUE),
    max.time = max(dd$time, na.rm = TRUE)
  )

  tibble(
    time = hz$est.grid,
    hazard = hz$haz.est,
    side = side_label
  )
}

hazard_df <- bind_rows(
  make_hazard_df(df0[[home_time_var]], df0[[home_status_var]], "Home"),
  make_hazard_df(df0[[away_time_var]], df0[[away_status_var]], "Away")
)

p_hazard <- ggplot(hazard_df, aes(x = time, y = hazard, linetype = side)) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "Minutes",
    y = "Estimated hazard",
    linetype = "Margin",
    title = "Smoothed hazard estimates by margin"
  ) +
  theme_minimal()

print(p_hazard)
save_plot(p_hazard, "eda_smoothed_hazard.png", width = 8, height = 5)
write.csv(hazard_df, file.path(output_dir, "eda_smoothed_hazard_values.csv"), row.names = FALSE)

# Weibull/PH diagnostic: log{-log S(t)} against log(t).
diag_df <- weibull_diag_df(sf_side)

p_weibull_diag <- ggplot(diag_df, aes(x = log_time, y = log_minus_log_survival)) +
  geom_line() +
  facet_wrap(~strata) +
  labs(
    x = "log(time)",
    y = "log(-log(S(t)))",
    title = "Weibull/PH diagnostic by margin"
  ) +
  theme_minimal()

print(p_weibull_diag)
save_plot(p_weibull_diag, "eda_weibull_ph_diag.png", width = 9, height = 5)
write.csv(diag_df, file.path(output_dir, "eda_weibull_ph_diag_values.csv"), row.names = FALSE)

# ----------------------------
# 8) Multicollinearity among rolling covariates
# ----------------------------
roll_num <- df0 %>%
  select(where(is.numeric)) %>%
  select(matches("_roll")) %>%
  select(-matches("^id$|^match_id$|^Unnamed|^index$"), everything())

if (ncol(roll_num) < 2) {
  stop("Fewer than two numeric rolling variables were found.")
}

zero_sd <- names(roll_num)[vapply(roll_num, function(x) sd(x, na.rm = TRUE) == 0, logical(1))]
if (length(zero_sd) > 0) roll_num <- roll_num %>% select(-all_of(zero_sd))

corr_roll <- cor(as.matrix(roll_num), use = "pairwise.complete.obs")
corr_roll[is.na(corr_roll)] <- 0

high_corr_idx <- which(abs(corr_roll) >= corr_threshold & upper.tri(corr_roll), arr.ind = TRUE)

eda_highcorr_pairs <- tibble(
  variable_1 = colnames(corr_roll)[high_corr_idx[, 1]],
  variable_2 = colnames(corr_roll)[high_corr_idx[, 2]],
  correlation = corr_roll[high_corr_idx]
) %>%
  arrange(desc(abs(correlation)))

eda_highcorr_pairs_top <- eda_highcorr_pairs %>%
  slice_head(n = 10) %>%
  mutate(
    variable_1_table = rename_for_table(variable_1),
    variable_2_table = rename_for_table(variable_2)
  )

write.csv(eda_highcorr_pairs, file.path(output_dir, "eda_highcorr_pairs_all.csv"), row.names = FALSE)
write.csv(eda_highcorr_pairs_top, file.path(output_dir, "eda_highcorr_pairs_top10.csv"), row.names = FALSE)

# Heatmap restricted to variables involved in the strongest high-correlation pairs.
vars_for_heatmap <- unique(c(eda_highcorr_pairs_top$variable_1, eda_highcorr_pairs_top$variable_2))
corr_heatmap <- corr_roll[vars_for_heatmap, vars_for_heatmap, drop = FALSE]

p_corr_heatmap <- ggcorrplot::ggcorrplot(
  corr_heatmap,
  lab = FALSE,
  type = "lower"
) +
  labs(title = "Correlation heatmap for strongly correlated rolling covariates") +
  theme_minimal()

print(p_corr_heatmap)
save_plot(p_corr_heatmap, "eda_corr_heatmap_roll_only.png", width = 9, height = 8)

cat("\n===== Top high-correlation rolling covariate pairs =====\n")
print(eda_highcorr_pairs_top)

# ----------------------------
# 9) Correlation between first-goal timing and number of goals
# ----------------------------
# The table in the thesis uses observed first-goal timings. Therefore, as in the EDA text,
# this is descriptive and does not correct for censoring.
goal_count_candidates <- c(home_goals_var, away_goals_var)
missing_goal_counts <- setdiff(goal_count_candidates, names(df0))

if (length(missing_goal_counts) > 0) {
  warning(
    "Goal-count columns not found: ", paste(missing_goal_counts, collapse = ", "),
    ". The goal-time correlation table will be skipped unless you rename the variables above."
  )
  eda_corr_goals_time <- tibble()
} else {
  df_goals <- df0 %>%
    mutate(
      total_goals = .data[[home_goals_var]] + .data[[away_goals_var]],
      overall_first_goal_minute = pmin(.data[[home_time_var]], .data[[away_time_var]], na.rm = TRUE)
    )

  eda_corr_goals_time <- bind_rows(
    cor_test_row(df_goals, "overall_first_goal_minute", "total_goals", "Overall First Goal Minute", "Total Goals"),
    cor_test_row(df_goals, home_time_var, home_goals_var, "Home First Goal Minute", "Home Goals"),
    cor_test_row(df_goals, away_time_var, away_goals_var, "Away First Goal Minute", "Away Goals"),
    cor_test_row(df_goals, home_time_var, away_goals_var, "Home First Goal Minute", "Away Goals"),
    cor_test_row(df_goals, away_time_var, home_goals_var, "Away First Goal Minute", "Home Goals")
  )

  write.csv(eda_corr_goals_time, file.path(output_dir, "eda_corr_goals_time.csv"), row.names = FALSE)

  cat("\n===== Correlation between first-goal timing and goals scored =====\n")
  print(eda_corr_goals_time)
}

# ----------------------------
# 10) Starting formation frequency comparison with Le Coz
# ----------------------------
formation_freq_home <- df0 %>%
  count(formation = home_formation_grouped, name = "home_n") %>%
  mutate(home_percent = 100 * home_n / sum(home_n))

formation_freq_away <- df0 %>%
  count(formation = away_formation_grouped, name = "away_n") %>%
  mutate(away_percent = 100 * away_n / sum(away_n))

eda_formations_comparison <- lecoz_formations %>%
  left_join(formation_freq_home, by = "formation") %>%
  left_join(formation_freq_away, by = "formation") %>%
  mutate(
    home_percent = tidy_percent(home_percent, 1),
    away_percent = tidy_percent(away_percent, 1),
    lecoz_percent = tidy_percent(lecoz_percent, 1),
    home_n = replace_na(home_n, 0L),
    away_n = replace_na(away_n, 0L)
  ) %>%
  select(formation, lecoz_percent, home_percent, away_percent, home_n, away_n)

write.csv(eda_formations_comparison, file.path(output_dir, "eda_formations_comparison_lecoz.csv"), row.names = FALSE)

cat("\n===== Formation frequency comparison =====\n")
print(eda_formations_comparison)

# ----------------------------
# 11) Compact EDA summary file
# ----------------------------
capture.output(
  {
    cat("EDA summary for time to first goal\n")
    cat("==================================\n\n")

    cat("Rows used:", nrow(df0), "\n\n")

    cat("Overall censoring summary\n")
    print(eda_censor_overall)

    cat("\nLeague censoring summary\n")
    print(eda_league_counts_censor)

    cat("\nLog-rank test, home vs away\n")
    print(logrank_summary)

    cat("\nKM median time by league\n")
    print(eda_median_by_league)

    cat("\nRaw Kendall tau on observed home/away times\n")
    print(eda_raw_kendall)

    cat("\nTop high-correlation rolling covariate pairs\n")
    print(eda_highcorr_pairs_top)

    cat("\nCorrelation between first-goal timing and goals scored\n")
    print(eda_corr_goals_time)

    cat("\nFormation frequency comparison\n")
    print(eda_formations_comparison)
  },
  file = file.path(output_dir, "eda_summary.txt")
)

log_msg("Done. EDA outputs saved to: ", output_dir)
