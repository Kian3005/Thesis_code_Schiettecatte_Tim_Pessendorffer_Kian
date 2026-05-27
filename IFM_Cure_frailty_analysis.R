library(tidyverse)
library(survival)
library(copula)
library(MASS) # Watch out if you want to use select from tidyverse!
library(broom)
library(dplyr)
library(CopulaCenR)
if (!require("flexsurv")) install.packages("flexsurv")
library(flexsurv)
library(smcure)
library(mvtnorm)
library(parfm)
library(flexsurvcure)



df <- read_rds("df_final_with_rolling_averages.rds")

# This code should only be runned once
# Check unique values for round
unique_rounds <- df |> 
  dplyr::select(Round) |> 
  distinct() 

# Change in round: regular season, Play-off I, Play-off II, Europe play-offs, Champions play-offs, Relegation play-offs, European competition play-off becomes all Jupiler Pro League
# European competition qualification play-offs — Semi-finals and European competition qualification play-offs — Finals as Eredivisie
# Relegation tie-breaker as Serie A
df <- df |>
  mutate(
    Round_grouped = case_when(
      str_detect(Round, "Regular season|Play-off I|Play-off II|Europe play-offs|Champions play-offs|Relegation play-offs|European competition play-off") ~ "Jupiler Pro League",
      str_detect(Round, "Eredivisie|European competition qualification play-offs — Semi-finals|European competition qualification play-offs — Finals") ~ "Eredivisie",
      str_detect(Round, "Serie A|Relegation tie-breaker") ~ "Serie A",
      str_detect(Round, "PL") ~ "Premier League",
      str_detect(Round, "Bundesliga") ~ "Bundesliga",
      str_detect(Round, "Ligue 1") ~ "Ligue 1",
      str_detect(Round, "LaLiga") ~ "La Liga",
      TRUE ~ NA_character_
    )
  )

# Save dataset with new Round_grouped variable
saveRDS(df, "df_final_with_rolling_averages.rds")
# Until here 

df <- read_rds("df_final_with_rolling_averages.rds")
# Check correlations number goals and first goal minute and also with home and away goals
cor.test(df$number_goals, df$first_goal_minute)
cor.test(df$home_first_goal_minute, df$Home_goals)
cor.test(df$away_first_goal_minute, df$Away_goals)
cor.test(df$home_first_goal_minute, df$Away_goals)
cor.test(df$away_first_goal_minute, df$Home_goals)

# Check frequency distributions of grouped formations, percentagewise
prop.table(table(df$home_formation_grouped))
prop.table(table(df$away_formation_grouped))

# General dependency analysis
scored_both <- subset(df, home_scored == 1 & away_scored == 1)

cor_pearson <- cor.test(scored_both$home_first_goal_minute, 
                        scored_both$away_first_goal_minute, method = "pearson")
cor_spearman <- cor.test(scored_both$home_first_goal_minute, 
                         scored_both$away_first_goal_minute, method = "spearman")

print(cor_pearson)
print(cor_spearman)


contingency_table <- table(df$home_scored, df$away_scored)
rownames(contingency_table) <- c("Home No Goal", "Home Goal")
colnames(contingency_table) <- c("Away No Goal", "Away Goal")

chi_test <- chisq.test(contingency_table)

print(contingency_table)
print(chi_test)

library(ggplot2)
ggplot(scored_both, aes(x = home_first_goal_minute, y = away_first_goal_minute)) +
  geom_jitter(alpha = 0.3, color = "darkblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Dependency time first goals",
       x = "Minute Home team",
       y = "Minute Away team") +
  theme_minimal()


base_vars <- c(
  "xG", "goals", "possesion", 
  "Number_passes", "Passing_accuracy", 
  "shots", "Shots_on_Target", 
  "Saves", "Corners"
)

times <- c("15", "30", "45", "60", "75", "90")
metrics <- c("xG", "PSxG", "Goal", "Saved", "Blocked", "OffTarget", "Woodwork")

time_vars <- outer(metrics, times, paste, sep="_") |> as.vector()

all_stat_names <- c(base_vars, time_vars)

suffix <- "_roll5" 

candidate_vars <- c(

  paste0("Home_", all_stat_names, "_for", suffix),      
  paste0("Home_", all_stat_names, "_against", suffix),  
  
  # --- AWAY TEAM HISTORIE ---
  paste0("Away_", all_stat_names, "_for", suffix),      
  paste0("Away_", all_stat_names, "_against", suffix),  
  
  "home_formation_grouped", 
  "away_formation_grouped", 
  "Round_grouped",
  "derby_type" 
)


df_model_home <- df |>
  dplyr::select(Home, home_first_goal_minute, home_scored, all_of(candidate_vars), match_id, season) |>
  drop_na() |> 
  mutate(across(where(is.numeric) & !all_of(c("home_first_goal_minute", "home_scored", "match_id")), scale))


# Comparing different baseline distributions without covariates based on AIC
dist <- c("weibull", "exponential", "gaussian", "logistic","lognormal", "loglogistic")
for (d in dist) {
  model <- survreg(
    formula = Surv(home_first_goal_minute, home_scored) ~ 1,
    data = df_model_home, 
    dist = d
  )
  cat("Distribution:", d, "AIC:", AIC(model), "\n")
}

# Also fit generalized gamma same way as done further in this script but now no covariates included
flexsurvreg(Surv(home_first_goal_minute, home_scored) ~ 1, data = df_model_home, dist = "gengamma")

# Search for best marginal model for Home
# Forward AIC
full_formula_string <- paste("Surv(home_first_goal_minute, home_scored) ~", 
                             paste(candidate_vars, collapse = " + "))
full_f <- as.formula(full_formula_string)
full <- survreg(
  formula = full_f,
  data = df_model_home, 
  dist = "weibull"
)

null <- survreg(
  Surv(home_first_goal_minute, home_scored) ~ 1, 
  data = df_model_home, 
  dist = "weibull"
)

AIC_forwards <- stepAIC(null, scope = list(lower = null, upper = full), direction = "forward")
final_model <- survreg(
  Surv(home_first_goal_minute, home_scored) ~ 
    Home_Number_passes_for_roll5 + 
    Away_Number_passes_for_roll5 + 
    Home_xG_for_roll5 + 
    Away_xG_against_roll5 + 
    Away_goals_for_roll5 + 
    Home_goals_for_roll5 + 
    Home_shots_for_roll5 + 
    Away_Blocked_30_for_roll5 + 
    Home_PSxG_15_for_roll5 + 
    Home_Woodwork_45_against_roll5 + 
    Away_Shots_on_Target_against_roll5 + 
    Away_Woodwork_45_against_roll5 + 
    home_formation_grouped + 
    Home_Passing_accuracy_for_roll5 + 
    Home_Blocked_90_for_roll5 + 
    Away_OffTarget_15_against_roll5 + 
    Away_Woodwork_60_for_roll5 + 
    Away_Blocked_15_against_roll5 + 
    Home_PSxG_15_against_roll5 + 
    Home_xG_15_against_roll5,
  data = df_model_home,
  dist = "weibull"
)

summary(final_model)

library(broom)
library(dplyr)




AIC_backwards <- stepAIC(full, scope = list(lower = null, upper = full), direction = "backward") # Higher AIC and more parameters

AIC_both <- stepAIC(full, scope = list(lower = null, upper = full), direction = "both") # Also higher and more terms

AIC_both_for <- stepAIC(null, scope = list(lower = null, upper = full), direction = "both") #same as just forward

# Only using total variables or total + first 15 and last 15 for example for xG and Goal because they look most important

base_candidates <- candidate_vars[!grepl("_(15|30|45|60|75|90)_", candidate_vars)]

time_candidates <- candidate_vars[grepl("_(15|30|45|60|75|90)_", candidate_vars)]

time_candidates_clean <- time_candidates[grepl("_(xG|Goal|PSxG)_", time_candidates)]

final_time_vars <- time_candidates_clean[grepl("_(15|90)_", time_candidates_clean)]

smart_candidates <- c(base_candidates, final_time_vars)

# Comparing different baseline distributions without covariates based on AIC
dist <- c("weibull", "exponential", "gaussian", "logistic","lognormal", "loglogistic")
for (d in dist) {
  model <- survreg(
    formula = Surv(away_first_goal_minute, away_scored) ~ 1,
    data = df_model_away, 
    dist = d
  )
  cat("Distribution:", d, "AIC:", AIC(model), "\n")
}

# Also fit generalized gamma same way as done further in this script but now no covariates included
flexsurvreg(Surv(away_first_goal_minute, away_scored) ~ 1, data = df_model_away, dist = "gengamma")




full_formula_string <- paste("Surv(home_first_goal_minute, home_scored) ~", 
                             paste(smart_candidates, collapse = " + "))

full_f <- as.formula(full_formula_string)

full <- survreg(
  formula = full_f,
  data = df_model_home, 
  dist = "weibull"
)

null <- survreg(
  Surv(home_first_goal_minute, home_scored) ~ 1, 
  data = df_model_home, 
  dist = "weibull"
)

AIC_forwards_home <- stepAIC(null, scope = list(lower = null, upper = full), direction = "forward")
summary(AIC_forwards_home)

AIC_both_home <- stepAIC(null, scope = list(lower = null, upper = full), direction = "both")
summary(AIC_both_home)

AIC_backwards_home <- stepAIC(full, scope = list(lower = null, upper = full), direction = "backward")
summary(AIC_backwards_home)

AIC_backwards_home <- stepAIC(full, scope = list(lower = null, upper = full), direction = "both")
summary(AIC_backwards_home)

# Based on p-values
impact_tabel <- drop1(AIC_both_home, test = "Chisq")
print(impact_tabel)
variable_to_remove <- rownames(impact_tabel)[which.max(impact_tabel$`Pr(>Chi)`)]
AIC_refined_home <- update(AIC_both_home, . ~ . - Home_PSxG_90_for_roll5)
summary(AIC_refined_home)


impact_tabel <- drop1(AIC_refined_home, test = "Chisq")
print(impact_tabel)
variable_to_remove <- rownames(impact_tabel)[which.max(impact_tabel$`Pr(>Chi)`)]
variable_to_remove
AIC_refined_home_2 <- update(AIC_refined_home, . ~ . - Away_Passing_accuracy_for_roll5)
summary(AIC_refined_home_2)

impact_tabel <- drop1(AIC_refined_home_2, test = "Chisq")
print(impact_tabel)
variable_to_remove <- rownames(impact_tabel)[which.max(impact_tabel$`Pr(>Chi)`)]
variable_to_remove
AIC_refined_home_3 <- update(AIC_refined_home_2, . ~ . - Away_Number_passes_against_roll5)
summary(AIC_refined_home_3)

impact_tabel <- drop1(AIC_refined_home_3, test = "Chisq")
print(impact_tabel)
variable_to_remove <- rownames(impact_tabel)[which.max(impact_tabel$`Pr(>Chi)`)]
variable_to_remove

summary(AIC_refined_home_3)
# Search best model away
df_model_away <- df |>
  dplyr::select(Away, away_first_goal_minute, away_scored, all_of(candidate_vars), season) |>
  drop_na() |> 
  mutate(across(where(is.numeric) & !all_of(c("away_first_goal_minute", "away_scored")), 
                ~ as.numeric(scale(.))))
base_candidates <- candidate_vars[!grepl("_(15|30|45|60|75|90)_", candidate_vars)]

time_candidates <- candidate_vars[grepl("_(15|30|45|60|75|90)_", candidate_vars)]

time_candidates_clean <- time_candidates[grepl("_(xG|Goal|PSxG)_", time_candidates)]

final_time_vars <- time_candidates_clean[grepl("_(15|90)_", time_candidates_clean)]

smart_candidates <- c(base_candidates, final_time_vars)

full_formula_string <- paste("Surv(away_first_goal_minute, away_scored) ~", 
                             paste(smart_candidates, collapse = " + "))
full_f <- as.formula(full_formula_string)
full <- survreg(
  formula = full_f,
  data = df_model_away, 
  dist = "weibull"
)

null <- survreg(
  Surv(away_first_goal_minute, away_scored) ~ 1, 
  data = df_model_away, 
  dist = "weibull"
)

AIC_forwards <- stepAIC(null, scope = list(lower = null, upper = full), direction = "forward")
summary(AIC_forwards)

AIC_forwards_both_away <- stepAIC(null, scope = list(lower = null, upper = full), direction = "both")
summary(AIC_forwards_both_away)

impact_tabel_away <- drop1(AIC_forwards_both_away, test = "Chisq")
print(impact_tabel_away)
AIC_away_refined_1 <- update(AIC_forwards_both_away, . ~ . - Home_Passing_accuracy_for_roll5)
drop1(AIC_away_refined_1, test = "Chisq")
AIC_away_refined_2 <- update(AIC_away_refined_1, . ~ . - Home_Corners_for_roll5)
drop1(AIC_away_refined_2, test = "Chisq")
AIC_away_refined_3 <- update(AIC_away_refined_2, . ~ . - Away_Saves_for_roll5)
drop1(AIC_away_refined_3, test = "Chisq")
AIC_away_refined_4 <- update(AIC_away_refined_3, . ~ . - Away_Goal_15_against_roll5)
drop1(AIC_away_refined_4, test = "Chisq")
AIC_away_refined_5 <- update(AIC_away_refined_4, . ~ . - Home_Shots_on_Target_for_roll5)
drop1(AIC_away_refined_5, test = "Chisq")
summary(AIC_away_refined_5)

AIC_backwards_away <- stepAIC(full, scope = list(lower = null, upper = full), direction = "backward")
summary(AIC_backwards_away)

AIC_backwards_both_away <- stepAIC(full, scope = list(lower = null, upper = full), direction = "both")
summary(AIC_backwards_both_away)






# Cure models 

data <- data.frame( Y = df_model_away$away_first_goal_minute, 
                    Delta = df_model_away$away_scored,
                    Away = df_model_away$Away,
                    X1 = df_model_away$Home_Number_passes_for_roll5, 
                    X2 = df_model_away$Away_Number_passes_for_roll5,
                    X3 = df_model_away$Home_xG_against_roll5,
                    X4 = df_model_away$Away_xG_for_roll5,
                    X5 = df_model_away$Away_goals_for_roll5,
                    X6 = df_model_away$Away_Shots_on_Target_for_roll5,
                    X7 = df_model_away$Away_Corners_for_roll5,
                    X8 = df_model_away$Away_PSxG_90_for_roll5,
                    X9 = df_model_away$Away_Passing_accuracy_against_roll5,
                    X10 = df_model_away$Away_goals_against_roll5,
                    
                    # Home Team variables
                    X11 = df_model_away$Home_shots_against_roll5,
                    X12 = df_model_away$Home_xG_for_roll5,
                    X13 = df_model_away$Home_possesion_against_roll5,
                    X14 = df_model_away$Home_xG_15_for_roll5,
                    X15 = df_model_away$Home_xG_90_for_roll5,
                    X16 = df_model_away$Home_PSxG_15_against_roll5,
                    X17 = df_model_away$Home_PSxG_90_against_roll5,
                    X18 = as.numeric(df_model_away$Round_grouped == "Eredivisie"),
                    X19 = as.numeric(df_model_away$Round_grouped == "Jupiler Pro League"),
                    X20 = as.numeric(df_model_away$Round_grouped == "La Liga"),
                    X21 = as.numeric(df_model_away$Round_grouped == "Ligue 1"),
                    X22 = as.numeric(df_model_away$Round_grouped == "Premier League"),
                    X23 = as.numeric(df_model_away$Round_grouped == "Serie A")
                    )

# Comparing different baseline distributions with covariates based on AIC, using survreg so not for cure models
dist <- c("weibull", "exponential", "gaussian", "logistic","lognormal", "loglogistic")
for (d in dist) {
  model <- survreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                   data = data, 
                   dist = d)
  cat("Distribution:", d, "AIC:", AIC(model), "\n")
}

# Fit generalized gamma
flexsurvreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
            data = data, 
            dist = "gengamma")


# Changing time for censored observations with no goal to 120 to create plateau
data_120 <- data %>%
  mutate(Y = ifelse(Delta == 0, 120, Y))

set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                          cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                          data = data_120, 
                          model = "ph", 
                          Var = TRUE)
#Delete X1 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)

#Delete X9 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X16 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3 + X4 + X5 + X6 + X7 + X8 + X10 + X11 + X12 + X13 + X14 + X15 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X5 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3 + X4 + X6 + X7 + X8 + X10 + X11 + X12 + X13 + X14 + X15 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X4 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3  + X6 + X7 + X8 + X10 + X11 + X12 + X13 + X14 + X15 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X13 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3 + X6 + X7 + X8 + X10 + X11 + X12 + X14 + X15 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X8 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X3 + X6 + X7 + X10 + X11 + X12 + X14 + X15 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X3 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X10 + X11 + X12 + X14 + X15 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X15 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7  + X10 + X11 + X12 + X14 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete X14 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X10 + X11 + X12 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete x7 from incidence
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7  + X10 + X11 + X12  + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete x17 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7  + X10 + X11 + X12 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete x10 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X11 + X12 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete x12 from latency
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X11  + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)
#Delete x17 from incidence
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X11  + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)

#Delete x14 from incidence
set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X11  + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X15 + X16 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph", 
                    Var = TRUE)

data_110 <- data %>%
  mutate(Y = ifelse(Delta == 0, 110, Y))

set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X11  + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X15 + X16 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_110, 
                    model = "ph", 
                    Var = TRUE)

data_100 <- data %>%
  mutate(Y = ifelse(Delta == 0, 100, Y))

set.seed(123)
Cure_away <- smcure(formula = Surv(Y, Delta) ~ X2 + X6 + X7 + X11  + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~ X1 + X2 + X3 + X4 + X5 + X6  + X8 + X9 + X10 + X11 + X12 + X13 + X15 + X16 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_100, 
                    model = "ph", 
                    Var = TRUE)

table(df_model_away$away_scored)
table(df_model_home$home_scored)



Cure_away_subset <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, 
                    data = data_120, 
                    model = "ph",
                    emmax = 50,
                    #eps = 0.1,
                    Var = T)



# Make data for home
data_home <- data.frame(Y = df_model_home$home_first_goal_minute, 
                    Delta = df_model_home$home_scored,
                    Home = df_model_home$Home,
                    X1 = df_model_home$Home_Number_passes_for_roll5, 
                    X2 = df_model_home$Away_Number_passes_for_roll5,
                    X3 = df_model_home$Home_xG_for_roll5,
                    X4 = df_model_home$Away_xG_against_roll5, 
                    X5 = df_model_home$Home_goals_for_roll5,
                    X6 = df_model_home$Away_goals_for_roll5,
                    X7 = df_model_home$Home_shots_for_roll5,
                    X8 = df_model_home$Home_PSxG_15_for_roll5,
                    X9 = df_model_home$Away_xG_for_roll5,
                    X10 = df_model_home$Away_Shots_on_Target_against_roll5, 
                    X11 = df_model_home$Home_PSxG_15_against_roll5, 
                    X12 = df_model_home$Home_xG_15_against_roll5, 
                    X13 = df_model_home$Home_Passing_accuracy_for_roll5, 
                    X14 = as.numeric(df_model_home$Round_grouped == "Eredivisie"), 
                    X15 = as.numeric(df_model_home$Round_grouped == "Jupiler Pro League"),
                    X16 = as.numeric(df_model_home$Round_grouped == "La Liga"),
                    X17 = as.numeric(df_model_home$Round_grouped == "Ligue 1"),
                    X18 = as.numeric(df_model_home$Round_grouped == "Premier League"),
                    X19 = as.numeric(df_model_home$Round_grouped == "Serie A"), 
                    X20 = as.numeric(df_model_home$home_formation_grouped == "3-5-2"), 
                    X21 = as.numeric(df_model_home$home_formation_grouped == "4-3-3"),
                    X22 = as.numeric(df_model_home$home_formation_grouped == "4-4-2"),
                    X23 = as.numeric(df_model_home$home_formation_grouped == "4-2-3-1"),
                    X24 = as.numeric(df_model_home$home_formation_grouped == "5-3-2"),
                    X25 = as.numeric(df_model_home$home_formation_grouped == "4-5-1"),
                    X26 = as.numeric(df_model_home$home_formation_grouped == "5-4-1")
                    )


# Comparing different baseline paramteric models for survreg based on AIC
model_weibull <- survreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, data = data_home, dist = "weibull")

model_exponential <- survreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, data = data_home, dist = "exponential")
model_loglogistic <- survreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, data = data_home, dist = "loglogistic")
model_lognormal <- survreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, data = data_home, dist = "lognormal")
model_gen_gamma <- flexsurvreg(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, data = data_home, dist = "gengamma")
aic_values <- c(Weibull = AIC(model_weibull), 
                Exponential = AIC(model_exponential), 
                LogLogistic = AIC(model_loglogistic), 
                LogNormal = AIC(model_lognormal))
print(aic_values)

                   
data_home_120 <- data_home %>%
  mutate(Y = ifelse(Delta == 0, 120, Y))

set.seed(123)
Cure_home <- smcure_debug(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                                 cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                                 data = data_home_120, 
                                 model = "ph",
                                 emmax = 50,
                                 #eps = 0.1,
                                 Var = T)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                          cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                          data = data_home_120, 
                          model = "ph",
                          emmax = 50,
                          #eps = 0.1,
                          Var = T)
summary(Cure_home_smcure)
# Delete X6 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X10 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X7 + X8 + X9 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X8 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X7 + X9 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X7 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X9 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X13 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X9 + X12 + X11 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X11 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X9 + X12 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X12 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X9 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X2 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X3 + X4 + X5 + X9 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)
# Delete X5 in latency part
set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X3 + X4 + X9 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_120, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)

data_home_110 <- data_home %>%
  mutate(Y = ifelse(Delta == 0, 110, Y))

set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X3 + X4 + X9 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_110, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)

data_home_100 <- data_home %>%
  mutate(Y = ifelse(Delta == 0, 100, Y))

set.seed(123)
Cure_home_smcure <- smcure(formula = Surv(Y, Delta) ~ X1 + X3 + X4 + X9 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           cureform = ~  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, 
                           data = data_home_100, 
                           model = "ph",
                           emmax = 50,
                           #eps = 0.1,
                           Var = T)



# Checking ph assumption via Schoenfeld residuals
model_home <- coxph(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23 + X24 + X25 + X26, data = data_home)
test_ph <- cox.zph(model_home)
print(test_ph)
# Plot significant variables schoenfeld residuals against time
test_ph_vis <- cox.zph(model_home, transform = "identity") # This is for visualizing, ensures equal distances between points of different times instead of based on the number of events
# Makes sure there is no gap between 45 and 46 minutes due to the handling of stoppage time goals
plot(test_ph_vis[18], main = "Schoenfeld Residuals: Premier League", ylab = "Schoenfeld residuals for Premier League")
abline(h = 0, col = "red", lty = 2)
?cox.zph
plot(test_ph_vis[21], main = "Schoenfeld Residuals: 4-3-3", ylab = "Schoenfeld residuals for 4-3-3")
abline(h = 0, col = "red", lty = 2)



model_away <- coxph(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + X21 + X22 + X23, data = data)
test_ph_away <- cox.zph(model_away)
print(test_ph_away)
test_ph_vis_away <- cox.zph(model_away, transform = "identity") # This is for visualizing, ensures equal distances between points of different times instead of based on the number of events
plot(test_ph_vis_away[7], lab = c(5, 5, 7), xlim = c(0, 90), ylab = "Schoenfeld residuals: Away corners for", main = "Schoenfeld residuals for Away corners for")
abline(h = 0, col = "red", lty = 2)
plot(test_ph_vis_away[9], lab = c(5, 5, 7), xlim = c(0, 90), ylab = "Schoenfeld residuals: Away passing accuracy against", main = "Schoenfeld residuals for Away passing accuracy against")
abline(h = 0, col = "red", lty = 2)



data_scored <- data_home %>% filter(Delta == 1)

model_scored <- coxph(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
                        X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + 
                        X21 + X22 + X23 + X24 + X25 + X26, 
                      data = data_scored)

test_ph_scored <- cox.zph(model_scored)
print(test_ph_scored)


data_scored <- data %>% filter(Delta == 1)

model_scored <- coxph(Surv(Y, Delta) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
                        X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 + 
                        X21 + X22 + X23 , 
                      data = data_scored)

test_ph_scored <- cox.zph(model_scored)
print(test_ph_scored)

plot(test_ph_scored[7])

# Diagnstics plots
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(ggplot2)

library(survminer)
library(survMisc)
data_plot <- data_home %>%
  mutate(League = case_when(
    X14 == 1 ~ "Eredivisie",
    X15 == 1 ~ "Jupiler Pro League",
    X16 == 1 ~ "La Liga",
    X17 == 1 ~ "Ligue 1",
    X18 == 1 ~ "Premier League",
    X19 == 1 ~ "Serie A",
    TRUE     ~ "Bundesliga"
  ))

fit_league <- survfit(Surv(Y, Delta) ~ League, data = data_plot)

summary_fit_league <- summary(fit_league)
df_fit_league <- data.frame(
  time = summary_fit_league$time,
  strata = gsub("League=", "", summary_fit_league$strata),
  cumhaz = summary_fit_league$cumhaz
)

df_wide_league <- df_fit_league %>%
  pivot_wider(
    names_from = strata,
    values_from = cumhaz,
    names_prefix = "H_"
  ) %>%
  arrange(time) %>%
  fill(starts_with("H_"), .direction = "down")

p_cumhaz <- ggsurvplot(
  fit_league,
  data = data_plot,
  fun = "cumhaz",
  conf.int = FALSE,
  legend.title = "Competition",
  xlab = "Time until first goal (min)",
  ylab = "Cumulative hazard H(t)",
  ggtheme = theme_minimal()
)
print(p_cumhaz)

p_logMinlog_linear <- df_fit_league %>%
  filter(cumhaz > 0) %>%
  ggplot(aes(x = time, y = log(cumhaz), color = strata)) +
  geom_step(linewidth = 0.8) +
  labs(
    x = "Time until first goal (min)",
    y = "log Cumulative Hazard log H(t)",
    title = "Log Cumulative Hazard vs Linear Time (PH Check)",
    color = "Competition"
  ) +
  scale_x_continuous(breaks = seq(0, 90,  by = 15)) + 
  theme_minimal() +
  theme(legend.position = "right")

print(p_logMinlog_linear)

diff_df <- df_wide_league %>%
  filter(`H_Bundesliga` > 0) %>%
  mutate(across(starts_with("H_"), ~log(.x), .names = "log_{.col}")) %>%
  mutate(
    diff_Eredivisie = `log_H_Eredivisie` - `log_H_Bundesliga`,
    diff_PL = `log_H_Premier League` - `log_H_Bundesliga`,
    diff_LaLiga = `log_H_La Liga` - `log_H_Bundesliga`
  ) %>%
  dplyr::select(time, starts_with("diff_")) %>%
  pivot_longer(cols = -time, names_to = "Competition", values_to = "logdiff")

p_logHazDiff <- ggplot(diff_df, aes(x = time, y = logdiff, color = Competition)) +
  geom_step() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Time (min)",
    y = "log(H_j) - log(H_ref)",
    title = "Log-hazard differences vs Bundesliga"
  ) +
  theme_minimal()
print(p_logHazDiff)

hazard_vs_reference <- df_wide_league %>%
  dplyr::select(time, starts_with("H_")) %>%
  tidyr::pivot_longer(
    cols = c(`H_Eredivisie`, `H_Jupiler Pro League`, `H_La Liga`, 
             `H_Ligue 1`, `H_Premier League`, `H_Serie A`),
    names_to = "Competition_Var",
    values_to = "H_j"
  ) %>%
  dplyr::mutate(
    Competition = gsub("H_", "", Competition_Var),
    H_ref = `H_Bundesliga` 
  )

p_compare_ref <- ggplot(hazard_vs_reference, aes(x = H_ref, y = H_j, color = Competition)) +
  geom_step(size = 0.8) +
  geom_abline(slope = 1, linetype = "dotted", color = "grey") + 
  labs(
    x = expression(hat(H)[Bundesliga](t)),
    y = expression(hat(H)[j](t)),
    title = "Proportionality Check: Cum. Hazard vs. Reference (Bundesliga)"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p_compare_ref)


data_tactics <- data_home %>%
  mutate(Formation = case_when(
    X20 == 1 ~ "3-5-2",
    X21 == 1 ~ "4-3-3",
    X22 == 1 ~ "4-4-2",
    X23 == 1 ~ "4-2-3-1",
    X24 == 1 ~ "5-3-2",
    X25 == 1 ~ "4-5-1",
    X26 == 1 ~ "5-4-1",
    TRUE     ~ "3-4-3"
  ))

fit_tactics <- survfit(Surv(Y, Delta) ~ Formation, data = data_tactics)

summary_fit_tactics <- summary(fit_tactics)
df_tactics <- data.frame(
  time = summary_fit_tactics$time,
  strata = gsub("Formation=", "", summary_fit_tactics$strata),
  cumhaz = summary_fit_tactics$cumhaz
)

p_tactics_ph <- df_tactics %>%
  filter(cumhaz > 0) %>%
  ggplot(aes(x = time, y = log(cumhaz), color = strata)) +
  geom_step(linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, 90, 15)) +
  labs(
    x = "Time to first goal (min)",
    y = "log Cumulatieve Hazard log H(t)",
    title = "PH Check: Starting formations (Home Model)",
    color = "Formation"
  ) +
  theme_minimal()

print(p_tactics_ph)


df_wide_tactics <- df_tactics %>%
  pivot_wider(
    names_from = strata,
    values_from = cumhaz,
    names_prefix = "H_"
  ) %>%
  arrange(time) %>%
  fill(starts_with("H_"), .direction = "down")

p_tactics_cumhaz <- ggsurvplot(
  fit_tactics,
  data = data_tactics,
  fun = "cumhaz",
  legend.title = "Formatie",
  xlab = "Tijd (min)",
  ylab = "Cumulatieve Hazard H(t)",
  ggtheme = theme_minimal()
)
print(p_tactics_cumhaz)

diff_df_tactics <- df_wide_tactics %>%
  filter(`H_3-4-3` > 0) %>%
  mutate(across(starts_with("H_"), ~log(.x), .names = "log_{.col}")) %>%
  mutate(
    diff_433 = `log_H_4-3-3` - `log_H_3-4-3`,
    diff_442 = `log_H_4-4-2` - `log_H_3-4-3`,
    diff_541 = `log_H_5-4-1` - `log_H_3-4-3`,
    diff_4231 = `log_H_4-2-3-1` - `log_H_3-4-3`
  ) %>% select(time, starts_with("diff_")) %>%
  pivot_longer(cols = -time, names_to = "Formatie", values_to = "logdiff")

p_tactics_diff <- ggplot(diff_df_tactics, aes(x = time, y = logdiff, color = Formatie)) +
  geom_step(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Tijd (min)",
    y = expression(log(H[j](t)) - log(H["3-4-3"](t))),
    title = "Log-hazard differences versus 3-4-3"
  ) +
  theme_minimal()
print(p_tactics_diff)

hazard_vs_ref_tactics <- df_wide_tactics %>%
  select(time, starts_with("H_")) %>%
  pivot_longer(
    cols = -c(time, `H_3-4-3`),
    names_to = "Formatie_Var",
    values_to = "H_j"
  ) %>%
  mutate(H_ref = `H_3-4-3`, Formatie = gsub("H_", "", Formatie_Var))

p_tactics_linear <- ggplot(hazard_vs_ref_tactics, aes(x = H_ref, y = H_j, color = Formatie)) +
  geom_step() +
  geom_abline(slope = 1, linetype = "dotted") +
  labs(
    x = expression(hat(H)["3-4-3"](t)),
    y = expression(hat(H)[j](t)),
    title = "Linearity Check: Formatie H_j vs H_3-4-3"
  ) +
  theme_minimal()
print(p_tactics_linear)



# Fitting copula based on the code from the CopulaCenR package
final_model_home <- survreg(
  Surv(home_first_goal_minute, home_scored) ~ 
    # Team Statistics (Roll 5)
    Home_Number_passes_for_roll5 +          # Home Number of Passes
    Away_Number_passes_for_roll5 +          # Away Number of Passes
    Home_xG_for_roll5 +                     # Home xG for
    Away_xG_against_roll5 +                 # Away xG Against
    Away_goals_for_roll5 +                  # Away Goals for
    Home_goals_for_roll5 +                  # Home Goals for
    Home_shots_for_roll5 +                  # Home Shots for
    Home_PSxG_15_for_roll5 +                # Home PSxG (First 15min) for
    Away_xG_for_roll5 +                     # Away xG for
    Away_Shots_on_Target_against_roll5 +     # Away SoT Against
    Home_PSxG_15_against_roll5 +            # Home PSxG (First 15min) Against
    Home_xG_15_against_roll5 +              # Home xG (First 15min) Against
    Home_Passing_accuracy_for_roll5 +       # Home Passing Accuracy for
    # Formations & Competition
    home_formation_grouped + 
    Round_grouped,
  data = df_model_home, 
  dist = "weibull"
)

vars_away_model <- c(
  # Away Team variabelen
  "Away_Number_passes_for_roll5",
  "Away_xG_for_roll5",
  "Away_goals_for_roll5",
  "Away_Shots_on_Target_for_roll5",
  "Away_Corners_for_roll5",
  "Away_PSxG_90_for_roll5",
  "Away_Passing_accuracy_against_roll5",
  "Away_goals_against_roll5",
  
  # Home Team variabelen
  "Home_xG_against_roll5",
  "Home_shots_against_roll5",
  "Home_Number_passes_for_roll5",
  "Home_xG_for_roll5",
  "Home_possesion_against_roll5",
  "Home_xG_15_for_roll5",
  "Home_xG_90_for_roll5",
  "Home_PSxG_15_against_roll5",
  "Home_PSxG_90_against_roll5"
)
away_formula <- as.formula(paste(
  "Surv(away_first_goal_minute, away_scored) ~", 
  paste(c(vars_away_model, "Round_grouped"), collapse = " + ")
))

final_model_away <- survreg(away_formula, data = df_model_away, dist = "weibull")

vars_home <- all.vars(formula(final_model_home))[-1] 
vars_away <- all.vars(formula(final_model_away))[-1]        
all_vars_needed <- unique(c(vars_home, vars_away))

# Based on the examples in the documentation of survreg, calculating the shape and scale parameters for the weibull distribution per game from survregs output
get_weibull_params <- function(model, data) {
  lp <- predict(model, newdata = data, type = "lp")
  scale_param <- exp(lp)
  shape_param <- 1 / model$scale
  list(scale = scale_param, shape = shape_param)
}


df_sync <- df %>%
  dplyr::select(
    home_first_goal_minute, home_scored,
    away_first_goal_minute, away_scored,
    all_of(all_vars_needed)
  ) %>%
  drop_na() %>%
  mutate(across(where(is.numeric) & !all_of(c(
    "home_first_goal_minute", "home_scored",
    "away_first_goal_minute", "away_scored"
  )), scale))

t1 <- df_sync$home_first_goal_minute
t2 <- df_sync$away_first_goal_minute
d1 <- df_sync$home_scored   # 1 = goal, 0 = censored
d2 <- df_sync$away_scored

params_h <- get_weibull_params(final_model_home, df_sync)
params_a <- get_weibull_params(final_model_away, df_sync)

# Generating S
u <- pweibull(t1, shape = params_h$shape, scale = params_h$scale, lower.tail = FALSE)
v <- pweibull(t2, shape = params_a$shape, scale = params_a$scale, lower.tail = FALSE)

# generating F
f1 <- dweibull(t1, shape = params_h$shape, scale = params_h$scale)
f2 <- dweibull(t2, shape = params_a$shape, scale = params_a$scale)

# Based on the likelihood construction in the code of CopulaCenR
CopulaCenR::rc_par_copula
CopulaCenR:::rc_copula_log_lik


rc_copula_log_lik_RAW <- function(eta, u1, u2, status1, status2, f1, f2, copula = "Clayton") {
  
  if (copula == "Clayton" && eta <= 0) return(1e10)
  if (copula == "Gumbel" && eta < 1) return(1e10)
  if (copula == "Joe" && eta < 1) return(1e10)
  
  if (copula == "Clayton") {
    C_val <- (u1^(-eta) + u2^(-eta) - 1)^(-1/eta)
    c_u1_val <- u1^(-eta - 1) * (u1^(-eta) + u2^(-eta) - 1)^(-1/eta - 1)
    c_u2_val <- u2^(-eta - 1) * (u1^(-eta) + u2^(-eta) - 1)^(-1/eta - 1)
    c_val <- (1 + eta) * (u1 * u2)^(-eta - 1) * (u1^(-eta) + u2^(-eta) - 1)^(-1/eta - 2)
  }
  
  if (copula == "Gumbel") {
    gh_F <- exp(-((-log(u1))^eta + (-log(u2))^eta)^(1/eta))
    C_val <- gh_F
    c_u1_val <- gh_F * ((-log(u1))^eta + (-log(u2))^eta)^(1/eta - 1) * (-log(u1))^(eta - 1)/u1
    c_u2_val <- gh_F * ((-log(u1))^eta + (-log(u2))^eta)^(1/eta - 1) * (-log(u2))^(eta - 1)/u2
    c_val <- ((-log(u1))^(eta - 1) * (-log(u2))^(eta - 1) * gh_F * ((-log(u1))^eta + (-log(u2))^eta)^(2/eta - 2))/(u1 * u2) - gh_F * ((-log(u1))^eta + (-log(u2))^eta)^(1/eta - 2) * (1/(u1 * u2)) * (1/eta - 1) * eta * (-log(u1))^(eta - 1) * (-log(u2))^(eta - 1)
  }
  
  if (copula == "Joe") {
    C_val <- 1 - ((1 - u1)^eta + (1 - u2)^eta - ((1 - u1)^eta) * ((1 - u2)^eta))^(1/eta)
    c_u1_val <- ((1 - u1)^eta + (1 - u2)^eta - ((1 - u1)^eta) * ((1 - u2)^eta))^(1/eta - 1) * ((1 - u1)^(eta - 1) - (1 - u1)^(eta - 1) * (1 - u2)^eta)
    c_u2_val <- ((1 - u1)^eta + (1 - u2)^eta - ((1 - u1)^eta) * ((1 - u2)^eta))^(1/eta - 1) * ((1 - u2)^(eta - 1) - (1 - u2)^(eta - 1) * (1 - u1)^eta)
    c_val <- ((1 - u1)^eta + (1 - u2)^eta - ((1 - u1)^eta) * ((1 - u2)^eta))^(1/eta - 2) * (1/eta - 1) * ((1 - u1)^(eta - 1) - (1 - u1)^(eta - 1) * (1 - u2)^eta) * ((1 - u2)^(eta - 1) - (1 - u2)^(eta - 1) * (1 - u1)^eta) * (-eta) + ((1 - u1)^eta + (1 - u2)^eta - ((1 - u1)^eta) * ((1 - u2)^eta))^(1/eta - 1) * (eta * ((1 - u1)^(eta - 1)) * ((1 - u2)^(eta - 1)))
  }
  
  if (copula == "Frank") {
    C_val <- (-1 / eta) * log(1 + ((exp(-eta * u1) - 1) * (exp(-eta * u2) - 1)) / (exp(-eta) - 1))
    
    c_u1_val <- ((exp(-eta * u2) - 1) * exp(-eta * u1)) / 
      ((exp(-eta) - 1) + (exp(-eta * u1) - 1) * (exp(-eta * u2) - 1))
    
    c_u2_val <- ((exp(-eta * u1) - 1) * exp(-eta * u2)) / 
      ((exp(-eta) - 1) + (exp(-eta * u1) - 1) * (exp(-eta * u2) - 1))
    
    c_val <- (-eta * (exp(-eta) - 1) * exp(-eta * (u1 + u2))) / 
      (((exp(-eta) - 1) + (exp(-eta * u1) - 1) * (exp(-eta * u2) - 1))^2)
  }
  
  term1 <- ifelse((status1 == 0) & (status2 == 0), C_val, 1)
  term1 <- log(abs(term1))
  
  term2 <- c_u1_val * f1
  term2 <- ifelse((status1 == 1) & (status2 == 0), term2, 1)
  term2 <- log(abs(term2))
  
  term3 <- c_u2_val * f2
  term3 <- ifelse((status1 == 0) & (status2 == 1), term3, 1)
  term3 <- log(abs(term3))
  
  term4 <- c_val * f1 * f2
  term4 <- ifelse((status1 == 1) & (status2 == 1), term4, 1)
  term4[term4 < 0] <- 1  
  term4 <- log(abs(term4))
  
  logL <- (-1) * sum(term1 + term2 + term3 + term4)
  return(logL)
}



fit_biv_weibull_copula_RAW <- function(u1, u2, status1, status2, f1, f2,
                                       family = c("Clayton", "Frank", "Gumbel", "Joe")) {
  
  family <- match.arg(tools::toTitleCase(tolower(family)), 
                      c("Clayton", "Frank", "Gumbel", "Joe"))
  
  start_val <- if(family %in% c("Gumbel", "Joe")) 1.1 else if(family == "Clayton") 0.1 else 0.5
  
  fit <- optim(par = start_val, 
               fn = rc_copula_log_lik_RAW, 
               u1 = u1, u2 = u2, 
               status1 = status1, status2 = status2, 
               f1 = f1, f2 = f2, 
               copula = family, 
               method = "BFGS", 
               hessian = TRUE,
               control = list(ndeps = rep(1e-6, 1)))
  
  eta_hat <- fit$par
  
  hess_inv <- try(solve(fit$hessian)[1,1], silent = TRUE)
  
  if (inherits(hess_inv, "try-error") || hess_inv < 0) { 
    se_eta <- NA
    z_val <- NA
  } else {
    var_eta <- hess_inv
    se_eta <- sqrt(var_eta)
    

    if (family %in% c("Gumbel", "Joe")) {
      z_val <- (eta_hat - 1) / se_eta   
    } else {
      z_val <- eta_hat / se_eta         
    }
  }
  
  # Kendall's Tau 
  tau <- switch(family,
                "Clayton"  = eta_hat / (eta_hat + 2),
                "Frank"    = {
                  if(abs(eta_hat) < 1e-5) { 0 } else {
                    debye_fun <- function(t) t / (exp(t) - 1)
                    D1 <- (1/eta_hat) * integrate(debye_fun, 1e-8, abs(eta_hat))$value
                    1 - 4/eta_hat * (1 - D1)
                  }
                },
                "Gumbel"   = 1 - 1/eta_hat,
                "Joe"      = 1 - 4 / (eta_hat * (eta_hat + 1) * (eta_hat + 3))
  )
  
  list(
    family = family,
    eta    = eta_hat,
    se_eta = se_eta,
    z_eta  = z_val,
    tau    = tau,
    loglik = -fit$value,
    AIC    = -2 * (-fit$value) + 2
  )
}


model_families <- c("Clayton", "Frank", "Gumbel", "Joe")

results_list_RAW <- lapply(model_families, function(fam) {
  tryCatch(fit_biv_weibull_copula_RAW(u1 = u, u2 = v, 
                                      status1 = d1, status2 = d2, 
                                      f1 = f1, f2 = f2, 
                                      family = fam),
           error = function(e) list(family = fam, error = "CRASHED"))
})

print(do.call(rbind, results_list_RAW))






# Frailty models, check variance for full marginal models
df_model_home <- df_model_home %>%
  mutate(
    Home = as.factor(Home),
    home_formation_grouped = as.factor(home_formation_grouped),
    Round_grouped = as.factor(Round_grouped),
    Home_s = as.factor(paste0(Home, season, sep = " "))
  )

# Comparison one frailty per season versus one frailty per team for home
Frailty_home <- parfm(
  Surv(home_first_goal_minute, home_scored) ~ 1,
  cluster = "Home", 
  data = df_model_home, 
  dist = "weibull", 
  frailty = "lognormal"
)
cat("Frailty per team (lognormal):\n", AIC(Frailty_home), "\n")

# Per season
Frailty_home <- parfm(
  Surv(home_first_goal_minute, home_scored) ~ 1,
  cluster = "Home_s", 
  data = df_model_home, 
  dist = "weibull", 
  frailty = "lognormal"
)
cat("Frailty per team (lognormal):\n", AIC(Frailty_home), "\n")


frailties <- c("none", "gamma", "ingau", "possta",
               "lognormal")

for (f in frailties) {
  model <- parfm(
    Surv(home_first_goal_minute, home_scored) ~ 1,
    cluster = "Home", 
    data = df_model_home, 
    dist = "weibull", 
    frailty = f
  )
  print(paste("Frailty:", f))
  print(model)
  print(paste("AIC:", AIC(model)))
}

parfm 

par_home_nc <- parfm(
  Surv(home_first_goal_minute, home_scored) ~ 1,
  cluster = "Home", 
  data = df_model_home, 
  dist = "weibull", 
  frailty = "loglogistic"
) # This code produces an error. R documentation shows that it accepts loglogistic as an input,
# but doesn't provide support for it later in the code of the function once it has to specify nFpar (number of Frailty parameters I presume)

Frailty_home <- parfm(
  Surv(home_first_goal_minute, home_scored) ~ 1,
  cluster = "Home", 
  data = df_model_home, 
  dist = "weibull", 
  frailty = "lognormal"
)
Frailty_home

Frailty_home <- parfm(
  Surv(home_first_goal_minute, home_scored) ~ 
    # Team Statistics (Roll 5)
    Home_Number_passes_for_roll5 +          # Home Number of Passes
    Away_Number_passes_for_roll5 +          # Away Number of Passes
    Home_xG_for_roll5 +                     # Home xG for
    Away_xG_against_roll5 +                 # Away xG Against
    Away_goals_for_roll5 +                  # Away Goals for
    Home_goals_for_roll5 +                  # Home Goals for
    Home_shots_for_roll5 +                  # Home Shots for
    Home_PSxG_15_for_roll5 +                # Home PSxG (First 15min) for
    Away_xG_for_roll5 +                     # Away xG for
    Away_Shots_on_Target_against_roll5 +     # Away SoT Against
    Home_PSxG_15_against_roll5 +            # Home PSxG (First 15min) Against
    Home_xG_15_against_roll5 +              # Home xG (First 15min) Against
    Home_Passing_accuracy_for_roll5 +       # Home Passing Accuracy for
    # Formations & Competition
    home_formation_grouped + 
    Round_grouped,
    cluster = "Home", 
    data = df_model_home, 
    dist = "weibull", 
    frailty = "lognormal", 
    method = "L-BFGS-B"
)


# Away 

df_model_away <- df_model_away %>%
  mutate(
    Away = as.factor(Away),
    Round_grouped = as.factor(Round_grouped),
    Away_s = as.factor(paste0(Away, season, sep = " "))
  )
# Comparison one frailty per season versus one frailty per team for home
Frailty_away <- parfm(
  Surv(away_first_goal_minute, away_scored) ~ 1,
  cluster = "Away", 
  data = df_model_away, 
  dist = "weibull", 
  frailty = "lognormal"
)
cat("Frailty per team (lognormal):\n", AIC(Frailty_away), "\n")

# Per season
Frailty_away <- parfm(
  Surv(away_first_goal_minute, away_scored) ~ 1,
  cluster = "Away_s", 
  data = df_model_away, 
  dist = "weibull", 
  frailty = "lognormal"
)
cat("Frailty per team (lognormal):\n", AIC(Frailty_away), "\n")

for (f in frailties) {
  model <- parfm(
    Surv(away_first_goal_minute, away_scored) ~ 1,
    cluster = "Away", 
    data = df_model_away, 
    dist = "weibull", 
    frailty = f
  )
  print(paste("Frailty:", f))
  print(model)
  print(paste("AIC:", AIC(model)))
}


vars_away_model <- c(
  # Away Team variabelen
  "Away_Number_passes_for_roll5",
  "Away_xG_for_roll5",
  "Away_goals_for_roll5",
  "Away_Shots_on_Target_for_roll5",
  "Away_Corners_for_roll5",
  "Away_PSxG_90_for_roll5",
  "Away_Passing_accuracy_against_roll5",
  "Away_goals_against_roll5",
  
  # Home Team variabelen
  "Home_xG_against_roll5",
  "Home_shots_against_roll5",
  "Home_Number_passes_for_roll5",
  "Home_xG_for_roll5",
  "Home_possesion_against_roll5",
  "Home_xG_15_for_roll5",
  "Home_xG_90_for_roll5",
  "Home_PSxG_15_against_roll5",
  "Home_PSxG_90_against_roll5"
)
away_formula <- as.formula(paste(
  "Surv(away_first_goal_minute, away_scored) ~", 
  paste(c(vars_away_model, "Round_grouped"), collapse = " + ")
))

Frailty_model_away <- parfm(away_formula, cluster = "Away", 
                          data = df_model_away, 
                          dist = "weibull", 
                          frailty = "lognormal", 
                          method = "L-BFGS-B")

Frailty_model_away <- parfm(Surv(away_first_goal_minute, away_scored) ~ 1, cluster = "Away", 
                            data = df_model_away, 
                            dist = "weibull", 
                            frailty = "lognormal", 
                            method = "L-BFGS-B")


