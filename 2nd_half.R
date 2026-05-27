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
df_2nd <- df |> 
  filter(first_goal_minute >= 46) |> 
  mutate(goal_second_home = (home_first_goal_minute - 45), 
         goal_second_away = away_first_goal_minute - 45)
df_2nd |> 
  filter(home_scored == 1) |> 
  group_by(goal_second_home) |> 
  summarise(n = n()) |> 
  ggplot(aes(x = goal_second_home, y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "Minute of goal second half", y = "Count of matches") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 45, by = 5)) +
  ggtitle("Distribution of the second goal scored by the home team in the second half")
names(df_2nd)

# Make first half covariates by adding all quarter statistics of first half
df_2nd_stats <- df_2nd |> 
  mutate(Home_xG_1st = Home_xG_15 + Home_xG_30 + Home_xG_45,
         Home_xG_for_roll_2nd = Home_xG_60_for_roll5 + Home_xG_75_for_roll5 + Home_xG_90_for_roll5,
         Home_PSxG_for_roll_2nd = Home_PSxG_60_for_roll5 + Home_PSxG_75_for_roll5 + Home_PSxG_90_for_roll5,
         Home_goals_for_roll_2nd = Home_Goal_60_for_roll5 + Home_Goal_75_for_roll5 + Home_Goal_90_for_roll5,
         Home_Saved_for_roll_2nd = Home_Saved_60_for_roll5 + Home_Saved_75_for_roll5 + Home_Saved_90_for_roll5,
         Home_Blocked_for_roll_2nd = Home_Blocked_60_for_roll5 + Home_Blocked_75_for_roll5 + Home_Blocked_90_for_roll5,
         Home_Off_Target_for_roll_2nd = Home_OffTarget_60_for_roll5 + Home_OffTarget_75_for_roll5 + Home_OffTarget_90_for_roll5,
         Home_Woodwork_for_roll_2nd = Home_Woodwork_60_for_roll5 + Home_Woodwork_75_for_roll5 + Home_Woodwork_90_for_roll5,
         Home_PSxG_1st = Home_PSxG_15+ Home_PSxG_30 + Home_PSxG_45,
         Home_Saved_1st = Home_Saved_15 + Home_Saved_30 + Home_Saved_45,
         Home_Blocked_1st = Home_Blocked_15 + Home_Blocked_30 + Home_Blocked_45,
         Home_Off_Target_1st = Home_OffTarget_15 + Home_OffTarget_30 + Home_OffTarget_45,
         Home_Woodwork_1st = Home_Woodwork_15 + Home_Woodwork_30 + Home_Woodwork_45,
         Away_xG_1st = Away_xG_15 + Away_xG_30 + Away_xG_45,
         Away_PSxG_1st = Away_PSxG_15 + Away_PSxG_30 + Away_PSxG_45,
         Away_Saved_1st = Away_Saved_15 + Away_Saved_30 + Away_Saved_45,
         Away_Blocked_1st = Away_Blocked_15 + Away_Blocked_30 + Away_Blocked_45,
         Away_Off_Target_1st = Away_OffTarget_15 + Away_OffTarget_30 + Away_OffTarget_45,
         Away_Woodwork_1st = Away_Woodwork_15 + Away_Woodwork_30 + Away_Woodwork_45,
         Away_xG_for_roll_2nd = Away_xG_60_for_roll5 + Away_xG_75_for_roll5 + Away_xG_90_for_roll5,
         Away_PSxG_for_roll_2nd = Away_PSxG_60_for_roll5 + Away_PSxG_75_for_roll5 + Away_PSxG_90_for_roll5,
         Away_goals_for_roll_2nd = Away_Goal_60_for_roll5 + Away_Goal_75_for_roll5 + Away_Goal_90_for_roll5,
         Away_Saved_for_roll_2nd = Away_Saved_60_for_roll5 + Away_Saved_75_for_roll5 + Away_Saved_90_for_roll5,
         Away_Blocked_for_roll_2nd = Away_Blocked_60_for_roll5 + Away_Blocked_75_for_roll5 + Away_Blocked_90_for_roll5,
         Away_Off_Target_for_roll_2nd = Away_OffTarget_60_for_roll5 + Away_OffTarget_75_for_roll5 + Away_OffTarget_90_for_roll5,
         Away_Woodwork_for_roll_2nd = Away_Woodwork_60_for_roll5 + Away_Woodwork_75_for_roll5 + Away_Woodwork_90_for_roll5)


fit_null_home <- survfit(Surv(goal_second_home, home_scored) ~  1, data = df_2nd_stats)
plot(fit_null_home, xlab = "Minute of goal second half", ylab = "Survival probability", main = "Kaplan-Meier Curve for Home Team Scoring in the Second Half", col = "blue")
summary(fit_null_home)
fit_null_away <- survfit(Surv(goal_second_away, away_scored) ~  1, data = df_2nd_stats)
lines(fit_null_away, xlab = "Minute of goal second half", ylab = "Survival probability", main = "Kaplan-Meier Curve for Away Team Scoring in the Second Half", col = "red")
summary(fit_null_away)

print(fit_null_home, print.rmean = TRUE) #45% censored
print(fit_null_away, print.rmean = TRUE) # 50.03% censored

survdiff(Surv(goal_second_home, home_scored) ~ Round_grouped + strata(home_formation_grouped), data = df_2nd_stats)
survdiff(Surv(goal_second_away, away_scored) ~ Round_grouped + strata(away_formation_grouped), data = df_2nd_stats)

# Scale variables 
base_vars <- c(
  "xG", "goals", "possesion", 
  "Number_passes", "Passing_accuracy", 
  "shots", "Shots_on_Target", 
  "Saves", "Corners"
)

times <- c("15", "30", "45", "60", "75", "90")
metrics <- c("xG", "PSxG", "Goal", "Saved", "Blocked", "OffTarget", "Woodwork")
metrics_subset_1st <- c("xG", "PSxG", "Saved", "Blocked", "Off_Target", "Woodwork")
metrics_subset_2nd <- c("xG", "PSxG", "goals", "Saved", "Blocked", "Off_Target", "Woodwork")
time_vars <- outer(metrics, times, paste, sep="_") |> as.vector()

all_stat_names <- c(base_vars, time_vars)

suffix <- "_roll5" 

candidate_vars <- c(
  
  paste0("Home_", base_vars, "_for", suffix),      
  paste0("Home_", base_vars, "_against", suffix),  
  
  # --- AWAY TEAM HISTORIE ---
  paste0("Away_", base_vars, "_for", suffix),      
  paste0("Away_", base_vars, "_against", suffix),
  paste0("Home_", metrics_subset_1st, "_1st"),
  paste0("Away_", metrics_subset_1st, "_1st"),
  paste0("Home_", metrics_subset_2nd, "_for_roll_2nd"),
  paste0("Away_", metrics_subset_2nd, "_for_roll_2nd"),
  
  "home_formation_grouped", 
  "away_formation_grouped", 
  "Round_grouped",
  "derby_type" 
)

df_2nd_stats <- df_2nd_stats |>
  dplyr::select(Away, Home, home_first_goal_minute, goal_second_home, goal_second_away, home_scored, away_scored, all_of(candidate_vars), match_id) |>
  drop_na() |> 
  mutate(across(where(is.numeric) & !all_of(c("home_first_goal_minute", "home_scored", "match_id", "goal_second_home", "goal_second_away", "away_scored")), scale))


full_formula_string <- paste("Surv(goal_second_home, home_scored) ~", 
                             paste(candidate_vars, collapse = " + "))
as.formula(full_formula_string)
n<-names(df_2nd_stats)
names <- n[42:67]
full_formula_string <- paste(full_formula_string, " + ",  paste(names, collapse = " + "))
formula <- as.formula(full_formula_string)


summary(df_2nd_stats$Home_xG_1st)
sd(df_2nd_stats$Home_xG_1st)
summary(df_2nd_stats$Home_PSxG_1st)
sd(df_2nd_stats$Home_PSxG_1st)
summary(df_2nd_stats$Away_xG_1st)
sd(df_2nd_stats$Away_xG_1st)

# Make cox with all covariates

cox_home_full <- coxph(formula = formula, data = df_2nd_stats)

cox_home_null <- coxph(Surv(goal_second_home, home_scored) ~ 1, data = df_2nd_stats)
AIC_forwards <- stepAIC(cox_home_null, scope = list(lower = cox_home_null, upper = cox_home_full), direction = "forward")
summary(AIC_forwards)

cox_final <- coxph(formula = Surv(goal_second_home, home_scored) ~ Home_xG_1st + 
        Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
        Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 + home_formation_grouped, data = df_2nd_stats)
summary(cox_final)
cox.zph(cox_final)

# Check PH for games with goal scored by home team in the second half
df_2nd_stats_home_scored <- df_2nd_stats |> 
  filter(home_scored == 1)

cox_home_scored <- coxph(formula = Surv(goal_second_home, home_scored) ~ Home_xG_1st + 
                           Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
                           Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 + home_formation_grouped, data = df_2nd_stats_home_scored)
cox.zph(cox_home_scored)



# Change censored time to 65
df_2nd_stats_65 <- df_2nd_stats |> 
  mutate(goal_second_home= ifelse(home_scored == 0, 65, goal_second_home), 
         goal_second_away= ifelse(away_scored == 0, 65, goal_second_away))


df_2nd_stats_65 <- df_2nd_stats_65 |> 
  mutate(home_formation_grouped_3_5_2 = ifelse(home_formation_grouped == "3-5-2", 1, 0), 
         home_formation_grouped_4_4_2 = ifelse(home_formation_grouped == "4-4-2", 1, 0), 
         home_formation_grouped_4_2_3_1 = ifelse(home_formation_grouped == "4-2-3-1", 1, 0), 
         home_formation_grouped_4_3_3 = ifelse(home_formation_grouped == "4-3-3", 1, 0), 
         home_formation_grouped_4_5_1 = ifelse(home_formation_grouped == "4-5-1", 1, 0), 
         home_formation_grouped_5_3_2 = ifelse(home_formation_grouped == "5-3-2", 1, 0), 
         home_formation_grouped_5_4_1 = ifelse(home_formation_grouped == "5-4-1", 1, 0))
        
set.seed(123)
cure_home <- smcure(formula = Surv(goal_second_home, home_scored) ~ Home_xG_1st + 
       Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
       Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 +
       home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
       home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
       home_formation_grouped_5_4_1 , 
       cureform = ~ Home_xG_1st + 
         Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
         Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 + home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
         home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
         home_formation_grouped_5_4_1, 
       data = df_2nd_stats_65, 
       model = "ph", Var = T)

set.seed(123)
cure_home <- smcure(formula = Surv(goal_second_home, home_scored) ~ Home_xG_1st + 
                      Home_Number_passes_for_roll5  + Away_Shots_on_Target_against_roll5 + 
                      Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 +
                      home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
                      home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
                      home_formation_grouped_5_4_1 , 
                    cureform = ~ Home_xG_1st + 
                      Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
                      Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 + home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
                      home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
                      home_formation_grouped_5_4_1, 
                    data = df_2nd_stats_65, 
                    model = "ph", Var = T)

set.seed(123)
cure_home <- smcure(formula = Surv(goal_second_home, home_scored) ~ Home_xG_1st + 
                      Home_Number_passes_for_roll5  + Away_Shots_on_Target_against_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 +
                      home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
                      home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
                      home_formation_grouped_5_4_1 , 
                    cureform = ~ Home_xG_1st + 
                      Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
                      Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 + home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
                      home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
                      home_formation_grouped_5_4_1, 
                    data = df_2nd_stats_65, 
                    model = "ph", Var = T)

set.seed(123)
cure_home <- smcure(formula = Surv(goal_second_home, home_scored) ~ Home_xG_1st + 
                      Home_Number_passes_for_roll5  + Away_Shots_on_Target_against_roll5 +
                      home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
                      home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
                      home_formation_grouped_5_4_1 , 
                    cureform = ~ Home_xG_1st + 
                      Home_Number_passes_for_roll5 + Away_goals_for_roll5 + Away_Shots_on_Target_against_roll5 + 
                      Home_xG_for_roll5 + Away_Saves_for_roll5 + Away_Number_passes_for_roll5 + home_formation_grouped_3_5_2 + home_formation_grouped_4_4_2 + home_formation_grouped_4_2_3_1 +
                      home_formation_grouped_4_3_3 + home_formation_grouped_4_5_1 + home_formation_grouped_5_3_2 + 
                      home_formation_grouped_5_4_1, 
                    data = df_2nd_stats_65, 
                    model = "ph", Var = T)



#Away
full_formula_string <- paste("Surv(goal_second_away, away_scored) ~", 
                             paste(candidate_vars, collapse = " + "))
as.formula(full_formula_string)
n<-names(df_2nd_stats)
names <- n[44:69]
full_formula_string <- paste(full_formula_string, " + ",  paste(names, collapse = " + "))
formula <- as.formula(full_formula_string)

cox_away_full <- coxph(formula = formula, data = df_2nd_stats)

cox_away_null <- coxph(Surv(goal_second_away, away_scored) ~ 1, data = df_2nd_stats)
AIC_forwards <- stepAIC(cox_away_null, scope = list(lower = cox_away_null, upper = cox_away_full), direction = "forward")
summary(AIC_forwards)

summary(coxph(formula = Surv(goal_second_away, away_scored) ~ 
                Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + 
                Away_xG_for_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
                Home_Blocked_1st, data = df_2nd_stats))
cox_away <- coxph(formula = Surv(goal_second_away, away_scored) ~ 
        Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
        Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + 
        Away_xG_for_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
        Home_Blocked_1st, data = df_2nd_stats)
summary(cox_away)
cox.zph(cox_away)

# Selects games where goal was scored by the away team and check PH
df_2nd_stats_away_scored <- df_2nd_stats |> 
  filter(away_scored == 1)

cox_away <- coxph(formula = Surv(goal_second_away, away_scored) ~ 
                    Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                    Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + 
                    Away_xG_for_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
                    Home_Blocked_1st, data = df_2nd_stats_away_scored)
summary(cox_away)
cox.zph(cox_away)


df_2nd_stats_65 <- df_2nd_stats |> 
  mutate(goal_second_home= ifelse(home_scored == 0, 65, goal_second_home), 
         goal_second_away= ifelse(away_scored == 0, 65, goal_second_away))

set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
       Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
       Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + 
       Away_xG_for_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
       Home_Blocked_1st, 
       cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
         Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + 
         Away_xG_for_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
         Home_Blocked_1st, 
       data = df_2nd_stats_65, model = "ph", Var = T)

set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
                      Home_Blocked_1st, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 + Home_xG_against_roll5 + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
                      Home_Blocked_1st, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 +  
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + Home_shots_for_roll5 + Away_Saved_1st + 
                      Home_Blocked_1st, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 +  
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5 + Home_shots_for_roll5 +
                      Home_Blocked_1st, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 +  
                      Away_Blocked_1st + Home_shots_for_roll5 +
                      Home_Blocked_1st, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 +  
                      Away_Blocked_1st  + Home_shots_for_roll5, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ 
                      Away_Number_passes_for_roll5 + Home_shots_for_roll5 , 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ away_Number_passes_for_roll5, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 + 
                      Home_shots_for_roll5 + Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)
set.seed(123)
cure_away <- smcure(formula = Surv(goal_second_away, away_scored) ~ Away_Number_passes_for_roll5, 
                    cureform = ~ Away_Number_passes_for_roll5 + Home_xG_against_roll5 + Away_Off_Target_1st + 
                      Away_Blocked_1st + Away_Passing_accuracy_against_roll5  + Away_xG_for_roll5 +
                      Away_Saved_1st + Home_Blocked_1st, 
                    data = df_2nd_stats_65, model = "ph", Var = T)