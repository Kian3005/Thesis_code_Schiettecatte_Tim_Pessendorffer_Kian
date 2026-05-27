# EDA thesis
library(tidyverse)
library(survival)
library(survminer)
df <- read_csv("lessen/Thesis/Processed_full_stats/Full_df/full_df_2.csv")
# Count plus percentage team id first goal
n_tot <- nrow(df)
counts <- df |> 
  count(team_id_first_goal)
counts |> 
  mutate(percentage = n / n_tot * 100)
# Fill minute first goal of 0-0 games with 90 (censored games, could have scored later if game continued)
df <- df |> 
  mutate(first_goal_minute =  ifelse(is.na(team_id_first_goal), 90, first_goal_minute))
# Plot survival curve of first goal minute
df$first_goal_minute <- as.numeric(df$first_goal_minute)

surv_curves <- function(df, league) {
  fit1 <- survfit(Surv(first_goal_minute, goal_scored) ~ 1, stype = 1, data = df)
  return(list(plot(fit1, conf.int = FALSE, main = league), 
  summary(fit1),
  print(fit1, print.rmean = TRUE),
  quantile(fit1, probs= c(0.10, 0.25, 0.50, 0.75, 0.90))
        )  
  )
}
surv_curves(df, "All leagues")

# Bundesliga
df_bund <- df[df$comp_season %in% c("Processed_full_stats/Bundesliga/Bundesliga_2021.csv",
                                    "Processed_full_stats/Bundesliga/Bundesliga_2122.csv", 
                                    "Processed_full_stats/Bundesliga/Bundesliga_2223.csv",
                                    "Processed_full_stats/Bundesliga/Bundesliga_2324.csv",
                                    "Processed_full_stats/Bundesliga/Bundesliga_2425.csv"), ]
surv_curves(df_bund, "Bundesliga")

# Serie A
df_serieA <- df[df$comp_season %in% c("Processed_full_stats/Serie A/Serie A_2021.csv",
                                    "Processed_full_stats/Serie A/Serie A_2122.csv", 
                                    "Processed_full_stats/Serie A/Serie A_2223.csv",
                                    "Processed_full_stats/Serie A/Serie A_2324.csv",
                                    "Processed_full_stats/Serie A/Serie A_2425.csv"), ]
surv_curves(df_serieA, "Serie A")

# Premier League
df_prem <- df[df$comp_season %in% c("Processed_full_stats/PL/PL_2021.csv",
                                    "Processed_full_stats/PL/PL_2122.csv", 
                                    "Processed_full_stats/PL/PL_2223.csv",
                                    "Processed_full_stats/PL/PL_2324.csv",
                                    "Processed_full_stats/PL/PL_2425.csv"), ]
surv_curves(df_prem, "Premier League")

# La Liga
df_laliga <- df[df$comp_season %in% c("Processed_full_stats/LaLiga/LaLiga_2021.csv",
                                    "Processed_full_stats/LaLiga/LaLiga_2122.csv",
                                    "Processed_full_stats/LaLiga/LaLiga_2223.csv",
                                    "Processed_full_stats/LaLiga/LaLiga_2324.csv",
                                    "Processed_full_stats/LaLiga/LaLiga_2425.csv"), ]
surv_curves(df_laliga, "La Liga")

# Ligue 1
df_ligue1 <- df[df$comp_season %in% c("Processed_full_stats/Ligue 1/Ligue 1_2021.csv",
                                    "Processed_full_stats/Ligue 1/Ligue 1_2122.csv",
                                    "Processed_full_stats/Ligue 1/Ligue 1_2223.csv",
                                    "Processed_full_stats/Ligue 1/Ligue/_2324.csv",
                                    "Processed_full_stats/Ligue 1/Ligue 1_2425.csv"), ]
surv_curves(df_ligue1, "Ligue 1")

# JPL
df_jpl <- df[df$comp_season %in% c( "Processed_full_stats/JPL/JPL_2021_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2021_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2122_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2122_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2223_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2223_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2324_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2324_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2425_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2425_regular_season.csv"), ]
surv_curves(df_jpl, "JPL")

# Eredivisie
df_ere <- df[df$comp_season %in% c("Processed_full_stats/Eredivisie/Eredivisie_2021.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2122.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2223.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2324.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2425.csv"), ]
surv_curves(df_ere, "Eredivisie")

# Table with calculated quantiles per competition
quantiles_all <- surv_curves(df, "All leagues")[[4]]$quantile
quantiles_bund <- surv_curves(df_bund, "Bundesliga")[[4]]$quantile
quantiles_serieA <- surv_curves(df_serieA, "Serie A")[[4]]$quantile
quantiles_prem <- surv_curves(df_prem, "Premier League")[[4]]$quantile
quantiles_laliga <- surv_curves(df_laliga, "La Liga")[[4]]$quantile
quantiles_ligue1 <- surv_curves(df_ligue1, "Ligue 1")[[4]]$quantile
quantiles_jpl <- surv_curves(df_jpl, "JPL")[[4]]$quantile
quantiles_ere <- surv_curves(df_ere, "Eredivisie")[[4]]$quantile
quantiles_df <- data.frame(
  Competition = c("All leagues", "Bundesliga", "Serie A", "Premier League", "La Liga", "Ligue 1", "JPL", "Eredivisie"),
  Q10 = c(quantiles_all[1], quantiles_bund[1], quantiles_serieA[1], quantiles_prem[1], quantiles_laliga[1], quantiles_ligue1[1], quantiles_jpl[1], quantiles_ere[1]),
  Q25 = c(quantiles_all[2], quantiles_bund[2], quantiles_serieA[2], quantiles_prem[2], quantiles_laliga[2], quantiles_ligue1[2], quantiles_jpl[2], quantiles_ere[2]),
  Q50 = c(quantiles_all[3], quantiles_bund[3], quantiles_serieA[3], quantiles_prem[3], quantiles_laliga[3], quantiles_ligue1[3], quantiles_jpl[3], quantiles_ere[3]),
  Q75 = c(quantiles_all[4], quantiles_bund[4], quantiles_serieA[4], quantiles_prem[4], quantiles_laliga[4], quantiles_ligue1[4], quantiles_jpl[4], quantiles_ere[4]),
  Q90 = c(quantiles_all[5], quantiles_bund[5], quantiles_serieA[5], quantiles_prem[5], quantiles_laliga[5], quantiles_ligue1[5], quantiles_jpl[5], quantiles_ere[5])
)


# Checking values of home_formation
print(df |> 
  count(home_formation, sort = TRUE)|> 
  mutate(percentage = n / n_tot * 100), n = 22)
# Some columns with missing home_formations (7), manual fix
df |> 
  filter(is.na(home_formation)) |> 
  select(Home, Away, link, home_formation)
# Override NA values for home_formation by correct formations. 
df <- df |> 
  mutate(home_formation = ifelse(link == "https://fbref.com/en/matches/6a2c2932/Clermont-Foot-Nantes-September-17-2023-Ligue-1", "3-4-3", home_formation),
         home_formation = ifelse(link == "https://fbref.com/en/matches/8ab9a941/Montpellier-Brest-November-10-2024-Ligue-1", "4-2-3-1", home_formation),
         home_formation = ifelse(link == "https://fbref.com/en/matches/5a7ff7ce/Reims-Auxerre-March-9-2025-Ligue-1", "4-2-3-1", home_formation),
         home_formation = ifelse(link == "https://fbref.com/en/matches/7055dbbe/Le-Havre-Saint-Etienne-March-9-2025-Ligue-1", "4-2-3-1", home_formation),
         home_formation = ifelse(link == "https://fbref.com/en/matches/01d100ee/Nantes-Strasbourg-March-9-2025-Ligue-1", "5-3-2", home_formation),
         home_formation = ifelse(link == "https://fbref.com/en/matches/4c2f6a7b/Wolverhampton-Wanderers-West-Ham-United-April-5-2021-Premier-League", "4-2-3-1", home_formation),
         home_formation = ifelse(link == "https://fbref.com/en/matches/daebaaa1/Burnley-Newcastle-United-April-11-2021-Premier-League", "4-4-2", home_formation)
    )
# Checking weird formation 4-2-4-0. Is just a 4-2-4.
df |> 
  filter(home_formation == "4-2-4-0") |> 
  select(Home, Away, link, home_formation)

df <- df |> 
  mutate(home_formation = ifelse(link %in% c("https://fbref.com/en/matches/e6a099a7/El-Clasico-Real-Madrid-Barcelona-March-20-2022-La-Liga", 
                                             "https://fbref.com/en/matches/fede7f6e/Manchester-United-Everton-March-9-2024-Premier-League"), "4-2-4", home_formation))
# Checking values of away_formation
print(df |> 
  count(away_formation, sort = TRUE) |> 
  mutate(percentage = n / n_tot * 100), n =22)
# Also 7 missings -> manual fix
df |> 
  filter(is.na(away_formation)) |> 
  select(Home, Away, link, away_formation)

df <- df |>
  mutate(away_formation = ifelse(link == "https://fbref.com/en/matches/6a2c2932/Clermont-Foot-Nantes-September-17-2023-Ligue-1", "4-2-3-1", away_formation),
         away_formation = ifelse(link == "https://fbref.com/en/matches/8ab9a941/Montpellier-Brest-November-10-2024-Ligue-1", "4-3-3", away_formation),
         away_formation = ifelse(link == "https://fbref.com/en/matches/5a7ff7ce/Reims-Auxerre-March-9-2025-Ligue-1", "5-4-1", away_formation),
         away_formation = ifelse(link == "https://fbref.com/en/matches/7055dbbe/Le-Havre-Saint-Etienne-March-9-2025-Ligue-1", "4-2-3-1", away_formation),
         away_formation = ifelse(link == "https://fbref.com/en/matches/01d100ee/Nantes-Strasbourg-March-9-2025-Ligue-1", "3-4-3", away_formation),
         away_formation = ifelse(link == "https://fbref.com/en/matches/4c2f6a7b/Wolverhampton-Wanderers-West-Ham-United-April-5-2021-Premier-League", "4-2-3-1", away_formation),
         away_formation = ifelse(link == "https://fbref.com/en/matches/daebaaa1/Burnley-Newcastle-United-April-11-2021-Premier-League", "5-3-2", away_formation)
  )

# Checking weird formation 4-2-4-0. Is just a 4-2-4.
df |> 
  filter(away_formation == "4-2-4-0") |> 
  select(Home, Away, link, away_formation)

df <- df |> 
  mutate(away_formation = ifelse(link == "https://fbref.com/en/matches/615eff06/Manchester-Derby-Manchester-City-Manchester-United-March-3-2024-Premier-League", "4-2-4", away_formation))

# Check for other missing values
colSums(is.na(df))

# Fix game with missings on home_goals and away goals
# Distribution Home_Number_passes
histogram <- function(df, var, binwidth) {
  var_name <- deparse(substitute(var))
  
  return(ggplot(df, aes(x = var)) +
    geom_histogram(binwidth = binwidth, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = paste0("Distribution of ", var_name),
        x = var_name,
        y = "Frequency") +
    theme_minimal()
  )
}

histogram(df, df$Home_Number_passes, 50)
# Distribution Away_Number_passes
histogram(df, df$Away_Number_passes, 50)
# Distribution Home_possesion
df$Home_possesion <- as.numeric(str_remove(df$Home_possesion, "%"))
histogram(df, df$Home_possesion, 5)
# Distribution Away_possesion
df$Away_possesion <- as.numeric(str_remove(df$Away_possesion, "%"))
histogram(df, df$Away_possesion, 5)
# Distribution Home_shots
histogram(df, df$Home_shots, 3)
# Distribution Away_shots
histogram(df, df$Away_shots, 3)
# Distribution Home_Shots_on_Target
histogram(df, df$Home_Shots_on_Target, 1)
# Distribution Away_Shots_on_Target
histogram(df, df$Away_Shots_on_Target, 1)
# Distribution Home_Corners
histogram(df, df$Home_Corners, 1)
# Distribution Away_Corners
histogram(df, df$Away_Corners, 1)
# Distribution Home_Passing_accuracy
df$Home_Passing_accuracy <- as.numeric(str_remove(df$Home_Passing_accuracy, "%"))
histogram(df, df$Home_Passing_accuracy, 5)
# Distribution Away_Passing_accuracy
df$Away_Passing_accuracy <- as.numeric(str_remove(df$Away_Passing_accuracy, "%"))
histogram(df, df$Away_Passing_accuracy, 5)
# Extract Home_goals and Away_goals from score
df <- df |> 
  mutate(Home_goals = as.numeric(str_extract(Score, "^[0-9]+")),
         Away_goals = as.numeric(str_extract(Score, "[0-9]+$"))) |> 
  relocate(Home_goals, Away_goals, .after = Score)

# Check for other missing values
colSums(is.na(df))

# Fix game with missings on home_goals and away goals. Replace Home_goals by 0 and Away_goals by 1
df <- df |> 
  mutate(Home_goals = ifelse(link == "https://fbref.com/en/matches/2a5f6e4f/Angers-Montpellier-October-29-2023-Ligue-1", 0, Home_goals),
         Away_goals = ifelse(link == "https://fbref.com/en/matches/2a5f6e4f/Angers-Montpellier-October-29-2023-Ligue-1", 1, Away_goals))

# Distribution Home_goals
histogram(df, df$Home_goals, 1)
# Distribution Away_goals
histogram(df, df$Away_goals, 1)
# Distribution of xG
histogram(df, df$xG, 1)
# Distribution of xG.1
histogram(df, df$xG.1, 1)

# Save cleaned dataframe
write_csv(df, "lessen/Thesis/Processed_full_stats/Full_df/full_df_3.csv")


### Further EDA
df <- read_csv("lessen/Thesis/Processed_full_stats/Full_df/full_df_3.csv")

# Bundesliga
df_bund <- df[df$comp_season %in% c("Processed_full_stats/Bundesliga/Bundesliga_2021.csv",
                                    "Processed_full_stats/Bundesliga/Bundesliga_2122.csv", 
                                    "Processed_full_stats/Bundesliga/Bundesliga_2223.csv",
                                    "Processed_full_stats/Bundesliga/Bundesliga_2324.csv",
                                    "Processed_full_stats/Bundesliga/Bundesliga_2425.csv"), ]

# Serie A
df_serieA <- df[df$comp_season %in% c("Processed_full_stats/Serie A/Serie A_2021.csv",
                                      "Processed_full_stats/Serie A/Serie A_2122.csv", 
                                      "Processed_full_stats/Serie A/Serie A_2223.csv",
                                      "Processed_full_stats/Serie A/Serie A_2324.csv",
                                      "Processed_full_stats/Serie A/Serie A_2425.csv"), ]

# Premier League
df_prem <- df[df$comp_season %in% c("Processed_full_stats/PL/PL_2021.csv",
                                    "Processed_full_stats/PL/PL_2122.csv", 
                                    "Processed_full_stats/PL/PL_2223.csv",
                                    "Processed_full_stats/PL/PL_2324.csv",
                                    "Processed_full_stats/PL/PL_2425.csv"), ]

# La Liga
df_laliga <- df[df$comp_season %in% c("Processed_full_stats/LaLiga/LaLiga_2021.csv",
                                      "Processed_full_stats/LaLiga/LaLiga_2122.csv",
                                      "Processed_full_stats/LaLiga/LaLiga_2223.csv",
                                      "Processed_full_stats/LaLiga/LaLiga_2324.csv",
                                      "Processed_full_stats/LaLiga/LaLiga_2425.csv"), ]

# Ligue 1
df_ligue1 <- df[df$comp_season %in% c("Processed_full_stats/Ligue 1/Ligue 1_2021.csv",
                                      "Processed_full_stats/Ligue 1/Ligue 1_2122.csv",
                                      "Processed_full_stats/Ligue 1/Ligue 1_2223.csv",
                                      "Processed_full_stats/Ligue 1/Ligue/_2324.csv",
                                      "Processed_full_stats/Ligue 1/Ligue 1_2425.csv"), ]

# JPL
df_jpl <- df[df$comp_season %in% c( "Processed_full_stats/JPL/JPL_2021_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2021_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2122_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2122_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2223_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2223_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2324_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2324_regular_season.csv",
                                    "Processed_full_stats/JPL/JPL_2425_playoffs.csv",
                                    "Processed_full_stats/JPL/JPL_2425_regular_season.csv"), ]

# Eredivisie
df_ere <- df[df$comp_season %in% c("Processed_full_stats/Eredivisie/Eredivisie_2021.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2122.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2223.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2324.csv",
                                   "Processed_full_stats/Eredivisie/Eredivisie_2425.csv"), ]

## Distribution home_number_passes per league
# Bundesliga
histogram(df_bund, df_bund$Home_Number_passes, 50)
# Serie A, looks a bit less than Bundesliga
histogram(df_serieA, df_serieA$Home_Number_passes, 50)
# Premier League
histogram(df_prem, df_prem$Home_Number_passes, 50)
# La Liga
histogram(df_laliga, df_laliga$Home_Number_passes, 50)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_Number_passes, 50)
# JPL
histogram(df_jpl, df_jpl$Home_Number_passes, 50)
# Eredivisie
histogram(df_ere, df_ere$Home_Number_passes, 50)

## Distribution away_number_passes per league
# Bundesliga
histogram(df_bund, df_bund$Away_Number_passes, 50)
# Serie A
histogram(df_serieA, df_serieA$Away_Number_passes, 50)
# Premier League
histogram(df_prem, df_prem$Away_Number_passes, 50)
# La Liga
histogram(df_laliga, df_laliga$Away_Number_passes, 50)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_Number_passes, 50)
# JPL
histogram(df_jpl, df_jpl$Away_Number_passes, 50)
# Eredivisie
histogram(df_ere, df_ere$Away_Number_passes, 50)

## Distribution home_possesion per league
# Bundesliga
histogram(df_bund, df_bund$Home_possesion, 5)
# Serie A
histogram(df_serieA, df_serieA$Home_possesion, 5)
# Premier League
histogram(df_prem, df_prem$Home_possesion, 5)
# La Liga
histogram(df_laliga, df_laliga$Home_possesion, 5)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_possesion, 5)
# JPL
histogram(df_jpl, df_jpl$Home_possesion, 5)
# Eredivisie
histogram(df_ere, df_ere$Home_possesion, 5)

## Distribution away_possesion per league
# Bundesliga
histogram(df_bund, df_bund$Away_possesion, 5)
# Serie A
histogram(df_serieA, df_serieA$Away_possesion, 5)
# Premier League
histogram(df_prem, df_prem$Away_possesion, 5)
# La Liga
histogram(df_laliga, df_laliga$Away_possesion, 5)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_possesion, 5)
# JPL
histogram(df_jpl, df_jpl$Away_possesion, 5)
# Eredivisie
histogram(df_ere, df_ere$Away_possesion, 5)

## Distribution home_shots per league
# Bundesliga
histogram(df_bund, df_bund$Home_shots, 3)
# Serie A
histogram(df_serieA, df_serieA$Home_shots, 3)
# Premier League
histogram(df_prem, df_prem$Home_shots, 3)
# La Liga
histogram(df_laliga, df_laliga$Home_shots, 3)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_shots, 3)
# JPL
histogram(df_jpl, df_jpl$Home_shots, 3)
# Eredivisie
histogram(df_ere, df_ere$Home_shots, 3)

## Distribution away_shots per league
# Bundesliga
histogram(df_bund, df_bund$Away_shots, 3)
# Serie A
histogram(df_serieA, df_serieA$Away_shots, 3)
# Premier League
histogram(df_prem, df_prem$Away_shots, 3)
# La Liga
histogram(df_laliga, df_laliga$Away_shots, 3)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_shots, 3)
# JPL
histogram(df_jpl, df_jpl$Away_shots, 3)
# Eredivisie
histogram(df_ere, df_ere$Away_shots, 3)

## Distribution home_shots_on_target per league
# Bundesliga
histogram(df_bund, df_bund$Home_Shots_on_Target, 1)
# Serie A
histogram(df_serieA, df_serieA$Home_Shots_on_Target, 1)
# Premier League
histogram(df_prem, df_prem$Home_Shots_on_Target, 1)
# La Liga
histogram(df_laliga, df_laliga$Home_Shots_on_Target, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_Shots_on_Target, 1)
# JPL
histogram(df_jpl, df_jpl$Home_Shots_on_Target, 1)
# Eredivisie
histogram(df_ere, df_ere$Home_Shots_on_Target, 1)

## Distribution away_shots_on_target per league
# Bundesliga
histogram(df_bund, df_bund$Away_Shots_on_Target, 1)
# Serie A
histogram(df_serieA, df_serieA$Away_Shots_on_Target, 1)
# Premier League
histogram(df_prem, df_prem$Away_Shots_on_Target, 1)
# La Liga
histogram(df_laliga, df_laliga$Away_Shots_on_Target, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_Shots_on_Target, 1)
# JPL
histogram(df_jpl, df_jpl$Away_Shots_on_Target, 1)
# Eredivisie
histogram(df_ere, df_ere$Away_Shots_on_Target, 1)

## Distribution home_corners per league
# Bundesliga
histogram(df_bund, df_bund$Home_Corners, 1)
# Serie A
histogram(df_serieA, df_serieA$Home_Corners, 1)
# Premier League
histogram(df_prem, df_prem$Home_Corners, 1)
# La Liga
histogram(df_laliga, df_laliga$Home_Corners, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_Corners, 1)
# JPL
histogram(df_jpl, df_jpl$Home_Corners, 1)
# Eredivisie
histogram(df_ere, df_ere$Home_Corners, 1)

## Distribution away_corners per league
# Bundesliga
histogram(df_bund, df_bund$Away_Corners, 1)
# Serie A
histogram(df_serieA, df_serieA$Away_Corners, 1)
# Premier League
histogram(df_prem, df_prem$Away_Corners, 1)
# La Liga
histogram(df_laliga, df_laliga$Away_Corners, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_Corners, 1)
# JPL
histogram(df_jpl, df_jpl$Away_Corners, 1)
# Eredivisie
histogram(df_ere, df_ere$Away_Corners, 1)

## Distribution home_passing_accuracy per league
# Bundesliga
histogram(df_bund, df_bund$Home_Passing_accuracy, 3)
# Serie A
histogram(df_serieA, df_serieA$Home_Passing_accuracy, 3)
# Premier League
histogram(df_prem, df_prem$Home_Passing_accuracy, 3)
# La Liga
histogram(df_laliga, df_laliga$Home_Passing_accuracy, 3)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_Passing_accuracy, 3)
# JPL
histogram(df_jpl, df_jpl$Home_Passing_accuracy, 3)
# Eredivisie
histogram(df_ere, df_ere$Home_Passing_accuracy, 3)

## Distribution away_passing_accuracy per league
# Bundesliga
histogram(df_bund, df_bund$Away_Passing_accuracy, 3)
# Serie A
histogram(df_serieA, df_serieA$Away_Passing_accuracy, 3)
# Premier League
histogram(df_prem, df_prem$Away_Passing_accuracy, 3)
# La Liga
histogram(df_laliga, df_laliga$Away_Passing_accuracy, 3)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_Passing_accuracy, 3)
# JPL
histogram(df_jpl, df_jpl$Away_Passing_accuracy, 3)
# Eredivisie
histogram(df_ere, df_ere$Away_Passing_accuracy, 3)

#### Cleaning of Home_goals and Away_goals because of penalties
df <- df |> 
  mutate(Home_goals = ifelse(link == "https://fbref.com/en/matches/7c44e3ed/Sparta-Rotterdam-Utrecht-June-4-2023-Eredivisie", 0, Home_goals),
         Away_goals = ifelse(link == "https://fbref.com/en/matches/7c44e3ed/Sparta-Rotterdam-Utrecht-June-4-2023-Eredivisie", 1, Away_goals))
write_csv(df, "lessen/Thesis/Processed_full_stats/Full_df/full_df_4.csv")

df <- read_csv("lessen/Thesis/Processed_full_stats/Full_df/full_df_4.csv")
## Distribution home_goals per league
# All leagues
histogram(df, df$Home_goals, 1)
# Bundesliga
histogram(df_bund, df_bund$Home_goals, 1)
# Serie A
histogram(df_serieA, df_serieA$Home_goals, 1)
# Premier League
histogram(df_prem, df_prem$Home_goals, 1)
# La Liga
histogram(df_laliga, df_laliga$Home_goals, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$Home_goals, 1)
# JPL
histogram(df_jpl, df_jpl$Home_goals, 1)
# Eredivisie
histogram(df_ere, df_ere$Home_goals, 1)

## Distribution away_goals per league
# All leagues
histogram(df, df$Away_goals, 1)
# Bundesliga
histogram(df_bund, df_bund$Away_goals, 1)
# Serie A
histogram(df_serieA, df_serieA$Away_goals, 1)
# Premier League
histogram(df_prem, df_prem$Away_goals, 1)
# La Liga
histogram(df_laliga, df_laliga$Away_goals, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$Away_goals, 1)
# JPL
histogram(df_jpl, df_jpl$Away_goals, 1)
# Eredivisie
histogram(df_ere, df_ere$Away_goals, 1)

## Distribution xG per league
# All leagues
histogram(df, df$xG, 1)
# Bundesliga
histogram(df_bund, df_bund$xG, 1)
# Serie A
histogram(df_serieA, df_serieA$xG, 1)
# Premier League
histogram(df_prem, df_prem$xG, 1)
# La Liga
histogram(df_laliga, df_laliga$xG, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$xG, 1)
# JPL
histogram(df_jpl, df_jpl$xG, 1)
# Eredivisie
histogram(df_ere, df_ere$xG, 1)

## Distribution xG.1 per league
# All leagues
histogram(df, df$xG.1, 1)
# Bundesliga
histogram(df_bund, df_bund$xG.1, 1)
# Serie A
histogram(df_serieA, df_serieA$xG.1, 1)
# Premier League
histogram(df_prem, df_prem$xG.1, 1)
# La Liga
histogram(df_laliga, df_laliga$xG.1, 1)
# Ligue 1
histogram(df_ligue1, df_ligue1$xG.1, 1)
# JPL
histogram(df_jpl, df_jpl$xG.1, 1)
# Eredivisie
histogram(df_ere, df_ere$xG.1, 1)


### Scatterplots target variabele (first_goal_minute) against covariates
scatter <- function(df, x_var, y_var) {
  x_var_name <- deparse(substitute(x_var))
  y_var_name <- deparse(substitute(y_var))
  
  return(ggplot(df, aes(x = x_var, y = y_var)) +
           geom_point(alpha = 0.5) +
           labs(title = paste0("Scatterplot of ", y_var_name, " vs ", x_var_name),
                x = x_var_name,
                y = y_var_name)
  )
}

## Scatterplot number_goals vs first_goal_minute
# First make a new variable number_goals
df <- df |>
  mutate(number_goals = Home_goals + Away_goals) |> 
  relocate(number_goals, .after = Away_goals)

# All leagues
scatter(df, df$number_goals, df$first_goal_minute)
# Bundesliga
scatter(df_bund, df_bund$number_goals, df_bund$first_goal_minute)
# Serie A
scatter(df_serieA, df_serieA$number_goals, df_serieA$first_goal_minute)
# Premier League
scatter(df_prem, df_prem$number_goals, df_prem$first_goal_minute)
# La Liga
scatter(df_laliga, df_laliga$number_goals, df_laliga$first_goal_minute)
# Ligue 1
scatter(df_ligue1, df_ligue1$number_goals, df_ligue1$first_goal_minute)
# JPL
scatter(df_jpl, df_jpl$number_goals, df_jpl$first_goal_minute)
# Eredivisie
scatter(df_ere, df_ere$number_goals, df_ere$first_goal_minute)

# correlation
cor(df$number_goals, df$first_goal_minute) # -0.59
cor(df_bund$number_goals, df_bund$first_goal_minute) # -0.6
cor(df_serieA$number_goals, df_serieA$first_goal_minute) # -0.59
cor(df_prem$number_goals, df_prem$first_goal_minute) #  -0.57
cor(df_laliga$number_goals, df_laliga$first_goal_minute) # -0.60
cor(df_ligue1$number_goals, df_ligue1$first_goal_minute) # -0.59
cor(df_jpl$number_goals, df_jpl$first_goal_minute) # -0.57
cor(df_ere$number_goals, df_ere$first_goal_minute) # -0.58

## Scatterplot xG vs first_goal_minute
# Make new variabele total_xg
df <- df |>
  mutate(total_xg = xG + xG.1) |> 
  relocate(total_xg, .after = xG.1)
# All leagues
scatter(df, df$total_xg, df$first_goal_minute)
# Bundesliga
scatter(df_bund, df_bund$total_xg, df_bund$first_goal_minute)
# Serie A
scatter(df_serieA, df_serieA$total_xg, df_serieA$first_goal_minute)
# Premier League
scatter(df_prem, df_prem$total_xg, df_prem$first_goal_minute)
# La Liga
scatter(df_laliga, df_laliga$total_xg, df_laliga$first_goal_minute)
# Ligue 1
scatter(df_ligue1, df_ligue1$total_xg, df_ligue1$first_goal_minute)
# JPL
scatter(df_jpl, df_jpl$total_xg, df_jpl$first_goal_minute)
# Eredivisie
scatter(df_ere, df_ere$total_xg, df_ere$first_goal_minute)

# correlations
cor(df$total_xg, df$first_goal_minute) # -0.32
cor(df_bund$total_xg, df_bund$first_goal_minute) # - 0.32 
cor(df_serieA$total_xg, df_serieA$first_goal_minute) # -0.32
cor(df_prem$total_xg, df_prem$first_goal_minute) # -0.30
cor(df_laliga$total_xg, df_laliga$first_goal_minute) # -0.35
cor(df_ligue1$total_xg, df_ligue1$first_goal_minute) # -0.32
cor(df_jpl$total_xg, df_jpl$first_goal_minute) # -0.32
cor(df_ere$total_xg, df_ere$first_goal_minute) # -0.28
# Less strong correlation compared to the actual number of goals. -> possible flaw in the calculation of xG?
# A difference of 0,27!
# Also a difference compared to correlation with PSxG?

## Extract PSxG_home, PSxG_away and total_PSxG
df <- df |>
  mutate(PSxG_home = Home_PSxG_15 + Home_PSxG_30 + Home_PSxG_45 + Home_PSxG_60 + Home_PSxG_75 + Home_PSxG_90,
         PSxG_away = Away_PSxG_15 + Away_PSxG_30 + Away_PSxG_45 + Away_PSxG_60 + Away_PSxG_75 + Away_PSxG_90,
         total_PSxG = PSxG_home + PSxG_away) |> 
  relocate(PSxG_home, .after = xG) |>
  relocate(PSxG_away, .after = xG.1) |>
  relocate(total_PSxG, .after = total_xg)

# All leagues
scatter(df, df$total_PSxG, df$first_goal_minute)
# Bundesliga
scatter(df_bund, df_bund$total_PSxG, df_bund$first_goal_minute)
# Serie A
scatter(df_serieA, df_serieA$total_PSxG, df_serieA$first_goal_minute)
# Premier League
scatter(df_prem, df_prem$total_PSxG, df_prem$first_goal_minute)
# La Liga
scatter(df_laliga, df_laliga$total_PSxG, df_laliga$first_goal_minute)
# Ligue 1
scatter(df_ligue1, df_ligue1$total_PSxG, df_ligue1$first_goal_minute)
# JPL
scatter(df_jpl, df_jpl$total_PSxG, df_jpl$first_goal_minute)
# Eredivisie
scatter(df_ere, df_ere$total_PSxG, df_ere$first_goal_minute)

# correlations
cor(df$total_PSxG, df$first_goal_minute) # -0.42
cor(df_bund$total_PSxG, df_bund$first_goal_minute) # -0.43
cor(df_serieA$total_PSxG, df_serieA$first_goal_minute) # -0.43
cor(df_prem$total_PSxG, df_prem$first_goal_minute) # -0.40
cor(df_laliga$total_PSxG, df_laliga$first_goal_minute) # -0.40
cor(df_ligue1$total_PSxG, df_ligue1$first_goal_minute) # -0.41
cor(df_jpl$total_PSxG, df_jpl$first_goal_minute) # -0.42
cor(df_ere$total_PSxG, df_ere$first_goal_minute) # -0.38
# Bit stronger corr than xG and closer to the corr with the actual number of goals.
# Still a difference of 0.17 but vetter than xG (0.10)

## Scatterplots total_xG VS number_goals
scatter(df, df$total_xg, df$number_goals)
# Bundesliga
scatter(df_bund, df_bund$total_xg, df_bund$number_goals)
# Serie A
scatter(df_serieA, df_serieA$total_xg, df_serieA$number_goals)
# Premier League
scatter(df_prem, df_prem$total_xg, df_prem$number_goals)
# La Liga
scatter(df_laliga, df_laliga$total_xg, df_laliga$number_goals)
# Ligue 1
scatter(df_ligue1, df_ligue1$total_xg, df_ligue1$number_goals)
# JPL
scatter(df_jpl, df_jpl$total_xg, df_jpl$number_goals)
# Eredivisie
scatter(df_ere, df_ere$total_xg, df_ere$number_goals)

# correlations
cor(df$total_xg, df$number_goals) # 0.56
cor(df_bund$total_xg, df_bund$number_goals) # 0.54
cor(df_serieA$total_xg, df_serieA$number_goals) # 0.55
cor(df_prem$total_xg, df_prem$number_goals) # 0.54
cor(df_laliga$total_xg, df_laliga$number_goals) # 0.59
cor(df_ligue1$total_xg, df_ligue1$number_goals) # 0.55
cor(df_jpl$total_xg, df_jpl$number_goals) # 0.52
cor(df_ere$total_xg, df_ere$number_goals) # 0.55
# Not that strong correlation. xG is not a perfect predictor of actual goals.

## Scatterplots total_PSxG VS number_goals
scatter(df, df$total_PSxG, df$number_goals)
# Bundesliga
scatter(df_bund, df_bund$total_PSxG, df_bund$number_goals)
# Serie A
scatter(df_serieA, df_serieA$total_PSxG, df_serieA$number_goals)
# Premier League
scatter(df_prem, df_prem$total_PSxG, df_prem$number_goals)
# La Liga
scatter(df_laliga, df_laliga$total_PSxG, df_laliga$number_goals)
# Ligue 1
scatter(df_ligue1, df_ligue1$total_PSxG, df_ligue1$number_goals)
# JPL
scatter(df_jpl, df_jpl$total_PSxG, df_jpl$number_goals)
# Eredivisie
scatter(df_ere, df_ere$total_PSxG, df_ere$number_goals)

# correlations
cor(df$total_PSxG, df$number_goals) # 0.72
cor(df_bund$total_PSxG, df_bund$number_goals) # 0.73
cor(df_serieA$total_PSxG, df_serieA$number_goals) # 0.72
cor(df_prem$total_PSxG, df_prem$number_goals) # 0.72
cor(df_laliga$total_PSxG, df_laliga$number_goals) # 0.75
cor(df_ligue1$total_PSxG, df_ligue1$number_goals) # 0.72
cor(df_jpl$total_PSxG, df_jpl$number_goals) # 0.72
cor(df_ere$total_PSxG, df_ere$number_goals) # 0.68 (lower bcs of impact outliers to the right)
# Strong correlation, better than xG. (+0.16)

## Corners versus first_goal_minute
# Make new variabele total_corners
df <- df |>
  mutate(total_corners = Home_Corners + Away_Corners) |> 
  relocate(total_corners, .after = Away_Corners)

# All leagues
scatter(df, df$total_corners, df$first_goal_minute)
# Bundesliga
scatter(df_bund, df_bund$total_corners, df_bund$first_goal_minute)
# Serie A
scatter(df_serieA, df_serieA$total_corners, df_serieA$first_goal_minute)
# Premier League
scatter(df_prem, df_prem$total_corners, df_prem$first_goal_minute)
# La Liga
scatter(df_laliga, df_laliga$total_corners, df_laliga$first_goal_minute)
# Ligue 1
scatter(df_ligue1, df_ligue1$total_corners, df_ligue1$first_goal_minute)
# JPL
scatter(df_jpl, df_jpl$total_corners, df_jpl$first_goal_minute)
# Eredivisie
scatter(df_ere, df_ere$total_corners, df_ere$first_goal_minute)

# correlations
cor(df$total_corners, df$first_goal_minute) # 0
cor(df_bund$total_corners, df_bund$first_goal_minute) # -0.04
cor(df_serieA$total_corners, df_serieA$first_goal_minute) # 0.02
cor(df_prem$total_corners, df_prem$first_goal_minute) # 0.03
cor(df_laliga$total_corners, df_laliga$first_goal_minute) # -0.02
cor(df_ligue1$total_corners, df_ligue1$first_goal_minute) # 0.03
cor(df_jpl$total_corners, df_jpl$first_goal_minute) # 0.02
cor(df_ere$total_corners, df_ere$first_goal_minute) # 0.03
# No correlation on total number of corners. 
# Maybe a correlation on difference in corners?
# Make new variabele diff_corners
df <- df |>
  mutate(diff_corners = Home_Corners - Away_Corners) |> 
  relocate(diff_corners, .after = total_corners)
# All leagues
scatter(df, df$diff_corners, df$first_goal_minute)
# Bundesliga
scatter(df_bund, df_bund$diff_corners, df_bund$first_goal_minute)
# Serie A
scatter(df_serieA, df_serieA$diff_corners, df_serieA$first_goal_minute)
# Premier League
scatter(df_prem, df_prem$diff_corners, df_prem$first_goal_minute)
# La Liga
scatter(df_laliga, df_laliga$diff_corners, df_laliga$first_goal_minute)
# Ligue 1
scatter(df_ligue1, df_ligue1$diff_corners, df_ligue1$first_goal_minute)
# JPL
scatter(df_jpl, df_jpl$diff_corners, df_jpl$first_goal_minute)
# Eredivisie
scatter(df_ere, df_ere$diff_corners, df_ere$first_goal_minute)
# correlations
cor(df$diff_corners, df$first_goal_minute) # 0.02
cor(df_bund$diff_corners, df_bund$first_goal_minute) # 0.03
cor(df_serieA$diff_corners, df_serieA$first_goal_minute) # 0.04
cor(df_prem$diff_corners, df_prem$first_goal_minute) # 0.02
cor(df_laliga$diff_corners, df_laliga$first_goal_minute) # 0.02
cor(df_ligue1$diff_corners, df_ligue1$first_goal_minute) # 0.03
cor(df_jpl$diff_corners, df_jpl$first_goal_minute) # 0.01
cor(df_ere$diff_corners, df_ere$first_goal_minute) # 0
# No correlation again. 
## What if we account for the team that scored the first goal (team_id_first_goal)?
df_h <- df |> 
  filter(team_id_first_goal == "Home")
df_a <- df |> 
  filter(team_id_first_goal == "Away")
# Home team scored first
cor(df_h$diff_corners, df_h$first_goal_minute) # 0.09
cor(df_a$diff_corners, df_a$first_goal_minute) # -0.08
# So small correlation
# correlation with home_corners and away_corners per team id.
cor(df_h$Home_Corners, df_h$first_goal_minute) # 0.07
cor(df_h$Away_Corners, df_h$first_goal_minute) # -0.08
cor(df_a$Home_Corners, df_a$first_goal_minute) # -0.06
cor(df_a$Away_Corners, df_a$first_goal_minute) # 0.07
# So small correlations again, a first goal of the home_team leads to a positive corr with home_corners and a negative correlation with away_corners. 
# the opposite holds for the away_team scoring first.

### Corners versus number_goals
# All leagues
scatter(df, df$total_corners, df$number_goals)
# Bundesliga
scatter(df_bund, df_bund$total_corners, df_bund$number_goals)
# Serie A
scatter(df_serieA, df_serieA$total_corners, df_serieA$number_goals)
# Premier League
scatter(df_prem, df_prem$total_corners, df_prem$number_goals)
# La Liga
scatter(df_laliga, df_laliga$total_corners, df_laliga$number_goals)
# Ligue 1
scatter(df_ligue1, df_ligue1$total_corners, df_ligue1$number_goals)
# JPL
scatter(df_jpl, df_jpl$total_corners, df_jpl$number_goals)
# Eredivisie
scatter(df_ere, df_ere$total_corners, df_ere$number_goals)

# correlations
cor(df$total_corners, df$number_goals) 
cor(df_bund$total_corners, df_bund$number_goals) 
cor(df_serieA$total_corners, df_serieA$number_goals) 
cor(df_prem$total_corners, df_prem$number_goals) 
cor(df_laliga$total_corners, df_laliga$number_goals) 
cor(df_ligue1$total_corners, df_ligue1$number_goals) 
cor(df_jpl$total_corners, df_jpl$number_goals) 
cor(df_ere$total_corners, df_ere$number_goals) 
# No correlations

## for home and away seperately?
# Home team
cor(df$Home_Corners, df$Home_goals) 
cor(df_bund$Home_Corners, df_bund$Home_goals) 
cor(df_serieA$Home_Corners, df_serieA$Home_goals) 
cor(df_prem$Home_Corners, df_prem$Home_goals)  # 0.08
cor(df_laliga$Home_Corners, df_laliga$Home_goals)
cor(df_ligue1$Home_Corners, df_ligue1$Home_goals) 
cor(df_jpl$Home_Corners, df_jpl$Home_goals) 
cor(df_ere$Home_Corners, df_ere$Home_goals) # 0.17
# Not strong correlations
# Away team
cor(df$Away_Corners, df$Away_goals) 
cor(df_bund$Away_Corners, df_bund$Away_goals) 
cor(df_serieA$Away_Corners, df_serieA$Away_goals) 
cor(df_prem$Away_Corners, df_prem$Away_goals) 
cor(df_laliga$Away_Corners, df_laliga$Away_goals) 
cor(df_ligue1$Away_Corners, df_ligue1$Away_goals) 
cor(df_jpl$Away_Corners, df_jpl$Away_goals) 
cor(df_ere$Away_Corners, df_ere$Away_goals) # 0.18
# Not strong correlations again

## Descriptive statistics of covariates
# Function to calculate descriptive statistics
descriptive_stats <- function(df, var_name) {
  df %>%
    summarise( # .data is used to call it as a string
      Mean = mean(.data[[var_name]], na.rm = TRUE),
      Median = median(.data[[var_name]], na.rm = TRUE),
      SD = sd(.data[[var_name]], na.rm = TRUE),
      Min = min(.data[[var_name]], na.rm = TRUE),
      Max = max(.data[[var_name]], na.rm = TRUE),
      N = sum(!is.na(.data[[var_name]])) 
    ) %>%
    mutate(Variable = var_name, .before = 1)
}

covariates <- c("Home_Number_passes", "Away_Number_passes", "Home_possesion", "Away_possesion",  
                 "Home_shots", "Away_shots", "Home_Shots_on_Target", "Away_Shots_on_Target",  
                 "Home_Corners", "Away_Corners", "Home_Passing_accuracy", "Away_Passing_accuracy",  
                 "Home_goals", "Away_goals", "xG", "xG.1", "PSxG_home", "PSxG_away", "total_PSxG",
                 "number_goals", "total_xg", "total_corners", "diff_corners", "Home_xG_15", "Home_xG_30",
                 "Home_xG_45", "Home_xG_60", "Home_xG_75", "Home_xG_90", "Away_xG_15", "Away_xG_30",
                 "Away_xG_45", "Away_xG_60", "Away_xG_75", "Away_xG_90")
descriptive_all <- map_dfr(covariates, ~ descriptive_stats(df, .x))


## Fill up missings attendance with 0 (no public allowed)
df <- df |>
  mutate(Attendance = ifelse(is.na(Attendance ), 0, Attendance ))
# Check reason missings save_accuracy
df |> 
  filter(is.na(Home_Saves_accuracy) | is.na(Away_Saves_accuracy)) |> 
  select(Home_Saves_accuracy, Away_Saves_accuracy, Home_Saves, Away_Saves) # When missing it's because there were no saves to be made

# Check missing saves 
df |> 
  filter(is.na(Home_Saves) | is.na(Away_Saves)) |> 
  select(link, Home_Saves, Away_Saves) # Missing because no shots on target

# Fill up missings with 0
df <- df |>
  mutate(Home_Saves = ifelse(is.na(Home_Saves), 0, Home_Saves),
         Away_Saves = ifelse(is.na(Away_Saves), 0, Away_Saves))

# Save final dataset 
write_csv(df, "lessen/Thesis/Processed_full_stats/Full_df/full_df_5.csv")