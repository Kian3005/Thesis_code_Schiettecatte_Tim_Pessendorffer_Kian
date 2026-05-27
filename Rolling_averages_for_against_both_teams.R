library(tidyverse)
library(slider)
# For the models we will use rolling averages of 4 last games instead of the last game only for the covariates. 
# This will smooth out the data a bit and make it less noisy. It will also enable us to take form more into account. 
# Rolling averages will be based on the date of the games and the Home/Away columns since games can be suspended and placed on a later week.

# Slider package: slide based on a date vector that will serve as an index. 
# Use of slide_index_mean() or slide_mean() where, .x is used for data and .i is used for the index (date here).
# Window is calculated based on the index.
# If we can arrange by data we don't need the index. 

# Calculation will be done per season and for each team.
# Extract season from comp_season, digits before .csv at the end of the string.
# Except for JPL_playoffs.csv and JPL_regular_season.csv, see fix missing values
df <- read_csv("lessen/Thesis/Processed_full_stats/Full_df/full_df_with_other_first_goal_time.csv")
df <- df |> 
  mutate(season = str_extract(comp_season, "\\d{4}(?=\\.csv)")) |> 
  relocate(season, .after = comp_season)
# Fix missing season values
df <- df |> 
  mutate(season = ifelse(is.na(season), str_extract(comp_season, "\\d{4}(?=_regular_season\\.csv|_playoffs\\.csv)"), season)) 

# Create a function to calculate rolling averages for a given team and season
# Try out
test_df <- df |>
  filter(season == "2021", Home == "Anderlecht" | Away == "Anderlecht") |>
  arrange(Date)

# Columns for rolling averages
numerics <- test_df |>
  select(where(is.numeric)) |> 
  select(-c(Wk, Attendance, diff_corners, number_goals, total_xg, total_PSxG, total_corners)) |> 
  colnames()

xg_rolling <- function(team){
  df <- df |> 
    filter(Home == team | Away == team)
  test_df_xg <- df |> 
    mutate(xg_team = ifelse(Home == team, xG, xG.1)) |> 
    relocate(xg_team, .before = xG.1)

  test_df_xg_rolling <- test_df_xg |>
    mutate(xg_rolling = slide_mean(xg_team, before = 4, after = -1)) |> 
    relocate(xg_rolling, .after = xg_team)

  test_df_xg_rolling_team <- test_df_xg_rolling |> 
    mutate(Home_xg_rolling = ifelse(Home == team, xg_rolling, NA), 
          Away_xg_rolling = ifelse(Home != team, xg_rolling, NA)) |> 
    relocate(c(Home_xg_rolling, Away_xg_rolling), .after = xg_rolling) |> 
    select(-c(xg_team, xg_rolling))
  return(test_df_xg_rolling_team)
}

# Test the function
df <- test_df
test_function_df <- xg_rolling("Anderlecht")

# Make df right for 2021
df <- read_csv("lessen/Thesis/Processed_full_stats/Full_df/full_df_with_other_first_goal_time.csv")
df <- df |> 
  mutate(season = str_extract(comp_season, "\\d{4}(?=\\.csv)")) |> 
  relocate(season, .after = comp_season)
# Fix missing season values
df <- df |> 
  mutate(season = ifelse(is.na(season), str_extract(comp_season, "\\d{4}(?=_regular_season\\.csv|_playoffs\\.csv)"), season))

df <- df |>
  filter(season == "2021") |>
  arrange(Date)

# Teams for in the map function
teams <- unique(df$Home)
# Add a match id to join the dataframes
df <- df |>
  mutate(match_id = row_number())

# Using map to do the function for all teams. The function now gives 2 rows for one game. 
# 
df_xgs <- map_dfr(teams, ~xg_rolling(.x))
df_xgs <- df_xgs |>
  group_by(match_id) |>
  mutate(
    Home_xg_rolling = sum(Home_xg_rolling, na.rm = TRUE),
    Away_xg_rolling = sum(Away_xg_rolling, na.rm = TRUE)
  ) |>
  ungroup() |>
  distinct(match_id, .keep_all = TRUE)




# Do all columns in once
df <- read_csv("lessen/Thesis/Processed_full_stats/Full_df/full_df_with_other_first_goal_time.csv") |>
  rename(
    Home_xG = xG,
    Away_xG = xG.1, 
    Home_PSxG = PSxG_home,
    Away_PSxG = PSxG_away
  ) |>
  mutate(Date = as.Date(Date)) |>
  arrange(Date) |>
  mutate(match_id = row_number())

df <- df |> 
  mutate(season = str_extract(comp_season, "\\d{4}(?=\\.csv)")) |> 
  relocate(season, .after = comp_season) |>
  mutate(season = ifelse(is.na(season), str_extract(comp_season, "\\d{4}(?=_regular_season\\.csv|_playoffs\\.csv)"), season))



numerics <- df |>
  select(where(is.numeric)) |>
  select(-c(match_id, Wk, Attendance, diff_corners, number_goals, total_xg, total_PSxG, total_corners, 
            contains("formation"), contains("id"), contains("minute"), contains("_rolling"))) |> 
  colnames()


vars_to_roll <- numerics[str_detect(numerics, "^Home_")] |>
  str_remove("^Home_")



prep_team_df <- function(df, is_home_team) {
  
  prefix_own <- if(is_home_team) "Home_" else "Away_"
  prefix_opp <- if(is_home_team) "Away_" else "Home_"
  
  df |>
    select(match_id, Date, season, 
           Team = if(is_home_team) "Home" else "Away", 
           all_of(paste0(prefix_own, vars_to_roll)),     
           all_of(paste0(prefix_opp, vars_to_roll))) |> 
    rename_with(~ str_remove(., prefix_own) |> paste0("_for"), starts_with(prefix_own)) |>
    rename_with(~ str_remove(., prefix_opp) |> paste0("_against"), starts_with(prefix_opp)) |>
    mutate(is_home = if(is_home_team) 1 else 0)
}


df_home_long <- prep_team_df(df, is_home_team = TRUE)
df_away_long <- prep_team_df(df, is_home_team = FALSE)

df_long <- bind_rows(df_home_long, df_away_long) |>
  arrange(Team, Date)


roll_fn <- function(x) {
  lag(slide_dbl(x, mean, .before = 4, .complete = FALSE), 1)
}


df_long_rolled <- df_long |>
  group_by(Team) |> 
  mutate(across(
    .cols = ends_with(c("_for", "_against")), 
    .fns = roll_fn,
    .names = "{.col}_roll5" 
  )) |>
  ungroup()



new_roll_cols <- colnames(df_long_rolled) |> str_subset("_roll5$")

home_feats <- df_long_rolled |>
  filter(is_home == 1) |>
  select(match_id, all_of(new_roll_cols)) |>
  rename_with(~ paste0("Home_", .), all_of(new_roll_cols))

away_feats <- df_long_rolled |>
  filter(is_home == 0) |>
  select(match_id, all_of(new_roll_cols)) |>
  rename_with(~ paste0("Away_", .), all_of(new_roll_cols))

df_final <- df |>
  left_join(home_feats, by = "match_id") |>
  left_join(away_feats, by = "match_id")



df_final |>
  filter(Home == "Anderlecht") |>
  select(Date, Home, Away, Home_xG_for_roll5, Home_shots_for_roll5) |>
  head(10)

# Print missing rows
missing_rows <- df_final |>
  filter(is.na(Home_xG_for_roll5) | is.na(Away_xG_against_roll5)) 

# Delete following columns containing _rolling
cols_to_remove <- colnames(df_final) |> str_subset("_rolling$")
df_final <- df_final |> select(-all_of(cols_to_remove))

# Save the final dataframe
saveRDS(df_final, "df_final_with_rolling_averages.rds")

# Check final df with no missings for all games on all _roll5 columns?
final_df <- readRDS("df_final_with_rolling_averages.rds")
cols_to_check <- colnames(final_df) |> str_subset("_roll5$")
missing_check <- final_df |>
  filter(if_any(all_of(cols_to_check), is.na))
# Count number of unique values for Home
unique_home <- final_df |> 
  select(Home) |> 
  distinct() |> 
  nrow()

