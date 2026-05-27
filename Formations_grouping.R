# Grouping the formations into 8 categories based on the paper of Le Coz et al.:
# A competing risk survival analysis of the impacts of team formation on goals in professional football
library(tidyverse)
full_df_6 <- read.csv("~/lessen/Thesis/Processed_full_stats/Full_df/full_df_6.csv")
df <- full_df_6
df |> colnames()
df |> 
  count(home_formation) |> 
  arrange(desc(n))
# Based on the above mentioned paper:
# we will group 3-5-2, 3-1-4-2, 3-5-1-1, 3-4-1-2 in 3-5-2.
# 3-4-3 and 3-3-3-1 in 3-4-3 (see screenshot why)
# 4-3-1-2, 4-4-2, 4-4-1-1, 4-2-2-2, 4-1-3-2, 4-1-2-1-2 in 4-4-2
# 4-5-1, 4-1-4-1, 4-3-2-1 in 4-5-1
# 5-3-2 stays 5-3-2
# 5-4-1, 3-2-4-1 in 5-4-1
# 4-2-3-1 stays 4-2-3-1 and also added 4-2-4 (looks like this formation is similar)
# 4-3-3 stays 4-3-3
df <- df |> 
  mutate(home_formation_grouped = case_when(
    home_formation %in% c("3-5-2", "3-5-1-1", "3-1-4-2", "3-4-1-2") ~ "3-5-2", 
    home_formation %in% c("3-4-3", "3-3-3-1") ~ "3-4-3",
    home_formation %in% c("4-3-1-2", "4-2-2-2", "4-4-2", "4-4-1-1", "4-1-3-2", "4-1-2-1-2") ~ "4-4-2", 
    home_formation %in% c("4-5-1", "4-3-2-1", "4-1-4-1") ~ "4-5-1", 
    home_formation %in% c("5-4-1", "3-2-4-1") ~ "5-4-1",
    home_formation %in% c("4-2-3-1", "4-2-4") ~ "4-2-3-1",
    home_formation %in% c("4-3-3") ~ "4-3-3",
    home_formation %in% c("3-2-4-1") ~ "3-2-4-1",
    home_formation %in% c("5-3-2") ~ "5-3-2"
  ))
df |> 
  count(home_formation_grouped) |> mutate(percentage = (n/nrow(df))*100) |> 
  arrange(desc(n))

# Same for the away team
df <- df |> 
  mutate(away_formation_grouped = case_when(
    away_formation %in% c("3-5-2", "3-5-1-1", "3-1-4-2", "3-4-1-2") ~ "3-5-2", 
    away_formation %in% c("3-4-3", "3-3-3-1") ~ "3-4-3",
    away_formation %in% c("4-3-1-2", "4-2-2-2", "4-4-2", "4-4-1-1", "4-1-3-2", "4-1-2-1-2") ~ "4-4-2", 
    away_formation %in% c("4-5-1", "4-3-2-1", "4-1-4-1") ~ "4-5-1", 
    away_formation %in% c("5-4-1", "3-2-4-1") ~ "5-4-1",
    away_formation %in% c("4-2-3-1", "4-2-4") ~ "4-2-3-1",
    away_formation %in% c("4-3-3") ~ "4-3-3",
    away_formation %in% c("3-2-4-1") ~ "3-2-4-1",
    away_formation %in% c("5-3-2") ~ "5-3-2"
  ))
df |> 
  count(away_formation_grouped) |> mutate(percentage = (n/nrow(df))*100) |> 
  arrange(desc(n))
write_csv(df, "~/lessen/Thesis/Processed_full_stats/Full_df/full_df_7.csv")
