### Load libraries -------------------------------------------------------------
# test comment

# Data wrangling libraries
library(dplyr)
library(tidyverse)
library(here)
library(zoo)

### Load datasets --------------------------------------------------------------

# Load land-use history data
hist <- read.delim(here("dat",
                        "derived",
                        "230720_histLong.csv"),
                   sep = ",")

# Load groundwater data
dat_big <- readRDS(here("dat",
                        "210215_bigData2_all.rds")) %>% 
  filter(Tiefenstufe2 == "0-30cm")


### Filter dataset -------------------------------------------------------------

dat <- dat_big %>%
  # Remove organic soils
  filter(is.na(BZE_Moor) | BZE_Moor == 0) %>% 
  # Remove sites with ley rotations
  filter(is.na(FB_L.Landnutzung_aktuell) | 
           FB_L.Landnutzung_aktuell != "Gruenland-Wechselw.") %>% 
  rename(sampled_soc = TOC_Stock)


### Wrangle history data -------------------------------------------------------

# Generate land use change variables
hist_reduced <- hist %>% 
  arrange(PointID, year) %>% 
  mutate(change = as.integer(factor(use))) %>% 
  group_by(PointID) %>% 
  mutate(change = change - lag(change),
         change = ifelse(is.na(change) | change == 0, 0, 1),
         change = cumsum(change)) %>% 
  mutate(use.changes = max(change)) %>% 
  mutate(init.use = first(na.omit(use)), 
         curr.use = last(na.omit(use))) %>% 
  mutate(not_ag = any(!grepl("^[AG]+$", use))) %>% 
  filter(not_ag == FALSE) %>% 
  mutate(history = 
           case_when((init.use == "A" & curr.use == "G" & use.changes == 1) 
                     ~ "A_to_G",
                     (init.use == "G" & curr.use == "A" & use.changes == 1) 
                     ~ "G_to_A",
                     (init.use == "A" & curr.use == "A" & use.changes == 0) 
                     ~ "only_A",
                     (init.use == "G" & curr.use == "G" & use.changes == 0) 
                     ~ "only_G",
                     TRUE ~ "flipflop")) %>%
  ungroup()

# Check number of sites in each land-use change category
hist_reduced %>%
  group_by(history) %>% 
  summarise(n = n_distinct(PointID)) %>% 
  ungroup()


# Establish dataset for permanent sites and those with only one land-use change
hist_one_change <- hist_reduced %>%
  filter((history == "A_to_G" & change == 1) | 
           (history == "G_to_A" & change == 1) |
           history == "only_A" | 
           history == "only_G") %>% 
  group_by(PointID) %>% 
  mutate(year_change = ifelse(history == "A_to_G" | history == "G_to_A",
                              min(year), NA),
         max_hist = max(year),
         time_since_change = max_hist - year_change) %>% 
  slice(1) %>% 
  ungroup()

### YOU ARE HERE###
# Collecting only unique values for change_dynamic removes some repeat GtoA_1
# Establish dataset for flipflop sites
hist_flipflop <- hist_reduced %>%
  filter(history == "flipflop") %>% 
  mutate(usechange = paste0(use, change)) %>% 
  group_by(PointID, usechange) %>% 
  mutate(duration = n(),
         change_direction = case_when(change == 0 ~ use,
                                      use == "A" ~ "G_to_A",
                                      use == "G" ~ "A_to_G"),
         year_changes = min(year)) %>% 
  group_by(PointID) %>% 
  mutate(change_years = list(unique(year_changes)),
         max_hist = max(year)) %>% 
  filter(year %in% unlist(change_years)) %>%
  ungroup()


### Combine filtered land history dataset with lab data and groundwater data ---

one_change <- hist_one_change %>% 
  inner_join(dat, by = "PointID") %>%
  # left_join(groundwater, by = "PointID") %>%
  select(PointID, history, year_change, max_hist, 
         time_since_change, sampled_soc) %>% 
  group_by(PointID) %>% 
  slice(1) %>% 
  ungroup()

flipflop <- hist_flipflop %>% 
  inner_join(dat, by = "PointID") %>%
  # left_join(groundwater, by = "PointID") %>%
  select(PointID, history, year_changes, change_direction, duration, 
         max_hist, sampled_soc) %>% 
  group_by(PointID) %>% 
  # slice(1) %>% 
  ungroup()


rm(list = setdiff(ls(), c("one_change", "flipflop")))

# # Test, remove after. Filter only sites that have undergone change
# onechange <- one_change |> 
#   filter(history != "only_A") |> 
#   filter(history != "only_G")
# n_distinct(onechange$PointID)
# n_distinct(flipflop$PointID)

# ff <- flipflop |> select(PointID) |> distinct()
# oc <- onechange |> select(PointID) |> distinct()

# both <- oc |> 
#   rbind(ff)

# test <- luc_summary |> 
#   filter(PointID %in% both$PointID)

### Write combined dataset to csv ----------------------------------------------
# write.csv(combined, 
#     file = "./dat/derived/fig_6_dat.csv", 
#     row.names = FALSE)
