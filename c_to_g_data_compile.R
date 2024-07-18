### Load libraries -------------------------------------------------------------

# Data wrangling libraries
library(tidyverse)
library(here)
library(readxl)
library(httr)


### Load datasets --------------------------------------------------------------

## Load land-use history data
hist <- read.delim(here("dat",
                        "230720_histLong.csv"),
                   sep = ",")

## Load lab and site data from Open Agrar
# https://doi.org/10.3220/DATA20200203151139

# URLs for lab and site data
lab_url <- "https://www.openagrar.de/servlets/MCRFileNodeServlet/openagrar_derivate_00028499/BZE_LW%20English%20Version/LABORATORY_DATA.xlsx"
site_url <- "https://www.openagrar.de/servlets/MCRFileNodeServlet/openagrar_derivate_00028499/BZE_LW%20English%20Version/SITE.xlsx"

# Create a temporary file
temp_file <- tempfile(fileext = ".xlsx")

# Download and read lab data file
GET(lab_url, write_disk(temp_file, overwrite = TRUE))
lab <- read_excel(temp_file, na = "NA") |>
  rename(depth = "Layer lower limit")

# Download and read site data file
GET(site_url, write_disk(temp_file, overwrite = TRUE))
site <- read_excel(temp_file, na = "NA") |>
  rename(gw_depth = "Groundwater class")

unlink(temp_file)


## List of organic soil PointIDs
organic_list <- site |>
  filter(BZE_peat == 1) |>
  pull(PointID)


## Generate sample lay rotation data
# These data do not represent actual site history of ley rotations and are only
# used for demonstration purposes
# ley_rot: 0 = no ley rotation, 1 = ley rotation
pointIDs <- lab |> distinct(PointID) |> pull(PointID)

# Generate second column with 10% 1s and 90% 0s (similar to actual data)
set.seed(123)
values <- sample(c(1,0), length(pointIDs), replace = TRUE, prob = c(0.10, 0.90))

# List of pointIDs with ley rotations
ley_list <- data.frame(PointID = pointIDs, ley_rot = values) |>
  filter(ley_rot == 1) |>
  pull(PointID)


### Filter dataset -------------------------------------------------------------

lab <- lab |>
  # Remove organic soils
  filter(!PointID %in% organic_list) |>
  # Remove sites with ley rotations
  filter(!PointID %in% ley_list)


### Aggregate lab data to three depth increments -------------------------------

lab_aggregated <- lab |>
  group_by(PointID) |>
  # Define depth increments
  mutate(depth_increment = case_when(depth <= 10 ~ "0-10cm",
                                  depth <= 30 ~ "10-30cm",
                                  depth <= 100 ~ "30-100cm")) |>
  filter(!is.na(depth_increment)) |>
  # Calculate horizon / soil layer thickness
  # ... for aggregation of e.g. bulk density
  group_by(PointID) |>
  arrange(PointID,
          depth) |>
  mutate(depth_diff = depth - lag(depth, default = 0)) |>
  # Aggregation
  group_by(PointID,
           depth_increment) |>
  summarise(
      # FSS is fine soil stock in Mg ha-1
      # Total organic carbon content
      # Unit: g kg-1 of fine soil < 2mm (mass per mass)
      TOC = weighted.mean(TOC,
                          w = FSS,
                          na.rm = TRUE),
      # Dry bulk density
      # Unit: g cm-3 (mass per volume)
      BD  = weighted.mean(BD_bulk,
                          w = depth_diff,
                          na.rm = TRUE),
      # Clay content
      # Unit: mass % of fine soil < 2mm
      clay = weighted.mean(Clay,
                          w = FSS,
                          na.rm = TRUE),
      # C:N ratio
      cn_ratio = weighted.mean((TOC/TN),
                          w = FSS,
                          na.rm = TRUE)) |>
  # Calculate organic carbon stock and stock per cm of soil depth
  mutate(TOC_stock =
      (TOC * BD * case_when(depth_increment == "0-10cm" ~ 10,
                          depth_increment == "10-30cm" ~ 20,
                          depth_increment == "30-100cm" ~ 70)) * 0.1) |>
  mutate(TOC_stock_cm =
      (TOC * BD * case_when(depth_increment == "0-10cm" ~ 10,
                          depth_increment == "10-30cm" ~ 20,
                          depth_increment == "30-100cm" ~ 70)) /
                          case_when(depth_increment == "0-10cm" ~ 10,
                              depth_increment == "10-30cm" ~ 20,
                              depth_increment == "30-100cm" ~ 70) * 0.1) |>
  ungroup()


### Wrangle history data -------------------------------------------------------

hist_reduced <- hist |> 
  # Generate land use change variables
  arrange(PointID, year) |>
  mutate(change = as.integer(factor(use))) |>
  group_by(PointID) |>
  mutate(change = change - lag(change),
          change = ifelse(is.na(change) | change == 0, 0, 1),
          change = cumsum(change)) |>
  mutate(use.changes = max(change)) |>
  mutate(init.use = first(na.omit(use)),
          curr.use = last(na.omit(use))) |>
  mutate(history = case_when((init.use == "A" & curr.use == "G" & use.changes == 1) ~ "C_to_G",
                              (init.use == "A" & curr.use == "A" & use.changes == 0) ~ "only_C",
                              (init.use == "G" & curr.use == "G" & use.changes == 0) ~ "only_G")) |>
  filter((history == "C_to_G" & change == 1) |
           history == "only_C" |
           history == "only_G") |>
  ungroup() |>
  # Carry forward only years since the first conversion of A to G
  mutate(min_year = min(year[history == "C_to_G"])) |>
  filter(year >= (min(min_year) - 1)) |>
  # Year of land-use change
  group_by(PointID) |>
  mutate(year_change = ifelse(history == "C_to_G", min(year), NA)) |>
  slice(1) |>
  ungroup()


### Combine filtered land history dataset with lab data and groundwater data ---

combined <- hist_reduced |> 
    inner_join(lab_aggregated, by = "PointID") |>
    left_join(site, by = "PointID") |>
    select(PointID, xcoord, ycoord, depth_increment, year_change, history, TOC,
            TOC_stock, TOC_stock_cm, clay, cn_ratio, gw_depth) |>
    ungroup()


### Write combined dataset to csv ----------------------------------------------
write.csv(combined,
    file = "./dat/derived/c_to_g_3_depths.csv",
    row.names = FALSE)
