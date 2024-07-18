### Load libraries -------------------------------------------------------------

# Data wrangling libraries
library(tidyverse)
library(here)
library(httr)
library(readxl)
library(sf)

### Load datasets --------------------------------------------------------------

## Load land-use history data
hist <- read.delim(here("dat",
                        "230720_histLong.csv"),
                   sep = ",")


## Load site data from Open Agrar
# https://doi.org/10.3220/DATA20200203151139

# URL for site data
site_url <- "https://www.openagrar.de/servlets/MCRFileNodeServlet/openagrar_derivate_00028499/BZE_LW%20English%20Version/SITE.xlsx"

# Create a temporary file
temp_file <- tempfile(fileext = ".xlsx")

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
pointIDs <- site |> distinct(PointID) |> pull(PointID)

# Generate second column with 10% 1s and 90% 0s (similar to actual data)
set.seed(123)
values <- sample(c(1,0), length(pointIDs), replace = TRUE, prob = c(0.10, 0.90))

# List of pointIDs with ley rotations
ley_list <- data.frame(PointID = pointIDs, ley_rot = values) |>
  filter(ley_rot == 1) |>
  pull(PointID)


### Filter dataset -------------------------------------------------------------

loc <- site |>
  # Remove organic soils
  filter(!PointID %in% organic_list) |>
  # Remove sites with ley rotations
  filter(!PointID %in% ley_list)


### Wrangle history data -------------------------------------------------------

# Generate land use change variables
hist_reduced <- hist |> 
  arrange(PointID, year) |> 
  mutate(change = as.integer(factor(use))) |> 
  group_by(PointID) |> 
  mutate(change = change - lag(change),
         change = ifelse(is.na(change) | change == 0, 0, 1),
         change = cumsum(change)) |> 
  mutate(use.changes = max(change)) |> 
  mutate(init.use = first(na.omit(use)), 
         curr.use = last(na.omit(use))) |> 
  mutate(not_ag = any(!grepl("^[AG]+$", use))) |> 
  filter(not_ag == FALSE) |> 
  mutate(history = 
           case_when((init.use == "A" & curr.use == "G" & use.changes == 1) 
                     ~ "A_to_G",
                     (init.use == "G" & curr.use == "A" & use.changes == 1) 
                     ~ "G_to_A",
                     (init.use == "A" & curr.use == "A" & use.changes == 0) 
                     ~ "only_A",
                     (init.use == "G" & curr.use == "G" & use.changes == 0) 
                     ~ "only_G",
                     TRUE ~ "flipflop")) |>
  slice(1) |> 
  ungroup()
  

# Check number of sites in each land-use change category
hist_reduced |>
  group_by(history) |> 
  summarise(n = n_distinct(PointID)) |> 
  ungroup()


### Combine datasets to filter for manuscript 1 sites --------------------------

bze_loc <- hist_reduced |> 
  inner_join(loc, by = "PointID") |>
  select(PointID, history, xcoord, ycoord) |>
  st_as_sf(coords = c("xcoord",
                      "ycoord"),
           crs = 25832)

## Remove all environmental variables except bze_loc
rm(list = setdiff(ls(), "bze_loc"))
