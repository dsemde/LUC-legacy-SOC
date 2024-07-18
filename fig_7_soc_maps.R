### Load libraries -------------------------------------------------------------

# Data wrangling libraries
library(tidyverse)
library(here)
library(httr)
library(readxl)

# Spatial libraries
library(sf)
library(ggspatial)

# Plotting libraries
library(ggplot2)
library(davesfaves)
library(ggpubr)
library(showtext)
library(cowplot)

### Load datasets and prepare data ---------------------------------------------

lu <- read.delim(here("dat",
                      "230720_histLong.csv"),
                 sep = ",") |> 
  
  # Reclassify
  # 1 = cropland
  # 0 = grassland
  # 9 = other
  mutate(use = case_when(use == "A" ~ 1,
                         use == "G" ~ 0,
                         TRUE ~ 9)) 

# Number, direction and duration of LUCs ####
# Sites with always agriculture
luc <- lu |> 
  group_by(PointID) |> 
  
  # History
  # 0 = non-agricultural history
  # 1 = always agriculture
  mutate(agr_hist = max(use),
         agr_hist = ifelse(agr_hist == 9, 0, 1)) |> 
  
  # Where always agriculture, identify land use change (LUC) direction
  # 0 = no change
  # 1 = grass to crop
  # -1 = crop to grass
  mutate(lu_agr = ifelse(use == 9, NA, use),
         lu_agr_diff = lu_agr - lag(lu_agr)) |> 
  
  # ... and calculate number of land use changes (LUCs)
  # 1 = no LUC
  # 2 = one LUC 
  # 9 = multiple LUCs
  mutate(tmp = ifelse(lu_agr_diff %in% c(1, -1), 1, 0),
         count_cum = cumsum(tmp),
         count = sum(tmp, na.rm = TRUE),
         count = ifelse(count >= 2, 9, count + 1)) |> 
  
  # Direction of first LUC
  # 1 = grass to crop
  # 2 = crop to grass
  mutate(direction = case_when(count_cum == 1 & lu_agr_diff == 1 ~ 1,
                               count_cum == 1 & lu_agr_diff == -1 ~ 2)) |> 

  # Separate years before and after 100 years before sampling
  # 1 = years within 100 years of sampling
  # 0 = years outside of 100 years of sampling
  mutate(year_100 = ifelse(year > max(year) - 100, 1, 0)) |> 
  ungroup()

luc_summary <- luc |> 
  group_by(PointID) |> 
  filter(count > 1) |> 
  filter(year >= year[!is.na(direction)]) |> 
  summarise(agr_hist = mean(agr_hist),
            count = mean(count),
            direction = mean(direction, na.rm = TRUE),
            # Number of years with new land use
            duration = sum(use == use[1], na.rm = TRUE), 
            pre_duration = sum(use == use[1] & year_100 == 0, na.rm = TRUE))

  
# ## Load SOC data
# soc_mineral <- readRDS(here("dat",
#                             "210215_bigData2_all.rds")) |>
#   filter(Tiefenstufe2 == "0-30cm") |>
#   # Remove sites with organic soils
#   filter(is.na(BZE_Moor) | BZE_Moor == 0) |>
#   # Remove sites with ley rotations
#   filter(is.na(FB_L.Landnutzung_aktuell) |
#            FB_L.Landnutzung_aktuell != "Gruenland-Wechselw.") |>
#   select(PointID,
#          lon_Grid,
#          lat_Grid,
#          TOC_Stock) |>
#   rename(soc_today = TOC_Stock)

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

## Combine sampling depths into study depths
soc_mineral <- lab |>
  group_by(PointID) |>
  # Define depth increments
  mutate(depth_increment = case_when(depth <= 30 ~ "0-30cm",
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
                          na.rm = TRUE)) |>
  # Calculate organic carbon stock and stock per cm of soil depth
  mutate(TOC_stock =
      (TOC * BD * case_when(depth_increment == "0-30cm" ~ 30,
                          depth_increment == "30-100cm" ~ 70) * 0.1)) |>
  ungroup() |>

  # Bind with site coordinates
  left_join(site, by = "PointID") |>
  filter(depth_increment == "0-30cm") |>
  filter(!is.na(TOC_stock)) |>
  select(PointID,
         xcoord,
         ycoord,
         TOC_stock) |>
  rename(soc_today = TOC_stock)


## Get sites with no LUC
no_luc <- luc |> 
  filter(count == 1 & agr_hist == 1) |>  
  group_by(PointID) |> 
  slice(1) |> 
  ungroup() |> 
  inner_join(soc_mineral) |> 
  select(PointID)


luc_summary <- luc_summary |> 
  # Only mineral soil
  inner_join(soc_mineral) |> 
  
  # Cumulative SOC change
  mutate(soc_change_rel = 
           case_when(direction == 1 ~ 
                       -33.57 * (1 - exp(-0.0152 * duration)) / 100,
                     direction == 2 ~  
                       47.26 * (1 - exp(-0.047 * duration)) / 100),
         pre_soc_change_rel = 
           case_when(direction == 1 ~ 
                       -33.57 * (1 - exp(-0.0152 * pre_duration)) / 100,
                     direction == 2 ~  
                       47.26 * (1 - exp(-0.047 * pre_duration)) / 100),
         soc_init = soc_today / (soc_change_rel + 1),
         soc_change_abs = soc_today - soc_init,
         pre_soc_change_abs = soc_init * pre_soc_change_rel, 
         soc_change_tot = soc_change_abs - pre_soc_change_abs) |> 
  
  # Current SOC change rates due to historic land-use change
  # Unit: Mg ha-1 yr-1
  # Analysis is restricted to sites with single LUC 
  mutate(soc_change_rel_previous_year = 
           case_when(direction == 1 ~ 
                       -33.57 * (1 - exp(-0.0152 * (duration - 1))) / 100,
                     direction == 2 ~  
                       47.26 * (1 - exp(-0.047 * (duration - 1))) / 100),
         soc_previous_year = soc_init * (1 + soc_change_rel_previous_year),
         soc_change_rate = soc_today - soc_previous_year,
         soc_change_rate = ifelse(count == 9, NA, soc_change_rate))


### Manuscript values ----------------------------------------------------------
#
# set.seed(123)
#
# ## Net SOC status of German agriclulture
# get_net_soc <- function(data, idx) {
#   dat <- data[idx, ]
#   sum(dat$soc_change_tot, na.rm = TRUE) * 6400
# }
#
# sum(luc_summary$soc_change_tot, na.rm = TRUE) * 6400
#
# # 95 % confidence interval
# boot.ci(boot(luc_summary, get_net_soc, R = 10000), type = "norm")$normal
#
#
# ## Net SOC change rate of German Agriculture
# # Separate numbers for increasing and decreasing sites
# get_soc_change_rate <- function(data, idx) {
#   dat <- data[idx, ]
#   mean(dat$soc_change_rate, na.rm = TRUE)
# }
#
# # Grass to crop
# luc_summary |>
#   filter(direction == 1) %>%
#   summarise(n = n(),
#             mean = mean(soc_change_rate, na.rm = TRUE),
#             ci = boot.ci(boot(., get_soc_change_rate, R = 10000),
#                          type = "norm")$normal)
#
# # Crop to grass
# luc_summary |>
#   filter(direction == 2) %>%
#   summarise(n = n(),
#             mean = mean(soc_change_rate, na.rm = TRUE),
#             ci = boot.ci(boot(., get_soc_change_rate, R = 10000),
#                          type = "norm")$normal)
#
# # Number of flip-flop sites
# luc_summary |>
#   filter(count == 9) |>
#   nrow()

### Load Germany map and geometry ----------------------------------------------

germany <- st_read(here("dat",
                        "ger",
                        "target_area.gpkg"))

coords <- site |>
  select(PointID, xcoord, ycoord)

fig_dat <- coords |>  
  right_join(luc_summary) |> 
  st_as_sf(coords = c("xcoord",
                      "ycoord"),
           crs = 25832)

no_luc <- no_luc |>
  left_join(coords) |> 
  st_as_sf(coords = c("xcoord",
                      "ycoord"),
           crs = 25832)

### SOC change maps ------------------------------------------------------------

## Load font for use in figures
font_add_google("Montserrat")
showtext_auto()

## Set plot text sizes
text_size <- 32
axis_label_size <- 36


## SOC change absolute
p1 <- ggplot() +
  ggtitle(bquote(bold("(a)") ~ "Cumulative SOC stock change over 100 years")) +
  geom_sf(data = germany,
          fill = NA) + 
  geom_sf(data = fig_dat,
          aes(col = soc_change_tot),
          shape = 16,
          size = 2) +
  geom_sf(data = no_luc,
          col = "#e6e6e6",
          shape = 16,
          size = 0.5) +
  scale_colour_dave_c("daves_faves",
                      direction = "reverse", 
                      mid = 0,
                      midcol = "#e6e6e6") +
  annotation_scale(pad_x = unit(0.75, "cm"),
                   pad_y = unit(0.5, "cm"),
                   style = "ticks",
                   text_cex = 2.5,
                   text_family = "Montserrat",
                   line_col = "#333333",
                   text_col = "#333333") +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_void() +
  annotate("text",
           x = -Inf,
           y = Inf,
           hjust = -1,
           vjust = 1.1,
           label = expression(SOC~stock~change~(Mg~ha^{-1})),
           family = "Montserrat",
           size = 8) +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(5, "pt"),
        axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.ticks = element_line(), 
        text = element_text(family = "Montserrat",
                            size = text_size),
        axis.text.y = element_text(margin = margin(l = -52, r = -15)),
        axis.text.x = element_text(margin = margin(t = 0, b = 0)),
        legend.position = "top",
        legend.margin = margin(t = 10, b = -20),
        legend.text = element_text(size = text_size))


## SOC change rate
p2 <- ggplot() +
  ggtitle(bquote(bold("(b)") ~ "Current SOC stock change rate")) +
  geom_sf(data = germany,
          fill = NA) + 
  geom_sf(data = fig_dat,
          aes(col = soc_change_rate),
          shape = 17,
          size = 2) +
  geom_sf(data = no_luc,
          col = "#e6e6e6",
          shape = 16,
          size = 0.5) +
  scale_colour_dave_c("daves_faves",
                      direction = "reverse", 
                      mid = 0,
                      midcol = "#e6e6e6") +
  annotation_north_arrow(pad_x = unit(0.5, "cm"),
                         pad_y = unit(0.75, "cm"),
                         location = "tr", 
                         style = north_arrow_minimal(text_size = 32,
                                                     text_col = "#333333",
                                                     line_col = "#333333",
                                                     fill = "#333333",
                                                     text_family = "Montserrat")) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_void() +
  annotate("text",
           x = -Inf,
           y = Inf,
           hjust = -0.65,
           vjust = 1.1,
           label = expression(SOC~stock~change~rate~(Mg~ha^{-1}~yr^{-1})),
           family = "Montserrat",
           size = 8) +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(5, "pt"),
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        text = element_text(family = "Montserrat",
                            size = text_size),
        axis.text.y = element_blank(),
        axis.text.x = element_text(margin = margin(t = 0, b = 0)),
        legend.position = "top",
        legend.margin = margin(t = 10, b = -20),
        legend.text = element_text(size = text_size))


## Combine maps
cowplot::set_null_device(cairo_pdf)
p3 <- cowplot::plot_grid(p1, p2)

annotate_figure(p3, 
                bottom = text_grob("Longitude", 
                                   family = "Montserrat",
                                   hjust = 0.45,
                                   size = axis_label_size),
                left = text_grob("Latitude",
                                 rot = 90,
                                 just = 0.5,
                                 family = "Montserrat",
                                 size = axis_label_size))

ggsave(filename = here("figures",
                       "fig_7_soc_change_maps_final.tiff"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "white",
       dpi = 300,
       compression = "lzw")

ggsave(filename = here("figures",
                       "fig_7_soc_change_maps_final.png"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "transparent",
       dpi = 300)