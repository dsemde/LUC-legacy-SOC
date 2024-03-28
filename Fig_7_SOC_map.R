### Load libraries -------------------------------------------------------------

# Data wrangling libraries
library(tidyverse)
library(here)
library(boot)

# Spatial libraries
library(sf)

# Plotting libraries
library(ggplot2)
library(davesfaves)
library(ggspatial)
library(ggpubr)

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

  
# SOC change ####
soc_mineral <- readRDS(here("dat",
                            "210215_bigData2_all.rds")) |> 
  filter(Tiefenstufe2 == "0-30cm") |> 
  # Remove sites with organic soils
  filter(is.na(BZE_Moor) | BZE_Moor == 0) |>
  # Remove sites with ley rotations
  filter(is.na(FB_L.Landnutzung_aktuell) | 
           FB_L.Landnutzung_aktuell != "Gruenland-Wechselw.") |> 
  select(PointID,
         lon_Grid,
         lat_Grid,
         TOC_Stock) |> 
  rename(soc_today = TOC_Stock)

# Get sites with no LUC
no_luc <- luc |> 
  filter(count == 1 & agr_hist == 1) |>  
  group_by(PointID) |> 
  slice(1) |> 
  ungroup() |> 
  inner_join(soc_mineral) |> 
  select(PointID)

# n_distinct(luc_summary$PointID[luc_summary$duration > 100]) # 115 sites

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

set.seed(123)

## Net SOC status of German agriclulture
get_net_soc <- function(data, idx) {
  dat <- data[idx, ]
  sum(dat$soc_change_tot, na.rm = TRUE) * 6400
}

sum(luc_summary$soc_change_tot, na.rm = TRUE) * 6400

# 95 % confidence interval
boot.ci(boot(luc_summary, get_net_soc, R = 10000), type = "norm")$normal


## Net SOC change rate of German Agriculture
# Separate numbers for increasing and decreasing sites
get_soc_change_rate <- function(data, idx) {
  dat <- data[idx, ]
  mean(dat$soc_change_rate, na.rm = TRUE)
}

# Grass to crop
luc_summary |> 
  filter(direction == 1) %>% 
  summarise(n = n(),
            mean = mean(soc_change_rate, na.rm = TRUE),
            ci = boot.ci(boot(., get_soc_change_rate, R = 10000),
                         type = "norm")$normal)

# Crop to grass
luc_summary |> 
  filter(direction == 2) %>% 
  summarise(n = n(),
            mean = mean(soc_change_rate, na.rm = TRUE),
            ci = boot.ci(boot(., get_soc_change_rate, R = 10000),
                         type = "norm")$normal)

# Number of flip-flop sites
luc_summary |> 
  filter(count == 9) |> 
  nrow()

### Load Germany map and geometry ----------------------------------------------

germany <- st_read(here("dat",
                        "ger",
                        "target_area.gpkg"))

coords <- readRDS(here("dat",
                       "raw",
                       "BZE_Standorte.rda")) |> 
  select(PointID, X_Koordinate_inUTM32, Y_Koordinate_inUTM32)

fig_dat <- coords |>  
  right_join(luc_summary) |> 
  st_as_sf(coords = c("X_Koordinate_inUTM32", 
                      "Y_Koordinate_inUTM32"),
           crs = 25832)

no_luc <- no_luc |>
  left_join(coords) |> 
  st_as_sf(coords = c("X_Koordinate_inUTM32", 
                      "Y_Koordinate_inUTM32"),
           crs = 25832)

### SOC change maps ------------------------------------------------------------

units_a <- expression(SOC~stock~change~(Mg~ha^{-1})) #nolint
units_b <- expression(SOC~stock~change~rate~(Mg~ha^{-1}~yr^{-1})) #nolint

## SOC change absolute
p <- ggplot() + 
  ggtitle(bquote(bold("A)") ~ "Cumulative SOC stock change over 100 years")) +
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
                   style = "ticks") +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_void() +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.15, "cm"),
        axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.ticks = element_line(), 
        text = element_text(family = "Calibri",
                            size = 12),
        plot.title = element_text(size = 12),
        axis.text = element_text(family = "Calibri",
                                 size = 8),
        axis.text.y = element_text(margin = margin(l = -15, r = 5)),
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        legend.position = "top",
        legend.margin = margin(t = 10),
        legend.text = element_text(size = 8))

p1 <- annotate_figure(p, text_grob(units_a,
                                   size = 8,
                                   family = "Calibri",
                                   hjust = 0.41,
                                   vjust = 7.5))

## SOC change rate
p <- ggplot() + 
  ggtitle(bquote(bold("B)") ~ "Current SOC stock change rate")) +
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
                         style = north_arrow_minimal) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_void() +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key.height = unit(0.15, "cm"),
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(), 
        text = element_text(family = "Calibri",
                            size = 12),
        plot.title = element_text(size = 12),
        axis.text = element_text(family = "Calibri",
                                 size = 8),
        axis.text.y = element_blank(),
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        legend.position = "top",
        legend.margin = margin(t = 10),
        legend.text = element_text(size = 8))

p2 <- annotate_figure(p, text_grob(units_b,
                                   size = 8,
                                   family = "Calibri",
                                   hjust = 0.45,
                                   vjust = 7.5))


## Combine maps
p3 <- ggarrange(p1, p2, ncol = 2, 
                widths = c(0.53, 0.49))

annotate_figure(p3, 
                bottom = text_grob("Longitude", 
                                   family = "Calibri",
                                   hjust = 0.25,
                                   size = 12),
                left = text_grob("Latitude", 
                                 rot = 90, 
                                 just = 0.5,
                                 family = "Calibri",
                                 size = 12))

ggsave(filename = here("figures",
                       "fig_7_soc_change_maps.png"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "white",
       dpi = 300)
