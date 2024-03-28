### Load libraries -------------------------------------------------------------

# Data wrangling libraries
library(tidyverse)
library(here)

# PLotting libraries
library(ggplot2)
library(colorspace)
library(ggpubr)
library(scales)

# Spatial libraries
library(sf)
library(giscoR)
library(rnaturalearth)
library(ggspatial)


### Load datasets and prepare data ---------------------------------------------

# crop_grass_hist <- read.delim(here("dat",
#                                    "cropland_grassland_history.csv"),
#                               sep = ",")

# loc <- readRDS(here("dat",
#                     "raw",
#                     "BZE_Standorte.rda")) |>
#   select(PointID, X_Koordinate_inUTM32, Y_Koordinate_inUTM32)

# bze_dat <- crop_grass_hist |> 
#   left_join(loc) |> 
#   st_as_sf(coords = c("X_Koordinate_inUTM32", 
#                       "Y_Koordinate_inUTM32"),
#            crs = 25832)

source(here("man1_site_locations.R"))


### Load Germany map -----------------------------------------------------------

# worldmap <- ne_countries(scale = "large", 
#                          type = "map_units",
#                          returnclass = "sf")

# germany <- worldmap[worldmap$name == "Germany", ]

germany <- st_read(here("dat",
                        "ger",
                        "target_area.gpkg"))


### Genarate land-use map ------------------------------------------------------

hist_lab <- c("only_A" = "Permanent cropland", 
              "only_G" = "Permanent grassland",
              "G_to_A" = "Grass to crop", 
              "A_to_G" = "Crop to grass",
              "flipflop" = "Multiple LUC")

hist_shape <- c("only_A" = 16, 
                "only_G" = 16, 
                "G_to_A" = 17, 
                "A_to_G" = 17,
                "flipflop" = 17)

hist_col <- c("only_A" = "#FFD079", 
              "only_G" = "#045275", 
              "G_to_A" = "#FFD079",
              "A_to_G" = "#045275",
              "flipflop" = "#089099")

## Permanent LU map 
permanent_use_dat <- bze_loc |> 
  filter(history == "only_A" | history == "only_G")

p1 <- ggplot() + 
  ggtitle(bquote(bold("A)") ~ "Sites under permanent land use")) +
  geom_sf(data = germany,
          fill = NA) + 
  geom_sf(data = permanent_use_dat,
          aes(shape = history, col = history),
          size = 2) +
  scale_shape_manual(labels = hist_lab, values = hist_shape) +
  scale_color_manual(labels = hist_lab, values = hist_col) +
  annotation_scale(pad_x = unit(0.75, "cm"),
                   pad_y = unit(0.5, "cm"),
                   style = "ticks") +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_void() +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.ticks = element_line(), 
        text = element_text(family = "Calibri",
                            size = 12),
        axis.text = element_text(family = "Calibri",
                                 size = 8),
        plot.title = element_text(size = 12),
        axis.text.y = element_text(margin = margin(l = -15, r = 5)),
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        legend.position = "top",
        legend.margin = margin(t = 10),
        legend.text = element_text(size = 8))


## Transitional LU map 
transitional_use_dat <- bze_loc |> 
  filter(history == "A_to_G" | history == "G_to_A" | history == "flipflop") |> 
  mutate(history = factor(history, levels = c("G_to_A", "A_to_G", "flipflop")))

p2 <- ggplot() + 
  ggtitle(bquote(bold("B)") ~ "Sites with a history of LUC")) +
  geom_sf(data = germany,
          fill = NA) + 
  geom_sf(data = transitional_use_dat,
          aes(shape = history, col = history),
          size = 2) +
  scale_shape_manual(labels = hist_lab, values = hist_shape) +
  scale_color_manual(labels = hist_lab, 
                     values = hist_col) +
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
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(), 
        text = element_text(family = "Calibri",
                            size = 12),
        axis.text = element_text(family = "Calibri",
                                 size = 8),
        plot.title = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        legend.position = "top",
        legend.margin = margin(t = 10),
        legend.text = element_text(size = 8))


## Combine maps
p3 <- ggarrange(p1, p2, ncol = 2, widths = c(0.53, 0.49))

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
                       "fig_2_lu_maps.tiff"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "white",
       dpi = 300)
