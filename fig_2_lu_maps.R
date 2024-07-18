### Load libraries -------------------------------------------------------------

## Data wrangling libraries
library(tidyverse)
library(here)

## PLotting libraries
library(ggplot2)
library(ggpubr)
library(showtext)
library(cowplot)

## Spatial libraries
library(sf)
library(ggspatial)


### Load datasets and prepare data ---------------------------------------------

## Load coordinates of BZE sites
source(here("site_locations.R"))


## Load Germany map
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

## Load font for use in figures
font_add_google("Montserrat")
showtext_auto()

## Set plot text sizes
text_size <- 32
axis_label_size <- 36

p1 <- ggplot() +
  ggtitle(bquote(bold("(a)") ~ "Sites under permanent land use")) +
  geom_sf(data = germany,
          fill = NA) +
  geom_sf(data = permanent_use_dat,
          aes(shape = history, col = history),
          size = 2) +
  scale_shape_manual(labels = hist_lab, values = hist_shape) +
  scale_color_manual(labels = hist_lab, values = hist_col) +
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
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.ticks = element_line(),
        text = element_text(family = "Montserrat",
                            size = text_size),
        axis.text.y = element_text(margin = margin(l = -40, r = -15)),
        axis.text.x = element_text(margin = margin(t = 0, b = 0)),
        legend.position = "top",
        legend.key.spacing.x =  unit(-205, "pt"),
        legend.margin = margin(l = 200),
        legend.text = element_text(margin = margin(l = -2),
                                   size = text_size))


## Transitional LU map
transitional_use_dat <- bze_loc |>
  filter(history == "A_to_G" | history == "G_to_A" | history == "flipflop") |>
  mutate(history = factor(history, levels = c("G_to_A", "A_to_G", "flipflop")))

p2 <- ggplot() +
  ggtitle(bquote(bold("(b)") ~ "Sites with a history of land-use change")) +
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
                         style = north_arrow_minimal(text_size = 32,
                                                     text_col = "#333333",
                                                     line_col = "#333333",
                                                     fill = "#333333",
                                                     text_family = "Montserrat")) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme_void() +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        text = element_text(family = "Montserrat",
                            size = text_size),
        axis.text.y = element_blank(),
        axis.text.x = element_text(margin = margin(t = 0, b = 0)),
        legend.position = "top",
        legend.key.spacing.x =  unit(-125, "pt"),
        legend.margin = margin(l = 150),
        legend.text = element_text(margin = margin(l = -2),
                                   size = text_size))


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
                       "fig_2_lu_maps_final.tiff"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "white",
       dpi = 300,
       compression = "lzw")

ggsave(filename = here("figures",
                       "fig_2_lu_maps_final.png"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "transparent",
       dpi = 300)