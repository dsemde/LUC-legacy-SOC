### Load libraries -------------------------------------------------------------

# Data wrangling libraries
library(tidyverse)
library(here)
library(radiant.data)
library(boot)

# Clustering libraries
library(forcats)
library(cartography)
library(cluster)

# Propensity score libraries
library(twang)

# Plotting libraries
library(ggplot2)
library(ggpubr)
library(scales)
library(gridExtra)

# Significance testing libraries
library(boot.pval)

# Parallel processing libraries
library(doParallel)
library(foreach)


# Set seed for reproducibility
set.seed(1337)

### Load dataset ---------------------------------------------------------------
source(here("A_to_G_PS_data_3_depths.R")) #nolint

dat <- combined  %>% 
  filter(history == "A_to_G" & depth_increment == "0-10cm")


### Generate year intervals-----------------------------------------------------

n <- 5
interval_list <- c()

for (i in 1:n) {
  groups <- getBreaks(dat$year_change, nclass = 2, method = "jenks") %>% 
    floor()

  dat <- dat  %>% 
    mutate(history = case_when(between(year_change, groups[2], groups[3]) ~ 
                                 "newest",
                               between(year_change, groups[1], groups[2]) ~ 
                                 paste0(groups[1], "-", groups[2])),
           .default = NA)
  
  dat <- dat %>% 
    filter(history == "newest")
  
  if (i < n) {
    interval_list[i] <- groups[1]
  } else {
    interval_list[i] <- groups[2]
    interval_list[i + 1] <- groups[3]
  }
}

### Assign year intervals to all data ------------------------------------------

conditions <- purrr::map(1:n, function(i) {
  between(combined$year_change, interval_list[i], interval_list[i + 1]) ~ 
    paste0(interval_list[i], "-", interval_list[i + 1])
})

years_cont <- purrr::map(1:n, function(i) {
  between(combined$year_change, interval_list[i], interval_list[i + 1]) ~ 
    max_yc - (floor((interval_list[i] + interval_list[i + 1]) / 2))
})

right <- purrr::map(1:n, function(i) {
  between(combined$year_change, interval_list[i], interval_list[i + 1]) ~ 
    max_yc - interval_list[i]
})

left <- purrr::map(1:n, function(i) {
  between(combined$year_change, interval_list[i], interval_list[i + 1]) ~ 
    (max_yc - 1) - interval_list[i + 1]
})

min_yc <- min(dat$year_change, na.rm = TRUE) - 1
max_yc <- max(dat$year_change, na.rm = TRUE) + 1

dat <- combined %>% 
  mutate(history = 
           case_when(history == "only_A" ~ 
                       "Permanent cropland",
                     history == "only_G" ~ 
                       "Permanent grassland",
                     !!!conditions)) %>% 
  mutate(years_cont = 
           case_when(history == "Permanent cropland" ~ 
                       max_yc,
                     history == "Permanent grassland" ~ 
                       min_yc,
                     !!!years_cont)) %>% 
  mutate(left = 
           case_when(history == "Permanent cropland" ~ 
                       1,
                     history == "Permanent grassland" ~ 
                       1,
                     !!!left)) %>% 
  mutate(right = 
           case_when(history == "Permanent cropland" ~ 
                       1,
                     history == "Permanent grassland" ~ 
                       1,
                     !!!right)) %>%
  mutate(history = factor(history),
         gw_depth = factor(gw_depth))


### Generate PS model ----------------------------------------------------------

depths <- unique(dat$depth_increment)

n_cores <- detectCores() - 1
cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

dat_w <- foreach(i = 1:length(depths), .combine = rbind, .packages = c("twang", "WeightIt")) %dopar% { #nolint
  subset_df <- as.data.frame(dat[dat$depth_increment == depths[i], ])

  set.seed(123)

  mnps_mod <- mnps(history ~ clay + gw_depth + cn_ratio, 
                   data = subset_df, 
                   estimand = "ATE",
                   verbose = FALSE,
                   stop.method = "ks.mean",
                   n.trees = 10000)

  # Convert ps scores to matrix
  ps_list <- mnps_mod$psList
  ps_matrix <- do.call(cbind, ps_list)

  subset_df$w <- get.weights(mnps_mod, stop.method = "ks.mean")
  subset_df
}

stopCluster(cl)

### Summary data ---------------------------------------------------------------
plot_dat <- dat_w %>% 
  group_by(years_cont, depth_increment) %>%
  summarise(n = n_distinct(PointID),
            mean = weighted.mean(TOC_stock_cm, w = w),
            ci = (weighted.sd(TOC_stock_cm, w = w) / sqrt(n)),
            left = max(left),
            right = max(right)) %>%
  ungroup() %>% 
  # Add facet variable for plotting
  mutate(facet = case_when(years_cont == max_yc ~ "Cropland", 
                           years_cont == min_yc ~ "Grassland",
                           TRUE ~ "Transitional"))


### Bootstrap CI determination -------------------------------------------------

# Weighted mean function to call inside bootstrapping
mean.fun <- function(data, idx) {
  dat <- data[idx, ]

  out <- dat %>% 
    summarize(mean = weighted.mean(TOC_stock_cm, w = w))

  return(as.numeric(out))
}

# Loop over year intervals and depths, performing bootstrapping for each
# Get a vector of unique years and depths
inter <- unique(dat_w$years_cont)

cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

# Loop over all possible interval-depth combinations
results <- foreach(i = 1:length(inter), .combine = rbind, .packages = c("boot", "dplyr")) %:% #nolint
  foreach(j = 1:length(depths), .combine = rbind) %dopar% { #nolint
  subset_df <- dat_w[dat_w$years_cont == inter[i] & 
                       dat_w$depth_increment == depths[j], ]
      
  # Calculate the bootstrapped confidence interval
  ci <- boot.ci(boot(subset_df, mean.fun, R = 10000), type = "norm")$normal

  # Add the results to the output data frame
  results <- data.frame(years_cont = inter[i],
                        depth_increment = depths[j],
                        ci_lower = ci[2],
                        ci_upper = ci[3])
  results
}

stopCluster(cl)

# Combine plot data and bootstrapped CI results
plot_dat <- plot_dat %>% 
  left_join(results, by = c("years_cont", "depth_increment"))

### Bootstrap p-value determination --------------------------------------------

# # Calculate the boostrapped pvalue between grassland and each LU interval
# cl <- makeCluster(n_cores, type = "PSOCK")
# registerDoParallel(cl)

# set.seed(123)

# grassdf <- dat_w %>% 
#     filter(history == "Permanent grassland")

# pvals <- foreach(i = 1:length(depths), .combine = rbind, .packages = c("boot", "boot.pval", "dplyr")) %:% #nolint 
#                 foreach(j = 1:length(inter), .combine = rbind) %dopar% { #nolint

#     if (n_distinct(grassdf$depth_increment) > 1) {
#         grass <- grassdf %>% 
#             filter(depth_increment == depths[i])

#         grassboot <- boot(grass, mean.fun, R = 10000)
#     } else if (grassboot$data$depth_increment[1] != grass$depth_increment) {
#         grass <- grassdf %>% 
#         filter(depth_increment == depths[i])

#         grassboot <- boot(grass, mean.fun, R = 10000)
#     }
        
#     subdf <- dat_w %>% 
#         filter(years_cont == inter[j]
#             & depth_increment == depths[i])
#     comparison <- boot(subdf, mean.fun, R = 10000)

#     results <- data.frame(years_cont = inter[j],
#                     depth_increment = depths[i],
#                     pvalue = boot.pval(comparison, theta_null = mean(grassboot$t0))) #nolint
#     results
# }

# stopCluster(cl)

# pvals <- pvals %>% 
#     arrange(years_cont, depth_increment)

### Plot -----------------------------------------------------------------------

# Create offset for lower two depths
plot_dat <- plot_dat %>% 
  mutate(new_left = case_when(depth_increment == "10-30cm" ~ left + 2,
                              depth_increment == "30-100cm" ~ left + 4,
                              TRUE ~ left)) %>% 
  mutate(new_right = case_when(depth_increment == "10-30cm" ~ right + 2,
                               depth_increment == "30-100cm" ~ right + 4,
                               TRUE ~ right))

# Adjust spacing for offset
plot_dat <- plot_dat %>% 
  mutate(new_left = ifelse(years_cont %in% unique(years_cont)[1:(length(unique(years_cont)) - 2)], #nolint
                           new_left, 
                           new_left + 6)) %>% 
  mutate(new_right = ifelse(years_cont %in% unique(years_cont)[1:(length(unique(years_cont)) - 3)], #nolint
                            new_right - 6, 
                            new_right))

# Set colours for depth increments
depth_colours <- c("0-10cm" = "#ffd079",
                   "10-30cm" = "#089099",
                   "30-100cm" = "#045275")

# Set y-limit for all facets
y_limit <- 4.50

# Select grassland limits for horizontal line
h_line <-  plot_dat %>% 
  filter(facet == "Grassland") %>% 
  select(depth_increment, mean)
h_line <- setNames(as.numeric(h_line$mean), 
                   as.character(h_line$depth_increment))

h_line_col <- "#292929"

## Facet 1 - Permanent cropland
p <- plot_dat %>% 
  filter(facet == "Cropland") %>% 
  ggplot(aes(x = facet, y = mean, fill = depth_increment)) + 
  geom_col(position = position_dodge(width = +0.3), 
           stat = "identity",
           width = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = +0.3),
                width = 0.3,
                size = 0.5,
                color = h_line_col) + # Width of the error bars
  geom_hline(yintercept = h_line[["0-10cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  geom_hline(yintercept = h_line[["10-30cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  geom_hline(yintercept = h_line[["30-100cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  xlab(element_blank()) +
  ylab(expression(SOC~density~(Mg~ha^{-1}~cm^{-1}))) + #nolint
  ylim(-0.1, y_limit) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(1, 0.93),
        axis.text = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.ticks = element_line(), 
        text = element_text(family = "Calibri")) +
  scale_x_discrete(labels = label_wrap(10)) +
  scale_fill_manual(values = depth_colours)

# Generate "n" annotation for number of sites per category
f1_annotation <- plot_dat %>% 
  filter(facet == "Cropland") %>%
  mutate(x = facet,
         y = rep(-0.08, 1),
         n = max(n)) %>% 
  slice(1)

# Add "n" annotation to plot
p1 <- p + annotate("text",
                   x = f1_annotation$x,
                   y = f1_annotation$y,
                   label = f1_annotation$n,
                   size = 3,
                   color = "#6a6b6a") 


## Facet 2 - Transitional sites

x_scale <- c(0, 9, 25, 40, 84, 234)
x_labels <- c("0", "10", "26", "41", "85", "230")

p <- plot_dat %>%
  filter(facet == "Transitional") %>% 
  ggplot(aes(x = years_cont - 1.5, y = mean, fill = depth_increment)) + 
  geom_rect(aes(xmin = new_left, xmax = new_right, ymin = 0, ymax = mean, 
                fill = depth_increment)) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = +5.5),
                width = 5,
                size = 0.5,
                color = h_line_col) + # Width of the error bars
  geom_hline(yintercept = h_line[["0-10cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  geom_hline(yintercept = h_line[["10-30cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  geom_hline(yintercept = h_line[["30-100cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  xlab(element_blank()) +
  ylab(element_blank()) + 
  ylim(-0.1, y_limit) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        text = element_text(family = "Calibri")) +
  scale_fill_manual(values = depth_colours) +
  scale_x_continuous(breaks = x_scale, labels = x_labels)

# Generate "n" annotation for number of sites per category
f2_annotation <- plot_dat %>% 
  filter(facet == "Transitional") %>%
  group_by(years_cont) %>% 
  mutate(x = years_cont - 1.25,
         y = rep(-0.08, length(years_cont)),
         n = max(n)) %>% 
  slice(1) %>% 
  ungroup()

# Add "n" annotation to plot
p2 <- p + annotate("text",
                   x = f2_annotation$x,
                   y = f2_annotation$y,
                   label = f2_annotation$n,
                   size = 3,
                   color = "#6a6b6a")


## Facet 3 - Permanent grassland
p <- plot_dat %>% 
  filter(facet == "Grassland") %>% 
  ggplot(aes(x = facet, y = mean, fill = depth_increment)) + 
  geom_col(position = position_dodge(width = +0.3), 
           stat = "identity",
           width = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = +0.3),
                width = 0.3,
                size = 0.5,
                color = h_line_col) + # Width of the error bars
  geom_hline(yintercept = h_line[["0-10cm"]], 
             linetype = "dashed", 
             size = 0.25, 
             color = h_line_col) +
  geom_hline(yintercept = h_line[["10-30cm"]], 
             linetype = "dashed", 
             size = 0.25,  
             color = h_line_col) +
  geom_hline(yintercept = h_line[["30-100cm"]], 
             linetype = "dashed", 
             size = 0.25,  
             color = h_line_col) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  ylim(-0.1, y_limit) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        text = element_text(family = "Calibri")) +
  scale_x_discrete(labels = label_wrap(10)) +
  scale_fill_manual(values = depth_colours)

# Generate "n" annotation for number of sites per category
f3_annotation <- plot_dat %>% 
  filter(facet == "Grassland") %>%
  mutate(x = facet,
         y = rep(-0.08, 1),
         n = max(n)) %>% 
  slice(1)

# Add "n" annotation to plot
p3 <- p + annotate("text",
                   x = f3_annotation$x,
                   y = f3_annotation$y,
                   label = f3_annotation$n,
                   size = 3,
                   color = "#6a6b6a") 


# Combine plots

f1_width <- .11
f2_width <- 0.8
f3_width <- 0.09

p4 <- ggarrange(p1, p2, p3,
                ncol = 3, 
                nrow = 1,
                widths = c(f1_width, f2_width, f3_width))

annotate_figure(p4, 
                bottom = text_grob("Years since conversion to grassland (from cropland)",#nolint
                                   family = "Calibri"))


### Save figure for publication ------------------------------------------------

ggsave(filename = here("figures",
                       "fig_5_crop_grass_luc_intervals.png"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "white",
       dpi = 300)
