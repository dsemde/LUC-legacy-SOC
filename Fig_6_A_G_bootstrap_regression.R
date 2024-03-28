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

# Parallel processing libraries
library(doParallel)
library(foreach)

# Regression libraries
library(minpack.lm)
library(nlstools)

# Weighted bootstrap regression tutorial
library(nls.multstart)
library(patchwork)
library(rTPC)
library(car)


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

dat_w <- foreach(i = 1:length(depths), .combine = rbind, .packages = c("twang")) %dopar% { #nolint
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

regression_dat <- dat_w %>% 
  group_by(years_cont, depth_increment) %>%
  summarise(n = n_distinct(PointID),
            mean = weighted.mean(TOC_stock_cm, w = w),
            sd = (weighted.sd(TOC_stock_cm, w = w)^2 / n),
            x = case_when(years_cont == 2007 ~ 400,
                          years_cont == 2018 ~ 0,
                          TRUE ~ years_cont)) %>%
  slice(1) %>% 
  ungroup() %>% 
  filter(depth_increment != "30-100cm") %>%
  group_by(years_cont) %>%
  mutate(temp_mean = (mean[depth_increment == "0-10cm"] + 
                        mean[depth_increment == "10-30cm"] * 2) / 3,
         sd = (sd[depth_increment == "0-10cm"] + 
                 sd[depth_increment == "10-30cm"] * 2) / 3) %>%
  ungroup() %>% 
  summarise(y = ((temp_mean - temp_mean[years_cont == 2018]) / 
                   temp_mean[years_cont == 2018]) * 100,
            x = x,
            sd = sd,
            n = n) %>% 
  group_by(x) %>%
  slice(1) %>% 
  ungroup() %>% 
  as.data.frame()


### Exponential Regression -----------------------------------------------------

c_start <- max(regression_dat$y)
k_start <- 2 * log(2) / c_start

fit_nlsLM <- minpack.lm::nlsLM(y ~ c * (1 - exp(-k * x)), #nolint
                               data = regression_dat,
                               start = list(c = c_start, k = k_start),
                               weights = 1 / sd)

coef <- coef(fit_nlsLM)

### Residual resampling bootstrap non-linear regression ------------------------

boot <- Boot(fit_nlsLM, method = "residual", R = 1000)

# Predict over new data
boot_preds <- boot$t %>% 
  as.data.frame() %>% 
  drop_na() %>% 
  mutate(iter = 1:n()) %>% #nolint
  group_by_all() %>% 
  do(data.frame(x = seq(min(regression_dat$x), 
                        max(regression_dat$x), 
                        length.out = 100))) %>% 
  ungroup() %>% 
  mutate(pred = predict(fit_nlsLM, newdata = ., type = "response"))

# Calculate bootstrapped confidence intervals
boot_conf_preds <- group_by(boot_preds, x) %>% 
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975),
            .groups = "drop")

ci_c_start <- min(regression_dat$y)
ci_k_start <- 2 * abs(log(2)) / abs(c_start)

fit_ci_lower <- minpack.lm::nlsLM(conf_lower ~ c * (1 - exp(-k * x)),
                                  data = boot_conf_preds,
                                  start = list(c = c_start, k = k_start))

fit_ci_upper <- minpack.lm::nlsLM(conf_upper ~ c * (1 - exp(-k * x)),
                                  data = boot_conf_preds,
                                  start = list(c = c_start, k = k_start))

ci_lower_coef <- coef(fit_ci_lower)
ci_upper_coef <- coef(fit_ci_upper)


### Locate asymptote and calculate CI cross ------------------------------------

# Return y for a given x
get_y <- function(x, coefs) {
  coefs <- unname(coefs)

  # y <- coefs[1] - coefs[2] * exp(-coefs[3] * x)
  y <- coefs[1] * (1 - exp(-coefs[2] * x))

  return(y)
}

# See how much change has occured after 20 years
get_y(20, coef)
get_y(20, ci_lower_coef)
get_y(20, ci_upper_coef)

# Get value of x when fit line is within 0.1 % yr-1 of the asymptote
# Poeplau et al. 2011
get_asymp_date <- function(coefs) {
  coefs <- unname(coefs)

  # Get asymptote approximation
  asymp_cutoff <- abs(coefs[1] * 0.001)

  # Find first instance when the line is within 0.01 % of the asymptote
  for (i in 1:max(regression_dat$x)) {
    asymp_date <- 
      ifelse((abs(get_y((i + 1), coefs)) - 
                abs(get_y(i, coefs))) <= asymp_cutoff, i + 1, 0)
    if (asymp_date != 0) break
  }

  return(asymp_date)
}

steadystate <- data.frame(x = get_asymp_date(coef), 
                          y = get_y(get_asymp_date(coef), coef))

x_ci <- data.frame(x = steadystate["x"], 
                   y = steadystate["y"], 
                   x_lower = unname(get_asymp_date(ci_lower_coef)), 
                   x_upper = unname(get_asymp_date(ci_upper_coef)))
y_ci <- data.frame(x = steadystate["x"], 
                   y = steadystate["y"],
                   y_upper = unname(get_y(x_ci["x_lower"], ci_lower_coef)), 
                   y_lower = unname(get_y(x_ci["x_upper"], ci_upper_coef)))



### Plotting -------------------------------------------------------------------

# Generate fit line
fit_line <- 
  stat_function(fun = function(x) {
    coef[1] * (1 - exp(-coef[2] * x))
  }, 
  colour = "#045275",
  alpha = 0.75)

ci_lower_fit <- 
  stat_function(fun = 
                function(x) {
                  ci_lower_coef[1] * (1 - exp(-ci_lower_coef[2] * x))
                }, 
                colour = "#662D3A",
                alpha = 0.5)

ci_upper_fit <- 
  stat_function(fun = 
                function(x) {
                  ci_upper_coef[1] * (1 - exp(-ci_upper_coef[2] * x))
                }, 
                colour = "#662D3A",
                alpha = 0.5)

# plot bootstrapped predictions
ggplot() +
  geom_line(aes(x, pred, group = iter), 
            boot_preds, 
            col = "#ffd079", 
            alpha = 0.05) +
  fit_line +
  ci_lower_fit + 
  ci_upper_fit +
  geom_point(data = steadystate, 
             aes(x, y),
             colour =  "#045275",
             size = 5, 
             shape = 18) +
  # geom_errorbar(data = y_ci, 
  #               aes(x = x, y = y, ymin = y_lower, ymax = y_upper),
  #               color = "#024B7A",
  #               width = 3) +
  geom_errorbar(data = x_ci, 
                aes(x = x, y = y, xmin = x_lower, xmax = x_upper),
                color = "#045275",
                width = 0.65) +
  geom_point(aes(x, y, size = log(n)), regression_dat, colour = "#089099") +
  annotate("text",
           x = 300, y = 5, label = "y == 47.26 (1 - e^{-0.04730*t})",
           parse = TRUE,
           size = 5,
           colour = "#4b4b4b",
           family = "Calibri") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        axis.ticks = element_line(), 
        text = element_text(family = "Calibri")) +
  labs(x = "Years since conversion to grassland (from cropland)",
       y = "Relative change in SOC stock (%)")


### Save figure for publication ------------------------------------------------

ggsave(filename = here("figures",
                       "fig_6_crop_grass_soc_regression.png"),
       width = 27,
       height = 18,
       units = "cm",
       bg = "transparent",
       dpi = 300)
