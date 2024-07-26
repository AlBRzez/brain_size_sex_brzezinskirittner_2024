library(dplyr)
library(readr)
library(here)
library(tidyr)
library(ggplot2)
library(cowplot)

## Code for creating subsamples

source(here("code", "analyses", "helpers.r"))
output_fold <- "samples"
# Get source data --------------------------------------------------------------
full_ukbb_data <- read_csv(here("data", "others", "full_ukbb_data.csv"))
brain_sizes <- read_csv(here("data", "fs_data", "brain_sizes.csv"))
tiv_dbm <- read_csv(here("data", "others", "dbm_qc_tiv.csv"))

# to make sure we don't have any missing values in the data we want to analyze
fs_data <- read_csv(here("data", "fs_data", "clean_fs_data.csv"))
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))

keep_columns <- 
  fs_guide |> 
  filter(ignore == 0) |> 
  pull(final_region)

fs_data <-
  fs_data |>
  select(subject_id, all_of(keep_columns)) |> 
  na.omit() 

## ----
info_or <- 
  left_join(full_ukbb_data, brain_sizes) |> 
  na.omit() |> 
  filter(subject_id %in% fs_data$subject_id)

#### Sanity check and explanations  ------
# I'm filtering this big df a little bit more, here's why
# At the end of the code, there are the "sizes" plots showing the removal effect
ggplot(info_or, aes(x = mask)) +
  geom_histogram(bins = 1000)
ggplot(info_or, aes(x = csf_wb)) +
  geom_histogram(bins = 1000)

# Cleaning ---------------------------------------------------------------------
info <- 
  info_or |>
  filter(mask < 2500000) 

# Excluding subjects with huge incongruency between scale_factor and fs_tiv
mod <- lm(fs_tiv ~ scale_fact, data = tiv_dbm)
tiv_dbm$rst <- scale(mod$residuals)[,1]

tiv_dbm <- #there's a plot with this exclusion's sanity check bellow
  tiv_dbm |>
  mutate(keep = ifelse(abs(rst) > 3, 0, 1)) 

short_info <- 
  info |>
  filter(euler_num > -217) |>
  left_join(tiv_dbm) |>
  filter(keep == 1) |>
  mutate(tiv = scale_fact * 1886539) |> 
  select(subject_id, sex, age_years, age_months, center, fs_tiv, tiv, 
         tbv, brain_seg, brain_seg_not_vent, euler_num, scale_fact)

write_csv(short_info, here("data", output_fold, "full_sample.csv"))


#### Plotting func -----
comp_plot <- function(df, title) {
  p1 <- 
    df |>
    pivot_longer(c(age_months, tiv)) |>
    ggplot(aes(x = value, color = sex, linetype = sex)) +
    facet_wrap(~name, scales = "free") +
    # facet_grid(sex~name, scales = "free_x") +
    geom_histogram(alpha = .3, position = "identity", fill = NA) +
    
    labs(title = paste("d = ", title),
         subtitle = paste(nrow(df), "rows")) +
    scale_color_manual(values = color_sex,
                       aesthetics = c("color", "fill")) +
    theme_light() +
    theme(
      legend.position = "none"
    )
  
  p2 <- 
    df |> 
    transmute(
      age_month_dif = age_months - lag(age_months),
      tiv_dif = tiv - lag(tiv),
      sex = sex
    ) |> 
    filter(sex == "Female") |>
    pivot_longer(-sex) |>
    ggplot(aes(x = value)) +
    facet_wrap(~name, scales = "free", ncol = 1) +
    geom_histogram(alpha = .5) +
    theme_light()
  
  plot_grid(p1, p2, rel_widths = c(.7, .3), nrow = 1)
}

# Matching process -------------------------------------------------------------

get_sample <- function(distance, female, male, var, mean_age, mean_var) {
  
  ## Function to create samples matched by age and another selected variable
  
  matched <- data.frame()
  d <- distance
  male <- 
    male |> 
    select(subject_id, age_months, !!sym(var)) |> 
    rename("var" = !!sym(var))
  
  female <- 
    female |> 
    select(subject_id, age_months, !!sym(var)) |> 
    rename("var" = !!sym(var))
  
  for (i in 1:nrow(male)) {
    ind <- which((abs(male$age_months[i]-female$age_months) < (d * mean_age)) & (abs(male$var[i]-female$var) < (d * mean_var)))
    
    if (sum(ind) > 0) {
      if (length(ind) == 1) {
        ind_matched <- ind
      } else {
        ind_matched <- sample(ind, 1)
      }
      mat <- rbind(male[i, ], female[ind_matched, ])
      mat$round  <- i
      matched <- rbind(matched, mat)
      female <- female[-ind_matched, ]
    }
    rm(ind)
  }
  print(dim(matched))
  names(matched) <- c("subject_id", "age_months", var, "round")
  return(matched)  
}

## by scale factor tiv ----------
mean_age <- mean(short_info$age_months)
mean_tiv <- mean(short_info$tiv)

short_info <- short_info |> arrange(sex, subject_id)
set.seed(123)
female <- short_info[short_info$sex == "Female", ]
female <- female[sample(nrow(female)), ]
male <- short_info[short_info$sex == "Male", ]
male <- male[sample(nrow(male)), ]

msf1 <- get_sample(.002, female, male, "tiv", mean_age, mean_tiv) # dfs with seed 123 
msf2 <- get_sample(.002, female, male, "tiv", mean_age, mean_tiv) # dfs with seed 123
msf3 <- get_sample(.002, female, male, "tiv", mean_age, mean_tiv) # dfs with seed 123
msf4 <- get_sample(.002, female, male, "tiv", mean_age, mean_tiv) # dfs with seed 123

msf1 <- msf1 |> left_join(short_info)
msf2 <- msf2 |> left_join(short_info)
msf3 <- msf3 |> left_join(short_info)
msf4 <- msf4 |> left_join(short_info)

comp_plot(msf1, ".002 v.1")
comp_plot(msf2, ".002 v.2")
comp_plot(msf3, ".002 v.3")
comp_plot(msf4, ".002 v.4")

final_matched <- msf3

write_csv(final_matched, here("data", output_fold,  "scalefact_matched.csv"))

## by fs tiv ----------
mean_fstiv <- mean(short_info$fs_tiv)
set.seed(123)
female <- short_info[short_info$sex == "Female", ]
female <- female[sample(nrow(female)), ]
male <- short_info[short_info$sex == "Male", ]
male <- male[sample(nrow(male)), ]

m1 <- get_sample(.002, female, male, "fs_tiv", mean_age, mean_fstiv) # dfs with seed 123 
m2 <- get_sample(.002, female, male, "fs_tiv", mean_age, mean_fstiv) # dfs with seed 123
m3 <- get_sample(.002, female, male, "fs_tiv", mean_age, mean_fstiv) # dfs with seed 123
m4 <- get_sample(.002, female, male, "fs_tiv", mean_age, mean_fstiv) # dfs with seed 123

m1 <- m1 |> left_join(short_info)
m2 <- m2 |> left_join(short_info)
m3 <- m3 |> left_join(short_info)
m4 <- m4 |> left_join(short_info)

comp_plot(m1, ".002 v.1")
comp_plot(m2, ".002 v.2")
comp_plot(m3, ".002 v.3")
comp_plot(m4, ".002 v.4")

tiv <- m3
tiv_mat <- tiv[1:nrow(final_matched),]

# Saving after inspecting with plots and checking sample size
write_csv(tiv_mat, here("data", output_fold, "fstiv_matched.csv"))

## Age matched ------
get_age_sample <- function(distance, female, male, mean_age) {
  matc <- data.frame()
  d <- distance
  for (i in 1:nrow(male)) {
    
    ind <- which((abs(male$age_months[i]-female$age_months) < (d * mean_age)))
    
    if (sum(ind) > 0) {
      if (length(ind) == 1) {
        ind_matched <- ind
      } else {
        ind_matched <- sample(ind, 1)
      }
      mat <- rbind(male[i, ], female[ind_matched, ])
      mat$round  <- i
      matc <- rbind(matc, mat)
      female <- female[-ind_matched, ]
    }
    rm(ind)
  }
  print(dim(matc))
  return(matc)
}

age <- get_age_sample(.002, female, male, mean_age)
age_mat <- age[1:nrow(final_matched),]

age_mat <- age_mat |> select(all_of(names(final_matched)))
comp_plot(age_mat, ".002 | only age")

write_csv(age_mat, here("data", output_fold, "age_matched.csv"))

## Not matched sample -----
nams <- names(final_matched)
nams <- nams[-4]

set.seed(18)
random <- 
  short_info |>
  group_by(sex) |>
  sample_n(nrow(final_matched)/2) |>
  ungroup() |>
  select(all_of(nams))
comp_plot(random, "random")

write_csv(random, here("data", output_fold,  "random_sample.csv"))

## Extremes sample -----
set.seed(36)
extremes1 <- 
  short_info |>
  filter(!subject_id %in% final_matched$subject_id) 

female <- extremes1[extremes1$sex == "Female", ]
female <- female[sample(nrow(female)), ]
male <- extremes1[extremes1$sex == "Male", ]
male <- male[sample(nrow(male)), ]

extremes <- get_age_sample(.15, female, male, mean_age)
extremes_sample <- extremes[1:nrow(final_matched),]
comp_plot(extremes_sample, "extremes")

write_csv(extremes_sample, here("data", output_fold,  "extreme_sample.csv"))

# Sizes ----------
# Code for sanity checks  
sample_distribution <- function(samp, samp_name, sizes) {
  df <- 
    inner_join(samp |> 
                 select(subject_id, sex), sizes) |> 
    pivot_longer(-c(subject_id, sex))
  
  ggplot(df, aes(x = value, color = sex, fill = sex)) +
    facet_wrap(~name, scales = "free_x") + 
    geom_histogram(aes(linetype = sex),
                   alpha = 0,
                   position = "identity") +
    labs(title = paste("Sizes:", samp_name),
         
         y = "Participants") +
    scale_color_manual(values = color_sex,
                       aesthetics = c("color", "fill")) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_text(color = "#595959")
    ) +
    guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                    alpha = c(1,1))))
}

sample_distribution(m3, "matched", brain_sizes)



info |>filter(csf_wb > 1000000)
info |>filter(mask > 2500000)

sample_distribution(info, "filtered", brain_sizes)
sample_distribution(info_or, "complete", brain_sizes)

## Scaling factor sanity check -----
ggplot(tiv_dbm, aes(x = fs_tiv, y = scale_fact)) +
  geom_point(alpha = .5, aes(color = factor(keep))) +
  geom_smooth(method = "lm") +
  scale_color_manual(values = c( "#ff4f4f", "black")) +
  theme_light() +
  theme(
    legend.position = "none"
  )

ggplot(info, aes(x = euler_num)) +
  geom_vline(xintercept = -217, color = "#ff4f4f") +
  geom_vline(xintercept = 0, color = "steelblue") +
  geom_bar(fill = "black", color = "black") +
  theme_light()
range(info$euler_num)
