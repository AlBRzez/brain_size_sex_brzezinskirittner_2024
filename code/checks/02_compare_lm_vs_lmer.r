#!/usr/bin/env Rscript
library(dplyr)
library(broom)
library(readr)
library(here)
library(tidyr)
library(ggplot2)

print("loaded packages")
source(here("code", "analyses", "helpers.r"))

## Code for comparing the model's estimates using lme and lm 
## The goal is to asses is center can be used as a fixed effect 

output_fold <- "trajectories"

# Get data ---------------------------------------------------------------------
samp_path <- here("data", "samples")

dfs <- get_samples(samp_path, fs = T)
list2env(dfs, envir = .GlobalEnv)
dfs_n <- names(dfs)

## For fs data
fs_data <- read_csv(here("data", "fs_data", "clean_fs_data.csv"))
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))
print("got data")

## lme results
lme_coef <- 
  read_csv(here("outputs", output_fold, "fs_coeff.csv")) |> 
  mutate(batch = "mixed") 

## Define functions ------------------------------------------------------------
get_mod <- function(reg, dem_df) {
  df <- inner_join(dem_df, reg) 
  
  # Remove tiv and regioanl outliers
  Q_tiv <- quantile(df$tiv, na.rm = TRUE)
  IQR_tiv <- Q_tiv[4] - Q_tiv[2]
  lower_tiv <- (Q_tiv[2] - 1.5 * IQR_tiv) 
  upper_tiv <- (Q_tiv[4] + 1.5 * IQR_tiv)
  
  Q_voi <- quantile(df$voi, na.rm = TRUE)
  IQR_voi <- Q_voi[4] - Q_voi[2]
  lower_voi <- (Q_voi[2] - 1.5 * IQR_voi) 
  upper_voi <- (Q_voi[4] + 1.5 * IQR_voi)
  
  df <- 
    df |> 
    filter(tiv > lower_tiv & tiv < upper_tiv) |> 
    filter(voi > lower_voi & voi < upper_voi) |> 
    mutate(log_euler = log(abs(euler_num) + 1))
  
  # + scale(log_euler) + (1 | center) 
  mod <- lm(scale(voi) ~ 1 + sex * scale(age_months) + center , data = df)
  
  return(list(mod))
}

mods <- c("lin_voi")

regs_guide <- 
  fs_guide |> 
  filter(ignore == 0)

all_data <- function(df, reg_number, dfs, mods, dfs_n, guide = NULL) { #added the "guide" for our ukbb data
  
  if(!is.null(guide)) {
    r <- guide$final_region[guide$reg == reg_number]
    print(r)
    reg_data <- df[,c("subject_id", r)]
    # names(reg_data) <- c("subject_id", "voi") # For ukbb arrow
  } else {
    reg_data <- df[,c(1, reg_number + 1)] # For ukbb arrow
  }
  names(reg_data) <- c("subject_id", "voi") # For ukbb arrow
  
  coef_tbl <- tibble()
  goodness_tbl <- tibble()
  # reg_data <- df[,c(1, reg_number + 1)] # For ukbb arrow
  
  coef_tbl_s <- tibble()
  goodness_tbl_s <- tibble()
  for(samp in 1:length(dfs)) {
    models <- get_mod(reg_data, dfs[[samp]])
    for(i in 1:length(models)) {
      coef_tbl_s <- tidy(models[[i]]) |> 
        mutate(df = dfs_n[samp], mod = mods[i])
      goodness_tbl_s <- glance(models[[i]]) |> 
        mutate(df = dfs_n[samp], mod = mods[i])
      coef_tbl <- bind_rows(coef_tbl, coef_tbl_s)
      goodness_tbl <- bind_rows(goodness_tbl, goodness_tbl_s)
    }
  }
  
  coef_tbl$reg <- reg_number
  goodness_tbl$reg <- reg_number
  
  return(list(data.frame(coef_tbl), data.frame(goodness_tbl))) 
}

print("declared functions")

## Running lm ------------------------------------------------------------------
n <- 1
out <- list()

for(r in sort(unique(regs_guide$reg))) {
  out[[n]] <- all_data(fs_data, r, dfs, mods, dfs_n, regs_guide)
  n <- n + 1
}

output <- purrr::transpose(out) %>% 
  purrr::map(bind_rows)

df_coefficients <- output[[1]]
df_goodness <- output[[2]]

got <- data.frame(table(df_goodness$reg))


write_csv(df_coefficients, here("outputs", output_fold, "fs_coeff_centerfixed.csv"))
write_csv(df_goodness, here("outputs", output_fold, "fs_goodness_centerfixed.csv"))

print("saved data")


# Comparison--------------------------------------------------------------------

final_df <- bind_rows(
  lme_coef |> 
    filter(mod == "lin_voi_regular" & effect == "fixed") |> 
    mutate(mod = "lin_voi") |> 
    select(-c(effect, group)),
  
  df_coefficients |> 
    mutate(batch = "fixed") |> 
    filter(!grepl("imaging", term)) 
) |> 
  left_join(terms_ref) |> 
  left_join(samples_ref)

final_2 <- 
  final_df |> 
  select(-c(term, mod, std.error, statistic, p.value)) |> 
  pivot_wider(names_from = batch, values_from = estimate) 

cors <- 
  final_2 |> 
  group_by(df, term_c) |> 
  summarise(c = cor(mixed, fixed),
            .groups = "drop")

ggplot(final_2, aes(x = mixed, y = fixed)) +
  facet_grid(df~term_c) +
  geom_vline(xintercept = 0, linetype = 2, color = "#595959") +
  geom_hline(yintercept = 0, linetype = 2, color = "#595959") +
  geom_abline(color = "#ff4f4f") +
  geom_text(data = cors, aes(label = paste0("r = ", round(c, 3))), x = -.3, y = .8) +
  geom_point(alpha = .6) +
  ggtitle("estimates comparison") +
  theme_light() +
  theme(
    aspect.ratio = 1
  )

