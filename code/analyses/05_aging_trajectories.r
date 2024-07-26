#!/usr/bin/env Rscript
library(dplyr)
library(broom.mixed)
library(readr)
library(parallel)
library(foreach)
library(here)
library(tidyr)
library(lme4)
library(lmerTest)

print("loaded packages")
source(here("code", "analyses", "helpers.r"))

output_fold <- "trajectories"
if(!dir.exists(here("outputs", output_fold))) {
  dir.create(here("outputs", output_fold))
}
# Get data ---------------------------------------------------------------------
samp_path <- here("data", "samples")

dfs <- get_samples(samp_path, fs = T)
list2env(dfs, envir = .GlobalEnv)
dfs_n <- names(dfs)

## For fs data
fs_data <- read_csv(here("data", "fs_data", "clean_fs_data.csv"))
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))
print("got data")

## Define functions ------------------------------------------------------------
get_mod <- function(reg, dem_df) {
  # Function to run different models for one specific region and one dataset
  
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
  
  models_list <- vector(mode = "list", length = 4)
  
  models_list[[1]] <- lmer(scale(voi) ~ 1 + sex * scale(age_months) + (1 | center) , data = df)
  models_list[[2]] <- lmer(scale(voi) ~ 1 + sex * (scale(age_months) + scale(I(age_months^2))) + (1 | center) , data = df)
  models_list[[3]] <- lmer(scale(voi) ~ 1 + sex * scale(age_months) + scale(log_euler) + (1 | center)  , data = df)
  models_list[[4]] <- lmer(scale(voi) ~ 1 + sex * (scale(age_months) + scale(I(age_months^2))) + (1 | center) , data = df)
  
  return(models_list)
}

mods <- c("lin_voi_regular", "sq_voi_regular", "lin_voi_euler", "sq_voi_euler")

regs_guide <- 
  fs_guide |> 
  filter(ignore == 0)

all_data <- function(df, reg_number, dfs, mods, dfs_n, guide = NULL) { # "guide" for our ukbb data

  # Function used to get regional data and iterate the models for each region
  
  if(!is.null(guide)) {
    r <- guide$final_region[guide$reg == reg_number]
    print(r)
    reg_data <- df[,c("subject_id", r)]
  } else {
    reg_data <- df[,c(1, reg_number + 1)] 
  }
  names(reg_data) <- c("subject_id", "voi") 
  
  coef_tbl <- tibble()
  goodness_tbl <- tibble()

  # gets the clean coefficients for all the models
  
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

# Test
# tic <- Sys.time()
# tst <- all_data(df = fs_data, reg_number = 208, dfs = dfs, mods = mods, dfs_n = dfs_n, guide = regs_guide)
# toc <- Sys.time()
# toc-tic
## Parallelization cluster -----------------------------------------------------
detectCores()
my_cluster <- makeForkCluster(
  nnodes = 4
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my_cluster)

## Iteration process -----------------------------------------------------------

combine <- function(x, ...) {  
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}


print("running models")
tic <- Sys.time()
output <- foreach(
  x = sort(unique(regs_guide$reg)),
  .errorhandling = "remove",
  .combine = "combine",
  .multicombine = TRUE,
  .verbose = TRUE,
  .inorder = TRUE
) %dopar% {
  all_data(fs_data, x, dfs, mods, dfs_n, regs_guide) 
  # all_data(fs_data, x, dfs, mods, dfs_n) 
}
toc <- Sys.time()

print(toc - tic)
print("got the models")
#models <- output[[3]]
#save(models, file = "mods_test.Rdata")

parallel::stopCluster(cl = my_cluster)


df_coefficients <- output[[1]]
df_goodness <- output[[2]]

got <- data.frame(table(df_goodness$reg))


write_csv(df_coefficients, here("outputs", output_fold, "fs_coeff.csv"))
write_csv(df_goodness, here("outputs", output_fold, "fs_goodness.csv"))

print("saved data")
