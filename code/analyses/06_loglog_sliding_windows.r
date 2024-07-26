#!/usr/bin/env Rscript
library(dplyr)
library(readr)
library(here)
library(tidyr)
library(ggplot2)

source(here("code", "analyses", "helpers.r"))
print("loaded packages")

output_fold <- "allometry"
if(!dir.exists(here("our_data", output_fold))) {
  dir.create(here("our_data", output_fold))
}
# Get data ---------------------------------------------------------------------
brain_sizes <- read_csv(here("data", "fs_data", "brain_sizes.csv"))

samp_path <- here("data", "samples")

dfs <- get_samples(samp_path, fs = T)
list2env(dfs, envir = .GlobalEnv)
dfs_n <- names(dfs)

## For fs data
fs_data <- read_csv(here("data", "fs_data", "clean_fs_data.csv"))
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))
print("got data")

keep_columns <- 
  fs_guide |> 
  filter(ignore == 0) |> 
  pull(final_region)

fs_data <-
  fs_data |>
  select(subject_id, all_of(keep_columns))

print("Got data")
#### ----------
# Make all the fs log-10
log_fs_data <- apply(fs_data, MARGIN = 2, FUN = log10)
log_fs_data <- as_tibble(log_fs_data)
log_fs_data$subject_id <- fs_data$subject_id

# Function to merge the log data with the demographics and filter
log_df_vol <- function(s_df, fs_df, outliers = NULL) {
  s_df$log_tbv <- log10(s_df$gm_cerebrum)
  s_df$log_tiv <- log10(s_df$tiv)
  out <- 
    inner_join(
      s_df, 
      fs_df
    ) 
  
  if(!is.null(outliers)){
    out <- 
      out |> 
      filter(!subject_id %in% outliers)
  }
  return(out)
}

new_tbv <- 
  brain_sizes |> 
  select(subject_id, gm_cerebrum)


full_log <- log_df_vol(full |> inner_join(new_tbv), log_fs_data)
matched_log <- log_df_vol(matched |> inner_join(new_tbv), log_fs_data)
extreme_log <- log_df_vol(extreme |> inner_join(new_tbv), log_fs_data)
random_log <- log_df_vol(random |> inner_join(new_tbv), log_fs_data)
age_mat_log <- log_df_vol(age_mat |> inner_join(new_tbv), log_fs_data)
fs_matc_log <- log_df_vol(fs_matc |> inner_join(new_tbv), log_fs_data)


# identify regions to run loglog
voi <- names(log_fs_data)
voi <- voi[-1]

get_est_mat <- function(df, f_sex, voi, win_size){#, m_names) {
  # last argument represents the total number of months
  # divided /2 tells how many months before and how many months after
  win <- win_size/2
  m_names <- vector()
  
  df_size <- tibble(
    months = numeric(),
    orig_n = numeric(),
    final_n_1d = numeric(),
    final_n_2d = numeric(),
  )
  
  j <- 1
  mat_c <- matrix(1, nrow = length(voi))
  
  sex_df <- 
    df |> 
    filter(sex == f_sex & (age_years > 49 & age_years < 76))
  for(month in ((min(sex_df$age_months) + win):(max(sex_df$age_months) - win))) {
    df_size[j, 1] <- month
    age_dat <- 
      sex_df |> 
      filter(age_months >= month - win & age_months <= month + win)
    df_size[j, 2] <- nrow(age_dat)
    
    
    # Outliers
    Q <- quantile(age_dat$tbv)
    IQR <- Q[4] - Q[2]
    lower <- (Q[2] - 1.5 * IQR) 
    upper <- (Q[4] + 1.5 * IQR)
    
    age_dat <- 
      age_dat |> 
      filter(tbv > lower & tbv < upper)
    
    df_size[j, 3] <- nrow(age_dat)
    
    v <- vector() 
    for (i in voi) {
      Qr <- quantile(age_dat |> pull(!!sym(i)))
      IQRr <- Qr[4] - Qr[2]
      lowerr <- (Qr[2] - 1.5 * IQRr) 
      upperr <- (Qr[4] + 1.5 * IQRr)
      
      reg_dat <- 
        age_dat |> 
        filter(!!sym(i) > lowerr & !!sym(i) < upperr)
      
      df_size[j, 4] <- nrow(reg_dat)
      
      mod <- lm(get(i) ~ log_tbv, data = reg_dat)
      v <- append(v, mod$coefficients[2])
      rm(mod)
    }
    mat_c <- cbind(mat_c, v)
    m_names <- c(m_names, paste0("month_", month))
    j <- j + 1
  }
  
  out_mat <- mat_c[,-1]
  out_mat <- as_tibble(out_mat)
  names(out_mat) <- m_names
  return(list(out_mat, df_size))
}

print("Declared functions")
print("Starting the itteration")

for(m in c(20, 50, 60, 76, 100)) {
# for(m in 60 ) {
  print(paste("M =", m))
  tic <- Sys.time()
  mat_full_fem <- get_est_mat(full_log, "Female", voi, m)
  mat_full_mal <- get_est_mat(full_log, "Male", voi, m)
  mat_matched_fem <- get_est_mat(matched_log, "Female", voi, m)
  mat_matched_mal <- get_est_mat(matched_log, "Male", voi, m)
  mat_extreme_fem <- get_est_mat(extreme_log, "Female", voi, m)
  mat_extreme_mal <- get_est_mat(extreme_log, "Male", voi, m)
  mat_agemat_fem <- get_est_mat(age_mat_log, "Female", voi, m)
  mat_agemat_mal <- get_est_mat(age_mat_log, "Male", voi, m)
  mat_random_fem <- get_est_mat(random_log, "Female", voi, m)
  mat_random_mal <- get_est_mat(random_log, "Male", voi, m)
  mat_fsmat_fem <- get_est_mat(fs_matc_log, "Female", voi, m)
  mat_fsmat_mal <- get_est_mat(fs_matc_log, "Male", voi, m)
  toc <- Sys.time()
  print(toc - tic)
  
  #### Saving --------------------------------------------------------------------
  save(mat_full_fem, mat_full_mal, mat_matched_fem, mat_matched_mal,
       mat_extreme_fem, mat_extreme_mal, mat_agemat_fem, mat_agemat_mal,
       mat_random_fem, mat_random_mal, mat_fsmat_fem, mat_fsmat_mal,
       file = here("outputs", output_fold, paste0("tbv_month_outliersrem_2d_",m,"m.RData")))
  
  rm(mat_full_fem, mat_full_mal, mat_matched_fem, mat_matched_mal,
     mat_extreme_fem, mat_extreme_mal, mat_agemat_fem, mat_agemat_mal,
     mat_random_fem, mat_random_mal, mat_fsmat_fem, mat_fsmat_mal)
}

voi_order <- tibble(voi = voi,
                    order = 1:length(voi))
write_csv(voi_order, here("outputs", output_fold, "voi_order.csv"))
