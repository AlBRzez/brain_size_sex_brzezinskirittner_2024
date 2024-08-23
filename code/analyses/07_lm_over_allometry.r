library(dplyr)
library(readr)
library(here)
library(tidyr)
library(broom)
library(rstatix)

source(here("code", "analyses", "helpers.r"))

output_folder <- "allometry"

# Get guide --------
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))
vol_guide <- 
  fs_guide |> 
  filter(ignore == 0) 

# voi_order <- read_csv(here("outputs", output_folder, "voi_order.csv"))
# 
# voi_order <- 
#   voi_order |> 
#   rename("final_region" = "voi") |> 
#   left_join(vol_guide)
# 
# reg_nums <- vol_guide$reg

voi_order <- vol_guide |> select(final_region, reg)
reg_nums <- voi_order$reg

# Functions -------
join_data <- function(female_obj, male_obj, reg_guide) {
  # Function for joining female and male data for a window
  fem <- 
    female_obj[[1]] |>
    mutate(
      sex = "Female",
      # reg = nums
    ) |> 
    left_join(reg_guide)
  
  mal <- 
    male_obj[[1]] |>
    mutate(
      sex = "Male",
      # reg = nums
    ) |> 
    left_join(reg_guide)
  
  out <- bind_rows(fem, mal) |> select(-final_region)
  return(out)
}

region_models <- function(allom_data, reg_num) {
  
  # Models for each window
  
  tmp <- 
    allom_data |>
    filter(reg == reg_num) |>
    select(-c(reg, df)) |>
    pivot_longer(-sex) |>
    mutate(month = as.numeric(gsub("month_", "", name)))
  
  ml <- lm(value ~ 1 + sex * scale(month), data = tmp)
  mq <- lm(value ~ 1 + sex * (scale(month) + scale(I(month^2))), data = tmp)
  
  coef <- 
    bind_rows(
      tidy(ml) |> 
        mutate(reg = reg_num, mod = "lin_int_regular"),
      tidy(mq) |> 
        mutate(reg = reg_num, mod = "sq_int_regular")
    )
  good <- 
    bind_rows(
      glance(ml) |> 
        mutate(reg = reg_num, mod = "lin_int_regular"),
      glance(mq) |> 
        mutate(reg = reg_num, mod = "sq_int_regular")
    )
  
  cohensd <- cohens_d(value ~ sex, data = tmp)
  return(list(coef, good, cohensd))
}

# Run models 
tic <- Sys.time()
complete_coefs <- tibble()
complete_good <- tibble()
complete_d <- tibble()
dfs_n <- dfs_n_c[1:5]

for(ws in c(20, 50, 60, 76, 100)) {
# for(ws in 60) {
  rdata_path <- paste0("outputs/", output_folder, "/tbv_month_outliersrem_2d_", ws,"m.RData")
  load(rdata_path)
  
  dfs <- vector(mode = "list", length = 5)
  dfs[[1]] <- join_data(mat_full_fem, mat_full_mal, voi_order) |> mutate(df = dfs_n_c[1])
  dfs[[5]] <- join_data(mat_matched_fem, mat_matched_mal, voi_order)|> mutate(df = dfs_n_c[2])
  dfs[[3]] <- join_data(mat_agemat_fem, mat_agemat_mal, voi_order)|> mutate(df = dfs_n_c[3])
  dfs[[4]] <- join_data(mat_random_fem, mat_random_mal, voi_order)|> mutate(df = dfs_n_c[4])
  dfs[[2]] <- join_data(mat_extreme_fem, mat_extreme_mal, voi_order)|> mutate(df = dfs_n_c[5])
  allometry_data <- bind_rows(dfs) |> mutate(window_size = ws)
  write_csv(allometry_data, paste0("outputs/", output_folder, "/complete_allometry_tbv_", ws, "months.csv"))
  
  coefficients <- tibble()
  goodness <- tibble()
  coh_d <- tibble()
  
  for(samp in 1:length(dfs)) {
    for(r in reg_nums) {
      reg <- region_models(dfs[[samp]], r)
      coefficients <- bind_rows(coefficients, 
                                reg[[1]] |> 
                                  mutate(df = dfs_n[[samp]]))
      goodness <- bind_rows(goodness, 
                            reg[[2]]|> 
                              mutate(df = dfs_n[[samp]]))    
      coh_d <- bind_rows(coh_d, 
                         reg[[3]]|> 
                           mutate(df = dfs_n[[samp]]))    
    }
  }
  
  complete_coefs <- bind_rows(complete_coefs, 
                              coefficients |> mutate(window_size = ws))
  complete_good <- bind_rows(complete_good, 
                             goodness |> mutate(window_size = ws))
  complete_d <- bind_rows(complete_d, 
                          coh_d |> mutate(window_size = ws))
}
toc <- Sys.time()

write_csv(complete_coefs, here("outputs", output_folder, "lm_over_betas_coef_tbv.csv"))
write_csv(complete_good, here("outputs", output_folder, "lm_over_betas_good_tbv.csv"))
write_csv(complete_d, here("outputs", output_folder, "lm_over_betas_cohens_tbv.csv"))

print(toc - tic)
