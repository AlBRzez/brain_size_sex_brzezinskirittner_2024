library(readr)
library(tidyr)
library(dplyr)
library(glue)
library(here)
library(purrr)

source(here("code", "analyses", "helpers.r"))

cerebra <- read_csv(here("code", "brain_plotting", "cerebra_plotting", "newLabels.csv"))
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))

# Load data --------------------------------------------------------------------
fs_coef_o <-
  read_csv(here("outputs", "trajectories", "fs_coeff.csv")) 

fs_allom_o <-
  read_csv(here("outputs", "allometry", "lm_over_betas_coef_tbv.csv")) 

# Data wrangling ---------------------------------------------------------------
fs_coef <- 
  fs_coef_o |>
  left_join(terms_ref) |> 
  filter(effect == "fixed" & df != "fs_matc") |> 
  group_by(df, mod, term_c) |> 
  nest() |>
  mutate(
    p_fdrc = map(data, function(df)  p.adjust(df$p.value, method = "fdr")),
  ) |>
  unnest(cols = c(data, p_fdrc)) |>
  ungroup() |> 
  mutate(
    p_corr_fdrc = ifelse(p_fdrc <= 0.05, 1, 0),
  ) 

fs_allom <- 
  fs_allom_o |>
  left_join(terms_ref) |> 
  filter(df != "fs_matc") |> 
  group_by(df, mod, term_c) |> 
  nest() |>
  mutate(
    p_fdrc = map(data, function(df)  p.adjust(df$p.value, method = "fdr")),
  ) |>
  unnest(cols = c(data, p_fdrc)) |>
  ungroup() |> 
  mutate(
    p_corr_fdrc = ifelse(p_fdrc <= 0.05, 1, 0),
  ) 

# Mapping dkt to cerebra -------------------------------------------------------
cerebra_r <- 
  cerebra |> 
  mutate(side = "right",
         name_use = gsub("right_", "", Name))
cerebra_l <- cerebra_r |> 
  mutate(side = "left",
         newlabels = newlabels + 51)

cerebra_use <- 
  bind_rows(cerebra_r, cerebra_l) |> 
  mutate(name_use = tolower(name_use),
         name_use = gsub("_", "", name_use))

vols <-
  fs_guide |> 
  filter(type == "volume") |> 
  transmute(
    fs_reg = reg,
    ignore = ignore,
    col_names = col_names,
    side = case_when(
      grepl("right", col_names) ~ "right",
      grepl("left", col_names) ~ "left",
      grepl("brain", col_names) ~ "brain"),
    name_use = gsub("volume_of_", "", col_names),
    name_use = gsub("_left|_right|_whole_brain", "", name_use),
    name_use = tolower(name_use),
    name_use = gsub("_", "", name_use),
  )

cerebra_mapping <- 
  left_join(vols, cerebra_use) |> 
  mutate(final_name = glue("{name_use}_{side}")) |> 
  filter(ignore == 0)
write_csv(cerebra_mapping, here("data", "others", "fs_cerebra_mapping.csv"))

# Arranging data for brainplots ------------------------------------------------
## Aging trajectories data ----
for_plotting_long <-
  fs_coef |>
  inner_join(cerebra_mapping, by = c("reg" = "fs_reg")) |>
  select(term_c, name_use, side, final_name, df, mod, estimate, newlabels)

comb <- expand.grid(
  df = unique(for_plotting_long$df),
  mod = unique(for_plotting_long$mod),
  term_c = unique(for_plotting_long$term_c),
  newlabels = 1:102
)

p_long <- 
  left_join(comb, for_plotting_long) |> 
  replace_na(list(estimate = 0)) 

write_csv(p_long, here("outputs", "plots", "brainplots", "brainplots_data", "trajectories_reg.csv"))
rm(for_plotting_long, p_long)

for_plotting_long <- #just significant regions
  fs_coef |>
  mutate(estimate = ifelse(p_corr_fdrc == 1, estimate, 0)) |>
  inner_join(cerebra_mapping, by = c("reg" = "fs_reg")) |>
  select(term_c, name_use, side, final_name, df, mod, estimate, newlabels)

p_long <- 
  left_join(comb, for_plotting_long) |> 
  replace_na(list(estimate = 0)) 

write_csv(p_long, here("outputs", "plots", "brainplots", "brainplots_data", "trajectories_thresh.csv"))
rm(for_plotting_long, comb, p_long)

## Data for lm over allometry --------
allom_male_est <- 
  fs_allom |>
  filter(mod == "lin_int_regular" & window_size == 60) |>
  select(df, reg, window_size, mod, term_c, estimate) |>
  pivot_wider(names_from = term_c, values_from = estimate) |>
  mutate(male_int = intercept + sex_male,
         male_slope = agem + sex_male_agem) |>
  pivot_longer(-c(df, reg, window_size, mod),
               names_to = "term_c", values_to = "estimate") |>
  inner_join(cerebra_mapping, by = c("reg" = "fs_reg"))


for_plotting_long <- # version to calculate the estimates for male
  allom_male_est |> 
  select(term_c, name_use, side, final_name, df, mod, estimate, newlabels)


# for_plotting_long <-
#   fs_allom |>
#   inner_join(cerebra_mapping, by = c("reg" = "fs_reg")) |>
#   select(term_c, name_use, side, final_name, df, mod, estimate, newlabels, window_size)


comb <- expand.grid(
  df = unique(for_plotting_long$df),
  mod = unique(for_plotting_long$mod),
  term_c = unique(for_plotting_long$term_c),
  # window_size = unique(for_plotting_long$window_size),
  newlabels = 1:102
)

p_long <- 
  left_join(comb, for_plotting_long) |> 
  replace_na(list(estimate = 0)) 

write_csv(p_long, here("outputs", "plots", "brainplots", "brainplots_data", "allom_60m.csv"))
rm(for_plotting_long, p_long)


for_plotting_long <- #just significant regions
  fs_allom |>
  filter(mod == "lin_int_regular" & window_size == 60) |> 
  mutate(estimate = ifelse(p_corr_fdrc == 1, estimate, 0)) |>
  inner_join(cerebra_mapping, by = c("reg" = "fs_reg")) |>
  select(term_c, name_use, side, final_name, df, mod, estimate, newlabels)

p_long <- 
  left_join(comb, for_plotting_long) |> 
  replace_na(list(estimate = 0)) 

write_csv(p_long, here("outputs", "plots", "brainplots", "brainplots_data", "allom_60m_thresh.csv"))
rm(for_plotting_long, p_long)

######
cerebra_use <- cerebra_use |> mutate(final_name = glue("{name_use}_{side}"))
write_csv(cerebra_use, here("code", "brain_plotting", "cerebra_plotting", "cerebra_reference.csv"))

