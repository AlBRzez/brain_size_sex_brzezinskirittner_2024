source("../ukbb_data_r/functions.R")
library(ggplot2)
library(here)
library(readr)
library(tidyr)
library(dplyr)
library(scales)
library(rstatix)
library(rvest)
library(glue)
library(purrr)
library(stringr)

source(here("code", "analyses", "helpers.r"))

## Code for modeling fluid intelligence scores

## Get data --------------------------------------------------------------------
## Get ukbb data
ukbb_path <- "/data/zeiyas/in_vivo_Datasets/UKbiobank/brzali"
ukbb <- open_dataset(paste0(ukbb_path, "/UKBB-tabular-processing/current_melt.arrow"), format = "ipc")
codings <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Codings.tsv"))
diction <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Data_Dictionary_Showcase.tsv"))

## Get subsamples
samp_path <- here("data", "samples")

dfs <- get_samples(samp_path, fs = F)
list2env(dfs, envir = .GlobalEnv)
dfs_n <- names(dfs)


## For fs data
fs_data <- read_csv(here("data", "fs_data", "clean_fs_data.csv"))
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))


cerebra_mapping <- read_csv(here("data", "others", "fs_cerebra_mapping.csv"))
# Fluid intelligence data extraction -------------------------------------------
fi_long <- get_ukbb_data(arrow = ukbb, diction = diction, 
                         fields_id = c(20016), 
                         subjs_id = full$subject_id,
                         instances = 2)

fi <- order_data(fi_long, diction, codings)
rm(fi_long)
gc()

regions <- fs_guide |> 
  filter(ignore == 0 & type == "volume") |> 
  pull(final_region)

volumes <- 
  fs_data |> 
  select(subject_id, all_of(regions))

final_data <- 
  inner_join(full, volumes) |> 
  inner_join(fi |> mutate(value = as.numeric(value)) |> 
               select(subject_id, fl_i_score = value))


# Estimates: fi ~ voi * sex + age ----------------------------------------------
estimates <- tibble()

for(df_u in c(2, 4)) {
  sm_df <- final_data |> filter(subject_id %in% dfs[[df_u]]$subject_id)
  sm_df$fl_i_score_z <- scale(sm_df$fl_i_score)
  
  dfname <- dfs_n[df_u]
  for(r in regions) {
    tmp_df <- sm_df |> select(fl_i_score_z, sex, voi = !!sym(r), age_months)
    mod_r <- lm(fl_i_score_z ~ 1 + sex * scale(voi) + scale(age_months), data = tmp_df)
    est <- broom::tidy(mod_r) |> mutate(final_region = r,
                                        df = dfname,
                                        mod = "fluid_int") 
    estimates <- bind_rows(estimates, est)
  }  
}

estimates_2 <- 
  estimates |> 
  mutate(term_c = case_when(
    term == "(Intercept)" ~ "intercept",
    term == "sexMale" ~ "sex_male",
    term == "scale(voi)" ~ "voi",
    term == "sexMale:scale(voi)" ~ "sex_male_voi",
    term == "scale(age_months)" ~ "agem"
  )) |> 
  left_join(fs_guide |> select(final_region, reg))

# Exploring estimates: plots and statistical testing
estimates_2 |>
  filter(term_c %in% c("voi", "sex_male_voi") & df == "matched") |>
  select(term_c, df, estimate, reg) |>
  pivot_wider(names_from = term_c, values_from = "estimate") |>
  ggplot(aes(x = voi, y = sex_male_voi)) +
  geom_point() +
  geom_abline(color = "#ff4f4f")  +
  theme_light()

ggplot(estimates_2, aes(x = estimate)) +
  facet_grid(df ~ term_c) +
  geom_vline(xintercept = 0, color = "#ff4f4f", linetype = 2) +
  geom_histogram(bins = 50) +
  theme_light()

ggplot(estimates_2, aes(y = estimate, x = df, group = df, color = df)) +
  facet_wrap(~ term_c) +
  geom_hline(yintercept = 0, color = "#595959", linetype = 2) +
  geom_boxplot() +
  scale_color_manual(values = c("#ee9b00", "#ae2012")) +
  theme_light()

# fdr correction
fdr_est <- 
  estimates_2 |>
  select(df, reg, term_c, mod,  final_region, estimate, p.value) |>
  group_by(df, term_c) |> 
  nest() |>
  mutate(
    p_fdrc = map(data, function(df)  p.adjust(df$p.value, method = "fdr")),
  ) |>
  unnest(cols = c(data, p_fdrc)) |>
  ungroup() |> 
  mutate(
    p_corr_fdrc = ifelse(p_fdrc <= 0.05, 1, 0),
  )

fdr_est |>
  count(df, term_c, p_corr_fdrc) |>
  filter(p_corr_fdrc == 1)  |>
  ggplot(aes(x = term_c, y = n)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = factor(df))) +
  geom_label(aes(label = n, color = factor(df)))

fdr_est |>
  select(-c(p.value, p_fdrc)) |>
  pivot_wider(names_from = df, values_from = c(estimate, p_corr_fdrc)) |> View()

fdr_est |>
  group_by(term_c) |>
  rstatix::cohens_d(estimate ~ df)

fdr_est |>
  select(-c(p.value, p_fdrc, p_corr_fdrc)) |>
  mutate(sign = ifelse(estimate < 0, "-", "+")) |>
  pivot_wider(names_from = df, values_from = c(estimate, sign)) |> View()

fdr_est |>
  select(-c(p.value, p_fdrc, p_corr_fdrc)) |>
  mutate(sign = ifelse(estimate < 0, "-", "+")) |>
  count(df, term_c, sign)

# saving estimates for brain viz 
for_plotting_long <-
  estimates_2 |>
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

write_csv(p_long, 
          here("outputs", "plots", "brainplots", "brainplots_data", "fluid_inteligence_estimates.csv"))
rm(for_plotting_long, p_long)

for_plotting_long <-
  fdr_est |>
  mutate(estimate = ifelse(p_corr_fdrc == 1, estimate, 0)) |> 
  inner_join(cerebra_mapping, by = c("reg" = "fs_reg")) |>
  select(term_c, name_use, side, final_name, df, mod, estimate, newlabels)

p_long <- 
  left_join(comb, for_plotting_long) |> 
  replace_na(list(estimate = 0)) 

write_csv(p_long, 
          here("outputs", "plots", "brainplots", "brainplots_data", "fluid_inteligence_estimates_thresh.csv"))

# Adding education in the model: fi ~ voi * sex + age + education --------------
#### Getting education data ----
sec_long <- get_ukbb_data(arrow = ukbb, diction = diction, 
                          fields_id = c(738, 6138, 845), subjs_id = full$subject_id,
                          instances = 0:2)

sec <- order_data(sec_long, diction, codings)

sec$rowid <- 1:nrow(sec)
sec$field_value <- as.numeric(sec$field_value)


qualif_relevel <- tibble(field_value = c(-3, -7, 1:6),
                         new_level = c(-3, -7, 7, 4, 3, 3, 5, 6),
                         new_meaning = c("Prefer not to answer", "None of the above",
                                         "College or University degree",
                                         "A levels", "Full secondary",
                                         "Full secondary", "NVQ/HND/HNC",
                                         "Other professional qualifications"))
qualif <-
  sec |> 
  filter(field_id == 6138) |> 
  select(-instance_id, rowid) |>
  left_join(qualif_relevel)

# Assign a "qualification level" to the years of study (for less than 18)
ed_relevel <- tibble(new_level = c(-3:2),
                     new_meaning = c("Prefer not to answer", "Never went",
                                     "Do not know", NA, "Some elementary",
                                     "Some secondary")
)
ed <- 
  sec |> 
  filter(field_id == 845) |> 
  mutate(
    new_level = case_when(
      field_value >= 5 & field_value <= 11 ~ 1,
      field_value >= 12 & field_value <= 17 ~ 2,
      field_value < 0 ~ field_value,
      T ~ 0
    )
  ) |> 
  left_join(ed_relevel)


edu_levels <-
  bind_rows(ed_relevel, qualif_relevel) |> 
  select(education_level = new_level, 
         education_meaning = new_meaning) |> 
  arrange(education_level) |> 
  distinct()

# Binds the qualifications with the education and gets the higher value for 
# each subject
education <- 
  bind_rows(qualif, ed) |> 
  select(subject_id, education_level = new_level, 
         education_meaning = new_meaning) |> 
  group_by(subject_id) |>
  filter(education_level == max(education_level)) |> 
  distinct() |> 
  ungroup()


education_data <- 
  left_join(final_data, education) |>
  mutate(education_level = factor(education_level))

#### Modeling ----
estimates_edu <- tibble()

for(df_u in c(2, 4)) {
  sm_df <- education_data |> filter(subject_id %in% dfs[[df_u]]$subject_id)
  sm_df$fl_i_score_z <- scale(sm_df$fl_i_score)
  
  dfname <- dfs_n[df_u]
  for(r in regions) {
    tmp_df <- sm_df |> select(fl_i_score_z, sex, voi = !!sym(r), age_months, education_level)
    mod_r <- lm(fl_i_score_z ~ 1 + sex * scale(voi) + scale(age_months) + education_level, data = tmp_df)
    est <- broom::tidy(mod_r) |> mutate(final_region = r,
                                        df = dfname,
                                        mod = "fluid_int") 
    estimates_edu <- bind_rows(estimates_edu, est)
  }  
}

glimpse(estimates_edu)

# Comparing the different model's estimates ------------------------------------
comb <- inner_join(estimates, estimates_edu, 
                   by = c("term", "final_region", "df", "mod"),
                   suffix = c("_reg", "_education"))

glimpse(comb)

ggplot(comb, aes(x = estimate_reg, y = estimate_education)) +
  facet_wrap(~term) +
  geom_point() +
  geom_abline(color = "#ff4f4f") +
  theme_light()

