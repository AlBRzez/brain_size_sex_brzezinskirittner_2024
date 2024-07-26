library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(tidyr)
library(scales)
library(ggtext)
library(stringr)
library(car)
library(glue)
library(purrr)


## Running some checks on the data and comparing estimates with and without 
## including euler number in the model

source(here("code", "analyses", "08_plotting_functions.r"))
colors_samples <- c(colors_samples, "Scale factor" = "#bb3e03")
samples_ref <- 
  samples_ref |> 
  add_row(df = "scale_factor", clean_sample = "Scale factor")
terms_ref <-
  terms_ref |> 
  add_row(term = "log(abs(euler_num) + 1)", term_c = "euler_number", clean_name = "Euler number")


#### Get and clean data --------------------------------------------------------
# references
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))
brain_sizes <- read_csv(here("data", "fs_data", "brain_sizes.csv"))

# Models
lm_coef_o <- 
  read_csv(here("outputs", "trajectories", "fs_coeff.csv"))

# Samples
samp_path <- here("data", "samples")

dfs <- get_samples(samp_path, fs = T)
df_n <- names(dfs)

# Checking brain sizes ---------------------------------------------------------

m_bz <- function(df) {
  out <- df |> 
    select(subject_id, sex, age_years, age_months, euler_num, scale_fact) |> 
    left_join(brain_sizes)
}

dfs_bz <- lapply(dfs, m_bz)
list2env(dfs_bz, envir = .GlobalEnv)

sample_distribution <- function(samp, samp_name) {
  df <- samp |> 
    # select(subject_id, sex) |> 
    pivot_longer(-c(subject_id, sex, age_years, age_months))
  
  p <- 
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
  return(p)
}

for(i in 1:length(df_n)) {
  print(sample_distribution(dfs_bz[[i]], df_n[i]))
}

# Correlation  -----------------------------------------------------------------
#### wrangle data ----------------

fs_c <- 
  lm_coef_o |> 
  filter(effect == "fixed") |> 
  left_join(terms_ref) |> 
  left_join(samples_ref)  |>
  left_join(fs_guide) |> 
  filter(ignore == 0) |> 
  rename("final_mod" = "mod") |> 
  group_by(df, final_mod, term_c) |> 
  nest() |>
  mutate(
    p_fdrc = map(data, function(df)  p.adjust(df$p.value, method = "fdr")),
  ) |>
  unnest(cols = c(data, p_fdrc)) |>
  ungroup() |> 
  mutate(
    p_corr_fdrc = ifelse(p_fdrc <= 0.05, 1, 0),
  )


models_guide <-
  tibble(final_mod = unique(fs_c$final_mod)) |> 
  mutate(
    type = ifelse(grepl("lin", final_mod), 
                  "y ~ 1 + age * sex + log(euler) + (1|site)",
                  "y ~ 1 + (age + age<sup>2</sup>) * sex + log(euler) + (1|site)"),
    type = ifelse(grepl("euler", final_mod), 
                  type,
                  gsub(" \\+ log\\(euler\\) ", "", type)),
    analysis = "aging trajectories"
  )
#########

arr_vdata <- function(df) {
  out <- df |> 
    left_join(models_guide, by = c("final_mod")) |> 
    mutate(
      clean_sample = factor(clean_sample, levels = c("Extreme sample",
                                                     "Not matched",
                                                     "Age matched",
                                                     "TIV and age matched",
                                                     "FS e-TIV and age matched"

      )),
      clean_name = factor(clean_name, levels = c(
        "Intercept",
        "Age",
        "Age (months)",
        "Sex (male)",
        "Sex (male) * age",
        "Sex (male) * age (months)",
        "Age<sup>2</sup>",
        "Age<sup>2</sup> (months)",
        "Sex (male) * age<sup>2</sup>",
        "Sex (male) * age<sup>2</sup> (months)",
        "TIV",
        "TBV",
        "Euler number"
      ))
    ) 
  return(out)
}




fs_c |> 
  mutate(
  # transmute(final_mod = final_mod,
    batch = str_split(final_mod, "_", simplify = T)[,3],
         mod = gsub("_regular|_euler", "", final_mod)) |> 
  select(df, term_c, batch, mod, reg, estimate) |> 
  pivot_wider(names_from = batch, values_from = estimate) |> 
  filter(term_c != "euler_number") |> 
  group_by(df, term_c, mod) |> 
  summarise(cor = cor(regular, euler)) |> 
  # pivot_wider(names_from = mod, values_from = cor) |> 
  ggplot(aes(x = term_c, y = df)) +
  facet_wrap(~mod) +
  geom_tile(aes(fill = cor)) +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0, 1)) +
  theme_light()



