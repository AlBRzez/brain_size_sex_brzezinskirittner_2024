color_guide <- c("#001219", "#005f73", "#0a9396", "#94d2bd", "#e9d8a6", 
                 "#ee9b00", "#ca6702", "#bb3e03", "#ae2012", "#9b2226")

terms_ref <- tibble(
  term = c("(Intercept)", 
           "sexMale", 
           "scale(age_years)",
           "sexMale:scale(age_years)", 
           "scale(I(age_years^2))",
           "sexMale:scale(I(age_years^2))", 
           "scale(tiv)", 
           "scale(tbv)", 
           "scale(month)",
           "sexMale:scale(month)", 
           "scale(I(month^2))",
           "sexMale:scale(I(month^2))", 
           "scale(age_months)",
           "sexMale:scale(age_months)", 
           "scale(I(age_months^2))",
           "sexMale:scale(I(age_months^2))", 
           "scale(log_euler)", 
           "log10(tbv)"),
  term_c = c("intercept", 
             "sex_male", 
             "age", 
             "sex_male_age", 
             "age_sq", 
             "sex_male_age_sq", 
             "tiv", 
             "tbv", 
             "agem", 
             "sex_male_agem", 
             "agem_sq", 
             "sex_male_agem_sq", 
             "agem", 
             "sex_male_agem", 
             "agem_sq", 
             "sex_male_agem_sq", 
             "euler", 
             "log10_tbv"),
  clean_name = c("Intercept (for female)", 
                 "Sex (male)", 
                 "Age", 
                 "Sex (male):age",
                 "Age<sup>2</sup>", 
                 "Sex (male):age<sup>2</sup>", 
                 "TIV", 
                 "TBV", 
                 "Age (months)", 
                 "Sex (male):age (months)",
                 "Age<sup>2</sup> (months)", 
                 "Sex (male):age<sup>2</sup> (months)",
                 "Age (months)", 
                 "Sex (male):age (months)",
                 "Age<sup>2</sup> (months)", 
                 "Sex (male):age<sup>2</sup> (months)",
                 "Euler number", 
                 "TBV"))

dfs_n_c <- c("full", "matched", "age_mat", "random", "extreme", "fs_matc")

samples_ref <- tibble(
  df = dfs_n_c,
  clean_sample = c(
    "Full sample", 
    "TIV and age matched", 
    "Age matched", 
    "Not matched",
    "Extreme sample",
    "FS e-TIV and age matched"))

colors_samples <- c(
  "TIV and age matched" = "#ee9b00", 
  "Age matched" = "#005f73", 
  "Not matched" = "#ae2012", 
  "Extreme sample" = "#001219",
  "Full sample" = "#595959",
  "FS e-TIV and age matched" = "#90a955"
  
)

samp_order <- c(
  "Full sample",
  "Extreme sample",
  "Not matched",
  "Age matched",
  "TIV and age matched",
  "FS e-TIV and age matched"
)
  
# color_sex <- c("Male" = "#26547C", "Female" = "#bb3e03")
color_sex <- c("Male" = "#F4632E", "Female" = "#6A559C")
# color_sex <- c("Male" = "#EF5423", "Female" = "#4C3681")

get_samples <- function(path, fs = FALSE) {
  full <- readr::read_csv(paste0(path, "/full_sample.csv"))
  matched <- readr::read_csv(paste0(path, "/scalefact_matched.csv"))
  age_mat <- readr::read_csv(paste0(path, "/age_matched.csv"))
  random <- readr::read_csv(paste0(path, "/random_sample.csv"))
  extreme <- readr::read_csv(paste0(path, "/extreme_sample.csv"))
  out <- list(
    "full" = full, 
    "matched" = matched, 
    "age_mat" = age_mat, 
    "random" = random, 
    "extreme" = extreme)
  
  if(fs) {
    fs_tiv <- readr::read_csv(paste0(path, "/fstiv_matched.csv"))
    out$fs_matc <- fs_tiv
  }
  
  return(out)
}

plots_theme <-  
  cowplot::theme_half_open() +
  # background_grid() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = NA),
    strip.text = ggtext::element_markdown(color = "#595959", face = "bold", 
                                          size = 14, hjust = 0),
    axis.ticks = ggplot2::element_blank()
  )


id_samples <- function(df) {
  df |> 
    dplyr::mutate(
      matched = ifelse(subject_id %in% matched$subject_id, 1, 0),
      age_mat = ifelse(subject_id %in% age_mat$subject_id, 1, 0),
      random = ifelse(subject_id %in% random$subject_id, 1, 0),
      extreme = ifelse(subject_id %in% extreme$subject_id, 1, 0),
      fs_matc = ifelse(subject_id %in% fs_matc$subject_id, 1, 0),
    )
}