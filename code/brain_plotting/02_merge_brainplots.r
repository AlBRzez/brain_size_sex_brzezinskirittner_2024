library(here)
library(dplyr)
library(ggplot2)
library(ggimage)
library(stringr)
library(glue)
library(ggtext)

source(here("code", "analyses", "helpers.r"))


terms_ref <- # for allometry
  terms_ref |>
  add_row(term = "int_males", term_c = "male_int",
          clean_name = "Intercept (for male)") |>
  add_row(term = "male_slope", term_c = "male_slope",
          clean_name = "Slope for male")


merge_brains <- function(fold, model, m_name, png_name, w = 14, h = 9) {
  
  folder <- paste0("./outputs/plots/brainplots/", fold)
  brains <- list.files(folder,
                       include.dirs  = T,
                       full.names = T)
  
  img_df <-
    tibble(img = brains) |> 
    mutate(
      name = str_split_fixed(img, "/", 6)[,6],
      name = gsub("\\.png", "", name),
      name = gsub(model, "", name),
      mod = model,
      df = str_extract(name, paste(samples_ref$df, collapse = "|")),
      name = gsub(paste(samples_ref$df, collapse = "|"), "", name),
      term_c = gsub("__", "", name)
    ) |> 
    left_join(samples_ref) |>
    left_join(terms_ref |> filter(!term %in% c("scale(month)",
                                               "sexMale:scale(month)",
                                               "scale(I(month^2))",
                                               "sexMale:scale(I(month^2))"))) |>
    left_join(
      stack(colors_samples),
      by = c("clean_sample" = "ind")
    ) |> 
    mutate(
      f_s = glue("<span style='color:{values}'>{clean_sample}</span>"),
      f_s = factor(f_s, levels = c(
        "<span style='color:#001219'>Extreme sample</span>",
        "<span style='color:#ae2012'>Not matched</span>",
        "<span style='color:#005f73'>Age matched</span>",
        "<span style='color:#ee9b00'>TIV and age matched</span>"
      )),
      clean_name = factor(clean_name, levels = c(
        "Intercept (for female)",
        "Intercept (for male)",
        "Age",
        "Age (months)",
        "Sex (male)",
        "Sex (male):age",
        "Sex (male):age (months)",
        "Age<sup>2</sup>",
        "Age<sup>2</sup> (months)",
        "Sex (male):age<sup>2</sup>",
        "Sex (male):age<sup>2</sup> (months)",
        "Slope for male",
        "TIV",
        "TBV"
      ))
    ) 
  
  
  if(sum(complete.cases(img_df)) != nrow(img_df)) {
    stop()
  }
  
  brains_plot <- 
    ggplot(img_df) + 
    geom_image(aes(x = 1, y = 1, image = img), size = 1) +
    facet_grid(f_s~clean_name, scales = "free", switch = "y") +
    theme_void() +
    labs(
      title = m_name
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      strip.text.x = element_textbox(color = "#595959", face = "bold", size = 13),
      strip.text.y = element_textbox(face = "bold", size = 16)
    )
  
  ggsave(here("outputs", "plots", "supplementary", paste0(png_name, ".png")),
         brains_plot, width = w, height = h, bg = "white")
  ggsave(here("outputs", "plots", "supplementary", paste0(png_name, ".svg")),
         brains_plot, width = w, height = h, bg = "white")
  
}


merge_brains("trajectories_reg", "lin_voi_regular", 
             "Sex differenciated aging trajectories", 
             "s2_complete_traj_brains")
# merge_brains("trajectories_reg_colorbar", "lin_voi_regular", 
#              "Sex differenciated aging trajectories", 
#              "s2_complete_traj_brains_cb")

merge_brains("trajectories_thresh", "lin_voi_regular", 
             "Sex differenciated aging trajectories - where p < 0.05", 
             "s3_complete_traj_brains_thresh")
# merge_brains("trajectories_reg_thresh", "lin_voi_regular", 
#              "Sex differenciated aging trajectories - where p < 0.05", 
#              "s3_complete_traj_brains_thresh_cb")

merge_brains("vertex_vol", "vol", 
             "Sex differenciated aging trajectories - volumetric vertexwise data", 
             "s7_vertexwise_volume")
merge_brains("vertex_area", "area", 
             "Sex differenciated aging trajectories - surface area vertexwise data", 
             "s8_vertexwise_area")
merge_brains("vertex_ct", "ct", 
             "Sex differenciated aging trajectories - cortical thickness vertexwise data", 
             "s9_vertexwise_ct")

merge_brains("allom_60m_reg", "lin_int_regular", "Allometry estimates", 
             "s11_complete_allometry", w = 17)
# merge_brains("allom_60m_reg_colorbar", "lin_int_regular", "Allometry estimates", 
#              "s11_complete_allometry_cb", w = 17)

merge_brains("fluid_inteligence_thresh", "fluid_int", 
             "Model over fluid inteligence - where p < 0.05", 
             "s13_fluid_inteligence_thresh")