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
  

fold1 <- "vertex_vol_colorbar"

folder <- paste0("./outputs/plots/brainplots/", fold1)

brains <- list.files(folder,
                     include.dirs  = T,
                     full.names = T)

model <- "vol"
# m_name <- "Model estimates: linear interaction model"
m_name <- "Sex differenciated aging trajectories - volumetric vertexwise data "
# m_name <- "Sex differenciated aging trajectories - where p < 0.05"
# m_name <- "Allometry estimates"

img_df <-
  tibble(img = brains) |> 
  mutate(
    name = str_split_fixed(img, "/", 6)[,6],
    # name = gsub(fold1, "", name),
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


# img_df <- img_df |> filter(!grepl("colorbar", img))

img_df

brains_plot <- 
  ggplot(img_df) + 
  geom_image(aes(x = 1, y = 1, image = img), size = 1) +
  # facet_grid(clean_name~f_s, scales = "free", switch = "y") +
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

# brains_plot

ggsave(here("outputs", "plots", "supplementary", "s7_vertexwise_vol.png"),
       brains_plot, width = 14, height = 9, bg = "white")
       # brains_plot, width = 17, height = 9, bg = "white")
ggsave(here("outputs", "plots", "supplementary", "s7_vertexwise_vol.svg"),
       brains_plot, width = 14, height = 9, bg = "white")
       # brains_plot, width = 17, height = 9, bg = "white")
       
       

  

