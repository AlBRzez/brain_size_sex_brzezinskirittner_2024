library(dplyr)
library(readr)
library(here)
library(tidyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(scales)
library(glue)
library(ggradar)
library(patchwork)

source(here("code", "analyses", "helpers.r"))
# Get atlases info --------
yeo <- 
  read_csv(here("data", "atlases", "Yeo_7_label.csv")) |> 
  mutate(Schaefer_7_1000_7 = 1:7,
         yeo_network = gsub("Attnention", "Attention", lebels),
         yeo_network_ab = c("VN", "SMN", "DAN", "VAN", "LBN", "FPN", "DMN"),
         yeo_network_lab = glue("{yeo_network} ({yeo_network_ab})"))

cyto <- read_csv(here("data", "atlases", "VonEconomo_label_Cyto.csv")) |> 
  mutate(VonEconomo_Cyto_lab = str_to_sentence(region_name), 
         VonEconomo_Cyto_lab = gsub("\\s*\\([^\\)]+\\)", "", VonEconomo_Cyto_lab),
         VonEconomo_Cyto_lab = glue("{VonEconomo_Cyto_lab} ({region_name_ab})"),
         VonEconomo_Cyto_lab = gsub("Secondarysensory", "Secondary sensory", VonEconomo_Cyto_lab),
         cyto_order = c(1, 5, 4, 2, 3, 7, 6)
  ) |> 
  rename("VonEconomo_Cyto" = "region_num")

laminar <- read_csv(here("data", "atlases", "VonEconomo_label_Laminar.csv")) |> 
  rename("VonEconomo_Laminar" = "region_num",
         "VonEconomo_Laminar_lab" = "region_name")

lobe <- read_csv(here("data", "atlases", "VonEconomo_label_lobe.csv")) |> 
  rename("VonEconomo_lobe" = "region_num",
         "VonEconomo_lobe_lab" = "region_name")

# Get and arrange lm data ------------------------------------------------------

atlases <- read_csv(here("data", "atlases", "all_atlases_surface.csv"))
area <- read_csv(here("outputs", "vertexwise", "b_area_civet_all_20.csv"))
volume <- read_csv(here("outputs", "vertexwise", "b_volume_civet_all_20.csv"))
ct <- read_csv(here("outputs", "vertexwise", "b_thickness_civet_all_20.csv"))

atlases <- 
  atlases |> 
  select(Schaefer_7_1000_7, VonEconomo_Laminar, VonEconomo_Cyto, VonEconomo_lobe)

# Identifying which rows to remove from the vertices df 
ind_remove <-   
  volume |> 
  transmute(
    ads = rowSums(across(where(is.numeric))),
    id = 1:nrow(area),
  ) |> 
  filter(ads == 0) |> 
  pull(id)


final_surface <- 
  bind_cols(volume[-ind_remove,] , atlases[-ind_remove,] ) |>
  mutate(id = row_number()) |> 
  pivot_longer(-c(id, Schaefer_7_1000_7, VonEconomo_Laminar, VonEconomo_Cyto, VonEconomo_lobe), 
               values_to = "volume") |> 
  left_join(
    bind_cols(area[-ind_remove,] , atlases[-ind_remove,] ) |>
      mutate(id = row_number()) |> 
      pivot_longer(-c(id, Schaefer_7_1000_7, VonEconomo_Laminar, VonEconomo_Cyto, VonEconomo_lobe), 
                   values_to = "area")
  ) |> 
  left_join(
    bind_cols(ct[-ind_remove,] , atlases[-ind_remove,] ) |>
      mutate(id = row_number()) |> 
      pivot_longer(-c(id, Schaefer_7_1000_7, VonEconomo_Laminar, VonEconomo_Cyto, VonEconomo_lobe), 
                   values_to = "ct")
  ) |> 
  mutate(name = gsub("age_matched", "agemat", name),
         name = gsub("age_sex", "agesex", name)) |>
  separate(name, into = c("df", "term_c"), sep = "_") |>
  mutate(df = gsub("agemat", "age_mat", df),
         term_c = gsub("agesex", "age_sex", term_c)) 

# Sex estimates ----------------------------------------------------------------
sex_estimates <- 
  final_surface |> 
  filter(term_c == "sex") |> # & df %in% c("matched", "random"))
  left_join(samples_ref) |> 
  mutate(clean_sample = factor(clean_sample, 
                               levels = c("Extreme sample",
                                          "Not matched",
                                          "Age matched",
                                          "TIV and age matched",
                                          "FS e-TIV and age matched")))

final_surface |>
  group_by(df, term_c) |>
  summarise(across(c(volume, area, ct), list(mean = mean, sd = sd), na.rm = TRUE)) |>
  filter(term_c == "sex") 

final_surface |>
  group_by(df, term_c) |>
  summarise(
    r_vol_ct = cor.test(volume, ct)$estimate,
    p_vol_ct = cor.test(volume, ct)$p.value,
    r_vol_ar = cor.test(volume, area)$estimate,
    p_vol_ar = cor.test(volume, area)$p.value,
    r_ar_ct = cor.test(ct, area)$estimate,
    p_ar_ct = cor.test(ct, area)$p.value,
  ) |> View()

## Sex estimates plot ----

p_sex_vol <- 
  ggplot(sex_estimates, 
         aes(x = factor(clean_sample), y = volume, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Volume",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_blank()
  ) 

p_sex_area <- 
  ggplot(sex_estimates, 
         aes(x = factor(clean_sample), y = area, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Surface area",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_blank()
  ) 

p_sex_ct <- 
  ggplot(sex_estimates, 
         aes(x = factor(clean_sample), y = ct, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Cortical thickness",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
  ) 

sex_estimates_plot <- plot_grid(p_sex_vol, p_sex_area, p_sex_ct, ncol = 1,
                                rel_heights = c(.75, .75, 1))

ggsave(here("outputs", "plots", "main_figures", "fig2_sex_samples.png"),
       sex_estimates_plot, width = 4, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig2_sex_samples.svg"),
       sex_estimates_plot, width = 4, height = 10, bg = "white", dpi = 600)

## Shifting the orientation
p_sex_vol_s <- 
  ggplot(sex_estimates, aes(x = factor(clean_sample), y = volume, 
                            color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Volume",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_text(size = 9)
  ) 

p_sex_area_s <- 
  ggplot(sex_estimates, aes(x = factor(clean_sample), y = area, 
                            color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Surface area",
    x = "",
    y = "",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_text(size = 9)
  ) 

p_sex_ct_s <- 
  ggplot(sex_estimates, aes(x = factor(clean_sample), y = ct, 
                            color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Cortical thickness",
    x = "",
    y = "",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_text(size = 9)
  ) 


sex_estimates_plot_h <- plot_grid(p_sex_vol_s, p_sex_area_s, p_sex_ct_s, nrow = 1)#,


ggsave(here("outputs", "plots", "main_figures", "fig2_sex_samples_v2.png"),
       sex_estimates_plot_h, width = 11, height = 4, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig2_sex_samples_v2.svg"),
       sex_estimates_plot_h, width = 11, height = 4, bg = "white", dpi = 600)

## For biological meaning ------------------------------------------------------
matched_sex <- 
  sex_estimates |> 
  filter(df == "matched") |> 
  left_join(yeo) |> 
  left_join(cyto) |> 
  left_join(laminar) |> 
  left_join(lobe) 

random_sex <- 
  sex_estimates |> 
  filter(df == "random") |> 
  left_join(yeo) |> 
  left_join(cyto) |> 
  left_join(laminar) |> 
  left_join(lobe) 

yeo_c <- tibble(Schaefer_7_1000_7 = 1:7,
                color = c("#781286", "#4682B4", "#00760E",  
                          "#C43AFA", "#DCF8A4", "#E69422", "#CD3E4E"))
cyto_c <- tibble(cyto_order = 1:7,
                 color = c("#781286", "#4682B4", "#00760E",  "#E69422",
                           "#DCF8A4", "#00afb9", "#C43AFA"))

laminar_c <- tibble(VonEconomo_Laminar = 1:4,
                    color = c("#99d98c", "#52b69a", "#168aad", "#1e6091"))

lobe_c <- tibble(VonEconomo_lobe = 1:7,
                 color = c("#781286", "#4682B4", "#00760E",  
                           "#C43AFA", "#DCF8A4", "#CD3E4E", "#E69422"))


# functions for boxplots for the different classifications
boxplot_func <- function(sm_df, l_order, l_names, subtitle, y, title) {
  
  ggplot(sm_df, aes(x = reorder(!!sym(l_names), !!sym(l_order)), y = !!sym(y), 
                    fill = color, group = !!sym(l_names))) +
    geom_boxplot(width = .3, show.legend = FALSE, alpha = .5) +
    geom_hline(yintercept = 0,
               color = "#595959", linetype = "dashed") +
    
    labs(
      subtitle = title,
      x = subtitle,
      y = "Estimates",
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    plots_theme +
    theme(
      plot.subtitle = element_text(size = 16, 
                                   face = "bold", color = "#595959"),
      legend.position = "none"
    ) 
}

clasifications <- function(df, l_order, l_names, subtitle, color_df,
                           orientation = c("vertical", "horizontal")) {
  tmp <- 
    df |> 
    drop_na(!!sym(l_names)) |> 
    left_join(color_df, by = l_order)
  
  v <- boxplot_func(tmp, l_order, l_names, subtitle, 
                    "volume", "Volume")
  a <- boxplot_func(tmp, l_order, l_names, subtitle, 
                    "area", "Surface Area")
  c <- boxplot_func(tmp, l_order, l_names, subtitle, 
                    "ct", "Cortical thickness")
  
  if(orientation[1] == "vertical") {
    v <- v + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank())
    a <- a + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank())
    plot_grid(v, a, c, ncol = 1, rel_heights = c(.75, .75, 1))
  } else {
    
    a <- a + theme(axis.title.y = element_blank())
    c <- c + theme(axis.title.y = element_blank())
    plot_grid(v, a, c, nrow = 1, rel_heights = c(1, .8, .8))
  }
  
}

## Sex estimates in the matched sample

yeo_sex_mv <- clasifications(matched_sex, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "vertical")
cyto_sex_mv <- clasifications(matched_sex, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "vertical")

ggsave(here("outputs", "plots", "main_figures", "fig2_matched_yeo_sex_v1.png"),
       yeo_sex_mv, width = 4, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig2_matched_cyto_sex_v1.png"),
       cyto_sex_mv, width = 4, height = 10, bg = "white", dpi = 600)

ggsave(here("outputs", "plots", "main_figures", "fig2_matched_yeo_sex_v1.svg"),
       yeo_sex_mv, width = 4, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig2_matched_cyto_sex_v1.svg"),
       cyto_sex_mv, width = 4, height = 10, bg = "white", dpi = 600)

yeo_sex_mh <- clasifications(matched_sex, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "horizontal")
cyto_sex_mh <- clasifications(matched_sex, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "horizontal")

ggsave(here("outputs", "plots", "main_figures", "fig2_matched_yeo_sex_v2.png"),
       yeo_sex_mh, width = 11, height = 4, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig2_matched_cyto_sex_v2.png"),
       cyto_sex_mh, width = 11, height = 4, bg = "white", dpi = 600)

ggsave(here("outputs", "plots", "main_figures", "fig2_matched_yeo_sex_v2.svg"),
       yeo_sex_mh, width = 11, height = 4, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig2_matched_cyto_sex_v2.svg"),
       cyto_sex_mh, width = 11, height = 4, bg = "white", dpi = 600)

## Sex estimates in the not-matched sample
yeo_sex_rv <- clasifications(random_sex, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "vertical")
cyto_sex_rv <- clasifications(random_sex, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "vertical")

yeo_sex_rh <- clasifications(random_sex, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "horizontal")
cyto_sex_rh <- clasifications(random_sex, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "horizontal")


## Correlations between metrics ------------------------------------------------
get_cors <- function(dat) {
  dat |>
    summarise(
      r_vol_ct = cor.test(volume, ct)$estimate,
      p_vol_ct = cor.test(volume, ct)$p.value,
      r_vol_ar = cor.test(volume, area)$estimate,
      p_vol_ar = cor.test(volume, area)$p.value,
      r_ar_ct = cor.test(ct, area)$estimate,
      p_ar_ct = cor.test(ct, area)$p.value,
    ) 
}

matched_sex |>
  drop_na(yeo_network) |>
  group_by(yeo_network) |>
  get_cors()

matched_sex |>
  drop_na(VonEconomo_Cyto_lab) |>
  group_by(VonEconomo_Cyto_lab) |>
  get_cors()


# Age estimates ----------------------------------------------------------------
age_estimates <- 
  final_surface |> 
  filter(term_c == "age") |> 
  left_join(samples_ref) |> 
  mutate(clean_sample = factor(clean_sample, 
                               levels = c("Extreme sample",
                                          "Not matched",
                                          "Age matched",
                                          "TIV and age matched",
                                          "FS e-TIV and age matched")))

## age estimates plot ----------------------------------------------------------
p_age_vol <- 
  ggplot(age_estimates, aes(x = factor(clean_sample), y = volume, 
                            color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = median(age_estimates$volume[age_estimates$df == "random"]),
             show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Volume",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_blank()
  ) 

p_age_area <- 
  ggplot(age_estimates, aes(x = factor(clean_sample), y = area, 
                            color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  geom_hline(yintercept = median(age_estimates$area[age_estimates$df == "random"]),
             show.legend = FALSE) +
  
  labs(
    subtitle = "Surface area",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_blank()
  ) 

p_age_ct <- 
  ggplot(age_estimates, aes(x = factor(clean_sample), y = ct, 
                            color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  geom_hline(yintercept = median(age_estimates$ct[age_estimates$df == "random"]),
             show.legend = FALSE) +
  
  labs(
    subtitle = "Cortical thickness",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 16, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
  ) 

age_estimates_plot <- plot_grid(p_age_vol, p_age_area, p_age_ct, ncol = 1,
                                rel_heights = c(.75, .75, 1))

age_estimates |>
  select(id, df, volume) |>
  pivot_wider(names_from = df, values_from = volume) |>
  ggplot(aes(x = random, y = matched)) +
  geom_point(size = .6, alpha = .4) +
  geom_abline(color = "#ff4f4f") +
  theme_light()


ggsave(here("outputs", "plots", "supplementary", "s4_age_samples.png"),
       age_estimates_plot, width = 4, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "supplementary", "s4_age_samples.svg"),
       age_estimates_plot, width = 4, height = 10, bg = "white", dpi = 600)

## For biological meaning --------
matched_age <- 
  age_estimates |> 
  filter(df == "matched") |> 
  left_join(yeo) |> 
  left_join(cyto) |> 
  left_join(laminar) |> 
  left_join(lobe) 

random_age <- 
  age_estimates |> 
  filter(df == "random") |> 
  left_join(yeo) |> 
  left_join(cyto) |> 
  left_join(laminar) |> 
  left_join(lobe) 


yeo_age_mv <- clasifications(matched_age, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "vertical")
cyto_age_mv <- clasifications(matched_age, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "vertical")

yeo_age_mh <- clasifications(matched_age, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "horizontal")
cyto_age_mh <- clasifications(matched_age, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "horizontal")

## age estimates in the not-matched sample
yeo_age_rv <- clasifications(random_age, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "vertical")
cyto_age_rv <- clasifications(random_age, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "vertical")

yeo_age_rh <- clasifications(random_age, "Schaefer_7_1000_7", "yeo_network_ab", "Functional networks", yeo_c, "horizontal")
cyto_age_rh <- clasifications(random_age, "cyto_order", "region_name_ab", "Cytoarchitectonic classes", cyto_c, "horizontal")

# Supplementary figures for biological meaning ---------------------------------

s_b1 <- plot_grid(yeo_sex_mv, cyto_sex_mv) + 
  plot_annotation(title = "Matched sample",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", color = "#ee9b00", hjust = .5)))
s_b2 <- plot_grid(yeo_sex_rv, cyto_sex_rv) + 
  plot_annotation(title = "Not matched sample",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", color = "#ae2012", hjust = .5)))
s_b3 <- plot_grid(s_b1, s_b2)
s_bf <- s_b3 + plot_annotation("Sex estimates",
                    theme = theme(plot.title = element_text(size = 18, face = "bold")))

ggsave(here("outputs", "plots", "supplementary", "s5_biol_sex_samples.png"),
       s_bf, width = 16, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "supplementary", "s5_biol_sex_samples.png"),
       s_bf, width = 16, height = 10, bg = "white", dpi = 600)


a_b1 <- plot_grid(yeo_age_mv, cyto_age_mv) + 
  plot_annotation(title = "Matched sample",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", color = "#ee9b00")))
a_b2 <- plot_grid(yeo_age_rv, cyto_age_rv) + 
  plot_annotation(title = "Not matched sample",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", color = "#ae2012")))
a_b3 <- plot_grid(a_b1, a_b2)
a_bf <- a_b3 + plot_annotation("Age estimates",
                               theme = theme(plot.title = element_text(size = 18, face = "bold")))

ggsave(here("outputs", "plots", "supplementary", "s6_biol_age_samples.png"),
       a_bf, width = 16, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "supplementary", "s6_biol_age_samples.png"),
       a_bf, width = 16, height = 10, bg = "white", dpi = 600)

# comparing between samples ----------------------------------------------------

sex_data <-   
  final_surface |>
  filter(term_c == "sex") |>  
  left_join(samples_ref)  

cyto_radar <- 
  sex_data |>
  filter(VonEconomo_Cyto != 0) |>
  group_by(clean_sample, VonEconomo_Cyto) |>
  summarise(
    vol = median(volume),
    area = median(area),
    ct = median(ct),
  ) |>
  
  ungroup() |> 
  pivot_longer(-c(clean_sample, VonEconomo_Cyto), 
               names_to = "metric", values_to = "val") |>
  left_join(cyto) |> 
  arrange(cyto_order) |> 
  select(val, metric, clean_sample, region_name_ab) |> 
  pivot_wider(names_from = region_name_ab, values_from = val) 


sch_radar <- 
  sex_data |>
  filter(Schaefer_7_1000_7 != 0) |>
  group_by(clean_sample, Schaefer_7_1000_7) |>
  summarise(
    vol = median(volume),
    area = median(area),
    ct = median(ct),
  ) |>
  ungroup(Schaefer_7_1000_7) |> 
  pivot_longer(-c(clean_sample, Schaefer_7_1000_7), 
               names_to = "metric", values_to = "val") |>
  left_join(yeo) |> 
  arrange(Schaefer_7_1000_7) |> 
  select(val, metric, clean_sample, yeo_network_ab) |> 
  pivot_wider(names_from = yeo_network_ab, values_from = val) 

# Normalizing values
cyto_radar2 <- cyto_radar[,-c(1,2)]/apply(abs(cyto_radar[,-c(1,2)]), 1, max)
cyto_radar2 <- cbind(cyto_radar[, 1:2], cyto_radar2)

sch_radar2 <- sch_radar[,-c(1,2)]/apply(abs(sch_radar[,-c(1,2)]), 1, max)
sch_radar2 <- cbind(sch_radar[, 1:2], sch_radar2)

vol_cyto <-
  cyto_radar2 |> 
  filter(metric == "vol") |> 
  select(-metric)

area_cyto <-
  cyto_radar2 |> 
  filter(metric == "area") |> 
  select(-metric)

ct_cyto <-
  cyto_radar2 |> 
  filter(metric == "ct") |> 
  select(-metric)

vol_sch <-
  sch_radar2 |> 
  filter(metric == "vol") |> 
  select(-metric) 

area_sch <-
  sch_radar2 |> 
  filter(metric == "area") |> 
  select(-metric) 

ct_sch <-
  sch_radar2 |> 
  filter(metric == "ct") |> 
  select(-metric) 



sex_radar <- function(dat, m1, m2)  {
  ggradar(
    dat,
    grid.min = m1, grid.max = m2, grid.mid = 0,
    values.radar = c(m1, 0, m2),
    group.colours = colors_samples[1:4],
    background.circle.colour = "#FFFFFF",
    legend.position = "bottom",
    base.size = 11,
    gridline.min.colour = "#cccfcf",
    gridline.mid.colour = "#191919",
    gridline.max.colour = "#cccfcf",
    legend.text.size = 10,
    group.line.width = 1,
    group.point.size = 4,
    gridline.min.linetype = 2,
    gridline.mid.linetype = 1,
    gridline.max.linetype = 2,
    axis.label.offset = 1.15,
    axis.label.size = 6,
    axis.labels = str_wrap(colnames(dat)[-1], 12),
    grid.label.size = 6,
    gridline.label.offset = 0,
    plot.legend = F
  ) #+
  # guides(color = guide_legend(override.aes = list(linetype = 0)))
}

min(vol_cyto[,-1])
max(vol_cyto[,-1])
radar_cyto_vol <- sex_radar(vol_cyto, -1, 1)
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_cyto_vol.png"),
       radar_cyto_vol, width = 5, height = 5, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_cyto_vol.svg"),
       radar_cyto_vol, width = 5, height = 5, dpi = 600, bg = "white")

min(ct_cyto[,-1])
max(ct_cyto[,-1])
radar_cyto_ct <- sex_radar(ct_cyto, -1, 1)
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_cyto_ct.png"),
       radar_cyto_ct, width = 5, height = 5, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_cyto_ct.svg"),
       radar_cyto_ct, width = 5, height = 5, dpi = 600, bg = "white")

min(area_cyto[,-1])
max(area_cyto[,-1])
radar_cyto_area <- sex_radar(area_cyto, -1, 1)
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_cyto_area.png"),
       radar_cyto_area, width = 5, height = 5, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_cyto_area.svg"),
       radar_cyto_area, width = 5, height = 5, dpi = 600, bg = "white")

min(vol_sch[,-1])
max(vol_sch[,-1])
radar_sch_vol <- sex_radar(vol_sch, -1, 1)
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_sch_vol.png"),
       radar_sch_vol, width = 5, height = 5, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_sch_vol.svg"),
       radar_sch_vol, width = 5, height = 5, dpi = 600, bg = "white")

min(ct_sch[,-1])
max(ct_sch[,-1])
radar_sch_ct <- sex_radar(ct_sch, -1, 1)
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_sch_ct.png"),
       radar_sch_ct, width = 5, height = 5, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_sch_ct.svg"),
       radar_sch_ct, width = 5, height = 5, dpi = 600, bg = "white")

min(area_sch[,-1])
max(area_sch[,-1])
radar_sch_area <- sex_radar(area_sch, -1, 1)
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_sch_area.png"),
       radar_sch_area, width = 5, height = 5, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "main_figures", "fig2_radar_sch_area.svg"),
       radar_sch_area, width = 5, height = 5, dpi = 600, bg = "white")

