library(dplyr)
library(readr)
library(ggplot2)
library(readr)
library(here)
library(scales)
library(cowplot)
library(tidyr)
library(purrr)
library(ggtext)
library(glue)
library(googlesheets4)

options(scipen = 999)

source(here("code", "analyses", "helpers.r"))

# Get samples ------------------------------------------------------------------
samp_path <- here("data", "samples")

dfs <- get_samples(samp_path, fs = T)
list2env(dfs, envir = .GlobalEnv)
df_n <- names(dfs)

brain_sizes <- read_csv(here("data", "fs_data", "brain_sizes.csv"))
# Samples descriptions ---------------------------------------------------------
get_description <- function(df) {
  df |> 
    group_by(sex) |> 
    summarise(
      mean_a = mean(age_years),
      sd_a = sd(age_years),
      mean_tiv = mean(tiv),
      sd_tiv = sd(tiv),
      n = n(),
      min_age = min(age_years),
      max_age = max(age_years),
    ) |> 
    mutate(
      lwr_a = mean_a - sd_a,
      upr_a = mean_a + sd_a,
      prc = n/sum(n)
    )
}

summary_stats <- map(dfs, get_description) |> bind_rows(.id = "df")
write_csv(summary_stats, here("outputs", "others", "summary_stats"))

nice_stats <- 
  summary_stats |> 
  transmute(
    df = df,
    Sex = sex,
    N = glue("{comma(n)} ({percent(prc)})"),
    `Age range` = glue("{min_age} - {max_age} years"),
    `Age (mean, sd)` = glue("{comma(mean_a, accuracy = .1)}, ({comma(sd_a, accuracy = .01)})"),
    `TIV (mean, sd)` = glue("{comma(mean_tiv, accuracy = .1)}, ({comma(sd_tiv, accuracy = .01)})"),
  )

write_csv(nice_stats, here("outputs", "others", "summary_stats"))

# Plots -----
# The following functions gives the distributions as histograms, first for age 
# and then for tiv
sample_distribution <- function(samp, samp_name) {
  ggplot(samp, aes(x = age_years, color = sex, fill = sex)) +
    geom_histogram(aes(linetype = sex),
                   # fill = NA,
                   alpha = 0,
                   binwidth = 1,
                   position = "identity") +
    labs(title = paste("Demographics:", samp_name),
         subtitle = paste(comma(nrow(samp)), "subjects"),
         x = "Age",
         y = "Participants") +
    scale_color_manual(values = color_sex,
                       aesthetics = c("color", "fill")) +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(labels = comma) +
    theme_light() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_text(color = "#595959")
    ) +
    guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                    alpha = c(1,1))))
}

full_sample_dist <- sample_distribution(get(samples_ref$df[1]), samples_ref$clean_sample[1])
sample_distribution(get(samples_ref$df[2]), samples_ref$clean_sample[2])
sample_distribution(get(samples_ref$df[3]), samples_ref$clean_sample[3])
sample_distribution(get(samples_ref$df[4]), samples_ref$clean_sample[4])
sample_distribution(get(samples_ref$df[5]), samples_ref$clean_sample[5])

# ggsave(here("outputs", "final_plots", "full_samp_age_dist.png"), full_sample_dist,
#        width = 6, height = 4)

sample_tiv_distribution <- function(samp, samp_name) {
  ggplot(samp, aes(x = tiv, color = sex, fill = sex)) +
    geom_histogram(aes(linetype = sex),
                   alpha = 0,
                   position = "identity",
                   bins = 50) +
    labs(title = paste("Demographics:", samp_name),
         subtitle = paste(comma(nrow(samp)), "subjects"),
         x = "Total intracranial volume (TIV)",
         y = "Participants") +
    scale_color_manual(values = color_sex,
                       aesthetics = c("color", "fill")) +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(labels = comma, limits = c(1000000, 2000000)) +
    plots_theme +
    guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                    alpha = c(1,1))))
}

full_sample_tiv_dist <- sample_tiv_distribution(get(samples_ref$df[1]), samples_ref$clean_sample[1])
matched_sample_tiv_dist <- sample_tiv_distribution(get(samples_ref$df[2]), samples_ref$clean_sample[2])


# ggsave(here("plots", "final_plots", "full_samp_tiv_dist.png"), full_sample_tiv_dist,
#        width = 6, height = 4)


#Composite plots for age and tiv ----
comp_plot <- function(samp, samp_name) {
  a <- ggplot(samp, aes(x = age_months, color = sex, fill = sex)) +
    geom_density(aes(y = after_stat(count) * 12,
                     linetype = sex, linewidth = sex),
                 alpha = 0,
                 position = "identity",
    ) +
    labs(
      subtitle = "Age distribution",
      x = "Age",
      y = "Participants"
    ) +
    scale_color_manual(values = color_sex,
                       aesthetics = c("color", "fill")) +
    scale_x_continuous(limits = c(530, 990), breaks = seq(540, 960, 60),
                       labels = seq(45, 80, 5)) +
    scale_y_continuous(labels = comma) +
    scale_linewidth_manual(values = c(1, 1.5)) +
    plots_theme +
    theme(legend.position = "none",
          plot.subtitle = element_text(face = "bold", size = 18, 
                                       color = "#595959"),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)
    ) 
  
  t <- ggplot(samp, aes(x = tiv, color = sex, fill = sex)) +
    geom_density(aes(y = after_stat(count) * 10000,
                     linetype = sex, linewidth = sex),
                 alpha = 0,
                 position = "identity",
    ) +
    labs(
      subtitle = "TIV distribution",
      x = "Total intracranial volume (cubic mm)",
      y = "Participants"
    ) +
    scale_color_manual(values = color_sex,
                       aesthetics = c("color", "fill")) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()), 
                       limits = c(1000000, 2000000)) +
    scale_y_continuous(labels = comma) +
    scale_linewidth_manual(values = c(1, 1.5)) +
    plots_theme +
    theme(
          plot.subtitle = element_text(face = "bold", size = 18, 
                                       color = "#595959"),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16)
    ) +
    guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                    alpha = c(1,1))))
  
  p <- plot_grid(a, t, nrow = 2)
  
  # now add the title
  title <- ggdraw() + 
    draw_label(
      samp_name,
      fontface = 'bold',
      x = 0,
      hjust = 0,
      size = 18,
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 50)
    ) 
  
  fp <- plot_grid(
    title, p,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  return(fp)
}


cp_f <- comp_plot(get(samples_ref$df[1]), samples_ref$clean_sample[1])
cp_m <- comp_plot(get(samples_ref$df[2]), samples_ref$clean_sample[2])
cp_am <- comp_plot(get(samples_ref$df[3]), samples_ref$clean_sample[3])
cp_nm <- comp_plot(get(samples_ref$df[4]), samples_ref$clean_sample[4])
cp_e <- comp_plot(get(samples_ref$df[5]), samples_ref$clean_sample[5])

ggsave(here("outputs", "plots", "others", "fig1_full_samp_dist.png"), cp_f,
       width = 7, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_nonmatch_samp_dist.png"), cp_nm,
       width = 7, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_matched_samp_dist.png"), cp_m,
       width = 7, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "others", "fig1_full_samp_dist.svg"), cp_f,
       width = 7, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_nonmatch_samp_dist.svg"), cp_nm,
       width = 7, height = 10, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_matched_samp_dist.svg"), cp_m,
       width = 7, height = 10, bg = "white", dpi = 600)

# full vs non-matched -----
# Plot to compare the age and tiv distribution in the full and not matched samples
fvsnm_age <- 
  full |> 
  mutate(sample = "Full sample") |> 
  bind_rows(
    random |> 
      mutate(sample = "Not matched") 
  ) |> 
  ggplot(aes(x = age_months, group = sample, color = sample)) +
  geom_density(
    aes(
      linetype = sample, linewidth = sample),
    alpha = 0, position = "identity"
  ) +
  labs(
    subtitle = "Age distribution",
    x = "Age",
    y = "Density"
  ) +
  scale_color_manual(values = colors_samples,
                     aesthetics = c("color", "fill")) +
  scale_x_continuous(limits = c(530, 990), breaks = seq(540, 960, 60),
                     labels = seq(45, 80, 5)) +
  scale_y_continuous(labels = percent) +
  scale_linewidth_manual(values = c(1, 1.5)) +
  plots_theme +
  theme(legend.position = "none",
        plot.subtitle = element_text(face = "bold", size = 12, color = "#595959")
  ) 

fvsnm_tiv <- 
  full |> 
  mutate(sample = "Full sample") |> 
  bind_rows(
    random |> 
      mutate(sample = "Not matched") 
  ) |> 
  ggplot(aes(x = tiv, group = sample, color = sample, fill = sample)) +
  
  geom_density(
    aes(#y = after_stat(ndensity), 
      linetype = sample, linewidth = sample),
    alpha = 0, position = "identity"
  ) +
  labs(
    subtitle = "TIV distribution",
    x = "Total intracranial volume (cubic mm)",
    y = "Density"
  ) +
  scale_color_manual(values = colors_samples,
                     aesthetics = c("color", "fill")) +
  scale_x_continuous(labels = comma, limits = c(1000000, 2000000)) +
  scale_y_continuous(labels = percent) +
  scale_linewidth_manual(values = c(1, 1.5)) +
  plots_theme +
  theme(plot.subtitle = element_text(face = "bold", size = 12, 
                                     color = "#595959"),
  ) +
  guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                  linewidth = c(1, 1),
                                                  alpha = c(1,1))))

fsnm <- plot_grid(fvsnm_age, fvsnm_tiv, nrow = 2)

# now add the title
title <- ggdraw() + 
  draw_label(
    "Samples distribution: Full vs Not matched samples",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 16,
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 50)
  ) 

f_fsnm <- plot_grid(
  title, fsnm,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

f_fsnm

ggsave(here("outputs", "plots", "supplementary", "s1_full_vs_notmatched_samp_dist.png"), 
       f_fsnm, width = 7, height = 10, bg = "white")
ggsave(here("outputs", "plots", "supplementary", "s1_full_vs_notmatched_samp_dist.svg"), 
       f_fsnm, width = 7, height = 10, bg = "white")

# plots for tiv vs tbv ------
tiv_tbv <- function(samp, samp_name) {
  df <- 
    full |> 
    select(subject_id, sex, age_years, age_months, 
           tiv, tbv = brain_seg_not_vent) |> 
    pivot_longer(c(tiv, tbv)) |> 
    mutate(name_b = toupper(name)) 
  
  ggplot(df, aes(x = age_months, y = value)) +
    facet_wrap(~sex, ncol = 1) +
    geom_point(alpha = .2, size = .5, shape = 20,
               aes(group = name, color = name)) +
    geom_smooth(method = "lm",
                aes(group = name_b, color = name_b)) +
    labs(
      title = samp_name,
      x = "Age",
      y = "Volume (mm<sup>3</sup>)") +
    scale_color_manual(values = c("TBV" = "#178276", "TIV" = "#2e294e",
                                  "tbv" = "#1fad9d",  "tiv" = "#594f96"),
                       aesthetics = c("color", "fill")) +
    scale_x_continuous(limits = c(530, 990), breaks = seq(540, 960, 60),
                       labels = seq(45, 80, 5)) +
    scale_y_continuous(limits = c(760000, 2000000), labels = comma,
                       expand = expansion(mult = c(0, 0.05))
    ) +
    plots_theme +
    theme(axis.title.y = element_markdown(),
          legend.position = "none",
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 12)
    ) +
    annotate("segment", x=-Inf, xend=990, y=-Inf, yend=-Inf,
             color = "#000000",  linewidth = 1.2) +
    guides(color = guide_legend(override.aes = list(linetype = c(1, 1),
                                                    linewidth = c(1, 1),
                                                    alpha = c(1,1))))
  
}

tvt_f <- tiv_tbv(get(samples_ref$df[1]), samples_ref$clean_sample[1]) 
tvt_m <- tiv_tbv(get(samples_ref$df[2]), samples_ref$clean_sample[2]) 
tvt_nm <- tiv_tbv(get(samples_ref$df[4]), samples_ref$clean_sample[4]) 

ggsave(here("outputs", "plots", "others", "fig1_full_tiv_tbv.png"), tvt_f,
       width = 6, height = 6, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_matched_tiv_tbv.png"), tvt_m,
       width = 6, height = 6, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_notmatched_tiv_tbv.png"), tvt_nm,
       width = 6, height = 6, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "others", "fig1_full_tiv_tbv.svg"), tvt_f,
       width = 6, height = 6, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_matched_tiv_tbv.svg"), tvt_m,
       width = 6, height = 6, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_notmatched_tiv_tbv.svg"), tvt_nm,
       width = 6, height = 6, bg = "white", dpi = 600)

# Brain sizes ------------------------------------------------------------------
# Visualization of different measures from fs segmentation 
sample_distribution <- function(samp, samp_name, sizes) {
  df <- inner_join(samp |> 
                     select(subject_id, sex), sizes) |> 
    pivot_longer(-c(subject_id, sex))
  
  ggplot(df, aes(x = value, color = sex, fill = sex)) +
    facet_wrap(~name, scales = "free_x") + 
    geom_histogram(aes(linetype = sex),
                   # fill = NA,
                   alpha = 0,
                   # binwidth = 1,
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
}


sample_distribution(full, "full", brain_sizes)
sample_distribution(random, "non matched sample", brain_sizes)
sample_distribution(matched, "matched sample", brain_sizes)
sample_distribution(fs_matc, "fs tiv matched sample", brain_sizes)
sample_distribution(extreme, "extreme sample", brain_sizes)

