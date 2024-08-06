## Set up ----------------------------------------------------------------------
library(here)
source(here("code", "analyses", "00_plotting_functions.r"))
library(ggcorrplot)
library(patchwork)
lm_data_folder <- "trajectories/"
allometry_data_folder <- "allometry/"

## Get data --------------------------------------------------------------------
lm_coef_o <- 
  read_csv(here("outputs", lm_data_folder, "fs_coeff.csv"))

allom_lm <- 
  read_csv(here("outputs", allometry_data_folder, "lm_over_betas_coef_tbv.csv")) 

## Wrangling -------------------------------------------------------------------
estimates_order <- c("Intercept (for female)",
                     "Age",
                     "Age (months)",
                     "Sex (male)",
                     "Sex (male):age",
                     "Sex (male):age (months)",
                     "Age<sup>2</sup>",
                     "Age<sup>2</sup> (months)",
                     "Sex (male):age<sup>2</sup>",
                     "Sex (male):age<sup>2</sup> (months)",
                     "TIV",
                     "TBV",
                     "Euler number")


lm_coef <- 
  lm_coef_o |> 
  filter(effect == "fixed") |> 
  filter(df != "fs_matc") |> 
  inner_join(terms_ref) |> 
  left_join(samples_ref)  |>
  left_join(fs_guide) |> 
  filter(ignore == 0 & type == "volume") |> 
  group_by(df, mod, term_c) |> 
  nest() |>
  mutate(
    p_fdrc = map(data, function(df)  p.adjust(df$p.value, method = "fdr")),
  ) |>
  unnest(cols = c(data, p_fdrc)) |>
  ungroup() |> 
  mutate(
    p_corr_fdrc = ifelse(p_fdrc <= 0.05, 1, 0),
  ) |> 
  rename("final_mod" = "mod")

allom <- 
  allom_lm |> 
  left_join(terms_ref) |> 
  left_join(samples_ref) |>
  left_join(fs_guide) |> 
  filter(ignore == 0) |> 
  group_by(df, mod, term_c, window_size) |> 
  nest() |>
  mutate(
    p_fdrc = map(data, function(df)  p.adjust(df$p.value, method = "fdr")),
  ) |>
  unnest(cols = c(data, p_fdrc)) |>
  ungroup() |> 
  mutate(
    p_corr_fdrc = ifelse(p_fdrc <= 0.05, 1, 0),
  ) |> 
  rename("final_mod" = "mod")

models_guide <- 
  bind_rows(
    tibble(mod = unique(lm_coef$final_mod)) |> 
      mutate(
        type = ifelse(grepl("lin", mod), 
                      "y ~ 1 + age * sex + log(euler) + (1|site)",
                      "y ~ 1 + (age + age<sup>2</sup>) * sex + log(euler) + (1|site)"),
        type = ifelse(grepl("euler", mod), 
                      type,
                      gsub(" \\+ log\\(euler\\) ", "", type)),
        analysis = "aging trajectories"
      ),
    expand.grid(mod = unique(allom$final_mod),
                window_size = 60) |>
      mutate(
        type = ifelse(grepl("lin", mod),
                      "y ~ 1 + age * sex",
                      "y ~ 1 + (age + age<sup>2</sup>) * sex"),
        analysis = "allometry")
  ) |> tibble()

# Fig 1 ------------------------------------------------------------------------
violin_fig1_v1 <- 
  fviolin(lm_coef, model = as.character(models_guide[1, 1]),
          subt = glue("{models_guide[1, 2]}: {models_guide[1, 3]}"),
          type = "fix", intercept = T, horizontal = F)

ggsave(here("outputs", "plots", "main_figures", "fig1_violin_v1_months.png"), 
       violin_fig1_v1, width = 3.7, height = 12, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig1_violin_v1_months.svg"), 
       violin_fig1_v1, width = 3.7, height = 12, bg = "white", dpi = 600)

# Fig 3 -----------------------------------------------------------------------
mod_df_allom <- get_plot_df(allom |> filter(window_size == 60), 
                            m = "lin_int_regular",
                            intercept = T)

violin_dat_allom <- 
  mod_df_allom |>
  mutate(
    clean_sample = factor(clean_sample, levels = rev(samples_ref$clean_sample)),
    clean_name = factor(clean_name, levels = estimates_order))

p_intercept <- 
  ggplot(violin_dat_allom |> filter(term_c == "intercept"),
         aes(x = factor(clean_sample), y = estimate, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 1,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Intercept (for female)",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(.49, 1.5)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 14, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title = element_blank()
    
  ) 

p_agem <- 
  ggplot(violin_dat_allom |> filter(term_c == "agem"),
         aes(x = factor(clean_sample), y = estimate, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  labs(
    subtitle = "Age (months)",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(-.2, .2)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 14, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title = element_blank()
  ) 

p_sex <- 
  ggplot(violin_dat_allom |> filter(term_c == "sex_male"),
         aes(x = factor(clean_sample), y = estimate, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Sex (male)",
    x = "",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(-.2, .2)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 14, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.title = element_blank()
  ) 

p_sex_agem <- 
  ggplot(violin_dat_allom |> filter(term_c == "sex_male_agem"),
         aes(x = factor(clean_sample), y = estimate, 
             color = clean_sample)) +
  geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
  geom_boxplot(width = .3, show.legend = FALSE) +
  geom_hline(yintercept = 0,
             color = "#595959", linetype = "dashed") +
  
  labs(
    subtitle = "Sex (male):age (months)",
    x = "Sample",
    y = "Estimates",
  ) +
  scale_color_manual(
    values = colors_samples,
    aesthetics = c("color", "fill")
  )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(limits = c(-.2, .2)) +
  plots_theme +
  theme(
    plot.subtitle = element_text(size = 14, 
                                 face = "bold", color = "#595959"),
    legend.position = "none",
    axis.title = element_blank()
  ) 


allom_viol <- plot_grid(p_intercept, p_agem, p_sex, p_sex_agem, ncol = 2,
                        rel_heights = c(.9, 1))

# now add the title
title <- ggdraw() + 
  draw_label(
    "Allometry models estimates",
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

f_allom_viol <- plot_grid(
  title, allom_viol,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave(here("outputs", "plots", "main_figures", "fig3_violins_v1.png"), 
       f_allom_viol, width = 8, height = 8, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig3_violins_v1.svg"), 
       f_allom_viol, width = 8, height = 8, bg = "white", dpi = 600)

#### Scatter ---------

scat_df <- 
  mod_df_allom |> 
  filter(df %in% c("random", "matched")) |> 
  select(clean_sample, term_c, estimate, reg, clean_name) |> 
  pivot_wider(names_from = clean_sample, values_from = estimate) 

allom_cor_df <- 
  scat_df |> 
  group_by(term_c) |> 
  summarise(r = cor(`TIV and age matched`, `Not matched`))


s_intercept <- 
  ggplot(scat_df |> filter(term_c == "intercept"),
         aes(x = `Not matched`, y = `TIV and age matched`)) +
  geom_hline(yintercept = 1, color = "#595959", linetype = 2) +
  geom_vline(xintercept = 1, color = "#595959", linetype = 2) +
  geom_point(alpha = .4) +
  geom_abline(color = "#ff4f4f") + 
  geom_richtext(data = allom_cor_df |> filter(term_c == "intercept"), 
                aes(label = glue("*r* = {comma(r, accuracy = .01)}")),
                x = .75, y = 1.25,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(
    subtitle = "Intercept (for female)",
    # x = ""
  ) +
  scale_x_continuous(limits = c(.5, 1.5)) +
  scale_y_continuous(limits = c(.5, 1.5)) +
  plots_theme +
  theme(aspect.ratio = 1,
        plot.subtitle = element_text(size = 14, 
                                     face = "bold", color = "#595959"))

s_agem <- 
  ggplot(scat_df |> filter(term_c == "agem"),
         aes(x = `Not matched`, y = `TIV and age matched`)) +
  geom_hline(yintercept = 0, color = "#595959", linetype = 2) +
  geom_vline(xintercept = 0, color = "#595959", linetype = 2) +
  geom_point(alpha = .4) +
  geom_abline(color = "#ff4f4f") +  
  geom_richtext(data = allom_cor_df |> filter(term_c == "agem"), 
                aes(label = glue("*r* = {comma(r, accuracy = .01)}")),
                x = -.1, y = .1,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(
    subtitle = "Age (months)",
    # x = "",
    # y = ""
  ) +
  scale_x_continuous(limits = c(-.2, .2)) +
  scale_y_continuous(limits = c(-.2, .2)) +
  plots_theme +
  theme(aspect.ratio = 1,
        plot.subtitle = element_text(size = 14, 
                                     face = "bold", color = "#595959"))

s_sex <- 
  ggplot(scat_df |> filter(term_c == "sex_male"),
         aes(x = `Not matched`, y = `TIV and age matched`)) +
  geom_hline(yintercept = 0, color = "#595959", linetype = 2) +
  geom_vline(xintercept = 0, color = "#595959", linetype = 2) +
  geom_point(alpha = .4) +
  geom_abline(color = "#ff4f4f") +  
  geom_richtext(data = allom_cor_df |> filter(term_c == "sex_male"), 
                aes(label = glue("*r* = {comma(r, accuracy = .01)}")),
                x = -.1, y = .1,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(
    subtitle = "Sex (male)",
  ) +
  scale_x_continuous(limits = c(-.2, .2)) +
  scale_y_continuous(limits = c(-.2, .2)) +
  plots_theme +
  theme(aspect.ratio = 1,
        plot.subtitle = element_text(size = 14, 
                                     face = "bold", color = "#595959"))

s_sex_agem <- 
  ggplot(scat_df |> filter(term_c == "sex_male_agem"),
         aes(x = `TIV and age matched`, y = `Not matched`)) +
  geom_hline(yintercept = 0, color = "#595959", linetype = 2) +
  geom_vline(xintercept = 0, color = "#595959", linetype = 2) +
  geom_point(alpha = .4) +
  geom_abline(color = "#ff4f4f") + 
  geom_richtext(data = allom_cor_df |> filter(term_c == "sex_male_agem"), 
                aes(label = glue("*r* = {comma(r, accuracy = .01)}")),
                x = -.1, y = .1,
                fill = NA, label.color = NA, # remove background and outline
                label.padding = grid::unit(rep(0, 4), "pt")) +
  labs(
    subtitle = "Sex (male):Age (months)",
    # y = ""
  ) +
  scale_x_continuous(limits = c(-.2, .2)) +
  scale_y_continuous(limits = c(-.2, .2)) +
  plots_theme +
  theme(aspect.ratio = 1,
        plot.subtitle = element_text(size = 14, 
                                     face = "bold", color = "#595959"))  

p2 <- plot_grid(s_intercept + theme(axis.title.x = element_blank()), 
                s_agem + theme(axis.title = element_blank()), 
                s_sex, 
                s_sex_agem + theme(axis.title.y = element_blank()), 
                ncol = 2)


ggsave(here("outputs", "plots", "main_figures", "fig3_scatter_v1.png"), 
       p2, width = 8, height = 8, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig3_scatter_v1.svg"), 
       p2, width = 8, height = 8, bg = "white", dpi = 600)

p3 <- plot_grid(s_intercept, 
                s_agem, #+ theme(axis.title.y = element_blank()), 
                s_sex, # + theme(axis.title.y = element_blank()), 
                s_sex_agem, # + theme(axis.title.y = element_blank()), 
                nrow = 1)


ggsave(here("outputs", "plots", "main_figures", "fig3_scatter_v2.png"), 
       p3, width = 16, height = 4, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig3_scatter_v2.svg"), 
       p3, width = 16, height = 4, bg = "white", dpi = 600)

#### Allom samp corrs ----------------------------------------------------------
heat_df <- 
  mod_df_allom |> 
  select(clean_sample, term_c, estimate, reg) 

make_corr_plot <- function(dat, term, cp_t) {
  
  tmp <- 
    dat |> 
    filter(term_c == term)  |> 
    select(-c(term_c)) |> 
    pivot_wider(names_from = clean_sample, values_from = estimate) |>
    arrange(reg) |> 
    select(-reg) |> 
    select(`Extreme sample`, 
           `Not matched`,
           `Age matched`,
           `TIV and age matched`) |> 
    as.matrix()
  
  corr_matrix <- cor(tmp)
  
  pc <- ggcorrplot(corr_matrix,
                   type = "lower",
                   show.diag = T,
                   outline.col = "white",
                   ggtheme = theme_void,
                   colors = c("#d53e4f", "#FFFFBF", "#3288bd"),
                   lab = T) +
    labs(
      title = cp_t,
    ) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 12)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 12)) +
    theme(
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(hjust = 1, size = 9),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(size = 14, 
                                face = "bold", color = "#595959")
    )
  
  return(pc)
}

cp_i <- make_corr_plot(heat_df, "intercept", "Intercept (for females)")
cp_s <- make_corr_plot(heat_df, "sex_male", "Sex (male)")
cp_a <- make_corr_plot(heat_df, "agem", "Age (months)")
cp_sa <- make_corr_plot(heat_df, "sex_male_agem", "Sex (male):age (months)")


ggarrange(cp_i, cp_a, cp_s, cp_sa, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")

est_cors_plot <- ggarrange(cp_i, cp_a, cp_s, cp_sa, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")#plot_grid(cp_i, cp_a, cp_s, cp_sa)
est_cors_plot <- est_cors_plot + plot_annotation(title = "Correlations of the model's estimates between samples",
                                                 theme = theme(plot.title = element_markdown(size = 16, face = "bold")))

ggsave(here("outputs", "plots", "supplementary", "s10_allom_est_corr.png"), 
       est_cors_plot, width = 6.5, height = 6.5, bg = "white")
ggsave(here("outputs", "plots", "supplementary", "s10_allom_est_corr.svg"), 
       est_cors_plot, width = 6.5, height = 6.5, bg = "white")

#####################
allom |>
  filter(df %in% c("matched", "random")) |>
  filter(final_mod == "lin_int_regular" & window_size == 60 & term_c == "intercept") |>
  mutate(sign = ifelse(estimate < 1, "-", "+")) |> 
  select(df, estimate, sign, Field) |> 
  pivot_wider(names_from = df, values_from = c(estimate, sign))  |>
  View()


