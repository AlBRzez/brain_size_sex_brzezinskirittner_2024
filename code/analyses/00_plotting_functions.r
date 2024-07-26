library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(scales)
library(ggtext)
library(stringr)
library(RColorBrewer)
library(broom)
library(cowplot)
library(car)
library(glue)
library(pheatmap)
library(ggpubr)
library(corrplot)

source(here("code", "analyses", "helpers.r"))
#### Get and clean data --------------------------------------------------------
# references
fs_guide <- read_csv(here("data", "fs_data", "fs_match_guide.csv"))


#### t tests and plotting functions --------------------------------------------
get_plot_df <- function(df, m, intercept = FALSE) {
  if(intercept) {
    out <-   
      df |> 
      filter(final_mod == m) |> 
      filter(df != "full") 
  } else {
    out <-   
      df |> 
      filter(final_mod == m) |> 
      filter(term_c != "intercept" & df != "full") 
  }
  out <- 
    out |> 
    group_by(clean_name) |> 
    mutate(h = max(estimate)) |> 
    ungroup() 
  return(out)
}

get_t <- function(data_frame) {
  out <- tibble()
  data_frame$df <- factor(data_frame$df, levels = c(
    "matched",
    "age_mat",
    "random",
    "extreme"
  ))
    
  for(i in unique(data_frame$term_c)) {
    df_tmp <- data_frame |> filter(term_c == i)
    df_t <- tidy(pairwise.t.test(df_tmp$estimate, df_tmp$df, 
                                 p.adjust.method = "fdr",
                                 paired = T
                                 )) |> 
      mutate(term_c = i)
    out <- bind_rows(out, df_t)
  }

  out <- 
    out |> 
    left_join(samples_ref, by = c("group1" = "df")) |> 
    left_join(samples_ref, by = c("group2" = "df")) |> 
    
    rename("start" = "clean_sample.x",
           "end" = "clean_sample.y") |> 
    left_join(terms_ref |> select(-term) |> distinct(), by = "term_c") |>
    left_join(data_frame |> select(term_c, h) |> distinct())  |>
    mutate(
      stars = case_when(
        p.value > 0.05 ~ " ",
        p.value < 0.05 & p.value > 0.01 ~ "*",
        p.value < 0.01 & p.value > 0.001 ~ "**",
        p.value < 0.001 ~ "***"
      )
    )
  
  return(out)
}

violins <- function(df, df_t, subt, sides = 1.1, central = 1.23, 
                    type = "free", horizontal) {
  
  orientation <- ifelse(horizontal == TRUE, 3, 1)
  
  if(type == "free") {
    df_t <- 
      df_t |> 
      mutate(y = case_when(
        start == "Extreme sample" & end == "Not matched" ~ (h * central),
        start == "Extreme sample" & end == "Age matched" ~ (h * central),
        start == "Not matched" & end == "Age matched" ~ (h * sides),
        start == "Extreme sample" & end == "TIV and age matched" ~ (h * central),
        start == "Not matched" & end == "TIV and age matched" ~ (h * sides),
        start == "Age matched" & end == "TIV and age matched" ~ (h * sides)
        
      ))
  } else {
    hm <- max(df_t$h)
    df_t <- 
      df_t |>
      mutate(y = case_when(
        start == "Extreme sample" & end == "Not matched" ~ hm + .4,
        start == "Extreme sample" & end == "Age matched" ~ hm + .4,
        start == "Not matched" & end == "Age matched" ~ hm + .2,
        start == "Extreme sample" & end == "TIV and age matched" ~ hm + .4,
        start == "Not matched" & end == "TIV and age matched" ~ hm + .2,
        start == "Age matched" & end == "TIV and age matched" ~ hm #+ .2
      ))
  }
  
  df_t <-
    df_t |>
    mutate(
      start = factor(start, levels = c("Extreme sample",
                                       "Not matched",
                                       "Age matched",
                                       "TIV and age matched",
                                       "FS e-TIV and age matched")),
      end = factor(end, levels = c("Extreme sample",
                                   "Not matched",
                                   "Age matched",
                                   "TIV and age matched",
                                   "FS e-TIV and age matched")),
    )
  
  
  
  df |> 
    mutate(
      clean_sample = factor(clean_sample, levels = samp_order),
      # clean_sample = factor(clean_sample, levels = c("Extreme sample",
      #                                                "Not matched",
      #                                                "Age matched",
      #                                                "TIV and age matched",
      #                                                "FS e-TIV and age matched")),
      clean_name = factor(clean_name, levels = c(
        "Intercept (for female)",
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
        "Euler number"
      ))
    ) |> 
    ggplot(aes(x = factor(clean_sample), y = estimate, color = clean_sample)) +
    geom_violin(scale = "width", aes(fill = clean_sample), alpha = .5) +
    geom_boxplot(width = .3, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "#595959", linetype = "dashed") +
    {
      if(type == "free") facet_wrap(~factor(clean_name), scales = "free", ncol = orientation) else facet_wrap(~factor(clean_name), ncol = orientation)
    } +
    geom_text(
      data = df_t,
      aes(label = stars, x = (as.numeric(end) - .5), y = y, color = end),
      show.legend = FALSE) +
    geom_segment(
      data = df_t,
      aes(x = start, xend = end, y = (y - .05), yend = (y - .05)),
      color = "#595959"
    ) +
    geom_segment(
      data = df_t,
      aes(x = start, xend = start, y = (y - .05), yend = (y - .1)),
      color = "#595959"
    ) +
    geom_segment(
      data = df_t,
      aes(x = end, xend = end, y = (y - .05), yend = (y - .1)),
      color = "#595959"
      
    ) +
    labs(
      title = "Model's estimates",
      subtitle = subt,
      x = "Sample",
      y = "Estimates",
      caption = "Paired t-tests: <b>\\*</b> = <i>p-value</i> < 0.05;<br><b>\\*\\*</b> = <i>p-value</i> < 0.01; <b>\\*\\*\\*</b> = <i>p-value</i> < 0.001"
    ) +
    scale_color_manual(
      values = colors_samples,
      aesthetics = c("color", "fill")
    )+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    plots_theme +
    theme(
      plot.caption = element_markdown(color = "#595959", size = 9),
      plot.subtitle = element_markdown(),
      # strip.text = element_text(size = 12),
      legend.position = "none",
      
    ) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,
             color = "#000000",  linewidth = 1.2) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,
             color = "#000000",  linewidth = 1.2)
  
}

fviolin <- function(df, batch_s, model, subt, sides = 1.1, central = 1.23, 
                    type = "free", intercept = FALSE, horizontal = TRUE) {
  mod_df <- get_plot_df(df, model, intercept = intercept)
  t_df <- get_t(mod_df)
  violins(df = mod_df, df_t = t_df, subt = subt, 
          type = type, horizontal = horizontal)
}



## Comparison/correlation plots ##############################
comp_scatter <- function(dat, models = c(), samp_df = c()) {
  if(length(samp_df) == 1) {
    tmp <- 
      dat |> 
      filter(final_mod %in% models &  df == samp_df) |> 
      select(reg, estimate, final_mod, clean_name) |>
      pivot_wider(names_from = final_mod, values_from = estimate) 
    xl <- glue("{models_guide$adj[models_guide$mod == models[1]]} ({samples_ref$clean_sample[samples_ref$df == samp_df[1]]})")
    yl <- glue("{models_guide$adj[models_guide$mod == models[2]]} ({samples_ref$clean_sample[samples_ref$df == samp_df[1]]})")
  } else {
    tmp <- 
      dat |> 
      filter(final_mod == models[[1]] &  df == samp_df[[1]]) |> 
      select(reg, estimate, final_mod, clean_name) |> 
      bind_rows(
        fs_vol |> 
          filter(final_mod == models[[2]] &  df == samp_df[[2]]) |> 
          select(reg, estimate, final_mod, clean_name) 
      ) |> 
      pivot_wider(names_from = final_mod, values_from = estimate) 
    
    xl <- glue("{models_guide$adj[models_guide$mod == models[1]]} ({samples_ref$clean_sample[samples_ref$df == samp_df[1]]})")
    yl <- glue("{models_guide$adj[models_guide$mod == models[2]]} ({samples_ref$clean_sample[samples_ref$df == samp_df[2]]})")
    if(grepl("PCP", yl)) {
      yl <- glue("PCP ({samples_ref$clean_sample[samples_ref$df == samp_df[2]]})")
    }
  }
  
  tmp <- 
    tmp |> 
    filter(!clean_name  %in% c("TIV", "Intercept")) |> 
    mutate(clean_name = factor(clean_name, levels = c(
      # "Intercept",
      "Age",
      "Sex (male)",
      "Sex (male):age",
      "Age<sup>2</sup>",
      "Sex (male):age<sup>2</sup>"
    )))
  
  correlations <- tibble(clean_name = character(), c = numeric())
  n <- 1
  for(e in unique(tmp$clean_name)) {
    sm <- tmp |> filter(clean_name == e)
    sm <- as.data.frame(sm)
    c <- cor.test(sm[,3], sm[,4])
    correlations[n, 1] <- e
    correlations[n, 2] <- c$estimate
    n <- n + 1
  }
  
  correlations <- correlations |> 
    mutate(clean_name = factor(clean_name, levels = c(
      "Intercept",
      "Age",
      "Sex (male)",
      "Sex (male):age",
      "Age<sup>2</sup>",
      "Sex (male):age<sup>2</sup>"
    )))
  
  cg <- tibble(clean_sample = names(colors_samples),
               col = colors_samples) |> left_join(samples_ref)
  
  p_a <- ggplot(tmp |> filter(clean_name == "Age"), 
                aes(x = !!sym(models[1]), y = !!sym(models[2]))) +
    geom_point(alpha = .7) +
    geom_abline(slope = 1, intercept = 0, color = "#ff4f4f") +
    geom_text(data = correlations|> filter(clean_name == "Age"),
              aes(x = -.4, y = .4, label = paste("r =", round(c, 3)))) +
    labs(
      title = "Age",
    ) +
    scale_x_continuous(limits = c(-.5, .5)) +
    scale_y_continuous(limits = c(-.5, .5)) +
    theme_light() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_textbox(color = "#595959"),
      plot.caption = element_markdown(color = "#595959"),
      plot.subtitle = element_markdown(),
      axis.title = element_blank(),
      aspect.ratio = 1
    )
  
  p_s <- ggplot(tmp |> filter(clean_name == "Sex (male)"), 
                aes(x = !!sym(models[1]), y = !!sym(models[2]))) +
    geom_point(alpha = .7) +
    geom_abline(slope = 1, intercept = 0, color = "#ff4f4f") +
    geom_text(data = correlations |> filter(clean_name == "Sex (male)"),
              aes(x = -.4, y = .4, label = paste("r =", round(c, 3)))) +
    labs(
      title = "Sex (male)",
    ) +
    scale_x_continuous(limits = c(-.5, .5)) +
    scale_y_continuous(limits = c(-.5, .5)) +
    theme_light() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_textbox(color = "#595959"),
      plot.caption = element_markdown(color = "#595959"),
      plot.subtitle = element_markdown(),
      axis.title = element_blank(),
      aspect.ratio = 1
    )
  
  p_sa <- ggplot(tmp |> filter(clean_name == "Sex (male):age"), 
                 aes(x = !!sym(models[1]), y = !!sym(models[2]))) +
    geom_point(alpha = .7) +
    geom_abline(slope = 1, intercept = 0, color = "#ff4f4f") +
    geom_text(data = correlations|> filter(clean_name == "Sex (male):age"),
              aes(x = -.075, y = .075, label = paste("r =", round(c, 3)))) +
    labs(
      title = "Sex (male):age",
    ) +
    scale_x_continuous(limits = c(-.1, .1)) +
    scale_y_continuous(limits = c(-.1, .1)) +
    theme_light() +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      strip.background = element_rect(fill = NA),
      strip.text = element_textbox(color = "#595959"),
      plot.caption = element_markdown(color = "#595959"),
      plot.subtitle = element_markdown(),
      axis.title = element_blank(),
      aspect.ratio = 1
    )
  
  b <- ifelse(grepl("voi", models[1]) & samp_df[1] == "matched",
              "Gold standard", xl)
  
  figure <- ggarrange(p_a, p_s, p_sa + font("x.text", size = 10),
                      ncol = 3, nrow = 1)
  annotate_figure(figure,
                  top = text_grob("Estimates comparison between models", 
                                  size = 14,
                                  hjust = 1.2),
                  # bottom = xl,
                  # bottom = "Gold standard",
                  bottom = b, 
                  left = text_grob(yl, rot = 90),
                  
  )
  
}

correlations_heatmap <- function(dat, samp_s) {
  
  tmp <- 
    dat |> 
    filter(df == samp_s) |> 
    mutate(mod = gsub("\\.y", "", mod)) |> 
    left_join(terms_ref) |> 
    left_join(models_guide) |> 
    mutate(clean_name = factor(clean_name, levels = rev(c(
      "Intercept (for female)",
      "Age",
      "Sex (male)",
      "Sex (male):age",
      "Age<sup>2</sup>",
      "Sex (male):age<sup>2</sup>"
    ))),
    adj = factor(adj, levels = c("Raw values",
                                 "Raw values - TIV as covariate",
                                 "Proportions",
                                 "Proportions - TIV as covariate",
                                 "Residuals",
                                 "Residuals - TIV as covariate",
                                 "Power-corrected proportions (PCP)",
                                 "Power-corrected proportions (PCP) - TIV as covariate"                
    ))
    ) 
  
  ggplot(tmp, aes(x = adj, y = clean_name, fill = c)) +
    geom_tile() +
    geom_text(aes(label = round(c, 3))) +
    scale_fill_gradient2(
      low = "#005f76",
      mid = "white",
      high = "#ae2012",
      midpoint = .5,
      limits = c(0, 1)
    )  +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    labs(
      title = glue("Correlation between the estimates in the matched sample and the adjusted volumes in the {tolower(samples_ref$clean_sample[samples_ref$df == samp_s])} sample"),
      x = "Method",
      y = "Estimate",
      fill = "Correlation"
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_markdown(),
      axis.ticks = element_blank()
    )
  
}


## Intramodel estimates balance ################################################


intramodel_comp <- function(df, samp_df, deltas_order_df) {
  tmp <-   
    df |>
    filter(df == samp_df) |>
    select(df, final_mod, reg, main_effects, interactions) |> 
    pivot_longer(c(main_effects, interactions)) |> 
    mutate(
      final = glue("{name}_{final_mod}"),
    ) |> 
    left_join(deltas_order_df)
  
  ggplot(tmp, aes(x = order, y = factor(reg))) +
    geom_tile(aes(fill = value)) +
    geom_vline(xintercept = seq(.5, 16.5, 2)) +
    geom_text(aes(label = round(value, 2)), size = 3) +
    labs(
      title = "Within model estimates comparisons",
      subtitle = paste("Sample:", samp_df),
      x = "",
      y = "Region",
      caption = "Main effects = age + age<sup>2</sup> | Interactions = age:sex + age<sup>2</sup>:sex"
    ) +
    scale_fill_gradient2(
      low = "#005f76",
      mid = "white",
      high = "#ae2012"
    ) +
    scale_x_continuous(labels = deltas_order_df$c_n,
                       breaks = deltas_order_df$order,
                       limits = c(.5, 16.5)) +
    theme_classic() +
    theme(
      plot.caption = element_markdown()
    )
}

# intramodel_comp(intramodel_delta, "matched", deltas_order_df)
# intramodel_comp(intramodel_delta, "unmatch", deltas_order_df)


## Intermodel comparison #######################################################

# matched_voi <- 
#   fs_deltas |> 
#   filter(final_mod %in% c("lin_voi_regular", "sq_voi_regular") & df == "matched") |> 
#   select(mod, y_type, batch, term_c, estimate, reg, type) |> 
#   pivot_wider(names_from = term_c, values_from = estimate) 
# 
# combinations <- 
#   expand_grid(
#     y_type = c("voi", "prop", "pcp", "res"),
#     type = c("lin_int", "sq_int"),
#     batch = c("regular", "tiv"),
#     sample = c("full", "unmatch", "age_mat", "random", "matched")
#   )
# 

comparing_models <- function(gold_std, est_data, combinations_df, r_row) {
  t <- combinations_df$type[r_row]
  b <- combinations_df$batch[r_row]
  s <- combinations_df$sample[r_row]
  y <- combinations_df$y_type[r_row]
  
  gold <- 
    gold_std |> 
    filter(type == t)
  
  tmp1 <- 
    fs_deltas |> 
    filter(y_type == y & type == t & batch == b & df == s) |> 
    select(mod, y_type, batch, term_c, estimate, reg, type) |> 
    pivot_wider(names_from = term_c, values_from = estimate)
  
  out <- tibble(.rows = 78)
  out$batch <- b
  out$sample <- s
  out$type <- t
  out$y_type <- y
  out$reg <- gold$reg
  
  # for all
  out$age <- gold$age - tmp1$age
  out$sex_male <- gold$sex_male - tmp1$sex_male
  out$sex_male_age <- gold$sex_male_age - tmp1$sex_male_age
  
  # for quadratic interaction models
  if(t == "sq_int") {
    out$age_sq <- gold$age_sq - tmp1$age_sq
    out$sex_male_age_sq <- gold$sex_male_age_sq - tmp1$sex_male_age_sq
  }
  
  return(out)
}


# final_comp <- tibble()
# for(i in 1:nrow(combinations)) {
#   o <- comparing_models(matched_voi, fs_deltas, combinations, i)
#   final_comp <- bind_rows(final_comp, o)
# }

clean_y_ref <- tibble(
  y_type = c("voi", "pcp", "res", "prop"),
  clean_y = c("Raw values", "PCP", "Residuals", "Proportions")
)

comp_viz <- function(comp_df, samp_df, batch_s, type_s, num = FALSE) {
  tmp <- 
    comp_df |> 
    filter(sample == samp_df & batch == batch_s & type == type_s) |>  
    pivot_longer(-c(batch, sample, type, y_type, reg), names_to = "term_c") |>
    filter(!is.na(value)) |> 
    left_join(clean_y_ref) |> 
    left_join(terms_ref)
  
  p <- 
    ggplot(tmp) +
    aes(x = clean_name, y = factor(reg)) +
    facet_wrap(~clean_y) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(
      low = "#005f76",
      mid = "white",
      high = "#ae2012",
      limits = c(-1, 1)
    ) +
    labs(
      title = "Inter-model comparison",
      # subtitle = paste("Raw values in the matched sample -", samp_df, "(", batch_s ,")"),
      subtitle = paste("Raw values in the matched sample -", 
                       samples_ref$clean_sample[samples_ref$df == samp_df], "sample"),
      x = "Estimate",
      y = "Regions",
      fill = "Diference",
      caption = "Difference in the model estimates of the age and TIV matched sample without adjustment 
      versus the estimates for the other sample with each one of the adjustment methods"
    ) +
    # scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_markdown(size = 7)
    )
  
  if(num) {
    p <- p + geom_text(aes(label = round(value, 2)), size = 3) 
  }
  return(p)
}

