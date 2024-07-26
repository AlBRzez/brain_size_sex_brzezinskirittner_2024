source("../ukbb_data_r/functions.R")
library(ggplot2)
library(here)
library(readr)
library(tidyr)
library(scales)
library(rstatix)

## Code for getting sociodemographic (education level and income) data from
## UKBB and make inter-samples comparissons

## Get data --------------------------------------------------------------------
ukbb_path <- "/data/zeiyas/in_vivo_Datasets/UKbiobank/brzali"
ukbb <- open_dataset(paste0(ukbb_path, "/UKBB-tabular-processing/current_melt.arrow"), format = "ipc")
codings <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Codings.tsv"))
diction <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Data_Dictionary_Showcase.tsv"))

source(here("code", "analyses", "helpers.r"))

dfs <- get_samples("data/samples/", fs = T)
list2env(dfs, envir = .GlobalEnv)
df_n <- names(dfs)
full_j <- full |> select(subject_id, age_years, age_months, sex, tiv)


# UKBB data extraction ---------------------------------------------------------
sec_long <- get_ukbb_data(arrow = ukbb, diction = diction, 
                          fields_id = c(738, 6138, 845), subjs_id = full$subject_id,
                          instances = 0:2)

sec <- order_data(sec_long, diction, codings)

sec_long2 <- get_ukbb_data(arrow = ukbb, diction = diction, 
                           fields_id = c(738, 6138, 845), subjs_id = full$subject_id,
                           instances = 2)

sec2 <- order_data(sec_long2, diction, codings)


sec$rowid <- 1:nrow(sec)
sec$field_value <- as.numeric(sec$field_value)
sec2$rowid <- 1:nrow(sec2)
sec2$field_value <- as.numeric(sec2$field_value)

# Wrangling data ---------------------------------------------------------------
#### reclasification of educational levels -------------------------------------

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
                                     "Some secondary"))

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

#### income variable -----------------------------------------------------------
income <- 
  sec2 |> 
  filter(field_id == 738) |> 
  select(subject_id, income_level = field_value, income_meaning = meaning)


#### Joining data --------------------------------------------------------------
final_sec <- left_join(education, income)

# Looking at differences -----------------------------------------------
# Checking if there are significant differences between the subjects in the
# matched sample and the ones that are not included

# These are the reported values
full_sec <- 
  left_join(full_j, final_sec) |> 
  id_samples() 

tabyl_education <- tabyl(full_sec, matched, education_meaning)
chisq.test(tabyl_education)

tabyl_income <- tabyl(full_sec, matched, income_meaning)
chisq.test(tabyl_income)

# other checks -----------------------------------------------------------------

# Sex in the matched sample
matched_sec <- left_join(matched, final_sec)
tmp_m <- tabyl(matched_sec, sex, education_meaning)
chisq.test(tmp_m)
ggplot(full_sec, aes(y = age_months, group = interaction(sex, matched), 
                     color =  interaction(sex, matched))) +
  geom_boxplot() +
  theme_light()

ggplot(full_sec, aes(y = age_months, group = matched, 
                     color =  factor(matched))) +
  facet_wrap(~sex) +
  geom_boxplot() +
  theme_light()


e <- full_sec |> filter(education_level > 2)

ggplot(e, aes(y = tiv, color = factor(education_level))) +
  geom_boxplot() +
  facet_wrap(~sex)

# Sex in the full sample
full_sec |> 
  filter(education_level > 0) |> 
  count(sex, education_level) |> 
  group_by(sex) |> 
  mutate(p = n/sum(n)) |> 
  ungroup() |> 
  left_join(edu_levels) |> 
  ggplot(aes(x = reorder(education_meaning, education_level), y = p, color =)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = sex)) +
  geom_text(aes(label = percent(p, accuracy = .1)), 
            position = position_dodge(width = .9), vjust = -.2, size = 3) +
  labs(title = "Educational level - maximum qualifications",
       x = "Educational level",
       y = "Percentage by sex") +
  scale_x_discrete(labels = label_wrap_gen(15)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = color_sex,
                     aesthetics = c("color", "fill")) +
  theme_light()


full_sec |> 
  filter(education_level > 0) |> 
  count(sex, matched, education_level) |> 
  group_by(sex, matched) |> 
  mutate(p = n/sum(n)) |> 
  ungroup() |> 
  left_join(edu_levels) |> 
  ggplot(aes(x = reorder(education_meaning, education_level), y = p)) +
  # facet_grid(~matched) +
  facet_grid(sex~matched) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = sex)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = color_sex,
                     aesthetics = c("color", "fill")) +
  theme_light()

full_sec |> 
  filter(income_level > 0) |> 
  count(sex, matched, income_level) |> 
  group_by(sex, matched) |> 
  mutate(p = n/sum(n),
         Matched = ifelse(matched == 1, "yes", "no"
         )) |> 
  ungroup() |> 
  left_join(diction |> 
              filter(FieldID == 738) |> 
              left_join(codings) |> 
              mutate(income_level = as.numeric(Value)) |> 
              select(income_level, income_meaning = Meaning) 
  ) |>
  ggplot(aes(x = reorder(income_meaning, income_level), y = p, color = sex)) +
  facet_grid(~Matched,  labeller = label_both) +
  # facet_grid(sex~matched) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = sex)) +
  geom_text(aes(label = percent(p, accuracy = .1)), 
            position = position_dodge(width = .9), vjust = -.2, size = 3) +
  labs(title = "Average total household income before tax",
       x = "Income",
       y = "Percentage inside group") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = color_sex,
                     aesthetics = c("color", "fill")) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(color = "#595959")
  )

full_sec |> 
  filter(education_level > 0) |> 
  count(sex, matched, education_level) |> 
  group_by(sex, matched) |> 
  mutate(p = n/sum(n),
         Matched = ifelse(matched == 1, "yes", "no"
         )) |> 
  ungroup() |> 
  left_join(edu_levels
  ) |>
  ggplot(aes(x = reorder(education_meaning, education_level), y = p, color = sex)) +
  facet_grid(~Matched,  labeller = label_both) +
  # facet_grid(sex~matched) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = sex)) +
  geom_text(aes(label = percent(p, accuracy = .1)), 
            position = position_dodge(width = .9), vjust = -.2, size = 3) +
  labs(title = "Educational level - maximum qualifications",
       x = "Educational level",
       y = "Percentage inside group") +
  scale_x_discrete(labels = label_wrap_gen(15)) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = color_sex,
                     aesthetics = c("color", "fill")) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(color = "#595959")
  )

ggplot() +
  geom_density(data = full, aes(x = age_months)) +
  geom_density(data = full_sec, aes(x = age_months, color = factor(matched))) +
  facet_wrap(~sex) +
  theme_light()

ggplot() +
  # geom_boxplot(data = full, aes(y = age_months )) +
  geom_boxplot(data = full_sec, aes(y = age_months, color = factor(matched))) +
  facet_wrap(~sex) +
  theme_light()
