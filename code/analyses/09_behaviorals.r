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
library(ggtext)
library(MetBrewer)
library(car)

## Code for getting health and cognitive (tests) data from
## UKBB and analyse sex differences in the different samples


## Get data --------------------------------------------------------------------
## Get ukbb data
ukbb_path <- "/data/zeiyas/in_vivo_Datasets/UKbiobank/brzali"
ukbb <- open_dataset(paste0(ukbb_path, "/UKBB-tabular-processing/current_melt.arrow"), format = "ipc")
codings <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Codings.tsv"))
diction <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Data_Dictionary_Showcase.tsv"))

source(here("code", "analyses", "helpers.r"))

dfs <- get_samples("data/samples/", fs = T)
list2env(dfs, envir = .GlobalEnv)
df_n <- names(dfs)
full_j <- full |> select(subject_id, age_years, age_months, sex, tiv)

# extraction -----
health_long <- get_ukbb_data(arrow = ukbb, diction = diction, 
                             fields_id = c(20002), 
                             subjs_id = full$subject_id,
                             instances = 2)

health <- order_data(health_long, diction, codings)
rm(health_long)

beh_long <-  get_ukbb_data(arrow = ukbb, diction = diction, 
                           fields_id = c(20116, 21001, 20016), 
                           subjs_id = full$subject_id,
                           instances = 2)
beh <- order_data(beh_long, diction, codings)
rm(beh_long)

url <- "https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100026"
url_list <- url %>% 
  read_html() %>% 
  html_nodes("table") %>% 
  html_table(fill = T) 

cognitive_fields_guide <- url_list[[1]] |> janitor::clean_names()

cog_long <- get_ukbb_data(arrow = ukbb, diction = diction, 
                          categories = cognitive_fields_guide$category_id,
                          subjs_id = full$subject_id,
                          instances = 2)

cognitive <- order_data(cog_long, diction, codings)
rm(cog_long)
gc()
# save(health, beh, cognitive, file = "our_data/data/behaviorals.Rdata")
# load(file = "our_data/data/behaviorals.Rdata")

# Health information -----------------------------------------------------------
glimpse(health)

health_guide <- 
  health |> 
  select(field_value, value) |> 
  distinct()

diseases <- c(1286, 1220, 1521, 1276, 1607, 1468, 1072, 1473, 1065, 1225, 
              1226, 1616, 1154, 1583, 1291, 1265, 1261, 1262, 1243, 1289,
              1222, 1223, 1263, 1287, 1081)

# All the diseases that might be relevant
dis <- 
  health |> 
  mutate(field_value = as.numeric(field_value)) |> 
  filter(field_value %in% diseases)

# grouping the four diseases that we are keeping
selected_dis <-
  dis |> 
  mutate(
    disease = case_when(
      value %in% c("diabetes", "type 1 diabetes", "type 2 diabetes", "diabetes insipidus") ~ "diabetes",
      value %in% c("hypertension", "essential hypertension") ~ "hypertension",
      value %in% c("stroke", "ischaemic stroke") ~ "stroke",
      value == "high cholesterol" ~ "high_cholesterol"
      
    )
  ) |> 
  filter(!is.na(disease)) |> 
  select(-c(field_id, instance_id, array_id, field_value, category, field, 
            value_type, coding, meaning, value)) |> 
  distinct() |> 
  mutate(a = 1) |> 
  pivot_wider(names_from = disease, values_from = a) 

# joining info with the subjects base info and samples
diseases_samp <- 
  full_j |> 
  left_join(
    selected_dis
  ) |> 
  id_samples() |> 
  replace_na(list(
    diabetes = 0,
    hypertension = 0,
    stroke = 0,
    high_cholesterol = 0
  ))

tyble_diabetes <- tabyl(diseases_samp, matched, diabetes)
chisq.test(tyble_diabetes)

tyble_hypertension <- tabyl(diseases_samp, matched, hypertension)
chisq.test(tyble_hypertension)

tyble_stroke <- tabyl(diseases_samp, matched, stroke)
chisq.test(tyble_stroke)

tyble_high_cholesterol <- tabyl(diseases_samp, matched, high_cholesterol)
chisq.test(tyble_high_cholesterol)

#### smoking and BMI -----------------------------------------------------------
smoke <- 
  full_j |> 
  left_join(
    beh |> 
      select(-c(instance_id, array_id, field_value, category, 
                value_type, coding, meaning)) |> 
      filter(field_id == 20116)
  ) |> 
  id_samples() 

tyble_smoke <- tabyl(smoke, matched, value)
chisq.test(tyble_smoke)


bmi <- 
  full_j |> 
  left_join(
    beh |> 
      select(-c(instance_id, array_id, field_value, category, 
                value_type, coding, meaning)) |> 
      filter(field_id == 21001) |> 
      mutate(bmi = as.numeric(value))
  ) |> 
  id_samples() 

t.test(bmi$bmi ~ bmi$matched)
cohens_d(bmi ~ matched, data = bmi)

# Cognitive tests --------------------------------------------------------------
cognitive <- 
  cognitive |> 
  left_join(cognitive_fields_guide |> 
              select(category = category_id, 
                     cat_desc = description)) 

samp_cog <-
  left_join(cognitive, full_j) |> 
  mutate(val = as.numeric(value)) |> 
  id_samples()

extended_cog_guide <- 
  cognitive |> 
  select(field_id, field, category, cat_desc) |> 
  distinct()

#### Tests with all the variables before cleaning ----
t_tests <- 
  samp_cog |> 
  drop_na(val) |> 
  group_by(category, field_id) |> 
  nest() |> 
  mutate(t_test = map(data, ~ t.test(.x$val ~ factor(.x$sex), var.eq = F, paired = F)),
         res = map(t_test, tidy)) |> 
  unnest(res) |>
  ungroup() |> 
  select(-c(data, t_test)) 

cohens_d_t <- 
  samp_cog |> 
  drop_na(val) |> 
  # filter(!field_id %in% c(4251, 4252, 4258, 6362, 6312)) |> 
  group_by(category, field_id) |> 
  nest() |> 
  mutate(#t_test = map(data, ~ t.test(.x$val ~ factor(.x$sex), var.eq = F, paired = F)),
    res = map(data, ~cohens_d(val ~ sex, data = .))) |> 
  unnest(res) |>
  ungroup() |> 
  select(-c(data)) 

full_cognitives <- 
  left_join(t_tests, cohens_d_t) |> 
  left_join(extended_cog_guide)


mat_cognitives <- 
  samp_cog |> 
  filter(matched == 1)

t_tests <- 
  mat_cognitives |> 
  drop_na(val) |> 
  # filter(!field_id %in% c(4251, 4252, 4258, 6362, 6312)) |> 
  group_by(category, field_id) |> 
  nest() |> 
  mutate(t_test = map(data, ~ t.test(.x$val ~ factor(.x$sex), var.eq = F, paired = F)),
         res = map(t_test, tidy)) |> 
  unnest(res) |>
  ungroup() |> 
  select(-c(data, t_test)) 

cohens_d_t <- 
  mat_cognitives |> 
  drop_na(val) |> 
  # filter(!field_id %in% c(4251, 4252, 4258, 6362, 6312)) |> 
  group_by(category, field_id) |> 
  nest() |> 
  mutate(#t_test = map(data, ~ t.test(.x$val ~ factor(.x$sex), var.eq = F, paired = F)),
    res = map(data, ~cohens_d(val ~ sex, data = .))) |> 
  unnest(res) |>
  ungroup() |> 
  select(-c(data)) 


matched_cognitives <- 
  left_join(t_tests, cohens_d_t) |> 
  left_join(extended_cog_guide)

table(full_cognitives$magnitude)
table(matched_cognitives$magnitude)



#### Clean cognitive values ---------------------------------------------------- 
keep_vars <- c("subject_id", "val", "category", "cat_desc",
               "field_id", "field")
#-- Category 501
cat_501_1 <- 
  samp_cog |> 
  filter(category == 501) |> 
  filter(field_id %in% c(6373, 6374)) |> 
  select(all_of(keep_vars))


cat_501_2 <- 
  cat_501_1 |> 
  select(subject_id, field_id, val) |> 
  pivot_wider(names_from = field_id, values_from = val) |> 
  transmute(
    subject_id = subject_id,
    val = `6373`/`6374`, 
    category = 501,
    cat_desc = "Matrix pattern completion",
    field_id = 30001,
    field = "Percent: solved / viewed puzzles"
  ) 

cat_501_3 <- 
  samp_cog |> 
  filter(field_id == 6333) |> 
  group_by(subject_id) |> 
  summarise(
    val = mean(val),
    category = 501,
    cat_desc = "Matrix pattern completion",
    field_id = 30002,
    field = "Mean Duration spent answering each puzzle"
  ) 


cat_501 <- 
  bind_rows(cat_501_1,
            cat_501_2,
            cat_501_3)

rm(cat_501_1,
   cat_501_2,
   cat_501_3)
gc()

#-- Category 502
cat_502_1 <- 
  samp_cog |> 
  filter(category == 502) |> 
  filter(field_id %in% c(23323, 23324)) |> 
  select(all_of(keep_vars))

cat_502_2 <- 
  cat_502_1 |> 
  select(subject_id, field_id, val) |> 
  pivot_wider(names_from = field_id, values_from = val) |> 
  transmute(
    subject_id = subject_id,
    val = `23324`/`23323`, 
    category = 502,
    cat_desc = "Symbol digit substitution",
    field_id = 30003,
    field = "Percent: solved / viewed matches"
  ) 

cat_502 <- 
  bind_rows(cat_502_1,
            cat_502_2)

rm(cat_502_1,
   cat_502_2)
gc()

#-- Category 503
cat_503_1 <- 
  samp_cog |> 
  filter(category == 503) |> 
  filter(field_id %in% c(21004, 6383)) |> 
  select(all_of(keep_vars))

cat_503_2 <- 
  cat_503_1 |> 
  select(subject_id, field_id, val) |> 
  pivot_wider(names_from = field_id, values_from = val) |> 
  transmute(
    subject_id = subject_id,
    val = `21004`/`6383`, 
    category = 503,
    cat_desc = "Tower rearranging",
    field_id = 30005,
    field = "Percent: solved / viewed puzzles"
  )

cat_503_3 <- 
  samp_cog |> 
  filter(field_id == 6313) |> 
  group_by(subject_id) |> 
  summarise(
    val = mean(val),
    category = 503,
    cat_desc = "Tower rearranging",
    field_id = 30006,
    field = "Mean Duration spent entering selection"
  ) 

cat_503 <- 
  bind_rows(cat_503_1,
            cat_503_2,
            cat_503_3)

rm(cat_503_1,
   cat_503_2,
   cat_503_3)
gc()

#-- Category 505
cat_505_1 <- 
  samp_cog |> 
  filter(category == 505) |> 
  filter(!field_id %in% 6770:6773) |> 
  select(all_of(keep_vars))

cat_505_2 <- 
  left_join(
    expand_grid(
      subject_id = unique(cat_501$subject_id),
      field_id = 6770:6771,
      array_id = 0:23),
    samp_cog |> 
      filter(field_id %in% 6770:6771) |> 
      select(subject_id, field_id, val, array_id) 
  ) |>
  mutate(val = ifelse(is.na(val), 0, val)) |>
  group_by(field_id, subject_id) |>
  summarise(val = mean(val),
            .groups = "drop") |>
  left_join(extended_cog_guide)

cat_505 <- 
  bind_rows(
    cat_505_1,
    cat_505_2
  )

rm(cat_505_1, cat_505_2)
gc()

#-- Category 100027
cat_1000027_ans <- 
  samp_cog |> 
  filter(category == 100027 & field_id != 4924) |> 
  mutate(val = ifelse(is.na(val), 0, 1))

cat_1000027_reg <- 
  samp_cog |> 
  #filter(category == 100027 & field_id != 4924) |> 
  filter(field_id %in% c(20016, 20128)) |>
  drop_na(val) |> 
  select(all_of(keep_vars))


#-- Category 100029
cat_100029 <- 
  samp_cog |> 
  filter(field_id == 4282) |> 
  select(all_of(keep_vars))

#-- Category 100030
cat_100030 <-
  samp_cog |> 
  filter(field_id %in% c(399, 400)) |> 
  drop_na(val) |> 
  mutate(field_id = case_when(
    array_id == 1 ~ as.numeric(glue("{field_id}.1")),
    array_id == 2 ~ as.numeric(glue("{field_id}.2")),
    array_id == 3 ~ as.numeric(glue("{field_id}.3"))
  )) |> 
  select(all_of(keep_vars))

#-- Category 100032 - the only thing we keep is "Mean time to correctly identify matches"
cat_100032 <- 
  samp_cog |> 
  filter(field_id == 20023) |> 
  select(all_of(keep_vars))

### ----
final_cog <- 
  bind_rows(
    cat_501,
    cat_502,
    cat_503,
    cat_505,
    cat_1000027_reg,
    cat_100029,
    cat_100030,
    cat_100032,
  ) |> 
  left_join(full_j) |>
  id_samples()

rm(cat_501,
   cat_502,
   cat_503,
   cat_505,
   cat_1000027_reg,
   cat_100029,
   cat_100030,
   cat_100032)

gc()

new_cog_guide <- 
  final_cog |>
  select(category:field) |>
  distinct()


behaviorals_guide_v3 <- 
  new_cog_guide |> 
  arrange(category, field_id) |> 
  mutate(origin = ifelse(field_id %in% extended_cog_guide$field_id, "regular", "new")) 

#write_csv(behaviorals_guide_v3, here("data", "others", "behaviorals_guide.csv"))


# hand modified prettyfied labels
behaviorals_labels <- read_csv(here("data", "others", "behaviorals_guide2_nice.csv"))

cog_tests <- function(dat) {
  tests <- 
    dat |> 
    drop_na(val) |> 
    group_by(category, field_id, field) |> 
    nest() |> 
    mutate(t_test = map(data, ~ t.test(.x$val ~ factor(.x$sex))),
           t = map(t_test, tidy),
           d = map(data, ~cohens_d(val ~ sex, data = .))) |> 
    unnest(t) |>
    unnest(d) |>
    ungroup() |> 
    select(-c(data, t_test)) 
  
  return(tests)
}

final_cog <- 
  left_join(final_cog, behaviorals_labels) |> 
  filter(keep == 1) |> 
  mutate(val2 = val, 
         val = val * reverse)

full_cognitives <- cog_tests(final_cog)
matched_cognitives <- cog_tests(final_cog |> filter(matched == 1))
random_cognitives <- cog_tests(final_cog |> filter(random == 1))

table(full_cognitives$magnitude)
table(matched_cognitives$magnitude)
table(random_cognitives$magnitude)


# Joining info about tests
cog_fin <- 
  random_cognitives |>
  select(category, field_id, p.value, effsize) |> 
  left_join(matched_cognitives |> 
              select(category, field_id, p.value, effsize), 
            by = c("category", "field_id"),
            suffix = c("_random", "_matched")) |>
  left_join(behaviorals_labels)

# write_csv(cog_fin, here("data", "others", "cognitive_results.csv"))

# First viz try
ggplot(cog_fin, aes(x = effsize_random, y = effsize_matched)) +
  geom_point() +
  geom_abline(color = "red") +
  geom_hline(yintercept = 0, linetype = 2, color = "#595959") +
  geom_vline(xintercept = 0, linetype = 2, color = "#595959") +
  ggrepel::geom_text_repel(aes(label = str_wrap(clean_field, 20)),
                           box.padding = .1) +
  plots_theme


# Getting correlations between performance and head size
cog_fin2 <- 
  cog_fin |> 
  group_by(cat_desc, field) |> 
  summarise(
    effsize_matched = mean(effsize_matched),
    effsize_random = mean(effsize_random)
  ) |>
  left_join(behaviorals_labels |> 
              select(-c(field_id, origin)) |> 
              filter(keep == 1) |> 
              distinct())

range(cog_fin2$effsize_matched)
range(cog_fin2$effsize_random)

cors <- 
  final_cog |>
  select(subject_id, val, category:field, tiv) |> 
  group_by(cat_desc, field) |>
  summarise(r = cor.test(val,  tiv, use = "complete.obs")$estimate,
            p = cor.test(val,  tiv, use = "complete.obs")$p.value,
            .groups = "drop") 

# Cognitive domains
doms <- tibble(category = c(503, 505, 501, 502, 100032, 100027, 100030, 100029),
               ev = c("exf", "exf", "nvr", "ps", "ps", "vnr", "vdm", "wm") )

fig4 <-
  ggplot(cog_fin2 |> right_join(cors),
         aes(x = effsize_random, y = effsize_matched)) +
  geom_point(aes(size = r, color = cat_desc)#, shape = ev)#, alpha = .7,
             ) +
  geom_abline(color = "red", alpha = .7) +
  geom_hline(yintercept = 0, linetype = 2, color = "#595959") +
  geom_vline(xintercept = 0, linetype = 2, color = "#595959") +
  labs(
    title = "Effect sizes of sex differences between the <span style='color:#ee9b00'>matched</span><br>
    and the <span style='color:#ae2012'>not matched</span> samples",
    subtitle = "Effect sizes of the difference between the scores of <span style='color:#bb3e03'>females</span> and <span style='color:#26547C'>males</span>",
    x = "Not matched sample",
    y = "Matched sample",
    color = "Test",
    # size = "Correlation with TIV",
  ) +
  # scale_color_viridis_d(guide = "legend", end = .9) +
  scale_color_manual(values=met.brewer("Tiepolo", 9)) +
  # scale_color_manual(values = nc) +
  ggrepel::geom_text_repel(aes(label = str_wrap(clean_field, 20),
                               color = cat_desc),
                           size = 3.7,
                           show.legend = F,
                           box.padding = .4,
                           point.padding = .2	
                            
                           ) +
  plots_theme +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown(),
    legend.title = element_text(),
    aspect.ratio = 1,
    legend.text = element_text(size = 10.5)
  ) +
  guides(size = "none",
         color = guide_legend(override.aes = list(alpha = c(1,1),
                                                  size = c(4, 4))))

fig4

ggsave(here("outputs", "plots", "main_figures", "fig4_cognitive.png"), fig4, 
       width = 9, height = 9, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig4_cognitive.svg"), fig4, 
       width = 9, height = 9, bg = "white", dpi = 600)


# helper to get part of the legend

fig4_h <- 
  ggplot(cog_fin2 |> right_join(cors),
       aes(x = effsize_random, y = effsize_matched)) +
  geom_point(aes(size = r, color = cat_desc)#, shape = ev)#, alpha = .7,
  ) +
  geom_abline(color = "red", alpha = .7) +
  geom_hline(yintercept = 0, linetype = 2, color = "#595959") +
  geom_vline(xintercept = 0, linetype = 2, color = "#595959") +
  labs(
    title = "Effect sizes of sex differences between the <span style='color:#ee9b00'>matched</span><br>
    and the <span style='color:#ae2012'>not matched</span> samples",
    subtitle = "Effect sizes of the difference between the scores of <span style='color:#bb3e03'>females</span> and <span style='color:#26547C'>males</span>",
    x = "Not matched sample",
    y = "Matched sample",
    color = "Test",
    size = "Correlation with TIV"
  ) +
  scale_color_manual(values=met.brewer("Tiepolo", 9)) +
  ggrepel::geom_text_repel(aes(label = str_wrap(clean_field, 20),
                               color = cat_desc),
                           size = 3.7,
                           show.legend = F,
                           box.padding = .4,
                           point.padding = .2	
                           
  ) +
  theme_light() +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown(),
    legend.title = element_text(),
    aspect.ratio = 1,
    legend.text = element_text(size = 10.5)
  ) 

ggsave(here("outputs", "plots", "main_figures", "fig4_helper.png"), fig4h, 
       width = 9, height = 9, bg = "white", dpi = 600)
ggsave(here("outputs", "plots", "main_figures", "fig4_helper.svg"), fig4h, 
       width = 9, height = 9, bg = "white", dpi = 600)

#### Checking participants ----
participants <- 
  final_cog |> 
  select(subject_id, category, cat_desc, sex:fs_matc) |> 
  distinct() |> 
  group_by(category, cat_desc, sex) |> 
  summarise(
    full = n(),
    matched = sum(matched), 
    age_mat = sum(age_mat), 
    random = sum(random), 
    extreme = sum(extreme), 
  ) 

part_by_test <- 
  participants |> 
  pivot_longer(-c(category, cat_desc, sex), names_to = "df") |> 
  left_join(samples_ref) |> 
  mutate(clean_sample = factor(clean_sample, levels = rev(samp_order))) |> 
  ggplot(aes(x = reorder(cat_desc, category), y = value)) +
  facet_wrap(~clean_sample, scales = "free_y") +
  geom_bar(aes(color = sex, fill = sex),
           stat = "identity", position = "dodge") +
  labs(
    title = "Participants with information for the cognitive tests",
    subtitle = "Number of <span style='color:#6A559C'>females</span> and <span style='color:#F4632E'>males</span> for each test",
    x = "",
    y = "Participants"
  ) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = color_sex, 
                     aesthetics = c("fill", "color")) +
  plots_theme +
  theme(
    plot.subtitle = element_markdown(),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)
  ) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,
           color = "#000000",  linewidth = 1.2) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,
           color = "#000000",  linewidth = 1.2)

ggsave(here("outputs", "plots", "supplementary", "s12_part_per_test.png"),
       part_by_test, width = 11, height = 7, dpi = 600, bg = "white")
ggsave(here("outputs", "plots", "supplementary", "s12_part_per_test.svg"),
       part_by_test, width = 11, height = 7, dpi = 600, bg = "white")


participants |> 
  pivot_longer(-c(category, cat_desc, sex), names_to = "df") |> 
  left_join(
    bind_rows(
      count(full, sex) |> mutate(df = "full"),
      count(matched, sex) |> mutate(df = "matched"),
      count(random, sex) |> mutate(df = "random"),
      count(age_mat, sex) |> mutate(df = "age_mat"),
      count(extreme, sex) |> mutate(df = "extreme"),
    )
  ) |> 
  ungroup() |> 
  mutate(prc = value / n) |> 
  View()