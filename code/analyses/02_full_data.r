library(dplyr)
library(readr)
library(here)
library(tidyr)
library(lubridate)
source("../ukbb_data_r/functions.R")

## This code extracts basic demographic information for all the UKBB sujbects
## for which he have fs data. It requieres ukbb extraction codes.

####
dbm_qc <- read_csv(here("data", "others", "info_for_qc.csv"))
brain_size <- read_csv(here("data", "fs_data", "brain_sizes.csv"))
euler_number <- read_csv(here("data", "fs_data", "euler_number.csv"))

## Get ukbb data
ukbb_path <- "/data/zeiyas/in_vivo_Datasets/UKbiobank/brzali"
ukbb <- open_dataset(paste0(ukbb_path, "/UKBB-tabular-processing/current_melt.arrow"), format = "ipc")
codings <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Codings.tsv"))
diction <- read_tsv(paste0(ukbb_path, "/UKBB-tabular-processing/Data_Dictionary_Showcase.tsv"))

# sanity check
# brain_size |> filter(is.na(tiv))

big_sample <- 
    dbm_qc |>
    filter(QC_vol == 1) |>
    select(subject_id) |>
    inner_join(
        brain_size |> select(subject_id)
    ) |>
    pull(subject_id)

ukbb_fields <- c(31,34, 52, 53, 54, 21003)

ukbb_long <- get_ukbb_data(arrow = ukbb, diction = diction, 
                             fields_id = ukbb_fields, 
                             subjs_id = big_sample,
                             instances = 2)

ukbb_data <- order_data(ukbb_long, diction, codings)
glimpse(ukbb_data)
unique(ukbb_data$field)

full_ukbb_data <- 
  ukbb_data |>
    select(subject_id, field, value) |> 
    pivot_wider(names_from = field, values_from = value) |> 
    clean_names() |>
    mutate(
      date_birth = as_date(paste(year_of_birth, month_of_birth, "15", sep = "-")),
      date_mri = as_date(date_of_attending_assessment_centre),
      age_months = interval(date_birth, date_mri) %/% months(1),
      age_when_attended_assessment_centre = as.numeric(age_when_attended_assessment_centre)
    ) |>
    rename(
      "age_years" = age_when_attended_assessment_centre,
      "center" = uk_biobank_assessment_centre
    ) 
  
glimpse(full_ukbb_data)

full_ukbb_data  <- 
  full_ukbb_data |>
  left_join(euler_number |> select(subject_id, euler_num = avg_num))

write_csv(full_ukbb_data, here("data", "samples", "full_ukbb_data.csv"))

