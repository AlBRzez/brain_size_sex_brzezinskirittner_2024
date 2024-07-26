library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(here)

## This code calculates different totals by tissue for each subject

# Get data ---------------------------------------------------------------------
fs_data <- read_csv(here("data", "fs_data", "clean_fs_data.csv"))

# Get totals -------------------------------------------------------------------
tissue_sums <- 
  fs_data |> 
  group_by(subject_id) |> 
  summarise(
    gm_cortex = cortex_left + cortex_right,
    gm_subcortical = sub_cort_gray,
    gm_cerebrum = gm_subcortical + gm_cortex,
    gm_cerebellum = cerebellum_cortex_left + 
      cerebellum_cortex_right,
    gm_total = gm_cerebrum + gm_cerebellum,
    gm_wb_ukbb = total_gray,
    wm_cerebrum = cerebral_white_matter_left +
      cerebral_white_matter_right,
    wm_cerebellum = cerebellum_white_matter_left +
      cerebellum_white_matter_right,
    wm_total = wm_cerebrum + wm_cerebellum,
    csf_wb = mask-brain_seg,
    ventricles = `3rd_ventricle` + `4th_ventricle` + `5th_ventricle` + ventricle_choroid + lateral_ventricle_left + lateral_ventricle_right,
    mask = mask,
    tbv = gm_total + wm_total,
    brain_stem = brain_stem,
    fs_tiv = tiv,
    brain_seg = brain_seg,
    brain_seg_not_vent, brain_seg_not_vent,
    .groups = "drop"
  ) 

# Save files -------------------------------------------------------------------
write_csv(tissue_sums, here("data", "fs_data", "brain_sizes.csv"))


