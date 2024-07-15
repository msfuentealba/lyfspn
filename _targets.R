library(targets)

tar_option_set(
  packages = c("tidyverse", "circlize", "ComplexHeatmap","RColorBrewer","ggh4x","ggpmisc","ggpubr") 
)

tar_source()

list(
  tar_target(name = file_epigenomics, command = file.path(here::here("data/epigenomics_data.csv")), format = "file"),
  tar_target(name = file_clinical, command = file.path(here::here("data/clinical_data.xlsx")), format = "file"),
  tar_target(name = file_proteomics, command = file.path(here::here("data/proteomics_data.xlsx")), format = "file"),
  tar_target(name = file_metabolomics, command = file.path(here::here("data/metabolomics_data.xlsx")), format = "file"),
  tar_target(name = file_lipidomics, command = file.path(here::here("data/lipidomics_data.xlsx")), format = "file"),
  tar_target(name = file_glycomics, command = file.path(here::here("data/glycomics_data.xlsx")), format = "file"),
  tar_target(name = file_cytomics, command = file.path(here::here("data/cytomics_data.xlsx")), format = "file"),
  
  tar_target(name = age_acceleration_dif, command = extract_age_acceleration_dif(file_epigenomics)),
  tar_target(name = baseline_clinical, command = extract_baseline_clinical(file_clinical)),
  tar_target(name = baseline_proteomics, command = extract_baseline_proteomics(file_proteomics)),
  tar_target(name = baseline_metabolomics, command = extract_baseline_metabolomics(file_metabolomics)),
  tar_target(name = baseline_lipidomics, command = extract_baseline_lipidomics(file_lipidomics)),
  tar_target(name = baseline_glycomics, command = extract_baseline_glycomics(file_glycomics)),
  tar_target(name = baseline_cytomics, command = extract_baseline_cytomics(file_cytomics)),
  
  tar_target(name = baseline_correlations, command = calculate_correlation_omics(age_acceleration_dif, baseline_clinical, baseline_proteomics, baseline_metabolomics, baseline_lipidomics, baseline_glycomics, baseline_cytomics)),
  tar_target(name = baseline_stats, command = stats_correlations_baseline(baseline_correlations)),
  tar_target(name = baseline_significant, command = extract_significant_baseline(baseline_stats)),
  tar_target(name = baseline_links, command = find_baseline_links(age_acceleration_dif, baseline_clinical, baseline_proteomics, baseline_metabolomics, baseline_lipidomics, baseline_glycomics, baseline_cytomics, baseline_stats, baseline_significant)),
  
  tar_target(name = circos_baseline, command = plot_circos_baseline(baseline_stats,baseline_significant,baseline_links)),
  tar_target(name = top_baseline, command = plot_top_baseline(age_acceleration_dif, baseline_clinical, baseline_proteomics, baseline_metabolomics, baseline_lipidomics, baseline_glycomics, baseline_cytomics, baseline_stats)),
  tar_target(name = print_plot, command = print_plot(top_baseline, "top_baseline_clinical.pdf", width=9.02, height=8.01))
)



