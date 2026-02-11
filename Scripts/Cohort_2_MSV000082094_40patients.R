setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# Load packages
library(tidyverse)
library(scales)
library(duckplyr)

metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/UC_MP_Emperor_Map.csv")
metadata_summary_2 <- metadata |> 
  ungroup() |> 
  group_by(Disease, current_5ASA, ASA_exposure) |> 
  summarize(count = n(), .groups = "drop")

feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/mzmine_version3_015_min_RT/version3_MSV000082094_iimn_fbmn_quant.csv")
feature_table <- feature_table |>
  rename_all(~gsub(" Peak area", "", .))
colnames(feature_table)

annotations <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/Version_3_10033b6efd11402fa83e4fedee5bab7e-merged_results_with_gnps.tsv")

annotations$`#Scan#` <- as.character(annotations$`#Scan#`)
annotations <- annotations |> 
  dplyr::rename(Feature = `#Scan#`)

# Data
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzXML")) |> 
  t() |>  
  as.data.frame() |>  
  rownames_to_column("filename") |> 
  dplyr::mutate(id = sub("_.*", "", filename)) |> 
  dplyr::select(filename, id, everything()) |> 
  dplyr::filter(!str_detect(filename, "Extract"))

# Filter for blank
data_blank <- data_transpose |>  
  dplyr::filter(str_detect(filename, "Blank")) |> 
  dplyr::select(-id)

blank_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                 Mean_blank = data_blank |> column_to_rownames("filename") |> colMeans(),
                                 SD_blank =  data_blank |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |> arrange(desc(Mean_blank))

# Filter for samples
data_sample <- data_transpose |> 
  dplyr::filter(!str_detect(filename, "Blank")) |> 
  dplyr::select(-id)

sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))

# remove features belonging to QCmix (Amitriptyline, Coumarin 314, Sulfamethazine, Sulfadimethoxine, Sulfamethizole, Sulfachloropyridazine)

feature_to_remove <- blank_feature_info |>
  left_join(sample_feature_info, by = "Feature") |>
  dplyr::filter(Mean_blank > 0) |>
  dplyr::mutate(sample_Blank = Mean_sample / Mean_blank) |>
  dplyr::filter(sample_Blank < 5 | is.na(sample_Blank)) |>
  distinct(Feature, .keep_all = TRUE)

# Remove any features listed in 'feature_to_remove'
data_clean <- data_transpose |>
  dplyr::select(-any_of(feature_to_remove$Feature))

data_5ASA <- data_clean |>
  dplyr::mutate(cholyl_5ASA = `2874` + `3207` + `3001`) |> 
  dplyr::mutate(X5ASA = `129`) |> 
  dplyr::mutate(dihydroxylated_5ASA = `3541` + `3414`) |>
  dplyr::mutate(monohydroxylated_5ASA = `3866` + `3793`) |>
  dplyr::mutate(all_bile_acids = cholyl_5ASA + dihydroxylated_5ASA + monohydroxylated_5ASA) |>
  dplyr::select(filename, id, X5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids)

plot_discordance <- function(df, x_var, x_label) {
  cor_res <- cor.test(
    df[[x_var]],
    df$X5ASA,
    method = "spearman",
    exact = FALSE)
  rho  <- unname(cor_res$estimate)
  pval <- cor_res$p.value
  df_discord <- df |>
    dplyr::mutate(
      rank_x   = rank(.data[[x_var]], ties.method = "average"),
      rank_ASA = rank(X5ASA, ties.method = "average"),
      rank_diff = rank_ASA - rank_x,
      discordant = abs(rank_diff) >
        quantile(abs(rank_diff), 0.8, na.rm = TRUE),
      x_plot = .data[[x_var]] + 1,
      y_plot = X5ASA + 1)
  ggplot(df_discord, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = discordant), size = 3, alpha = 1) +
    scale_color_manual(values = c("FALSE" = "grey30", "TRUE" = "firebrick")) +
    scale_x_log10(expand = expansion(mult = c(0.15, 0.05))) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.05))) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      x = paste0(x_label, " (raw + 1; log10)"),
      y = "5-ASA (raw + 1; log10)",
      subtitle = sprintf(
        "Spearman ρ = %.2f, p = %.2g\nRed = top 20%% discordant by rank difference",
        rho, pval)) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none")
}

plot_cholyl <- plot_discordance(
  data_5ASA,
  x_var   = "cholyl_5ASA",
  x_label = "cholyl-5-ASA")
plot_cholyl

plot_dihydroxy <- plot_discordance(
  data_5ASA,
  x_var   = "dihydroxylated_5ASA",
  x_label = "dihydroxylated 5-ASA")
plot_dihydroxy

plot_monohydroxy <- plot_discordance(
  data_5ASA,
  x_var   = "monohydroxylated_5ASA",
  x_label = "monohydroxylated 5-ASA")
plot_monohydroxy

plot_all_bile <- plot_discordance(
  data_5ASA,
  x_var   = "all_bile_acids",
  x_label = "All bile acid–5-ASA conjugates")
plot_all_bile

ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/plot_cholyl.pdf", plot = plot_cholyl, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/plot_dihydroxy.pdf", plot = plot_dihydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/plot_monohydroxy.pdf", plot = plot_monohydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/plot_all_bile.pdf", plot = plot_all_bile, width = 3, height = 5, dpi = 900)

# After you build merged_df_UC with the join:
merged_df_UC_clean <- data_5ASA |>
  dplyr::left_join(metadata, by = c("id" = "vial_name")) |>
  dplyr::select(id, X5ASA, cholyl_5ASA, Disease, current_5ASA, ASA_exposure, everything()) |>
  ungroup() |> dplyr::filter(!id == "Blank") |> 
  dplyr::filter(!id == "FIT147E3CAL")

#write_csv(merged_df_UC_clean, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/MSV000082094_raw_data.csv")

count_X5ASA_positive <- merged_df_UC_clean |>
  dplyr::filter(current_5ASA == "1") |>
  nrow()
count_X5ASA_positive

count_X5ASA_exposure <- merged_df_UC_clean |>
  dplyr::filter(current_5ASA == "0", ASA_exposure == "1") |>
  nrow()
count_X5ASA_exposure

count_cholylX5ASA_double_negative <- merged_df_UC_clean |>
  dplyr::filter(current_5ASA == "0", ASA_exposure == "0") |> 
  nrow()
count_cholylX5ASA_double_negative

df_clean_subset <- merged_df_UC_clean |>
  dplyr::mutate(cholyl_5ASA = as.numeric(cholyl_5ASA), X5ASA = as.numeric(X5ASA)) |>
  dplyr::filter(!if_any(where(is.character), ~ str_detect(.x, regex("Blank", ignore_case = TRUE))), current_5ASA == 0, ASA_exposure == 0)

# Then continue with your existing cleaning / colors
merged_df_UC <- merged_df_UC_clean |> 
  dplyr::mutate(
    cholyl_5ASA      = if_else(cholyl_5ASA == 0, 1, cholyl_5ASA),
    current_5ASA    = as.factor(current_5ASA),
    ASA_exposure    = as.factor(ASA_exposure),
    point_color = case_when(
      current_5ASA == "1" & ASA_exposure == "1" ~ "#27AAE1",
      ASA_exposure == "1"                      ~ "#EC008C",
      TRUE                                     ~ "grey")) |> 
  dplyr::filter(!is.na(current_5ASA)) |> 
  dplyr::select(id, cholyl_5ASA, current_5ASA, ASA_exposure, everything())


plot_UC <- ggplot(merged_df_UC, aes(x = current_5ASA, y = cholyl_5ASA)) +
  geom_boxplot(aes(fill = current_5ASA), alpha = 1, outlier.shape = NA) +
  geom_jitter(aes(color = point_color), width = 0.2, size = 6, alpha = 0.7, shape = 16) +
  scale_color_identity() +
  scale_y_log10(breaks = 10^(0:5), labels = parse(text = paste0("10^", 0:5))) +
  labs(y = "Peak area (log scale)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 18, face = "bold"),
    legend.position = "none")

plot_UC

#ggsave("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/MSV00082094_CA_5ASA_version3.pdf", plot_UC, width = 3, height = 7, dpi = 900)

df2 <- merged_df_UC |>
  dplyr::mutate(cholyl_5ASA = if_else(cholyl_5ASA == 1, 0, cholyl_5ASA)) |> 
  dplyr::mutate(presence = if_else(cholyl_5ASA > 0, "yes", "no"))

diag_lvls <- levels(factor(df2$current_5ASA))
pairs     <- combn(diag_lvls, 2, simplify = FALSE)

fisher_res_1 <- map_dfr(pairs, function(lvls) {
  sub <- df2 |> filter(current_5ASA %in% lvls)
  
  # IMPORTANT: fix the table orientation ONCE
  tbl <- table(
    current_5ASA = factor(sub$current_5ASA, levels = lvls),
    presence     = factor(sub$presence,    levels = c("yes", "no"))
  )
  
  ft <- fisher.test(tbl)
  
  OR      <- unname(ft$estimate)      # odds of presence in lvls[1] vs lvls[2]
  OR_inv  <- 1 / OR                   # odds of presence in lvls[2] vs lvls[1]
  CI_low  <- ft$conf.int[1]
  CI_high <- ft$conf.int[2]
  CI_low_inv  <- 1 / CI_high          # flipped CI
  CI_high_inv <- 1 / CI_low
  
  tibble(
    group1        = lvls[1],
    group2        = lvls[2],
    OR            = OR,
    OR_inv        = OR_inv,
    CI_low        = CI_low,
    CI_high       = CI_high,
    CI_low_inv    = CI_low_inv,
    CI_high_inv   = CI_high_inv,
    p             = ft$p.value
  )
})
fisher_res_1

# plotting 5-ASA
merged_df_UC <- merged_df_UC_clean |> 
  dplyr::mutate(
    X5ASA  = if_else(X5ASA == 0, 1, X5ASA),
    current_5ASA  = as.factor(current_5ASA),
    ASA_exposure  = as.factor(ASA_exposure),
    point_color = case_when(
      current_5ASA == "1" & ASA_exposure == "1" ~ "#27AAE1",
      ASA_exposure == "1"                      ~ "#EC008C",
      TRUE                                     ~ "grey")) |> 
  dplyr::filter(!is.na(current_5ASA)) |> 
  dplyr::select(id, X5ASA, current_5ASA, ASA_exposure, everything())

plot_UC <- ggplot(merged_df_UC, aes(x = current_5ASA, y = X5ASA)) +
  geom_boxplot(aes(fill = current_5ASA), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = point_color), width = 0.2, size = 6, alpha = 1, shape = 16) +
  scale_color_identity() +
  scale_y_log10(breaks = 10^(0:5), labels = parse(text = paste0("10^", 0:5))) +
  labs(y = "Peak area (log scale)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 18, face = "bold"),
    legend.position = "none")

plot_UC

#ggsave("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/MSV00082094_5ASA_3.pdf", plot_UC_5ASA, width = 3, height = 7, dpi = 900)

#write_csv(merged_df_UC, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/IBD_cohort_2.csv")








