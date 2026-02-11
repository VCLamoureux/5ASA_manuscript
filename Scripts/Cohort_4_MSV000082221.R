setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# Load packages
library(tidyverse)
library(data.table)
library(pheatmap)
library(ggbeeswarm)
library(broom)

metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/Metadata_complete_MSV000082221.csv")
metadata_1 <- metadata %>%
  dplyr::mutate(plcl_id = anonymized_name)%>%
  dplyr::mutate(patient = str_extract(sample_name, "(?<=\\.)[^.]+(?=\\.[^.]+$)")) |> 
  dplyr::select(plcl_id, sample_name, patient, use_5asa, everything())

n_patients <- metadata_1 %>%
  dplyr::distinct(patient) %>%
  nrow()
n_patients

patients_per_disease <- metadata_1 %>%
  dplyr::distinct(patient, disease) %>%  
  dplyr::count(disease, name = "n_patients")
patients_per_disease

feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/Version2_MSV000082221_iimn_fbmn_quant.csv")
feature_table <- feature_table %>%
  dplyr::rename_all(~gsub(" Peak area", "", .))
colnames(feature_table)

annotations <- fread("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/4f5257efa6b1452cb9e72b4b0d333fbc-merged_results_with_gnps.tsv")
annotations$`#Scan#` <- as.character(annotations$`#Scan#`)
annotations <- annotations |> 
  dplyr::rename(Feature = `#Scan#`)

info_feature <- feature_table[,1:3]
colnames(info_feature) <- c("Feature", "mz", "RT")
info_feature$Feature <- as.character(info_feature$Feature)
info_feature_name <- left_join(info_feature, annotations, by = "Feature")

# Data
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzML")) |> 
  t() |>  
  as.data.frame() |>  
  rownames_to_column("filename")

data_transpose_1 <- data_transpose %>%
  dplyr::mutate(plcl_base = stringr::str_extract(filename, "(?<=PLCL)\\d+(?:\\.\\d+)?"),
    plcl_id = dplyr::case_when(
      plcl_base == "17712" & stringr::str_detect(filename, "^PLCL17712_B") ~ "17712.B",
      plcl_base == "17712" & stringr::str_detect(filename, "^PLCL17712_R") ~ "17712.A", TRUE ~ plcl_base)) %>%
  dplyr::mutate(X5ASA = rowSums(dplyr::across(c(`67`, `78`, `137`)), na.rm = TRUE),
    cholyl_5ASA = rowSums(dplyr::across(c(`655`, `720`)), na.rm = TRUE),
    dihydroxylated_5ASA = rowSums(dplyr::across(c(`786`, `788`)), na.rm = TRUE),
    monohydroxylated_5ASA = rowSums(dplyr::across(c(`813`)), na.rm = TRUE),
    all_bile_acids = rowSums(
      dplyr::across(c(cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA)), na.rm = TRUE)) %>%
  dplyr::select(plcl_id, X5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids)

merged_df <- data_transpose_1 %>%
  left_join(metadata_1, by = "plcl_id") |> 
  dplyr::select(use_5asa, medication_complianc, treatment_response, disease, X5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids, has_surgery, everything()) |> 
  dplyr::filter(!is.na(sample_name))

#write_csv(merged_df, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/MSV000082221_raw_data.csv")


# do some stats
n_use5asa_pos <- merged_df %>%
  dplyr::filter(use_5asa == "y", is.finite(X5ASA), X5ASA > 0) %>%
  dplyr::distinct(patient, .keep_all = TRUE) |> 
  dplyr::filter(disease == "cd")
n_use5asa_pos

n_use5asa <- merged_df %>%
  dplyr::filter(use_5asa == "y") %>%
  dplyr::distinct(patient)
n_use5asa

# recode
merged_df_just_all_bas <- merged_df |>
  dplyr::filter(X5ASA > 0) |> 
  dplyr::mutate(has_5asa_conjugate = if_else(
      (cholyl_5ASA > 0| dihydroxylated_5ASA > 0 | monohydroxylated_5ASA > 0), "present", "absent"), has_5asa_conjugate = factor(has_5asa_conjugate, levels = c("absent", "present"))) |> 
  dplyr::select(has_5asa_conjugate, everything())

merged_df_just_all_bas_1 <- merged_df_just_all_bas |> 
  dplyr::mutate(
    treatment_response = case_when(
      treatment_response %in% c("response", "responder") ~ "response",
      treatment_response %in% c("loss_response", "non_response", "non-responder") ~ "loss_response", TRUE ~ NA_character_),
    treatment_response = factor(
      treatment_response,
      levels = c("loss_response", "response")))

merged_df_just_all_bas_2 <- merged_df_just_all_bas_1 |> 
  dplyr::mutate(has_surgery = case_when(
      has_surgery %in% c("y", "yes", "Y", "Yes") ~ "yes",
      has_surgery %in% c("n", "no", "N", "No") ~ "no",
      TRUE ~ NA_character_), has_surgery = factor(has_surgery, levels = c("no", "yes")))

plot_discordance <- function(df, x_var, x_label) {
  cor_res <- cor.test(
    df[[x_var]],
    df$X5ASA,
    method = "spearman",
    exact = FALSE)
  rho  <- unname(cor_res$estimate)
  pval <- cor_res$p.value
  df_discord <- df %>%
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
  merged_df_just_all_bas_2,
  x_var   = "cholyl_5ASA",
  x_label = "cholyl-5-ASA")
plot_cholyl

plot_dihydroxy <- plot_discordance(
  merged_df_just_all_bas_2,
  x_var   = "dihydroxylated_5ASA",
  x_label = "dihydroxylated 5-ASA")
plot_dihydroxy

plot_monohydroxy <- plot_discordance(
  merged_df_just_all_bas_2,
  x_var   = "monohydroxylated_5ASA",
  x_label = "monohydroxylated 5-ASA")
plot_monohydroxy

plot_all_bile <- plot_discordance(
  merged_df_just_all_bas_2,
  x_var   = "all_bile_acids",
  x_label = "All bile acid–5-ASA conjugates")
plot_all_bile

#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/plot_cholyl.pdf", plot = plot_cholyl, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/plot_dihydroxy.pdf", plot = plot_dihydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/plot_monohydroxy.pdf", plot = plot_monohydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/plot_all_bile.pdf", plot = plot_all_bile, width = 3, height = 5, dpi = 900)

# After you build merged_df_UC with the join:
# Then continue with your existing cleaning / colors
merged_df_color <- merged_df |> 
  dplyr::mutate(
    cholyl_5ASA = if_else(cholyl_5ASA == 0, 1, cholyl_5ASA),
    X5ASA_1 = if_else(X5ASA == 0, 1, X5ASA),
    use_5asa = factor(use_5asa),
    point_color = dplyr::case_when(
      use_5asa == "y" ~ "#27AAE1",
      use_5asa == "n" ~ "#EC008C", TRUE ~ "grey")) |> 
  dplyr::select(point_color, everything()) |>
  dplyr::mutate(disease = factor(disease, levels = c("uc", "cd")))

n_patients_uc <- merged_df_color %>%
  dplyr::filter(disease == "uc") |> 
  dplyr::distinct(patient) %>%
  nrow()
n_patients_uc

n_sample_uc <- merged_df_color %>%
  dplyr::filter(disease == "uc") |> 
  nrow()
n_sample_uc

plot_UC <- ggplot(merged_df_color, aes(x = disease, y = cholyl_5ASA)) +
  geom_boxplot(fill = "grey85", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = point_color), width = 0.2, size  = 6, alpha = 1, shape = 16) +
  scale_color_identity() +
  scale_y_log10(limits = c(1, 1e5), breaks = 10^(0:5), labels = parse(text = paste0("10^", 0:5))) +
  labs(x = "disease", y = "Peak area (log scale)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 18, face = "bold"), legend.position = "none")

plot_UC
#ggsave("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/MSV000082221_5ASA_CA.pdf", plot_UC, width = 3, height = 7, dpi = 900)

# for 5-ASA 
plot_UC_5ASA <- ggplot(merged_df_color, aes(x = disease, y = X5ASA_1)) +
  geom_boxplot(fill = "grey85", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = point_color), width = 0.2, size  = 6, alpha = 1, shape = 16) +
  scale_color_identity() +
  scale_y_log10(limits = c(1, 1e5), breaks = 10^(0:5), labels = parse(text = paste0("10^", 0:5))) +
  labs(x = "disease", y = "Peak area (log scale)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14),
        axis.title  = element_text(size = 18, face = "bold"), legend.position = "none")

plot_UC_5ASA
#ggsave("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/MSV000082221_5ASA.pdf", plot_UC_5ASA, width = 3, height = 7, dpi = 900)


#write_csv(merged_df_color, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082221/MSV000082221_5ASA_bile_acid_data.csv")





