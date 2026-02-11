setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# Load packages

library(tidyverse)
library(data.table)
library(scales)

metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/metadata (1).csv")
metadata_red <- metadata %>%
  dplyr::select(sample_ID, Cohort, ALT, AST, Alkaline_Phosphatase)
medication <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/medication_updated.csv")
medication$sample_ID <- gsub("-", "_", medication$sample_ID)
metadata <- metadata_red %>%
  left_join(medication, by = c("sample_ID" = "sample_ID"))

feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/MZmine_output/fbmn_pediatric_iimn_fbmn_quant.csv")
feature_table <- feature_table %>%
  rename_all(~gsub(" Peak area", "", .))
colnames(feature_table)

annotations <- fread("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/99cebd129c6748e3b67a8834589b9faa-merged_results_with_gnps.tsv")
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

data_sixmix <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "sixmix", filename))
data_sixmix$filename <- gsub(".mzML", "", data_sixmix$filename)

sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix |>  column_to_rownames("filename") |> colMeans(),
                                  SD_sixmix =  data_sixmix |>  column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) |>
  dplyr::filter(Mean_sixmix > 0) |> arrange(desc(Mean_sixmix))
sixmix_feature_info$Feature <- as.character(sixmix_feature_info$Feature)
sixmix_feature_info_name <- left_join(sixmix_feature_info, info_feature_name, by = "Feature" ) |> 
  dplyr::select("Feature", "Mean_sixmix", "SD_sixmix", "CV_sixmix", "mz", "RT", "Compound_Name")

# Filter for blank
data_blank <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "Blank", filename))

blank_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                 Mean_blank = data_blank |> column_to_rownames("filename") |> colMeans(),
                                 SD_blank =  data_blank |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |> arrange(desc(Mean_blank))

blank_feature_info_name <- left_join(blank_feature_info, info_feature_name, by = "Feature" )

# Filter for samples
data_sample <- data_transpose |> 
  dplyr::filter(!str_detect(filename, "Blank|sixmix|pool"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

# Filter for Pool
data_pool <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "pool", filename))
data_pool$filename <-gsub(".mzML", "", data_pool$filename)

pool_feature_info <- data.frame(Feature = colnames(data_pool)[-1],
                                Mean_pool = data_pool |> column_to_rownames("filename") |> colMeans(),
                                SD_pool =  data_pool |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_pool = SD_pool/Mean_pool) |>
  dplyr::filter(Mean_pool > 0) |> arrange(desc(Mean_pool))
pool_feature_info_name <- left_join(pool_feature_info, info_feature_name, by = "Feature" )


# remove features belonging to QCmix (Amitriptyline, Coumarin 314, Sulfamethazine, Sulfadimethoxine, Sulfamethizole, Sulfachloropyridazine)
specified_features <- c("64986", "93883", "34687", "46935", "35250", "38229")

feature_to_remove <- blank_feature_info |>  
  left_join(pool_feature_info) |> 
  dplyr::filter(Mean_blank > 0) |> 
  dplyr::mutate(pool_Blank = Mean_pool/Mean_blank) |> 
  dplyr::filter(pool_Blank < 5 | is.na(pool_Blank)) |> 
  dplyr::bind_rows(blank_feature_info |> 
                     dplyr::filter(Feature %in% specified_features)) |> 
  dplyr::distinct(Feature, .keep_all = TRUE)

# Remove any features listed in 'feature_to_remove'
data_clean <- data_transpose |> 
  dplyr::select(-c(feature_to_remove$Feature))

# Features to be removed samples/QCmix < 5
feature_to_remove_qcmix <- sixmix_feature_info |> 
  left_join(sample_feature_info) |> 
  dplyr::filter(Mean_sixmix > 0) |> 
  dplyr::mutate(sample_Mix = Mean_sample/Mean_sixmix) |> 
  dplyr::filter(sample_Mix < 5 | is.na(sample_Mix)) |> 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with QCmix removal, pool, blanks
data_clean2 <- data_clean |> 
  dplyr::select(-c(feature_to_remove_qcmix$Feature)) |> 
  dplyr::filter(!str_detect(filename, "pool|Blank|sixmix"))

data_5ASA <- data_clean2 |> 
  dplyr::mutate(cholyl_5ASA = `91304` + `67925` + `75945` + `82989`) |> 
  dplyr::mutate(ASA_pyruvate = `6985`) |>
  dplyr::mutate(X5ASA = `5927` + `16165`+ `17230` + `14253`+ `10363` + `9393`+ `11043`+ `5239`) |> 
  dplyr::mutate(acetyl_5ASA = `30801`) |>
  dplyr::mutate(dihydroxylated_5ASA = `108607` + `111062`) |>
  dplyr::mutate(monohydroxylated_5ASA = `127284`) |>
  dplyr::mutate(all_bile_acids = cholyl_5ASA + dihydroxylated_5ASA + monohydroxylated_5ASA) |>
  dplyr::select(filename, X5ASA, acetyl_5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids) |> dplyr::filter(X5ASA > 1.3e5)

data_aggre <- data_clean2 |> 
  dplyr::mutate(cholyl_5ASA = `91304` + `67925` + `75945` + `82989`) |> 
  dplyr::mutate(ASA_pyruvate = `6985`) |>
  dplyr::mutate(X5ASA = `5927` + `16165`+ `17230` + `14253`+ `10363` + `9393`+ `11043`+ `5239`) |> 
  dplyr::mutate(acetyl_5ASA = `30801`) |>
  dplyr::mutate(dihydroxylated_5ASA = `108607` + `111062`) |>
  dplyr::mutate(monohydroxylated_5ASA = `127284`) |>
  dplyr::mutate(all_bile_acids = cholyl_5ASA + dihydroxylated_5ASA + monohydroxylated_5ASA) |>
  dplyr::select(filename, X5ASA, ASA_pyruvate, acetyl_5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids)

plot_discordance <- function(df, x_var, x_label) {
  cor_res <- cor.test(
    df[[x_var]],
    df$X5ASA,
    method = "spearman",
    exact  = FALSE)
  rho  <- unname(cor_res$estimate)
  pval <- cor_res$p.value
  df_discord <- df %>%
    dplyr::mutate(
      rank_x    = rank(.data[[x_var]], ties.method = "average"),
      rank_ASA  = rank(X5ASA,          ties.method = "average"),
      rank_diff = rank_ASA - rank_x,
      discordant = abs(rank_diff) >
        quantile(abs(rank_diff), 0.8, na.rm = TRUE),
      x_plot = .data[[x_var]] + 1,
      y_plot = X5ASA + 1)
  ymax <- max(df_discord$y_plot, na.rm = TRUE)
  y_upper <- if (ymax >= 5e7) 1e8 else NA_real_
  ggplot(df_discord, aes(x = x_plot, y = y_plot)) +
    geom_point(aes(color = discordant), size = 3, alpha = 1) +
    scale_color_manual(values = c("FALSE" = "grey30", "TRUE" = "firebrick")) +
    scale_x_log10(
      breaks  = scales::breaks_log(n = 6),
      expand  = expansion(mult = c(0.15, 0.05))) +
    scale_y_log10(
      limits  = c(NA, y_upper),
      breaks  = scales::breaks_log(n = 6),
      expand  = expansion(mult = c(0.05, 0.05))) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      x = paste0(x_label, " (peak area + 1; log\u2081\u2080)"),
      y = "5-ASA (peak area + 1; log\u2081\u2080)",
      subtitle = sprintf(
        "Spearman \u03c1 = %.2f, p = %.2g\nRed = top 20%% discordant by rank difference",
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

#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/plot_cholyl.pdf", plot = plot_cholyl, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/plot_dihydroxy.pdf", plot = plot_dihydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/plot_monohydroxy.pdf", plot = plot_monohydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/plot_all_bile.pdf", plot = plot_all_bile, width = 3, height = 5, dpi = 900)

# Merge with metadata using the matching filename columns

merged_df_UC <- data_aggre |> 
  dplyr::mutate(Status = case_when(str_detect(filename, "HC") ~ "healthy", str_detect(filename, "PT") ~ "IBD", TRUE ~ NA_character_)) |> 
  dplyr::select(Status, filename, everything())
merged_df_UC$filename <- gsub(".mzML", "", merged_df_UC$filename)

merged_df_UC_metadata <- merged_df_UC |> 
  left_join(metadata, by = c("filename" = "sample_ID")) |> 
  dplyr::filter(!is.na(Status))

# apply a threshold as this data was contaminated with PEG causing carry-over between samples - threshold based on metadata as well
merged_df_UC_metadata_91304 <- merged_df_UC_metadata |> 
  dplyr::mutate(dplyr::across(c(X5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids), ~ if_else(.x < 10000, 0, .x))) |>
  dplyr::mutate(dplyr::across(c(X5ASA), ~ if_else(.x < 1.3e5, 0, .x))) |>
  dplyr::mutate(log_cholyl_5ASA = log(cholyl_5ASA + 1)) |>
  dplyr::mutate(log_5ASA = log(X5ASA + 1)) |>
  dplyr::select(log_cholyl_5ASA, filename, Status, Medications, ALT, AST, Alkaline_Phosphatase, everything())

#write_csv(merged_df_UC_metadata_91304, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/Pediatric_raw_data.csv")

# Plot with jitter points colored by the custom point_color column
plot_UC <- ggplot(merged_df_UC_metadata_91304, aes(x = Status, y = log_cholyl_5ASA)) +
  geom_boxplot(aes(fill = Status), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Status), width = 0.2, size = 6, alpha = 1, shape = 16) +
  scale_fill_discrete() +
  scale_color_discrete() +
  scale_y_continuous(limits = c(0, 20)) +
  labs(y = "Log (Peak area)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title  = element_text(size = 18, face = "bold"), plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.position = "none")

plot_UC
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/Pediatric_log_5-ASA.pdf", plot_UC, width = 2, height = 7, dpi = 900)

df2 <- merged_df_UC_metadata_91304 %>%
  dplyr::mutate(presence = if_else(log_cholyl_5ASA > 0, "yes", "no"))
diag_lvls <- levels(factor(df2$Status))
pairs     <- combn(diag_lvls, 2, simplify = FALSE)
fisher_res <- map_dfr(pairs, function(lvls) {
  sub <- df2 %>% filter(Status %in% lvls)
  tbl <- table(
    factor(sub$Status, levels = lvls),
    factor(sub$presence, levels = c("yes","no"))
  )
  ft  <- fisher.test(tbl)
  tibble(
    group1 = lvls[1],
    group2 = lvls[2],
    p      = ft$p.value
  )
})
fisher_res <- fisher_res %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "none"))
fisher_res

