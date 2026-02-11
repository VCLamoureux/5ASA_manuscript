setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# load packages
library(tidyverse)
library(readxl)
library(data.table)
library(vegan)
library(caret)
library(mixOmics)
library(rstatix)
library(broom)
library(ggpubr)
library(patchwork)
library(viridis)
library(pheatmap)
library(duckplyr)
library(scales)

# load metadata
metadata <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/Metadata_MSV000084775.csv") |>
  dplyr::mutate(Metabolomics_FileName_Run1 = Metabolomics_FileName_Run1)
metadata$Metabolomics_FileName_Run1 <- gsub(".mzXML", "", metadata$Metabolomics_FileName_Run1)
metadata_summary   <- metadata |> dplyr::count(Diagnosis, name = "count")
metadata_summary_2 <- metadata |> dplyr::count(Diagnosis, `Current_5ASA`, ASA_Exposure, name = "count")

# load tables
feature_table <- readr::read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/MZmine_output/fbmn_iimn_fbmn_quant.csv") |>
  dplyr::rename_with(~ gsub(" Peak area", "", .x, fixed = TRUE))

drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/MZmine_output/fbmn_quant_mzmine")
colnames(drift)[1] <- "Feature"

annotations <- data.table::fread("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/Annotation_MSV000084775.tsv") |>
  as.data.frame() |>
  dplyr::rename(Feature = `#Scan#`) |>
  dplyr::mutate(Feature = as.character(Feature))

data_transpose <- feature_table |>
  tibble::column_to_rownames("row ID") |>
  dplyr::select(tidyselect::contains(".mzML")) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("filename")
data_transpose$filename <- gsub(".mzML", "", data_transpose$filename)

info_feature <- feature_table[, 1:3]
info_feature <- info_feature |> 
  dplyr::rename(Feature = `row ID`)
info_feature$Feature <- as.character(info_feature$Feature)
colnames(info_feature) <- c("Feature", "mz", "RT")

info_feature_name <- info_feature |> 
  left_join(annotations, by = "Feature")

# quality control RT drift
feature_row <- drift |>  
  dplyr::filter(Feature == 817)
rt_columns <- feature_row |>  
  dplyr::select(contains("RT"))
rt_data <- tidyr::gather(rt_columns, key = "filename", value = "RT")
rt_data$filename <- gsub(".mzML Feature RT", "", rt_data$filename)
rt_data_ordered <- rt_data

# Arrange the data by 'order'
rt_data_ordered <- rt_data_ordered |> 
  dplyr::filter(!str_detect(filename, "Blk"))
rt_stats <- rt_data |> 
  summarise(Mean_RT = mean(RT, na.rm = TRUE), 
            SD_RT = sd(RT, na.rm = TRUE), 
            Variance_RT = var(RT, na.rm = TRUE))

ggplot(rt_data_ordered, aes(x = filename, y = RT)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "Retention Time (RT)", title = "Retention Time Drift for Feature 9183")


# quality control peak areas drift
feature_row_peakarea <- drift |>  
  dplyr::filter(Feature == 760)

peakarea_columns <- feature_row_peakarea |>  
  dplyr::select(contains("Peak area"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data$filename <- gsub(".mzML Peak area", "", peakarea_data$filename)

# Merge the two dfs by 'filename'
peakarea_data_ordered <- peakarea_data |> 
  dplyr::filter(!str_detect(filename, "Blk"))

peakarea_stats <- peakarea_data |> 
  summarise(Mean_peakarea = mean(Peak, na.rm = TRUE), 
            SD_peakarea = sd(Peak, na.rm = TRUE), 
            Variance_peakarea = var(Peak, na.rm = TRUE))

ggplot(peakarea_data_ordered, aes(x = filename, y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "peakarea", title = "peakarea Drift for Feature 9183")

tic_df <- data_transpose %>%
  dplyr::mutate(TIC = rowSums(across(-filename), na.rm = TRUE)) |> 
  dplyr::select(filename, TIC, everything()) |> 
  dplyr::filter(!str_detect(filename, "Wash|Blk|Extract"))

p_tic <- ggplot(tic_df, aes(x = filename, y = TIC)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line   = element_line(color = "black")) +
  labs(x = "Sample filename", y = "Total Ion Current (TIC)", title = "TIC Across Samples")
p_tic

data_sixmix <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "Std_Mix", filename)) |> 
  dplyr::filter(!filename == "Std_Mix_03_RA2_01_50152")
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
  dplyr::filter(str_detect(pattern = "Blk", filename))

blank_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                 Mean_blank = data_blank |> column_to_rownames("filename") |> colMeans(),
                                 SD_blank =  data_blank |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |> arrange(desc(Mean_blank))
blank_feature_info_name <- left_join(blank_feature_info, info_feature_name, by = "Feature" )

# Filter for samples
data_sample <- data_transpose |> 
  dplyr::filter(!str_detect(filename, "Blk|Wash|Extrac|Std"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

# remove features belonging to QCmix (Amitriptyline, Coumarin 314, Sulfamethazine, Sulfadimethoxine, Sulfamethizole, Sulfachloropyridazine)
specified_features <- c("1131", "760", "1374", "972", "773", "817")
feature_to_remove <- blank_feature_info |>  
  left_join(sample_feature_info) |> 
  dplyr::filter(Mean_blank > 0) |> 
  dplyr::mutate(sample_Blank = Mean_sample/Mean_blank) |> 
  dplyr::filter(sample_Blank < 5 | is.na(sample_Blank)) |> 
  dplyr::bind_rows(blank_feature_info |> dplyr::filter(Feature %in% specified_features)) |> 
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
  dplyr::filter(!str_detect(filename, "Wash|Std|Extraction|Blk"))

data_5ASA <- data_clean2 %>%
  dplyr::mutate(cholyl_5ASA = `1323` + `1423`) |> 
  dplyr::mutate(X5ASA = `52` + `121`) |> 
  dplyr::mutate(dihydroxylated_5ASA = `1527` + `1618` +`1547`) |>
  dplyr::mutate(monohydroxylated_5ASA = `1767` + `1730`) |>
  dplyr::mutate(all_bile_acids = cholyl_5ASA + dihydroxylated_5ASA + monohydroxylated_5ASA) |>
  dplyr::select(filename, X5ASA, cholyl_5ASA, dihydroxylated_5ASA, monohydroxylated_5ASA, all_bile_acids)

plot_discordance <- function(df, x_var, x_label) {
  
  # Spearman correlation (RAW values)
  cor_res <- cor.test(
    df[[x_var]],
    df$X5ASA,
    method = "spearman",
    exact = FALSE
  )
  
  rho  <- unname(cor_res$estimate)
  pval <- cor_res$p.value
  
  # Rank-based discordance (RAW values)
  df_discord <- df %>%
    mutate(
      rank_x   = rank(.data[[x_var]], ties.method = "average"),
      rank_ASA = rank(X5ASA, ties.method = "average"),
      rank_diff = rank_ASA - rank_x,
      discordant = abs(rank_diff) >
        quantile(abs(rank_diff), 0.8, na.rm = TRUE),
      
      # 🔹 PSEUDOCOUNT FOR PLOTTING ONLY
      x_plot = .data[[x_var]] + 1,
      y_plot = X5ASA + 1
    )
  
  # Plot with pseudocount
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
        rho, pval
      )
    ) +
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

#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/plot_cholyl.pdf", plot = plot_cholyl, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/plot_dihydroxy.pdf", plot = plot_dihydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/plot_monohydroxy.pdf", plot = plot_monohydroxy, width = 3, height = 5, dpi = 900)
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/plot_all_bile.pdf", plot = plot_all_bile, width = 3, height = 5, dpi = 900)

merged_df_all_clean <- data_5ASA %>%
  left_join(metadata, by = c("filename" = "Metabolomics_FileName_Run1")) %>%
  dplyr::filter(Diagnosis %in% c("Healthy_control", "UC", "CD")) %>%
  dplyr::mutate(Diagnosis = factor(Diagnosis, levels = c("Healthy_control", "UC", "CD")))

merged_healthy <- data_5ASA %>%
  left_join(metadata, by = c("filename" = "Metabolomics_FileName_Run1")) |> 
  dplyr::filter(Diagnosis == "Healthy_control")
merged_df_UC <- data_5ASA %>%
  left_join(metadata, by = c("filename" = "Metabolomics_FileName_Run1")) |> 
  dplyr::filter(Diagnosis == "UC")
merged_df_CD <- data_5ASA %>%
  left_join(metadata, by = c("filename" = "Metabolomics_FileName_Run1")) |> 
  dplyr::filter(Diagnosis == "CD")

merged_healthy %>% 
  summarize(n_nonzero = sum(X5ASA != 0, na.rm = TRUE), total = n())
merged_df_UC %>% 
  summarize(n_nonzero = sum(X5ASA != 0, na.rm = TRUE), total = n())
merged_df_CD %>% 
  summarize(n_nonzero = sum(X5ASA != 0, na.rm = TRUE), total  = n())
merged_df_all_clean %>% 
  summarize(n_nonzero = sum(X5ASA != 0, na.rm = TRUE), total = n())

#write_csv(merged_df_UC, "MSV000084775_merged_df_UC.csv")

# Perform pairwise t-test on Feature 1323 by Diagnosis with BH adjustment
pairwise_results_UC <- pairwise.wilcox.test(merged_df_all_clean$X5ASA, merged_df_all_clean$Diagnosis, p.adjust.method = "none")
pairwise_results_UC$p.value
# Create a comparisons list from the unique factor levels in Diagnosis
comparisons_list_UC <- combn(levels(as.factor(merged_df_all_clean$Diagnosis)), 2, simplify = FALSE)

merged_df_UC_1 <- merged_df_all_clean |> 
  dplyr::mutate(X5ASA = if_else(X5ASA == 0, 1, X5ASA), Current_5ASA  = as.factor(Current_5ASA), ASA_Exposure = as.factor(ASA_Exposure), point_color = case_when(Current_5ASA == "1" & ASA_Exposure == "1" ~ "#27AAE1", ASA_Exposure == "1" ~ "#EC008C", Current_5ASA == "Missing" & ASA_Exposure == "Missing"  ~ "#662d91", Current_5ASA == "0" & ASA_Exposure == "0" ~ "grey", TRUE ~ "#1d9c78"))

merged_df_UC_filtered <- merged_df_UC_1 |> 
  dplyr::filter(X5ASA > 0)
#replace by 1.000000e-06 to see how many for each categories 
merged_df_UC_filtered_2 <- merged_df_UC_filtered |> 
  dplyr::select(filename, X5ASA, cholyl_5ASA, Diagnosis, Current_5ASA, ASA_Exposure, point_color)

#write_csv(merged_df_UC_filtered_2, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/MSV000084775_merged_df_UC_filtered_5ASA.csv")

metadata_summary_3 <- merged_df_UC_filtered |> 
  ungroup() %>% 
  group_by(Diagnosis, `Current_5ASA`, ASA_Exposure) |>  
  summarize(count = n(), .groups = "drop")

# Plot with jitter points colored by the custom point_color column
plot_UC <- ggplot(merged_df_UC_filtered_2, aes(x = Diagnosis, y = X5ASA)) +
  geom_boxplot(aes(fill = Diagnosis), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = point_color), width = 0.2, height = 0.1, size = 6, alpha = 0.7, shape = 16) +
  scale_color_identity() +
  scale_y_log10(limits = c(1, 1e6), breaks = c(1, 10, 100, 1000, 10000, 100000, 1e6), labels = c(expression(10^0),expression(10^1),expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
  labs(y = "Peak area (log scale)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title  = element_text(size = 18, face = "bold"), plot.title  = element_text(size = 20, face = "bold", hjust = 0.5), legend.position = "none")
plot_UC

# Save as high-resolution PDF
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/MSV00084775_V5-5ASA.pdf", plot_UC, width = 3, height = 7, dpi = 00)

# calculate Fisher's exact test
df2 <- merged_df_all_clean %>%
  dplyr::mutate(presence = if_else(cholyl_5ASA > 0, "yes", "no"))
diag_lvls <- levels(factor(df2$Diagnosis))
pairs     <- combn(diag_lvls, 2, simplify = FALSE)
fisher_res <- map_dfr(pairs, function(lvls) {
  sub <- df2 %>% filter(Diagnosis %in% lvls)
  tbl <- table(
    factor(sub$Diagnosis, levels = lvls),
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

#write_csv(merged_df_UC_filtered_2, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000084775/MSV000084775_sixmix_feature_info_name.csv")


