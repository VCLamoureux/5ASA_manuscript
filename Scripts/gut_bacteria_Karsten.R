# set working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# import libraries
library(duckplyr)
library(tidyverse)
library(scales)
library(igraph)
library(ggraph)

# import tables
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/hCom_new_members_Feb2025/mzmine_library_output_2/fbmn_2_iimn_fbmn_quant.csv")
feature_table <- feature_table |> 
  dplyr::rename_with(~gsub(" Peak area", "", .))
colnames(feature_table)

drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/hCom_new_members_Feb2025/mzmine_library_output_2/fbmn_2_quant_mzmine")
colnames(drift)[1] <- "Feature"

randomized_sequence <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/hCom_new_members_Feb2025/randomized_sequence_New_members_feb2025.csv")
colnames(randomized_sequence) <- as.character(unlist(randomized_sequence[1, ]))
randomized_sequence <- randomized_sequence[-1, ]
colnames(randomized_sequence)[2] <- "filename"
randomized_sequence$order <- as.numeric(randomized_sequence$order)

annotation <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/hCom_new_members_Feb2025/mzmine_library_output_2/ec12e706a21d4d178983c9380174fb82-merged_results_with_gnps.tsv")
annotation <- annotation |> 
  dplyr::rename(Feature = `#Scan#`)
annotation$Feature <- as.character(annotation$Feature)

info_feature <- feature_table[, 1:3]
info_feature <- info_feature |> 
  dplyr::rename(Feature = `row ID`)
info_feature$Feature <- as.character(info_feature$Feature)
colnames(info_feature) <- c("Feature", "mz", "RT")

info_feature_name <- info_feature |> 
  left_join(annotation, by = "Feature")

data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzML")) |> 
  t() |>  
  as.data.frame() |> 
  rownames_to_column("filename")

# quality control RT drift
feature_row <- drift |>  
  dplyr::filter(Feature == 37054)
rt_columns <- feature_row |>  
  dplyr::select(contains("RT"))
rt_data <- tidyr::gather(rt_columns, key = "filename", value = "RT")
rt_data$filename <- gsub(".mzML Feature RT", "", rt_data$filename)
rt_data_ordered <- rt_data |> 
  left_join(randomized_sequence, by = "filename") |> 
  dplyr::filter(!filename == "5ASA_CA")

rt_data_ordered$order <- as.numeric(rt_data_ordered$order)

# Arrange the data by 'order'
rt_data_ordered <- rt_data_ordered |> 
  arrange(order) |> 
  dplyr::filter(!(str_detect(filename, "Blank|Asn_UDCA|UDCA")))
rt_stats <- rt_data |> 
  summarise(Mean_RT = mean(RT, na.rm = TRUE), 
            SD_RT = sd(RT, na.rm = TRUE), 
            Variance_RT = var(RT, na.rm = TRUE))

ggplot(rt_data_ordered, aes(x = reorder(filename, order), y = RT)) +  # reorder based on 'order' column
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "Retention Time (RT)", title = "Retention Time Drift for Feature 9183")

# quality control peak areas drift
feature_row_peakarea <- drift |>  
  dplyr::filter(Feature == 37054)

peakarea_columns <- feature_row_peakarea |>  
  dplyr::select(contains("Peak area"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data$filename <- gsub(".mzML Peak area", "", peakarea_data$filename)

# Merge the two dfs by 'filename'
peakarea_data_ordered <- peakarea_data |> 
  left_join(randomized_sequence, by = "filename") |> 
  arrange(order) |> 
  dplyr::filter(!str_detect(filename, "sixmix"))

peakarea_stats <- peakarea_data |> 
  summarise(Mean_peakarea = mean(Peak, na.rm = TRUE), 
            SD_peakarea = sd(Peak, na.rm = TRUE), 
            Variance_peakarea = var(Peak, na.rm = TRUE))

ggplot(peakarea_data_ordered, aes(x = reorder(filename, order), y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "peakarea", title = "peakarea Drift for Feature 9183")

# Filter for QCmix
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
  dplyr::filter(!str_detect(filename, "sixmix|Blank|pool|5ASA_CA"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

# remove features belonging to QCmix (Amitriptyline, Coumarin 314, Sulfamethazine, Sulfadimethoxine, Sulfamethizole, Sulfachloropyridazine)
specified_features <- c("42928", "47335", "26350", "37054", "26731", "30058")

feature_to_remove <- blank_feature_info |>  
  left_join(pool_feature_info) |> 
  dplyr::filter(Mean_blank > 0) |> 
  dplyr::mutate(pool_Blank = Mean_pool/Mean_blank) |> 
  dplyr::filter(pool_Blank < 5 | is.na(pool_Blank)) |> 
  dplyr::bind_rows(blank_feature_info |> dplyr::filter(Feature %in% specified_features)) |> 
  dplyr::distinct(Feature, .keep_all = TRUE)

# Remove any features listed in 'feature_to_remove'
data_clean <- data_transpose |> 
  dplyr::select(-c(feature_to_remove$Feature))
  
# Features to be removed samples/QCmix < 5
feature_to_remove_qcmix <- sixmix_feature_info |> 
  left_join(pool_feature_info) |> 
  dplyr::filter(Mean_sixmix > 0) |> 
  dplyr::mutate(pool_Mix = Mean_pool/Mean_sixmix) |> 
  dplyr::filter(pool_Mix < 5 | is.na(pool_Mix)) |> 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with QCmix removal, pool, blanks
data_clean2 <- data_clean |> 
  dplyr::select(-c(feature_to_remove_qcmix$Feature)) |> 
  dplyr::filter(!str_detect(filename, "sixmix|pool|Blank")) |> 
  dplyr::filter(!filename %in% c("5ASA_CA.mzML"))

data_clean_transpose <- data_clean2 |>
  column_to_rownames("filename") |>
  t() |> 
  as.data.frame() |> 
  rownames_to_column("Feature")

merged <- left_join(data_clean_transpose, annotation, by = "Feature")
merged_filtered <- merged |> 
  dplyr::select(Compound_Name, Feature, contains("mzML"))

# RT of the standard 8.90 min feature 47081
# 48736 LCA-5ASA RT 9.40
# 47996 Dihydroxy-5ASA RT 9.16
# 47081 CA-5ASA RT 8.90

# 48762 LCA - 5ASA delta 117 RT 9.40 ###ISF###
# 47066 CA - 5ASA delta 117 RT 8.90 ###ISF###

# TCA 42993 7.49 min
# TLCA 47399 8.99 min 
# TCDCA 45431 8.35 min
# TDCA 45847 8.55 min

merged_filtered_red <- merged_filtered |> 
  dplyr::filter(Feature %in% c("48736", "47996", "47081", "42993", "47399", "45431", "45847", "44983")) |>
  dplyr::select(-Compound_Name) |> 
  pivot_longer(cols = -Feature, names_to = "sample", values_to = "intensity") 

merged_formatted <- merged_filtered_red |> 
  dplyr::mutate(condition = str_remove(sample, "_[123]\\.mzML$"))

merged_taurine <- merged_formatted |> 
  dplyr::filter(str_detect(sample, "taurine")) |> 
  group_by(Feature, condition) |> 
  dplyr::mutate(zero_count = sum(intensity == 0), intensity = if_else(zero_count >= 2, 0, intensity)) |> 
  ungroup() |> 
  dplyr::select(-zero_count) 

merged_taurine_condition <- merged_taurine |> 
  dplyr::mutate(Time_Point = case_when(grepl("T0", sample) ~ "T0", 
                                       grepl("T72", sample) ~ "T72", TRUE ~ NA_character_),
                condition = case_when(
                  grepl("^Cs", sample) ~ "Clostridium sporogenes",
                  grepl("^Ecc", sample) ~ "Enterobacter cloacae subsp. cloacae",
                  grepl("^Ef", sample) ~ "Enterococcus faecalis",
                  grepl("^Lpp", sample) ~ "Lactobacillus plantarum subsp plantarum",
                  grepl("^Pc", sample) ~ "Phocaeicola coprocola",
                  grepl("^Pd", sample) ~ "Phocaeicola dorei",
                  grepl("^Rf", sample) ~ "Ruminococcus flavefaciens",
                  grepl("^Vd", sample) ~ "Veillonella dispar",
                  TRUE ~ NA_character_))

# plot

df_line_plot <- merged_taurine_condition |>
  dplyr::filter(Time_Point %in% c("T0", "T72")) |>
  dplyr::mutate(
    replicate  = str_extract(sample, "(?<=_)\\d(?=\\.mzML$)"),
    intensity  = if_else(intensity == 0, 1e-6, intensity),
    Time_Point = factor(Time_Point, levels = c("T0", "T72"))) |>
  dplyr::filter(Feature == "47081") |>
  group_by(condition, Time_Point) |>
  summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    sd_intensity   = sd(intensity,   na.rm = TRUE),
    .groups        = "drop")

df_line_plot_5ASA_CA <- ggplot(df_line_plot, aes(
  x = Time_Point,
  y = mean_intensity,
  color = condition,
  group = condition
)) +
  geom_line(size = 0.8) +
  geom_point(size = 6) +
  geom_errorbar(
    aes(
      ymin = mean_intensity - sd_intensity,
      ymax = mean_intensity + sd_intensity
    ),
    width = 0.1,
    size  = 0.8,
    alpha = 0.7
  ) +
  scale_y_continuous(limits = c(0, 3E6), 
                     labels = scales::label_scientific(digits = 2)) +
  scale_color_manual(values = c(
    "Clostridium sporogenes"               = "#1B9E77",
    "Enterobacter cloacae subsp. cloacae"   = "#D95F02",
    "Enterococcus faecalis"                = "#F23DE9",
    "Lactobacillus plantarum subsp plantarum" = "#7FCFF7",
    "Phocaeicola coprocola"                = "#EDC948",
    "Phocaeicola dorei"                    = "#CA9A96",
    "Ruminococcus flavefaciens"            = "#FF9DA6",
    "Veillonella dispar"                   = "#27AAE1"
  )) +
  labs(
    x     = "Time Point",
    y     = "Mean Peak Area",
    color = "Species"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position   = "none"
  )

df_line_plot_5ASA_CA
#ggsave("TCDCA.pdf", df_line_plot_5ASA_CA, width = 3, height = 4, dpi = 900)
getwd()
# repeat with the Feature listed above 


## Plot the different bile acids
features_to_include <- annotation |>
  dplyr::filter(LibraryName %in% c("BILELIB19.mgf", "GNPS-BILE-ACID-MODIFICATIONS.mgf")) |>
  pull(Feature) |>
  unique()

merged_all_bas <- merged |> 
  dplyr::select(Compound_Name, Feature, contains("mzML")) |> 
  dplyr::select(-Compound_Name) |>
  dplyr::filter(Feature %in% features_to_include) |> 
  pivot_longer(cols = -Feature, names_to = "sample", values_to = "intensity")

merged_all_bas_extracted <- merged_all_bas |> 
  extract(col = sample, into = c("condition", "Time_Point", "replicate"), regex  = "^(.*)_(T[0-9]+)_(.*)$", remove = FALSE)
merged_all_bas_extracted$Time_Point <- as.numeric(gsub("T", "", merged_all_bas_extracted$Time_Point))
merged_all_bas_extracted$replicate <- as.numeric(gsub(".mzML", "", merged_all_bas_extracted$replicate))

merged_all_bas_formatted <- merged_all_bas_extracted |> 
  dplyr::mutate(Time_Point = as.numeric(Time_Point), replicate  = as.factor(replicate), condition  = as.factor(condition)) |> 
  dplyr::filter(str_detect(sample, "taurine")) |> 
  dplyr::mutate(condition = case_when(
                  grepl("^Cs", sample) ~ "Clostridium sporogenes",
                  grepl("^Ecc", sample) ~ "Enterobacter cloacae subsp. cloacae",
                  grepl("^Ef", sample) ~ "Enterococcus faecalis",
                  grepl("^Lpp", sample) ~ "Lactobacillus plantarum subsp plantarum",
                  grepl("^Pc", sample) ~ "Phocaeicola coprocola",
                  grepl("^Pd", sample) ~ "Phocaeicola dorei",
                  grepl("^Rf", sample) ~ "Ruminococcus flavefaciens",
                  grepl("^Vd", sample) ~ "Veillonella dispar",
                  TRUE ~ NA_character_))
merged_all_bas_formatted$Feature <- as.character(merged_all_bas_formatted$Feature)

annotation_red <- annotation |>
  dplyr::select(Feature, LibraryName, Compound_Name, SpecMZ, SpectrumID)

## Enterococcus faecalis
merged_all_bas_formatted_e_f <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Enterococcus faecalis") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_e_f <- merged_all_bas_formatted_e_f |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## lactobacillus plantarum
merged_all_bas_formatted_lp <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Lactobacillus plantarum subsp plantarum") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_lp <- merged_all_bas_formatted_lp  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## C. sporogenes
merged_all_bas_formatted_cs <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Clostridium sporogenes") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_cs <- merged_all_bas_formatted_cs  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## E. cloacae
merged_all_bas_formatted_ec <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Enterobacter cloacae subsp. cloacae") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_ec <- merged_all_bas_formatted_ec  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## P. dorei
merged_all_bas_formatted_pd <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Phocaeicola dorei") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_pd <- merged_all_bas_formatted_pd  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## V. dispar
merged_all_bas_formatted_vd <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Veillonella dispar") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_vd <- merged_all_bas_formatted_vd  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## Ruminococcus flavefaciens
merged_all_bas_formatted_rf <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Ruminococcus flavefaciens") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_rf <- merged_all_bas_formatted_rf  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()

## Phocaeicola coprocola
merged_all_bas_formatted_pc <- merged_all_bas_formatted |> 
  dplyr::filter(condition == "Phocaeicola coprocola") |> 
  left_join(annotation_red, by = "Feature")

df_plot_mean_pc <- merged_all_bas_formatted_pc  |> 
  group_by(Feature, condition, Time_Point) |> 
  summarise(mean_intensity = mean(intensity, na.rm = TRUE)) |> 
  ungroup()


## Phocaeicola coprocola
df_increase_pc <- df_plot_mean_pc %>%
  group_by(Feature, condition) %>%
  pivot_wider(
    names_from   = Time_Point,
    names_prefix = "T",
    values_from  = mean_intensity
  ) %>%
  dplyr::filter(T0 > 0, T72 >= 3 * T0) %>%
  ungroup()

df_final_pc <- df_plot_mean_pc |>
  semi_join(df_increase_pc, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(
    hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))

## V. dispar
df_increase_vd <- df_plot_mean_vd |>
  group_by(Feature) |>
  pivot_wider(names_from = "Time_Point", 
              names_prefix = "T",
              values_from = "mean_intensity") |>
  dplyr::filter(!is.na(T0), !is.na(T72)) |>
  dplyr::filter(T72 >= 3 * T0, T72 > 0)

df_final_vd <- df_plot_mean_vd |>
  semi_join(df_increase_vd, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(
    hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))

## Ruminococcus flavefaciens
df_increase_rf <- df_plot_mean_rf |>
  group_by(Feature) |>
  pivot_wider(names_from = "Time_Point", 
              names_prefix = "T",
              values_from = "mean_intensity") |>
  dplyr::filter(!is.na(T0), !is.na(T72)) |>
  dplyr::filter(T72 >= 3 * T0, T72 > 0)

df_final_rf <- df_plot_mean_rf |>
  semi_join(df_increase_rf, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(
    hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))

## P. dorei
df_increase_pd <- df_plot_mean_pd |>
  group_by(Feature) |>
  pivot_wider(names_from = "Time_Point", 
              names_prefix = "T",
              values_from = "mean_intensity") |>
  dplyr::filter(!is.na(T0), !is.na(T72)) |>
  dplyr::filter(T72 >= 3 * T0, T72 > 0)

df_final_pd <- df_plot_mean_pd |>
  semi_join(df_increase_pd, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(
    hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))

## lactobacillis plantarum
df_increase_lp <- df_plot_mean_lp |>
  group_by(Feature) |>  
  pivot_wider(names_from = "Time_Point", 
              names_prefix = "T",
              values_from = "mean_intensity") |>
  dplyr::filter(!is.na(T0), !is.na(T72)) |>
  dplyr::filter(T72 >= 3 * T0, T72 > 0)

df_final_lp <- df_plot_mean_lp |>
  semi_join(df_increase_lp, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(
    hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))

## Enterococcus Faecalis
df_increase_ef <- df_plot_mean_e_f |>
  group_by(Feature) |>
  pivot_wider(names_from = "Time_Point",
              names_prefix  = "T",
              values_from = "mean_intensity") |>
  dplyr::filter(!is.na(T0), !is.na(T72)) |>
  dplyr::filter(T72 >= 3 * T0, T72 > 0) |> left_join(annotation_red, by = "Feature") 


df_final_ef <- df_plot_mean_e_f |>
  semi_join(df_increase_ef, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(
    hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))

## C. sporogenes
df_increase_cs <- df_plot_mean_cs |>
  group_by(Feature) |>
  pivot_wider(names_from = "Time_Point", 
              names_prefix = "T",
              values_from = "mean_intensity") |>
  dplyr::filter(!is.na(T0), !is.na(T72)) |>
  dplyr::filter(T72 >= 3 * T0, T72 > 0)

df_final_cs <- df_plot_mean_cs |>
  semi_join(df_increase_cs, by = "Feature") |> 
  left_join(annotation_red, by = "Feature") |> 
  dplyr::filter(Time_Point == "T72") |> 
  dplyr::filter(!str_detect(Compound_Name, "18.01|36.02|54.03|88.99|72.04")) |>
  dplyr::mutate(delta_mass = if_else(
    str_detect(Compound_Name, "delta mass"),
    str_extract(Compound_Name, "(?<=delta mass )[-+]?[0-9]*\\.?[0-9]+"),
    Compound_Name
  )) |>
  dplyr::mutate(hydroxylation = case_when(
      str_detect(Compound_Name, regex("DCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("HDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("CDCA", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("gMCA", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("lithocholic acid", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("taurocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hyocholic acid", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("chenodeoxycholic acid", ignore_case = TRUE)) ~ "di",
      str_detect(Compound_Name, regex("monohydroxylated", ignore_case = TRUE)) ~ "mono",
      str_detect(Compound_Name, regex("\\bCA\\b", ignore_case = TRUE)) &
        !str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ "tri",
      str_detect(Compound_Name, regex("hydroxylated", ignore_case = TRUE)) ~ 
        tolower(str_match(Compound_Name, "(?i)(mono|di|tri|tetra|penta)hydroxylated")[,2]),
      TRUE ~ NA_character_))



# Separate the data for both Bile acid library
## E. feacalis network
df_bile <- df_final_ef |>
  dplyr::filter(LibraryName == "BILELIB19.mgf")

df_gnps <- df_final_ef |>
  dplyr::filter(LibraryName == "GNPS-BILE-ACID-MODIFICATIONS.mgf")

bilelib_compounds <- unique(df_bile$delta_mass)
gnps_compounds    <- unique(df_gnps$delta_mass)

center_node <- "E. faecalis"
group_bile  <- "BILELIB19"
group_gnps  <- "GNPS"

edges_center <- data.frame(
  from = c(center_node, center_node),
  to   = c(group_bile, group_gnps),
  stringsAsFactors = FALSE)

edges_bile <- data.frame(
  from = rep(group_bile, length(bilelib_compounds)),
  to   = bilelib_compounds,
  stringsAsFactors = FALSE)

edges_gnps <- data.frame(
  from = rep(group_gnps, length(gnps_compounds)),
  to   = gnps_compounds,
  stringsAsFactors = FALSE)

edges <- dplyr::bind_rows(edges_center, edges_bile, edges_gnps)

g <- graph_from_data_frame(edges, directed = TRUE)

ggraph(g, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(color = "gray40") +
  geom_node_point(aes(filter = (name == center_node)), size = 6, color = "darkred") +  # emphasize center
  geom_node_point(aes(filter = (name %in% c(group_bile, group_gnps))), size = 4, color = "darkorange") +
  geom_node_point(aes(filter = (!name %in% c(center_node, group_bile, group_gnps))), size = 3, color = "steelblue") +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 size = 3, 
                 color = "black") +
  theme_void() +
  ggtitle("E. faecalis network")

#ggsave("Radial_Network_ef.pdf", plot = p, width = 11, height = 10, dpi = 900)
#getwd()
# repeat for the other microbe by changing df_final_ef by any other microbes

df_bile <- df_final_lp |> dplyr::filter(LibraryName == "BILELIB19.mgf")
#df_bile <- df_final_rf |> dplyr::filter(LibraryName == "BILELIB19.mgf")
#df_bile <- df_final_pc |> dplyr::filter(LibraryName == "BILELIB19.mgf")
#df_bile <- df_final_vd |> dplyr::filter(LibraryName == "BILELIB19.mgf")

df_gnps <- df_final_lp |> dplyr::filter(LibraryName == "GNPS-BILE-ACID-MODIFICATIONS.mgf")
#df_gnps <- df_final_rf |> dplyr::filter(LibraryName == "GNPS-BILE-ACID-MODIFICATIONS.mgf")
#df_gnps <- df_final_pc |> dplyr::filter(LibraryName == "GNPS-BILE-ACID-MODIFICATIONS.mgf")
#df_gnps <- df_final_vd |> dplyr::filter(LibraryName == "GNPS-BILE-ACID-MODIFICATIONS.mgf")

bilelib_compounds <- unique(df_bile$delta_mass)
gnps_compounds    <- unique(df_gnps$delta_mass)

center_node <- "L. plantarum"
group_bile  <- "BILELIB19"
group_gnps  <- "GNPS"

edges_center <- data.frame(
  from = c(center_node, center_node),
  to   = c(group_bile, group_gnps),
  stringsAsFactors = FALSE)

edges_bile <- data.frame(
  from = rep(group_bile, length(bilelib_compounds)),
  to   = bilelib_compounds,
  stringsAsFactors = FALSE)

edges_gnps <- data.frame(
  from = rep(group_gnps, length(gnps_compounds)),
  to   = gnps_compounds,
  stringsAsFactors = FALSE)

edges <- dplyr::bind_rows(edges_center, edges_bile, edges_gnps)
g <- graph_from_data_frame(edges, directed = TRUE)

ggraph(g, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(color = "gray40") +
  geom_node_point(aes(filter = (name == center_node)), size = 6, color = "darkred") +  # emphasize center
  geom_node_point(aes(filter = (name %in% c(group_bile, group_gnps))), size = 4, color = "darkorange") +
  geom_node_point(aes(filter = (!name %in% c(center_node, group_bile, group_gnps))), size = 3, color = "steelblue") +
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 size = 3, 
                 color = "black") +
  theme_void() +
  ggtitle("L. plantarum network")
#ggsave("Radial_Network_ef.pdf", plot = p, width = 11, height = 10, dpi = 900)
#getwd()



## 1) Filter feature 47081  -------------------------------------------------
## 1) Filter feature 47081  -------------------------------------------------
df_47081 <- merged_taurine_condition %>%
  dplyr::filter(
    Feature == "47081",
    Time_Point %in% c("T72"),
    !is.na(intensity)) %>%
  dplyr::mutate(
    intensity  = if_else(is.na(intensity), 0, intensity),
    Time_Point = factor(Time_Point, levels = c("T72"))) |> 
  dplyr::mutate(
    x_group = interaction(condition, Time_Point, sep = "_"),
    x_label = paste0(condition, "\n", Time_Point))
  
write_csv(df_47081, "OneDrive-UniversityofCalifornia,SanDiegoHealth/df_47081_all_conditions_T72.csv")

## 2) BUILD SAME COLOR PALETTE AS p_173005  ---------------------------------
n_cond <- length(unique(df_47081$condition))

custom_colors        <- scales::hue_pal(l = 70, c = 60)(n_cond)
custom_colors_darker <- scales::hue_pal(l = 40, c = 80)(n_cond)

# assign names so ggplot matches correctly
names(custom_colors)        <- levels(factor(df_47081$condition))
names(custom_colors_darker) <- levels(factor(df_47081$condition))

p_47081_grouped <- ggplot(
  df_47081,
  aes(x = condition, y = intensity, fill = condition, color = condition)) +
  geom_boxplot(
    aes(group = interaction(condition, Time_Point)),
    color    = "black",
    fill     = NA,         # remove white overlay
    size     = 0.5,
    alpha    = 0.6,
    width    = 0.55,
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(group = interaction(condition, Time_Point)),
    shape    = 19,
    size     = 4,
    alpha    = 1,
    stroke   = NA,
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width  = 0.8
    )
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors_darker) +
  
  # --- FIXED Y-AXIS ---
  scale_y_continuous(
    limits = c(0, 3e6),   # adjust as needed
    breaks = c(1e6, 2e6, 3e6),
    labels = scales::label_scientific(digits = 1)
  ) +
  
  labs(
    x = "Condition",
    y = "Raw peak area"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 1, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line    = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks   = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position = "none"
  )

p_47081_grouped


#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Main_project/hCOM1_SynCOm/hCom_new_members_Feb2025/CA_5ASA_new_members_t0.pdf", p_47081_grouped, width = 4, height = 7, dpi = 900)
#getwd()

