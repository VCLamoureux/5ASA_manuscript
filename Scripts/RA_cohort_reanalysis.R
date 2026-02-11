setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# import packages
library(tidyverse)
library(pheatmap)
library(reshape2)
library(ggbeeswarm)
library(ggrepel)
library(duckplyr)
# Import tables
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/MZmine_3_01292025/fbmn_iimn_fbmn_quant.csv")
feature_table <- feature_table |> 
  dplyr::rename_with(~gsub(" Peak area", "", .))
colnames(feature_table)

# Import tables for RT and peak area drift of IS
drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/MZmine_3_01292025/fbmn_quant_mzmine.csv")
colnames(drift)[1] <- "Feature"

# import ReDU formatted metadata
metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/ReDU_Diet_Nov19_Template.csv")
metadata <- metadata |> 
  dplyr::select(-MassiveID) 

# import clinical metadata and format them
clinical_metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/clinical_data_20patients_updated_may20_final_response_plus_diet_score.csv")
clinical_metadata_renamed_1 <- clinical_metadata |> 
  dplyr::mutate(Patient = paste0("F", Patient))
clinical_metadata_renamed_1$Patient <- gsub("FD010d-0", "FD010d0", clinical_metadata_renamed_1$Patient)
clinical_metadata_renamed_1$Patient <- gsub("FD001D0", "FD001d0", clinical_metadata_renamed_1$Patient)
clinical_metadata_renamed_plasma <- clinical_metadata |> 
  dplyr::mutate(Patient = paste0("P", Patient))
clinical_metadata_renamed_plasma$Patient <- gsub("PD010d-0", "PD010d0", clinical_metadata_renamed_plasma$Patient)
clinical_metadata_renamed_plasma$Patient <- gsub("PD001D0", "PD001d0", clinical_metadata_renamed_plasma$Patient)
clinical_metadata_renamed_2 <- clinical_metadata_renamed_1 |> 
  dplyr::mutate(Patient = str_replace(Patient, "d0", "_T2"),
                Patient = str_replace(Patient, "d14", "_T3"))
clinical_metadata_renamed_plasma_2 <- clinical_metadata_renamed_plasma |> 
  dplyr::mutate(Patient = str_replace(Patient, "d0", "_T2"),
                Patient = str_replace(Patient, "d14", "_T3"))
clinical_metadata_renamed_2 <- clinical_metadata_renamed_2 |> 
  dplyr::rename(filename = Patient)
clinical_metadata_renamed_plasma_2 <- clinical_metadata_renamed_plasma_2 |> 
  dplyr::rename(filename = Patient)

# import annotation 
annotation <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/Annotation_Marta_RA_29_January_2025.tsv")
annotation$`#Scan#` <- as.character(annotation$`#Scan#`)
annotation <- annotation |> 
  dplyr::rename(Feature = `#Scan#`)
annotation$Feature <- as.character(annotation$Feature)

# get info features
info_feature <- feature_table[, 1:3]
info_feature <- info_feature |> 
  dplyr::rename(Feature = `row ID`)

info_feature$Feature <- as.character(info_feature$Feature)
colnames(info_feature) <- c("Feature", "mz", "RT")

info_feature_name <- info_feature |> 
 left_join(annotation, by = "Feature")

data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzXML")) |> 
  t() |>  
  as.data.frame() |> 
  rownames_to_column("filename")

# Check retention drift on internal standard Sulfamethizole
feature_row <- drift |>  
  dplyr::filter(Feature == 2394)
rt_columns <- feature_row |> 
  dplyr::select(contains("Feature RT"))
rt_data <- tidyr::gather(rt_columns, key = "filename", value = "RT") 

rt_data_cleaned <- rt_data |> 
  dplyr::filter(filename != 0) |> 
  dplyr::filter(!str_detect(filename, "Blank"))
rt_data_cleaned$RT <- as.numeric(rt_data_cleaned$RT)

rt_stats <- rt_data_cleaned |>
  summarise(Mean_RT = mean(RT, na.rm = TRUE), SD_RT = sd(RT, na.rm = TRUE), Variance_RT = var(RT, na.rm = TRUE))
ggplot(rt_data_cleaned, aes(x = filename, y = RT)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "Retention Time (RT)", title = "Retention Time Drift for Feature 9233")

# Check for intensity drift on internal standard Sulfamethizole
feature_row_peakarea <- drift |> 
  dplyr::filter(Feature == 2394)

peakarea_columns <- feature_row_peakarea |> 
  dplyr::select(contains("Peak area"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data_cleaned <- peakarea_data |> 
  #dplyr::filter(filename != 0) |> 
  dplyr::filter(!str_detect(filename, "Blank"))

peakarea_stats <- peakarea_data_cleaned |>
  summarise(Mean_peakarea = mean(Peak, na.rm = TRUE), SD_peakarea = sd(Peak, na.rm = TRUE), Variance_peakarea = var(Peak, na.rm = TRUE))

ggplot(peakarea_data_cleaned, aes(x = filename, y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "peakarea", title = "peakarea Drift for Feature 9233")

data_sixmix <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "6Mix", filename))
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
  dplyr::filter(str_detect(filename, "PD|almond|black|bread|chia|FD|flaxseeds|ginger|miso|oats|sesame|tahini|turmeric|vinegar"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

specified_features <- c("5552", "2324", "4471","2394", "2872", "7140")

feature_to_remove <- blank_feature_info %>% 
  left_join(sample_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>%
  dplyr::mutate(sample_Blank = Mean_sample / Mean_blank) %>%
  dplyr::filter(sample_Blank < 5 | is.na(sample_Blank)) %>%
  dplyr::bind_rows(blank_feature_info %>% dplyr::filter(Feature %in% specified_features)) %>%
  dplyr::distinct(Feature, .keep_all = TRUE)  

data_clean <- data_transpose |>
  dplyr::select(-dplyr::any_of(feature_to_remove$Feature)) |>
  dplyr::filter(!str_detect(filename, "6Mix|Blank"))

data_clean_transpose <- data_clean |>
  tibble::column_to_rownames("filename") |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Feature") |>
  dplyr::mutate(Feature = as.character(Feature))

merged <- dplyr::left_join(data_clean_transpose, annotation, by = "Feature")

# keep ID columns + only sample columns with _T1/_T2/_T3 (which are the mzXML files)
merged_mzml_mean_all <- merged |>
  dplyr::select(dplyr::any_of(c("Feature", "Compound_Name", "LibraryName")), dplyr::matches("(_T1|_T2|_T3)\\.mzXML$", ignore.case = TRUE))
mzml_columns <- grep("(_T1|_T2|_T3)\\.mzXML$", names(merged_mzml_mean_all), value = TRUE,ignore.case = TRUE)

df_melted <- reshape2::melt(
  merged_mzml_mean_all,
  id.vars = c("Feature", "Compound_Name", "LibraryName"),
  measure.vars = mzml_columns,
  variable.name = "Sample_Name",
  value.name = "Peak_Area") |>
  dplyr::mutate(Peak_Area = as.numeric(Peak_Area))

df_melted_feces <- df_melted |>
  dplyr::filter(str_detect(Sample_Name, "FD"))
df_melted_feces$Sample_Name <- gsub(".mzXML", "", df_melted_feces$Sample_Name)
df_melted_feces_clinical <- df_melted_feces |>
  dplyr::left_join(clinical_metadata_renamed_2, by = c("Sample_Name" = "filename"))

name_map <- c(
  "6052"  = "Cholyl_5ASA_1",
  "7049"  = "Cholyl_5ASA_2",
  "6885"  = "Cholyl_5ASA_3",
  "8986"  = "Lithocholyl_5ASA_1",
  "9248"  = "Lithocholyl_5ASA_2",
  "8182"  = "Deoxycholyl_5ASA_1",
  "7891"  = "Deoxycholyl_5ASA_1",
  "7044"  = "Deoxycholyl_5ASA_2",
  "7405"  = "Deoxycholyl_5ASA_3")

feature_ids <- names(name_map)

df_melted_conjugate <- df_melted_feces_clinical |>
  dplyr::filter(Feature %in% feature_ids) |>
  dplyr::filter(!stringr::str_detect(Sample_Name, "_T1$|_T2$")) |>
  dplyr::mutate(
    Peak_Area = as.numeric(Peak_Area),
    sulfasalazine = dplyr::if_else(Sample_Name %in% c("FD003_T3", "FD028_T3"), "user", "non-user"),
    Feature = dplyr::recode(as.character(Feature), !!!name_map, .default = as.character(Feature)),
    BA_hydroxylate = stringr::str_remove(Feature, "_\\d+$"))

df_melted_conjugate_sum <- df_melted_conjugate |>
  dplyr::group_by(Sample_Name, BA_hydroxylate) |>
  dplyr::summarise(Peak_Area_sum = sum(Peak_Area, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(df_melted_conjugate |> dplyr::select(Sample_Name, sulfasalazine, Pain_50_improvement, DMARD) |>
      dplyr::distinct(), by = "Sample_Name") |>
  dplyr::mutate(Peak_Area_Log_sum = log10(Peak_Area_sum + 1))

df_plot <- df_melted_conjugate_sum |>
  dplyr::mutate(
    sulfasalazine = factor(sulfasalazine, levels = c("non-user", "user")),
    BA_hydroxylate = factor(BA_hydroxylate),
    Peak_Area_Log_sum = log10(as.numeric(Peak_Area_sum) + 1))

p_facet <- ggplot(df_plot,
  aes(x = sulfasalazine, y = Peak_Area_Log_sum, color = sulfasalazine)) +
  geom_jitter(width = 0.15, size = 9, alpha = 1, stroke = NA, na.rm = TRUE) +
  facet_wrap(~ BA_hydroxylate, scales = "free_y", ncol = 4, strip.position = "bottom") +
  scale_color_manual(values = c("non-user" = "gray50", "user" = "#d1495b")) +
  labs(x = "Sulfasalazine", y = "log10(summed peak area + 1)", color = NULL) +
  coord_cartesian(ylim = c(0, 8)) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.placement = "outside")

p_facet

# export the formatted table, if needed
#write_csv(df_melted_conjugate_sum, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/df_melted_conjugate_sum_RA_COHORT.csv")
#ggsave(filename = "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/drug_metabolites_5asa.pdf", plot = p_facet, device = "pdf", width = 10, height = 7, units = "in", dpi = 900)
#getwd()
