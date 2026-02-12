# set working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# import libraries
library(duckplyr)
library(tidyverse)
library(scales)
library(igraph)
library(ggraph)

# import tables
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/mzmine_output/bifidos_iimn_fbmn_quant.csv")
feature_table <- feature_table |> 
  dplyr::rename_with(~gsub(" Peak area", "", .))
colnames(feature_table)

drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/mzmine_output/bifidos_quant_mzmine.csv")
colnames(drift)[1] <- "Feature"

randomized_sequence <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/randomized_sequence_Bifidos_5ASA_CA_12122025.csv")
colnames(randomized_sequence) <- as.character(unlist(randomized_sequence[1, ]))
randomized_sequence <- randomized_sequence[-1, ]
colnames(randomized_sequence)[2] <- "filename"
randomized_sequence$order <- as.numeric(randomized_sequence$order)

annotation <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/a876e775080a4cabacf44578463e8bbe-merged_results_with_gnps.tsv")
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
  dplyr::filter(Feature == 21022)
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
  dplyr::filter(!str_detect(filename, "5-ASA_CA_std"))
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
  dplyr::filter(Feature == 21022)

peakarea_columns <- feature_row_peakarea |>  
  dplyr::select(contains("Peak area"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data$filename <- gsub(".mzML Peak area", "", peakarea_data$filename)

# Merge the two dfs by 'filename'
peakarea_data_ordered <- peakarea_data |> 
  left_join(randomized_sequence, by = "filename") |> 
  arrange(order) |> 
  dplyr::filter(!str_detect(filename, "sixmix|pool|Blank|5-ASA_CA_std"))

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
  dplyr::filter(str_detect(filename, "Blank")) |>
  dplyr::filter(!filename %in% c("Blank_1_20251215170622.mzML","Blank_7.mzML"))

blank_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                 Mean_blank = data_blank |> column_to_rownames("filename") |> colMeans(),
                                 SD_blank =  data_blank |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |> arrange(desc(Mean_blank))

blank_feature_info_name <- left_join(blank_feature_info, info_feature_name, by = "Feature" )

# Filter for samples
data_sample <- data_transpose |> 
  dplyr::filter(!str_detect(filename, "sixmix|Blank|pool|5-ASA_CA"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

# remove features belonging to QCmix (Amitriptyline, Coumarin 314, Sulfamethazine, Sulfadimethoxine, Sulfamethizole, Sulfachloropyridazine)
specified_features <- c("26765", "10933", "28376", "21022", "11139", "13699")

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
  dplyr::filter(!filename %in% c("5-ASA_CA_std.mzML")) |> 
  dplyr::filter(str_detect(filename, "subs"))

metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/metadata_bifidos.csv")

data_clean_3 <- data_clean2 |> 
  dplyr::select(filename, '28355') |> 
  left_join(metadata, by = "filename")

df_bile <- data_clean_3 %>%
  dplyr::transmute(
    filename,
    condition = as.character(condition),
    time_h    = time_h,
    value     = as.numeric(`28355`)
  ) %>%
  tidyr::drop_na(condition, value)

# Put control first
control_name <- "medium"
other_conds  <- setdiff(sort(unique(df_bile$condition)), control_name)

df_bile_1 <- df_bile %>%
  dplyr::mutate(
    condition = factor(condition, levels = c(control_name, other_conds))) |> 
  dplyr::filter(time_h == 72)

#write_csv(df_bile_1, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/bifidos_feature28355_monoculture_data.csv")


# Colors
n_cond <- nlevels(df_bile_1$condition)
custom_colors        <- scales::hue_pal(l = 70, c = 60)(n_cond)
custom_colors_darker <- scales::hue_pal(l = 40, c = 80)(n_cond)
names(custom_colors)        <- levels(df_bile_1$condition)
names(custom_colors_darker) <- levels(df_bile_1$condition)

# Plot
p_monoculture <- ggplot(
  df_bile_1,
  aes(x = condition, y = value, fill = condition, color = condition)
) +
  geom_boxplot(
    size  = 0.5,
    alpha = 0.6,
    width = 0.55,
    outlier.shape = NA
  ) +
  geom_jitter(
    shape    = 19,
    size     = 5,
    alpha    = 1,
    stroke   = NA,
    position = position_jitter(width = 0.2)
  ) +
  scale_fill_manual(values  = custom_colors) +
  scale_color_manual(values = custom_colors_darker) +
  labs(
    x = "Condition",
    y = "Feature 28355 (intensity)"
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

p_monoculture
#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Bifidos_5ASA_BAs_cultures/bifidos_28355_monoculture_boxplot.pdf", plot = p_monoculture, width = 4, height = 7, dpi = 900)

















