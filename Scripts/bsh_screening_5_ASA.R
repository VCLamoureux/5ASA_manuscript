# set working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# import libraries
library(tidyverse)
library(pheatmap)
library(duckplyr)

# import tables
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/mzmine/bsh_engineered_iimn_fbmn_quant.csv")
feature_table <- feature_table |>
  dplyr::rename_with(~ gsub(" Peak area$", "", .x)) |>
  dplyr::rename_with(~ str_replace_all(.x, c("bu1" = "bv1", "bu2" = "bv2")))

drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/mzmine/bsh_engineered_quant_mzmine.csv")
colnames(drift)[1] <- "Feature"

randomized_sequence <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/randomized_bsh_engineered_20122025.csv")
colnames(randomized_sequence) <- as.character(unlist(randomized_sequence[1, ]))
randomized_sequence <- randomized_sequence[-1, ]
colnames(randomized_sequence)[2] <- "filename"
randomized_sequence$order <- as.numeric(randomized_sequence$order)

randomized_sequence <- randomized_sequence  |>
  dplyr::mutate(filename = filename |> str_replace_all(c("bu1" = "bv1", "bu2" = "bv2")))

annotation <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/94afdef1461e4bbcbbe6e56b7509ac6f-merged_results_with_gnps.tsv")
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
  dplyr::filter(Feature == 25129)
rt_columns <- feature_row |>  
  dplyr::select(contains("RT"))
rt_data <- tidyr::gather(rt_columns, key = "filename", value = "RT")
rt_data$filename <- gsub(".mzML Feature RT", "", rt_data$filename)
rt_data_ordered <- rt_data |> 
  left_join(randomized_sequence, by = "filename")

rt_data_ordered$order <- as.numeric(rt_data_ordered$order)

# Arrange the data by 'order'
rt_data_ordered <- rt_data_ordered |> 
  arrange(order)
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
  dplyr::filter(Feature == 25129)

peakarea_columns <- feature_row_peakarea |>  
  dplyr::select(contains("Peak area"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data$filename <- gsub(".mzML Peak area", "", peakarea_data$filename)

# Merge the two dfs by 'filename'
peakarea_data_ordered <- peakarea_data |> 
  left_join(randomized_sequence, by = "filename") |> 
  arrange(order) |> 
  dplyr::filter(!str_detect(filename, "sixmix|pool|Blank"))

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
specified_features <- c("29037", "32891", "14175", "25129", "14439", "17531")

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
  dplyr::filter(!str_detect(filename, "sixmix|pool|Blank"))

metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/metadata_bsh.csv")
metadata <- metadata  |>
  dplyr::mutate(filename = filename |> str_replace_all(c("bu1" = "bv1", "bu2" = "bv2")))

data_clean_3 <- data_clean2 |> 
  dplyr::select(filename, '32781', '34496', '35808', '30694', '31534', '33336', '29287', '31801', '33344', '31518') |> 
  left_join(metadata, by = "filename") |> 
  dplyr::filter(!str_detect(filename,"ls"))

feature_cols <- c(
  "32781","34496","35808","30694","31534",
  "33336","29287","31801","33344", "31518")

data_clean_3 <- data_clean_3 %>%
  group_by(condition, time) %>%
  dplyr::mutate(across(all_of(feature_cols),
      ~ {
        n_zero <- sum(.x == 0, na.rm = TRUE)
        if (n_zero >= 2) {
          rep(0, length(.x))
        } else {
          .x
        }
      })) %>%ungroup()

df_hm <- data_clean_3 |>
  dplyr::rename(
    CA_5ASA        = `32781`,
    DCA_5ASA       = `34496`,
    LCA_5ASA       = `35808`, 
    His_DCA        = `30694`,
    Cadaverine_DCA = `31534`,
    Leu_ile_CA     = `33336`,
    TCA            = `29287`,
    TDCA           = `31801`,
    TLCA           = `33344`,
    GABA_CA        = `31518`)

# Rename df
df_hm_main <- df_hm

# define order for the pheatmap
met_order <- c("TLCA","TDCA","TCA","LCA_5ASA","DCA_5ASA","CA_5ASA", "His_DCA","Cadaverine_DCA","Leu_ile_CA", "GABA_CA")
met_order <- met_order[met_order %in% colnames(df_hm)]
df_hm_sum <- df_hm_main |>
  dplyr::mutate(time = as.numeric(time)) |>
  group_by(condition, time) |>
  summarise(across(all_of(met_order), ~mean(.x, na.rm = TRUE)), .groups = "drop") |>
  arrange(time, condition) |>
  dplyr::mutate(sample = paste(condition, paste0("T", time), sep = " | ")) |>
  dplyr::select(sample, all_of(met_order))

mat <- df_hm_sum |>
  column_to_rownames("sample") |>
  as.matrix() |> 
  t()
log_mat <- log(1 + mat)

# export for source data table
log_mat <- as.data.frame(log_mat)
log_mat_df <- log_mat |> 
  rownames_to_column("Condition")
#write_csv(log_mat_df, "OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/bsh_heatmap_log_matrix.csv")

df_long <- df_hm |>
  dplyr::mutate(time = as.numeric(time)) |>
  pivot_longer(cols = all_of(met_order), names_to = "metabolite", values_to = "value") |>
  dplyr::filter(time %in% c(0, 48))

df_levene <- df_long |>
  dplyr::group_by(condition, metabolite) |>
  dplyr::summarise(p_levene = {
      d <- dplyr::tibble(
        value = log1p(value),
        time  = factor(time)
      )
      
      # Need ≥ 2 observations per group
      if (all(table(d$time) >= 2)) {
        tryCatch(
          car::leveneTest(value ~ time, data = d, center = median)[["Pr(>F)"]][1],
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      }
    },.groups = "drop")

df_levene

df_stats <- df_long %>%
  group_by(condition, metabolite) %>%
  summarise(
    v0  = value[time == 0],
    v48 = value[time == 48],
    p_value = {
      x <- log1p(v0)
      y <- log1p(v48)
      if (length(x) >= 2 && length(y) >= 2) {
        t.test(x, y, var.equal = TRUE)$p.value
      } else NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(condition) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    star  = case_when(
      is.na(p_adj)  ~ "",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "")) %>%
  ungroup()

df_stats <- df_long |>
  group_by(condition, metabolite) |>
  summarise(
    v0  = list(value[time == 0]),
    v48 = list(value[time == 48]),
    p_value = {
      v0  <- value[time == 0]
      v48 <- value[time == 48]
      if (all(v0 == 0) && any(v48 > 0)) {0.001
      } else if (length(v0) >= 2 && length(v48) >= 2) {
        tryCatch(t.test(v0, v48)$p.value, error = function(e) NA_real_)
      } else {NA_real_}}, .groups = "drop") |>
  dplyr::mutate(star = case_when(is.na(p_value)  ~ "", p_value < 0.001 ~ "***", p_value < 0.01  ~ "**", p_value < 0.05  ~ "*", TRUE ~ ""))

hm_cols <- colnames(log_mat)
hm_rows <- rownames(log_mat)
star_mat <- matrix("", nrow = nrow(log_mat), ncol = ncol(log_mat), dimnames = dimnames(log_mat))

for (j in seq_len(ncol(star_mat))) {
  col_name <- hm_cols[j]
  if (!grepl("T48$", col_name)) next
  cond <- sub(" \\| T48$", "", col_name)
  for (i in seq_len(nrow(star_mat))) {
    met <- hm_rows[i]
    s <- df_stats |>
      filter(condition == cond, metabolite == met) |>
      pull(star)
    if (length(s) == 1) star_mat[i, j] <- s}}

# Define custom gradient colors and plot heatmap
gradient_colors <- colorRampPalette(c("#FFFFFF", "#C7D6F0", "#EBB0A6"))(30)
bsh_plot <- pheatmap(
  log_mat,
  color = gradient_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 90,
  fontsize = 10,
  cellwidth = 18,
  cellheight = 14,
  border_color = "#D3D3D3",
  legend = TRUE,
  display_numbers = star_mat,
  number_color = "black",
  fontsize_number = 12)

#ggsave("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/bsh_engineered/bsh_full_2.pdf", plot = bsh_plot, width = 10, height = 10, dpi = 900)












