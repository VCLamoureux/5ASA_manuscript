# set working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/cCom/")

# import libraries
library(duckplyr)
library(readxl)
library(tidyverse)
library(pheatmap)
library(mixOmics)
library(vegan)
library(caret)
library(ggpubr)
library(stringr)
library(car)
library(scales)
library(rstatix)
library(ggh4x)
library(igraph)
library(ggraph)

# import tables
feature_table <- read_csv("cCom_iimn_fbmn_quant.csv")
feature_table <- feature_table |> 
  dplyr::rename_with(~gsub(" Peak area", "", .))
colnames(feature_table)

randomized_sequence <- read_csv("cCOM_randomized_aequence_24_august_2025.csv")
colnames(randomized_sequence) <- as.character(unlist(randomized_sequence[1, ]))
randomized_sequence <- randomized_sequence[-1, ]
colnames(randomized_sequence)[2] <- "filename"
randomized_sequence$order <- as.numeric(randomized_sequence$order)

annotation <- read_tsv("babec751a2c541c6844dc602c19bf379-merged_results_with_gnps.tsv")
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
  rownames_to_column("filename") #|> dplyr::select('filename', '174170', '179381', '173450')


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
  dplyr::filter(str_detect(pattern = "pool|ID_", filename))
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
  dplyr::filter(!str_detect(filename, "sixmix|Blank|pool|ID_"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

# remove features belonging to QCmix (Amitriptyline, Coumarin 314, Sulfamethazine, Sulfadimethoxine, Sulfamethizole, Sulfachloropyridazine)
specified_features <- c("146177", "174380", "72728", "114132", "73907", "85329")

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
  dplyr::filter(Feature %in% c("182394", "124943", "143724", "188843", "77301", "156611", "173005", "174170", "182805")) |>
  dplyr::select(-Compound_Name) |> 
  pivot_longer(cols = -Feature, names_to = "sample", values_to = "intensity") 

merged_formatted <- merged_filtered_red |> 
  dplyr::mutate(condition = str_remove(sample, "_[123]\\.mzML$"))

merged_taurine <- merged_formatted |> 
  dplyr::filter(str_detect(sample, "drug")) |> 
  group_by(Feature, condition) |> 
  dplyr::mutate(zero_count = sum(intensity == 0), intensity = if_else(zero_count >= 2, 0, intensity)) |> 
  ungroup() |> 
  dplyr::select(-zero_count) 

merged_taurine_condition <- merged_taurine |> 
  dplyr::mutate(Time_Point = case_when(grepl("T0", sample) ~ "T0", 
                                       grepl("T72", sample) ~ "T72", TRUE ~ NA_character_),
                condition = case_when(
                  grepl("CDHMUM8", sample) ~ "Bacteroides_caccae_CDHMUM8",
                  grepl("HI_110723", sample) ~ "Bacteroides_fragilis_HI_110723_9",
                  grepl("A48", sample) ~ "Bacteroides_luhongzhouii_A48",
                  grepl("CD5", sample) ~ "Bacteroides_thetaiotaomicron_CD5",
                  grepl("CD15AnD", sample) ~ "Bacteroides_uniformis_CD15AnD",
                  grepl("A15UNF", sample) ~ "Bacteroides_vulgatus_A15UNF",
                  grepl("BruPT-FE", sample) ~ "BruPT-FE",
                  grepl("CrfA", sample) ~ "Clostridium_innocuum_CrfA",
                  grepl("CD11", sample) ~ "Clostridium_symbiosum_CD11",
                  grepl("CD3", sample) ~ "Collinsella_aerofaciens_CD3",
                  grepl("A14CFS", sample) ~ "Enterocloster_bolteae_A14CFS",
                  grepl("CD4", sample) ~ "Erysipelatoclostridium_ramosum_CD4",
                  grepl("A04A", sample) ~ "Parabacteroides_distasonis_A04A",
                  grepl("RCM", sample) ~ "RCM",
                  grepl("A2", sample) ~ "Ruminococcus_gnavus_A2",
                  grepl("cCOM", sample) ~ "cCOM",
                  grepl("BHIS", sample) ~ "BHIS",
                  TRUE ~ NA_character_))

# plot
target_feature <- 173005

df_line_plot <- merged_taurine_condition |>
  dplyr::filter(Time_Point %in% c("T0","T72"),
                Feature == target_feature,
                !is.na(condition)) |>
  dplyr::mutate(
    replicate  = condition,
    intensity  = dplyr::if_else(is.na(intensity) | intensity == 0, 0, intensity),
    Time_Point = factor(Time_Point, levels = c("T0","T72"))) |>
  dplyr::group_by(condition, Time_Point) |>
  dplyr::summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    sd_intensity   = sd(intensity,   na.rm = TRUE),
    .groups        = "drop")


cols_microbiomeMASST <- c(
  # Bacteroides (blues)
  "Bacteroides_caccae_CDHMUM8"            = "#1F77B4",
  "Bacteroides_fragilis_HI_110723_9"      = "#2E86C1",
  "Bacteroides_luhongzhouii_A48"          = "#5DADE2",
  "Bacteroides_thetaiotaomicron_CD5"      = "#2874A6",
  "Bacteroides_uniformis_CD15AnD"         = "#7FB3D5",
  "Bacteroides_vulgatus_A15UNF"           = "#154360",
  
  # Clostridia (oranges)
  "Clostridium_innocuum_CrfA"             = "#E67E22",
  "Clostridium_symbiosum_CD11"            = "#D35400",
  
  # Actinobacteria (purples/reds)
  "Collinsella_aerofaciens_CD3"           = "#8E44AD",
  "Enterocloster_bolteae_A14CFS"          = "#C0392B",
  "Erysipelatoclostridium_ramosum_CD4"    = "#A04000",
  
  # Others
  "Parabacteroides_distasonis_A04A"       = "#27AE60", 
  "Ruminococcus_gnavus_A2"                = "#17A589",
  "BruPT-FE"                               = "#E91E63",
  
  # Media/controls (neutrals)
  "RCM"                                    = "#34495E",
  "cCOM"                                   = "#7F8C8D",
  "BHIS"                                   = "#95A5A6")

df_173005 <- merged_taurine_condition %>%
  dplyr::filter(Feature == 173005) %>%
  dplyr::mutate(intensity = ifelse(is.na(intensity), 0, intensity)) %>%
  dplyr::filter(Time_Point == "T72") |> 
  dplyr::filter(!str_detect(condition, "luhongzhouii"))

#write_csv(df_173005, "CA_5ASA_crohn_cCOM_culturing.csv")
getwd()

controls    <- c("BHIS", "BruPT-FE", "RCM")
bacteroides <- sort(unique(df_173005$condition[grepl("^Bacteroides_", df_173005$condition)]))
ccom        <- "cCOM"
other       <- sort(setdiff(unique(df_173005$condition),
                            c(controls, ccom, bacteroides)))

desired_order <- c(controls, ccom, bacteroides, other)

df_173005 <- df_173005 %>%
  dplyr::mutate(
    condition = factor(condition, levels = desired_order),
    intensity_log10 = log10(intensity + 1),
    presence = factor(if_else(intensity > 0, "yes", "no"), levels = c("yes", "no")))

df_stats <- df_173005 %>%
  group_by(condition) %>%
  summarize(
    mean = mean(intensity_log10, na.rm = TRUE),
    sd   = sd(intensity_log10,   na.rm = TRUE),
    n    = n(),
    .groups = "drop")

df_stats
ref_cond <- "BHIS"

# Create presence/absence
df_173005 <- df_173005 %>%
  dplyr::mutate(presence_raw = factor(if_else(intensity > 0, "yes", "no"), levels = c("yes", "no")))

conditions  <- levels(df_173005$condition)
other_conds <- setdiff(conditions, ref_cond)

n_cond <- length(levels(df_173005$condition))
custom_colors        <- hue_pal(l = 70, c = 60)(n_cond) 
custom_colors_darker <- hue_pal(l = 40, c = 80)(n_cond)  

y_min <- -5e4
p_173005 <- ggplot(df_173005,
                   aes(x = condition, y = intensity, fill = condition)) +
  geom_boxplot(color = "black", size = 0.5, alpha = 0.5,
               width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(color = condition),
              shape = 19,
              position = position_jitter(0.15),
              size = 5,
              alpha = 1,
              stroke = NA) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors_darker) +
  scale_y_continuous(
    limits = c(y_min, 1.6e6), 
    breaks = c(0, 2.5e5, 5e5, 1e6, 1.5e6),
    labels = scales::label_scientific(digits = 2)) +
  labs(
    x = "Condition",
    y = "Raw peak area") +
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
    legend.position = "none")

p_173005

#ggsave("CA_5ASA_crohn_cCOM_T72.pdf", p_173005, width = 5, height = 7, dpi = 900)
#getwd()




