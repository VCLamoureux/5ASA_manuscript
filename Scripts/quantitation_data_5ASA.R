# set the working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/")

# load libraries
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggbreak)
library(duckplyr)

# import table 
quant_table_Lreuteri <- read_csv("Quantification_L_reuteri_CA_5ASA.csv")
df_plot_Lreu <- quant_table_Lreuteri %>%
  dplyr::mutate(condition = str_remove(sample, "_\\d+$")) |> 
  dplyr::mutate(presence = if_else(concentration_um > 0, "yes", "no"))

desired_order <- c("SterileControlSeries1_LDM4_Gluc_5ASA", "SterileControlSeries2_LDM4_Gluc_5ASA_CA", "MicrobialControlSeries1_LDM4_Lreuteri6475_Gluc_5ASA", "LreuteriTestSeries_LDM4_Lreuteri6475_Gluc_5ASA_CA")

df_plot_Lreu_2 <- df_plot_Lreu %>%
  dplyr::mutate(condition = factor(condition, levels = desired_order))
#write_csv(df_plot_Lreu_2, "L_reuteri_5ASA_CA_quant_processed.csv")

df_plot_Lreu_3 <- ggplot(df_plot_Lreu_2, aes(x = condition, y = concentration_um, fill = condition)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = condition), width = 0.2, size = 4, alpha = 0.8, shape = 16) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  scale_color_brewer(palette = "Set2", guide = FALSE) +
  scale_y_continuous() +
  labs(x = NULL, y = "Cholyl-5-ASA (µM)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y     = element_text(size = 14),
    axis.title      = element_text(size = 18, face = "bold"),
    legend.position = "none")

df_plot_Lreu_3
#ggsave("Quant_L_reuteri.pdf", df_plot_3, width = 2, height = 8, dpi = 900)


# import quantification table and format
quant_table <- read_csv("Quantification_IBD_cohort_human.csv")
quant_table_1 <- quant_table |> 
  dplyr::mutate(concentration_3 = if_else(concentration_3 == 0, 0, concentration_3), status = factor(status, levels = c("age matched control", "5-ASA treated")))

cv_conc4 <- quant_table_1 %>%
  dplyr::filter(status == "5-ASA treated") %>%
  summarize(
    mean_2 = mean(concentration_3, na.rm = TRUE),
    sd_2   = sd(concentration_3,   na.rm = TRUE),
    cv2    = sd_2 / mean_2 * 100)
cv_conc4

eps <- 1e-3

qt1_log <- quant_table_1 %>% dplyr::mutate(conc_floor = if_else(concentration_3 <= 0, eps, concentration_3))
#write_csv(qt1_log, "IBD_5ASA_CA_cohort4.csv")

plot_UC_conc2 <- ggplot(qt1_log, aes(x = status, y = conc_floor)) +
  geom_boxplot(aes(fill = status), alpha = 0.5, outlier.shape = NA, width = 0.6) +
  geom_point(aes(color = status), position = position_jitter(width = 0.2, height = 0), size = 5, alpha = 0.85, stroke = NA) +
  scale_fill_discrete(guide = FALSE) +
  scale_color_discrete(guide = FALSE) +
  labs(x = NULL, y = "Cholyl-5-ASA (µM, log10 scale)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 18, face = "bold"),
    legend.position = "none") +
  scale_y_log10(
    breaks = c(eps, 0.01, 0.1, 1, 10, 50),
    labels = c("0", "0.01", "0.1", "1", "10", "50"))
plot_UC_conc2
#ggsave("quant_plot.pdf", plot_UC_conc2, width = 3, height = 7, dpi = 900)


df2 <- quant_table_1 %>% dplyr::mutate(presence = if_else(concentration_3 > 0, "yes", "no"))
diag_lvls <- levels(factor(df2$status))
pairs     <- combn(diag_lvls, 2, simplify = FALSE)

# Fisher function
fisher_res <- map_dfr(pairs, function(lvls) {
  sub <- df2 %>% filter(status %in% lvls)
  tbl <- table(
    factor(sub$status, levels = lvls),
    factor(sub$presence, levels = c("yes","no")))
  ft  <- fisher.test(tbl)
  tibble(
    group1 = lvls[1],
    group2 = lvls[2],
    p      = ft$p.value
  )
})
fisher_res

quant_table_monoculture_syncom <- read_csv("dysbiotic_community_quantitation - Sheet1.csv")
df_bile <- quant_table_monoculture_syncom %>%
  dplyr::filter(!is.na(um_concentration)) %>%
  dplyr::filter(!str_detect(Bacteria, regex("healthy|dysbiotic|blank", ignore_case = TRUE))) %>% 
  dplyr::filter(!str_detect(Sample, regex("24 hrs", ignore_case = TRUE))) %>%
  dplyr::mutate(Bacteria = factor(Bacteria))
#write_csv(df_bile, "dysbiotic_quant.csv")

# Define the control group
control_name <- "control"

# Extract all bacteria except the control
other_bacteria <- setdiff(levels(df_bile$Bacteria), control_name)
df_bile$Bacteria <- factor(df_bile$Bacteria, levels = c(control_name, other_bacteria))  
n_cond <- length(levels(df_bile$Bacteria))

custom_colors        <- scales::hue_pal(l = 70, c = 60)(n_cond)
custom_colors_darker <- scales::hue_pal(l = 40, c = 80)(n_cond)
names(custom_colors)        <- levels(df_bile$Bacteria)
names(custom_colors_darker) <- levels(df_bile$Bacteria)

p_monoculture <- ggplot(df_bile,
  aes(x = Bacteria, y = um_concentration, fill = Bacteria, color = Bacteria)) +
  geom_boxplot(
    color        = "black",
    fill         = NA,
    size         = 0.5,
    alpha        = 0.6,
    width        = 0.55,
    outlier.shape = NA) +
  geom_jitter(
    shape    = 19,
    size     = 4,
    alpha    = 1,
    stroke   = NA,
    position = position_jitter(width = 0.2)) + scale_fill_manual(values  = custom_colors) + scale_color_manual(values = custom_colors_darker) + labs(
    x = "Bacterial strain",
    y = "Cholyl-5-ASA (µM)") +
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

p_monoculture

#ggsave("Thomas_dysbiotic_community.pdf", p_monoculture, width = 4, height = 5, dpi = 900)
getwd()

# import quantification table for microbial monoculture / synthetic communities
Defined_complex <- quant_table_monoculture_syncom %>%
  dplyr::filter(str_detect(Bacteria, regex("control|healthy|dysbiotic", ignore_case = TRUE))) |> 
  dplyr::mutate(Bacteria = str_remove(Bacteria, regex("\\s+replicate.*$", ignore_case = TRUE)))

pal5 <- c(
  "healthy stool 1"  = "#1B9E77",
  "healthy stool 2"  = "#D95F02",
  "healthy stool 3"  = "#7570B3",
  "Healthy SynCom"   = "#E7298A",
  "Dysbiotic SynCom" = "#66A61E", 
  "control" = "grey")

df_plot <- Defined_complex %>%
  dplyr::mutate(point_color = pal5[Bacteria]) |> 
  dplyr::mutate(Bacteria = factor(Bacteria, levels = c("control","healthy stool 1","healthy stool 2","healthy stool 3","Healthy SynCom", "Dysbiotic SynCom")))
#write_csv(df_plot, "human_donors_consortia.csv")

cv_result <- df_plot %>% 
  dplyr::filter(Bacteria %in% c("healthy stool 1","healthy stool 2","healthy stool 3")) %>%
  group_by(Bacteria) %>%
  summarize(mean_conc = mean(um_concentration, na.rm = TRUE)) %>%
  summarize(cv_percent = sd(mean_conc) / mean(mean_conc) * 100)

df_plot_1 <- ggplot(df_plot, aes(x = Bacteria, y = um_concentration)) +
  geom_boxplot(aes(fill = Bacteria), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Bacteria), width = 0.2, size = 6, alpha = 1, shape = 16) +
  scale_fill_manual(values = pal5) +
  scale_color_manual(values = pal5) +
  scale_y_continuous(limits = c(0, 0.40)) +
  labs(x = NULL, y = "Concentration (µM)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y     = element_text(size = 14),
    axis.title      = element_text(size = 18, face = "bold"),
    legend.position = "none")
df_plot_1
#ggsave("Community_quant_monoculture_donor.pdf", df_plot_1, width = 4, height = 6, dpi = 900)
getwd()

# Read the Enterococcus strain level quantification table
strain_level <- read_csv("20250916_PDorrestein_5-ASA-CA_EnterococcusStudy_Summary_1.csv", show_col_types = FALSE)

strain_data <- strain_level %>%
  dplyr::filter(!is.na(condition)) %>%
  dplyr::mutate(
    Bacteria = str_remove(condition, "^_"),
    Bacteria = str_replace(Bacteria, "_T\\d+$", "")) %>%
  dplyr::mutate(Bacteria = case_when(
      str_detect(condition, "^_BlkBlk")   ~ "BlkBlk",
      str_detect(condition, "^_sterile")  ~ "sterile", TRUE ~ Bacteria))

ordered_levels <- c("BlkBlk","sterile", sort(setdiff(unique(strain_data$Bacteria), c("BlkBlk", "sterile"))))
strain_data_1 <- strain_data %>%
  dplyr::mutate(Bacteria = factor(Bacteria, levels = ordered_levels))
#write_csv(strain_data, "Processed_5ASA_CA_Enterococcus_strain_data.csv")

pal <- colorRampPalette(brewer.pal(8, "Set2"))(length(ordered_levels))
names(pal) <- ordered_levels

# plot the monocultures
df_plot_strains <- ggplot(strain_data_1, aes(x = Bacteria, y = um_concentration)) +
  geom_boxplot(aes(fill = Bacteria), alpha = 1, outlier.shape = NA) +
  geom_jitter(aes(color = Bacteria), width = 0.2, size = 5, alpha = 1, shape = 16, stroke = NA) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(x = NULL, y = "Concentration (µM)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y     = element_text(size = 14),
    axis.title      = element_text(size = 18, face = "bold"),
    legend.position = "none")

df_plot_strains
#ggsave("strain_quant.pdf", df_plot_strains, width = 5, height = 7, dpi = 900)
getwd()

# read in the quant data for 5-ASA
ASA_quant <- readr::read_csv("Quantitation_5ASA_10P.csv", show_col_types = FALSE, na = c("N/A", "NA", ""))
df_meta <- ASA_quant %>%
  dplyr::transmute(sample_name = stringr::str_replace(sample_name, "_DF20_", ""), status = factor(status, levels = c("age matched control", "5-ASA treated")), ASA_concentration_um_3)

# check for cv
cv_conc4 <- df_meta %>%
  dplyr::filter(status == "5-ASA treated") %>%
  summarize(
    mean_2 = mean(ASA_concentration_um_3, na.rm = TRUE),
    sd_2   = sd(ASA_concentration_um_3,   na.rm = TRUE),
    cv2    = sd_2 / mean_2 * 100)
cv_conc4

eps <- 1e-3
qt1_log <- df_meta %>%
  dplyr::mutate(conc_floor = dplyr::if_else(is.na(ASA_concentration_um_3) | ASA_concentration_um_3 <= 0, eps, ASA_concentration_um_3))

# plot 5-ASA concentrations in a prospective cohort
plot_UC_conc2 <- ggplot(qt1_log, aes(x = status, y = conc_floor)) +
  geom_boxplot(aes(fill = status), alpha = 0.5, outlier.shape = NA, width = 0.6) +
  geom_point(aes(color = status), position = position_jitter(width = 0.2, height = 0), size = 5, alpha = 0.85, stroke = NA) +
  scale_fill_discrete(guide = FALSE) +
  scale_color_discrete(guide = FALSE) +
  labs(x = NULL, y = "5-ASA (µM, log10 scale)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 18, face = "bold"),
    legend.position = "none") +
  scale_y_log10()

plot_UC_conc2
#ggsave("quant_plot_5ASA.pdf", plot_UC_conc2, width = 3, height = 7, dpi = 900)


