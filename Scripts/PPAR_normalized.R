# set the working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/")

library(tidyverse)
library(dplyr)
library(rstatix)
library(ggplot2)
library(car)

ppar_df <- readr::read_csv("Data_ppar_y.csv", show_col_types = FALSE)
ppar_df_clean <- ppar_df |>
  dplyr::mutate(group = dplyr::case_when(
      str_detect(condition, "^control")   ~ "Control",
      str_detect(condition, "^5ASA")      ~ "5-ASA",
      str_detect(condition, "^CA_5ASA")   ~ "Cholyl–5-ASA", TRUE ~ NA_character_)) |>
  dplyr::filter(!is.na(group))

# Order: Control, 5-ASA, Cholyl–5-ASA
ppar_df_clean <- ppar_df_clean |>
  dplyr::mutate(group = factor(group, levels = c("Control", "5-ASA", "Cholyl–5-ASA"))) |>
  dplyr::filter(is.na(concentration_mM) | concentration_mM != 0.5)

mean_control <- ppar_df_clean |>
  dplyr::filter(group == "Control") |>
  summarise(m = mean(rlu, na.rm = TRUE)) |>
  pull(m)

ppar_df_norm <- ppar_df_clean |>
  dplyr::mutate(rlu_norm = rlu / mean_control)
#write_csv(ppar_df_norm, "PPAR_normalized_woRosiglitazone_5mM.csv")

ppar_summary <- ppar_df_norm |>
  group_by(group) |>
  summarise(mean_norm = mean(rlu_norm, na.rm = TRUE), sem_norm  = sd(rlu_norm, na.rm = TRUE) / sqrt(n()), .groups   = "drop")

# Define custom colors for groups
ppar_colors <- c("Control" = "black", "5-ASA" = "#ffb1ae", "Cholyl–5-ASA" = "#ffe6e5")
plot_ppar_adj <- ggplot(ppar_summary, aes(x = group, y = mean_norm, fill = group)) +
  geom_col(color = "black", width = 0.6, alpha = 0.95) +
  geom_errorbar(aes(ymin = mean_norm - sem_norm, ymax = mean_norm + sem_norm), width = 0.2, size = 0.9) +
  geom_jitter(data  = ppar_df_norm, aes(x = group, y = rlu_norm, color = group), inherit.aes = FALSE, width = 0.12, size = 6, alpha = 0.9, stroke = NA) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 6)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = ppar_colors, guide = "none") +
  scale_color_manual(values = ppar_colors, guide = "none") +
  labs(x = NULL, y = "PPAR-γ activity\nControl-normalized (RLU, Control = 1)") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title  = element_text(size = 18, face = "bold"))

plot_ppar_adj
#ggsave("PPAR_05mM.pdf", plot = plot_ppar_adj, width = 3, height = 7, dpi = 900)

#Subset to main groups (no Rosiglitazone if you want just the 3)
ppar_main <- ppar_df_norm |>
  dplyr::filter(group %in% c("Control", "5-ASA", "Cholyl–5-ASA")) |>
  droplevels()

# Shapiro test for normality
shapiro_results <- ppar_main |>
  group_by(group) |> 
  shapiro_test(rlu_norm)
shapiro_results

# Levene test for homogeneity of variances
levene_results <- leveneTest(rlu_norm ~ group, data = ppar_main)
levene_results

# if assumptions are met, then proceed with ANOVA
anova_main <- aov(rlu_norm ~ group, data = ppar_main)
summary(anova_main)
tukey_results <- TukeyHSD(anova_main)
tukey_results


