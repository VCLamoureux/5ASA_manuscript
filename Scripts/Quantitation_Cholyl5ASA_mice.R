# set the working directory
setwd("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/5ASA_drug_project/")

# load libraries
library(tidyverse)
library(duckplyr)
library(ggbreak)
library(scales)

# import table 
quant_table_ileum_flush <- read_csv("ileum_flush_CA5ASA.csv")
quant_table_ileum_flush_formatted <- quant_table_ileum_flush %>%
  dplyr::mutate(across(where(is.character), ~ ifelse(. == "N/A", "0", .))) %>%
  dplyr::filter(!str_detect(`Sample Name`, "Blk")) %>%
  dplyr::mutate(Condition = case_when(
      str_detect(`Sample Name`, regex("PBS", ignore_case = TRUE)) ~ "Control",
      str_detect(`Sample Name`, regex("5ASA-CA|CA5ASA|Cholyl-5ASA", ignore_case = TRUE)) ~ "Cholyl-5ASA",
      str_detect(`Sample Name`, regex("5ASA", ignore_case = TRUE)) ~ "5ASA", TRUE ~ "Unknown"),
    `Calculated Concentration` = as.numeric(`Calculated Concentration`), MW = case_when(
      Condition == "5ASA"        ~ 153.135,
      Condition == "Cholyl-5ASA" ~ 543.7010, TRUE ~ NA_real_), conc_uM = case_when(!is.na(MW) ~ `Calculated Concentration` / MW, TRUE ~ 0)) |> 
  dplyr::mutate(body_part = str_extract(`Sample Name`, "(Flush|Tissue)_[A-Za-z]+"))

desired_order <- c("Tissue: Liver", "Flush: Ileum", "Tissue: Ileum", "Flush: Colon", "Tissue: Colon")

# to plot cholyl-5-ASA in the 5ASA group, uncomment the next line and comment the following one
df_cholyl <- quant_table_ileum_flush_formatted %>%
  #dplyr::filter(Condition == "5ASA", !is.na(body_part))
  dplyr::filter(Condition == "Cholyl-5ASA", !is.na(body_part))

df_cholyl_1 <- df_cholyl %>%
  dplyr::mutate(body_part_clean = body_part %>%
      str_replace("Tissue_", "Tissue: ") %>%
      str_replace("Flush_",  "Flush: ") %>%
      factor(levels = desired_order))
#write_csv(df_cholyl_1, "Cholyl5ASA_mice_conc_uM.csv")

plot_cholyl <- ggplot(
  df_cholyl_1,
  aes(x = body_part_clean, y = conc_uM)) +
  geom_boxplot(
    fill         = "white",
    colour       = "black",
    alpha        = 0.8,
    outlier.shape = NA,
    width        = 0.6) +
  geom_jitter(
    colour = "#27AAE1",
    width  = 0.15,
    size   = 4,
    alpha  = 0.9,
    stroke = NA) +
  scale_y_continuous(
    trans  = "log10",
    labels = scales::label_number(accuracy = 0.01)) +
  labs(x = "Disease Location", y = "Cholyl–5-ASA\nLog (Peak area)") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 18, face = "bold"),
    plot.title  = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none")

plot_cholyl
#ggsave("Cholyl5ASA_levels_in_CA5ASA_treated_mice_ileum_flush_vq.pdf", plot = plot_cholyl, width = 3, height = 5, dpi = 900)













