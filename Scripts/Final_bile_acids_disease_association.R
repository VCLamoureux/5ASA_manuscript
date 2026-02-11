# Set your working directory
setwd("/Users/vincentlamoureux/OneDrive - University of California, San Diego Health/")

# These packages are part of CRAN repository
#install.packages("data.table", dependencies = TRUE)
#install.packages("tidyverse", dependencies = TRUE)
#install.packages("pheatmap", dependencies = TRUE)

## Load the packages required for the analysis
library(data.table)
library(tidyverse)
library(pheatmap)
library(circlize)
library(ComplexUpset)
library(tibble)

# Specify the folder path - it should be the folder inside the working directory
folder_path <- "Postdoc_UCSD/Postdoc_projects/Main_project/Nature_protocols/Bile_acids_class/"

## Download/Import the ReDU metadata file - it should be in the working directory folder and NOT be in the sub-folder with the csv files from the Fast Search

# Define the filename for the ReDU metadata
processed_redu_metadata <- "all_sampleinformation.tsv"
options(timeout = max(600L, getOption("timeout")))

# Check if the pre-processed metadata file exists in the working directory
if (!file.exists(file.path(getwd(), processed_redu_metadata))) {
  redu_url <- "https://redu.gnps2.org/dump"
  download.file(redu_url, file.path(getwd(), processed_redu_metadata), mode = "wb", quiet = TRUE)
}

redu_metadata <- data.table::fread(file.path(getwd(), processed_redu_metadata))
# Path to the single FASST results file
fasst_file <- file.path(folder_path, "FASST_all_Candidate.csv")

# Read the file and create a Compound column from NAME
molecules_interest <- readr::read_csv(fasst_file, show_col_types = FALSE) %>%
  dplyr::mutate(Compound = NAME)

# Filter by cosine and matching peaks
molecules_interest_filtered <- molecules_interest %>%
  dplyr::filter(Cosine >= 0.7, `Matching Peaks` >= 4)

# Prepare the data tables for merging
## Create a function to extract the desired segment from the USI column
USI_key <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":", fixed = TRUE))
  paste(parts[2], parts[length(parts)], sep = ":")
}

# Apply the function to each row of the USI column in the molecules_interest
molecules_interest_filtered$USI <- vapply(molecules_interest_filtered$USI, USI_key, FUN.VALUE = character(1))

# Prepare the ReDU metadata USI column for merging with FASST output table
redu_metadata$USI <- vapply(redu_metadata$USI, USI_key, FUN.VALUE = character(1))

# Merge the ReDU metadata table and the FASST MASST output table
ReDU_MASST <- dplyr::left_join(molecules_interest_filtered, redu_metadata, by = "USI")
unique_accessions <- ReDU_MASST %>%
  #dplyr::filter(!is.na(ATTRIBUTE_DatasetAccession)) %>%
  dplyr::distinct(ATTRIBUTE_DatasetAccession)

dataset_prefix_list <- unique_accessions %>%
  dplyr::mutate(prefix = dplyr::case_when(
      str_starts(ATTRIBUTE_DatasetAccession, "ST")   ~ "ST",
      str_starts(ATTRIBUTE_DatasetAccession, "MSV")  ~ "MSV",
      str_starts(ATTRIBUTE_DatasetAccession, "MTBL") ~ "MTBL", TRUE ~ "Other")) %>%
  dplyr::group_by(prefix) %>%
  dplyr::summarise(n_unique = n(), datasets = toString(ATTRIBUTE_DatasetAccession), .groups = "drop")

dataset_prefix_list

# Standardize the body parts and Health Status
ReDU_MASST_standardize <- ReDU_MASST |>
  dplyr::mutate(UBERONBodyPartName = str_replace_all(UBERONBodyPartName,
      "skin of trunk|skin of pes|head or neck skin|axilla skin|skin of manus|arm skin|skin of leg","skin"),
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, "blood plasma|blood serum", "blood"),
    HealthStatus = str_replace(HealthStatus, "Chronic Illness", "chronic illness"),
    HealthStatus = str_replace(HealthStatus, "Healthy", "healthy"))

# Separate humans and rodents from the merged data table
df_humans <- ReDU_MASST_standardize |>
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

df_humans_unique <- df_humans |>
  dplyr::distinct(filename, .keep_all = TRUE)

analyze_counts <- function(df, column_interest) {
  df_body_parts <- df |> dplyr::distinct(dplyr::across(dplyr::all_of(column_interest)))
  
  df_counts <- df |>
    dplyr::count(dplyr::across(dplyr::all_of(column_interest)), name = "Counts_fastMASST")
  
  compounds <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(column_interest))) |>
    dplyr::summarise(
      Compounds = dplyr::n_distinct(Compound),
      CompoundsList = toString(unique(Compound)),
      .groups = "drop"
    )
  
  df_body_parts |>
    dplyr::left_join(df_counts, by = column_interest) |>
    dplyr::left_join(compounds, by = column_interest)
}

# Get a glimpse of the number of counts per organ
body_counts_humans <- analyze_counts(df_humans, "UBERONBodyPartName")
head(body_counts_humans)

# Create a function to pivot the table for data visualization
# Create a data frame for humans, retaining both columns:
df_humans <- ReDU_MASST_standardize |>
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

df_humans_remv_missing <- df_humans |>
  dplyr::mutate(DOIDCommonName = dplyr::case_when(DOIDCommonName == "missing value" & HealthStatus != "missing value" ~ HealthStatus, TRUE ~ DOIDCommonName)) |>
  dplyr::filter(DOIDCommonName != "missing value")

# Pivot table function remains the same
prepare_pivot_table <- function(df, column_interest, compound_col) {
  grouped_df <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(compound_col, column_interest)))) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop")
  grouped_df %>%
    tidyr::pivot_wider(
      names_from = dplyr::all_of(compound_col),
      values_from = Count,
      values_fill = list(Count = 0))
}

# Define the variables based on your research question
## Here we are interesting in organ distribution in humans and rodents of the molecule of interest
variable <- "UBERONBodyPartName"
pivot_table_humans <- prepare_pivot_table(df_humans, variable, "Compound")

# Prepare the table to be compatible with pheatmap package
humans_molecules_counts_by_bodypart <- pivot_table_humans |>
  dplyr::arrange(.data[[variable]]) |>
  tibble::column_to_rownames(variable)

# Convert all columns to numeric for the humans df
humans_molecules_counts_by_bodypart <- humans_molecules_counts_by_bodypart |>
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

# Define your chosen colors
colors_version <- c("#FFFFFF", "#C7D6F0", "#EBB0A6")
# Creating the gradient function
color_gradient <- colorRampPalette(colors_version)
# Generate 30 discrete colors from this gradient
gradient_colors <- color_gradient(30)

# The users can log scale or not the data
log_humans_molecules_counts_by_bodypart <- log2(1 + humans_molecules_counts_by_bodypart)

# Organ distribution in humans
## Use heatmap for data visualization or organ distribution - humans
### If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE
Organ_humans <- pheatmap(
  log_humans_molecules_counts_by_bodypart,
  color = gradient_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  angle_col = 90,
  main = "Organ distribution in humans",
  fontsize = 10,
  cellwidth = 10,
  cellheight = 10,
  treeheight_row = 100,
  fontsize_row = 12,
  fontsize_col = 1,
  legend_fontsize = 10,
  border_color = NA)
Organ_humans
#ggsave("Organ_distribution_in_humans.pdf", plot = Organ_humans, width = 10, height = 10, dpi = 900)
#getwd()

# Organ distribution in rodents
## Use heatmap for data visualization or organ distribution - rodents
### If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE

#ggsave("Organ_distribution_in_rodents.pdf", plot = Organ_rodents, width = 10, height = 10, dpi = 90#getwd()

# Health phenotype association
## Filter for human information in the ReDU metadata
df_redu_humans <- redu_metadata |> dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

# Humans - filter for lifestage, DOIDCommonName, HealthStatus, BiologicalSex
## For LifeStage
human_ReDU_LifeStage <- df_redu_humans |>
  dplyr::count(LifeStage) |>
  dplyr::rename(LifeStage_counts = n, LifeStage = LifeStage)
human_ReDU_LifeStage$LifeStage_counts <- as.numeric(human_ReDU_LifeStage$LifeStage_counts)
# For DOIDCommonName
human_ReDU_DOIDCommonName <- df_redu_humans |>
  dplyr::count(DOIDCommonName) |>
  dplyr::rename(DOIDCommonName_counts = n, DOIDCommonName = DOIDCommonName)
human_ReDU_DOIDCommonName$DOIDCommonName_counts <- as.numeric(human_ReDU_DOIDCommonName$DOIDCommonName_counts)
# For HealthStatus
human_ReDU_HealthStatus <- df_redu_humans |>
  dplyr::count(HealthStatus) |>
  dplyr::rename(HealthStatus_counts = n, HealthStatus = HealthStatus)
human_ReDU_HealthStatus$HealthStatus_counts <- as.numeric(human_ReDU_HealthStatus$HealthStatus_counts)
# For BiologicalSex
human_ReDU_BiologicalSex <- df_redu_humans |>
  dplyr::count(BiologicalSex) |>
  dplyr::rename(BiologicalSex_counts = n, BiologicalSex = BiologicalSex)
human_ReDU_BiologicalSex$BiologicalSex_counts <- as.numeric(human_ReDU_BiologicalSex$BiologicalSex_counts)

# Normalization of the FASST search output
## Grouping and counting occurrences
### (Optional) remove the # symbol to activate the line and apply a minimum of 3 counts to be include before normalization
#### Here we are interested in finding disease association
grouped_df_humans <- df_humans |>
  dplyr::group_by(Compound, DOIDCommonName) |>
  dplyr::summarise(Count = dplyr::n(), .groups = "drop")
#dplyr::filter(Count >= 3)

grouped_df_humans_pivot_table <- grouped_df_humans |>
  tidyr::pivot_wider(names_from = Compound, values_from = Count, values_fill = list(Count = 0))

# Merging
merged_DOID_humans <- dplyr::left_join(grouped_df_humans_pivot_table, human_ReDU_DOIDCommonName, by = "DOIDCommonName")

# new from 20nov2024 because if all values are 0, the column is removed
merged_DOID_humans <- merged_DOID_humans |>
  dplyr::filter(DOIDCommonName != "missing value") |>
  dplyr::select(DOIDCommonName, dplyr::where(~ !all(. == 0))) |>
  dplyr::filter(rowSums(dplyr::across(dplyr::where(is.numeric))) != 0)

merged_DOID_humans_transpose <- merged_DOID_humans |>
  tibble::column_to_rownames("DOIDCommonName") |>
  dplyr::select(-DOIDCommonName_counts) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("Bile_acids") |>
  dplyr::filter(!str_detect(Bile_acids, "-18\\.00|-18\\.01|-36\\.02|36\\.02|36\\.01|54\\.02|54\\.03|88\\.99|72\\.04")) |>
  dplyr::filter(Bile_acids != "NA") |>
  tibble::column_to_rownames("Bile_acids")

df <- merged_DOID_humans_transpose

# Convert to binary (presence/absence) format
df_binary <- df %>%
  tibble::rownames_to_column("Bile_acids") %>%
  dplyr::mutate(dplyr::across(-Bile_acids, ~ dplyr::if_else(.x > 0, 1L, 0L))) %>%
  dplyr::filter(rowSums(dplyr::across(-Bile_acids)) >= 2)

upset_plot <- ComplexUpset::upset(df_binary, intersect = names(df_binary)[-1], name = "Bile Acids Intersections", width_ratio = 0.2, min_size = 2) +
  labs(title = "UpSet Plot of Bile Acids Across Diseases", x = "Diseases", y = "Intersection Size") +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
upset_plot
#ggsave("UpSet_plot.pdf", upset_plot, width = 22, height = 13, dpi = 900)

compounds_kept <- df_binary$Bile_acids
n_spectra_df_binary <- df_humans %>%
  dplyr::filter(Compound %in% compounds_kept)
n_spectra_df_binary
diseases_kept <- setdiff(colnames(df_binary), "Bile_acids")
n_spectra_df_binary_disease <- df_humans %>%
  dplyr::filter(Compound %in% df_binary$Bile_acids, DOIDCommonName %in% diseases_kept) #|> nrow()
n_spectra_df_binary_disease

#write_csv(n_spectra_df_binary_disease, "Postdoc_UCSD/Postdoc_projects/Main_project/Nature_protocols/Bile_acids_class/Bile_acids_disease_association_filtered_counts.csv")
# A tibble: 23,656 × 52

diseases_target <- c(
  #"inflammatory bowel disease",
  #"rheumatoid arthritis"
  #"ulcerative colitis",
  #"Crohn's disease",
  "diabetes mellitus",
  "obesity")

# Identify all diseases in df_binary
all_diseases <- names(df_binary)[-1]
# All other diseases
other_diseases <- setdiff(all_diseases, diseases_target)
# Bile acids EXCLUSIVELY found in the four diseases
exclusive_4 <- df_binary %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(diseases_target), ~ .x == 1L)) %>%
  dplyr::filter(dplyr::if_all(dplyr::all_of(other_diseases), ~ .x == 0L))
# Extract names
exclusive_bile_acids <- exclusive_4$Bile_acids
exclusive_bile_acids

# Count
length(exclusive_bile_acids)
exclusive_bile_acids <- exclusive_4$Bile_acids
df_humans_filtered_6 <- df_humans_remv_missing %>%
  dplyr::filter(Compound %in% exclusive_bile_acids)

# Save as high-resolution figure
#ggsave("dot_plot_fixed_color_not_exact_y.pdf", plot, width = 12, height = 6, dpi = 300)

if (exists("plot")) print(plot)

# Remove the "Bile_acids" column and convert the remaining data to a matrix.
binary_matrix <- as.matrix(df_binary[, -1])

# Ensure the matrix has non-NA column names.
if (is.null(colnames(binary_matrix)) || any(is.na(colnames(binary_matrix)))) {
  colnames(binary_matrix) <- paste0("Disease_", seq_len(ncol(binary_matrix)))
}

# Remove unwanted diseases from the binary matrix.
unwanted <- c("Chagas disease", "dental caries")
binary_matrix_filtered <- binary_matrix[, !(colnames(binary_matrix) %in% unwanted), drop = FALSE]

# Compute the co-occurrence matrix from the filtered binary matrix.
co_occurrence_filtered <- t(binary_matrix_filtered) %*% binary_matrix_filtered

# Assign proper row and column names.
rownames(co_occurrence_filtered) <- colnames(binary_matrix_filtered)
colnames(co_occurrence_filtered) <- colnames(binary_matrix_filtered)

# Remove self-matches and duplicate links:
diag(co_occurrence_filtered) <- 0
co_occurrence_filtered[upper.tri(co_occurrence_filtered)] <- 0

# Create the chord (circos) plot.
#pdf("Circos_plot_no_duplicates_filtered_5.pdf", width = 10, height = 10, onefile = TRUE)
#circlize::chordDiagram(co_occurrence_filtered, transparency = 0.6, annotationTrack = c("name", "grid"))
#dev.off()

#getwd()
