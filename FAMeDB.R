# This script was adapted from HADEG.R available at https://github.com/jarojasva/HADEG


# Set the path to the directory with the proteinortho file and the FAMeDB pathways table
# setwd("Path_to_your_directory_files")

library(stringr)
library(dplyr)
library(readr)
library(scales)
library(svglite)
library(ggplot2)
library(reshape2)

# Read the FAMeDB_pathways.csv file
all_pathways <- read.csv("FAMeDB_pathways_2026_01_19.csv", stringsAsFactors = FALSE)

# Ensure that there are no leading or trailing whitespaces in the Protein_ID column
all_pathways$Protein_ID <- trimws(all_pathways$Protein_ID)

# Read the orthology results file with the first row as header
df_ortho <- read.table("Results_FAMeDB.proteinortho.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "", fill = TRUE)

# Rename the first three columns as ColumnA, ColumnB, and ColumnC
colnames(df_ortho)[1:3] <- c("ColumnA", "ColumnB", "ColumnC")

# Sort all columns of df_ortho alphabetically based on their names.
df_ortho <- df_ortho[, sort(colnames(df_ortho))]

# Adjust the header to have the 'FAMeDB_2026_01_19.faa' column at the beginning
df_ortho <- df_ortho[, c(grep("FAMeDB_2026_01_19.faa", colnames(df_ortho)),
                         setdiff(seq_along(df_ortho), grep("FAMeDB_2026_01_19.faa", colnames(df_ortho))))]

# Remove specific columns
columns_to_remove <- c("ColumnA", "ColumnB", "ColumnC")
df_selected <- df_ortho[, !(colnames(df_ortho) %in% columns_to_remove)]

## Remove lines containing "*" in the "FAMeDB_2026_01_19.faa" column
df_selected <- df_selected[!grepl("\\*", df_selected$FAMeDB_2026_01_19.faa), ]

## Save the selected data in a new file
write.table(df_selected, file = "0_table_FAMeDB_raw_hits.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Split the 'FAMeDB_2026_01_19.faa' column into two new columns
split_columns <- str_split_fixed(df_selected$FAMeDB_2026_01_19.faa, "\\|", 2)

# Assign names to the new columns
colnames(split_columns) <- c("Protein_ID", "Protein_desc")

# Combine the new columns with the existing dataframe
df_selected <- cbind(df_selected, split_columns)

# Remove the original column
df_selected <- df_selected[, !grepl("FAMeDB_2026_01_19.faa", colnames(df_selected))]

# Remove "Protein_desc" column
df_selected <- df_selected[, !grepl("Protein_desc", colnames(df_selected))]

# Ensure that there are no leading or trailing whitespaces in the Protein_ID column
df_selected$Protein_ID <- trimws(df_selected$Protein_ID)

# Select Protein_ID, Gene and Pathway columns from all_pathways
pathways_subset <- all_pathways %>%
  select(Protein_ID, Gene, Pathway)

# Inner join of dataframes by "Protein_ID"
df_selected <- df_selected %>%
  inner_join(pathways_subset, by = "Protein_ID")

# Move the columns "Protein_ID", "Gene" and "Pathway" as the first three columns
df_selected <- df_selected[, c("Protein_ID", "Gene", "Pathway",
                               setdiff(names(df_selected), c("Protein_ID", "Gene", "Pathway")))]

# Order by the column "Pathway"
df_selected <- df_selected[order(df_selected$Pathway), ]

# Save the result in a new file
write.table(df_selected, file = "1_table_FAMeDB_codes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Create a new table with the first three columns and the corresponding headers
df_counts <- data.frame(Protein_ID = df_selected$Protein_ID, Gene = df_selected$Gene, Pathway = df_selected$Pathway)

# Function to count the number of codes in each cell
count_codes <- function(x) {
  ifelse(grepl("\\*", x), 0, str_count(x, ",") + 1)
}

# Apply the function to all remaining columns and keep the headers
df_counts <- cbind(df_counts, sapply(df_selected[, 4:ncol(df_selected)], count_codes))

# Save the result in a new file
write.table(df_counts, file = "2_table_FAMeDB_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Remove "Gene" and "Pathway" columns
columns_to_remove <- c("Gene", "Pathway")
df_counts2 <- df_counts[, !(colnames(df_counts) %in% columns_to_remove)]

# Merge df_counts with all_pathways based on Protein_ID
df_counts_merged <- merge(df_counts2, all_pathways, by = "Protein_ID", all.x = TRUE)

# Drop columns "Mechanism", "central_pathway", "code_compound", and "code_pathway"
df_counts_merged <- df_counts_merged[, !(colnames(df_counts_merged) %in% c("Mechanism", "central_pathway", "code_compound", "code_pathway"))]

# Reorder the columns
df_counts_merged <- df_counts_merged[, c("Compound", "Superpathway", "Pathway", "Protein_ID", "Gene", setdiff(names(df_counts_merged), c("Compound", "Superpathway", "Pathway", "Protein_ID", "Gene")))]

# Organize df_counts_merged by Gene alphabetically
df_counts_merged <- df_counts_merged %>%
  arrange(Gene)

# Group by Gene and sum the counts across all columns ending with ".faa"
df_counts_aggregated <- df_counts_merged %>%
  group_by(Gene) %>%
  summarise(across(ends_with(".faa"), \(x) sum(x, na.rm = TRUE))) %>%
  distinct()

# Create a summary data frame for all_pathways (not used as all entries are distinct - gene names are all different)
#all_pathways_summary <- all_pathways %>%
  #select(-Protein_ID) %>%
  #distinct()

# Merge df_counts_aggregated with all_pathways_summary based on Gene
df_counts_aggregated_2 <- merge(df_counts_aggregated, all_pathways, by = "Gene", all.x = TRUE)

# Remove duplicates based on Gene
df_counts_aggregated_2 <- df_counts_aggregated_2 %>% distinct(Gene, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_2 <- df_counts_aggregated_2[, c("Compound", "Superpathway", "Pathway", "Protein_ID", "Gene", setdiff(names(df_counts_aggregated_2), c("Compound", "Superpathway", "Pathway", "Protein_ID", "Gene")))]

# Drop columns "central_pathway", "code_compound", and "code_pathway"
df_counts_aggregated_2 <- df_counts_aggregated_2[, !(colnames(df_counts_aggregated_2) %in% c("central_pathway", "code_compound", "code_pathway"))]

# Save the result in a new file
write.table(df_counts_aggregated_2, file = "3_table_FAMeDB_final.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Get all unique names in the "Compound" column
unique_compounds <- unique(df_counts_aggregated_2$Compound)

# Iterate over the list of unique compounds
for (compound in unique_compounds) {
  # Filter the DataFrame to get only rows with the current compound
  compound_df <- df_counts_aggregated_2[df_counts_aggregated_2$Compound == compound, ]
  
  # Create the file name with the desired format and save the DataFrame to a separate file
  file_name <- paste0("4_", gsub("/", "_", compound), ".tsv")
  write.table(compound_df, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}



# Define the base color palette (use enough colors for unique compounds)
base_colors_bubble <- c("#D62828", "#FDDE9F", "#197835", "snow", "#FDEFD5", 
                        "#D9F0BF", "#89ABBE", "#DE75AE", "#CE7C71", "#19974F",
                        "#C2A3CF", "#955540", "#6FB078", "#00616D", "#4291C2",
                        "#00619A") 



# Check if you have enough colors
if (length(unique_compounds) > length(base_colors_bubble)) {
  stop("Need more colors in the base_colors_bubble vector for all unique compounds.")
}

# Create the named color vector
# Maps each unique compound name to a specific color from the base palette
compound_colors_map <- setNames(
  base_colors_bubble[1:length(unique_compounds)],
  sort(unique_compounds) # Sort to ensure reproducible mapping
)

# Example (you would see: "CompoundA" = "#d62828", "CompoundB" = "#FDDE7A", etc.)
# print(compound_colors_map)

## BASED ON GENES

# Group by Gene and sum the counts across all columns ending with ".faa"
df_counts_aggregated_by_gene <- df_counts_merged %>%
  group_by(Gene) %>%
  summarise(across(ends_with(".faa"), \(x) sum(x, na.rm = TRUE))) %>%
  distinct()


# Create a summary data frame for all_pathways
all_pathways_summary <- all_pathways %>%
  select(-Protein_ID) %>%
  distinct()

# Merge df_counts_aggregated with all_pathways_summary based on Gene
df_counts_aggregated_by_gene_2 <- merge(df_counts_aggregated_by_gene, all_pathways_summary, by = "Gene", all.x = TRUE)

# Remove duplicates based on Gene
df_counts_aggregated_by_gene_2 <- df_counts_aggregated_by_gene_2 %>% distinct(Gene, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_by_gene_2 <- df_counts_aggregated_by_gene_2[, c("Mechanism", "Compound", "Superpathway", "Pathway", "Gene", setdiff(names(df_counts_aggregated_by_gene_2), c("Mechanism", "Compound", "Superpathway", "Pathway", "Gene")))]

# Drop columns "central_pathway", "code_compound", and "code_pathway"
df_counts_aggregated_by_gene_2 <- df_counts_aggregated_by_gene_2[, !(colnames(df_counts_aggregated_by_gene_2) %in% c("central_pathway", "code_compound", "code_pathway"))]

# Save the result in a new file
write.table(df_counts_aggregated_by_gene_2, file = "3_table_FAMeDB_final_by_gene.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Get all unique names in the "Compound" column
unique_compounds <- unique(df_counts_aggregated_by_gene_2$Compound)

# Iterate over the list of unique compounds
for (compound in unique_compounds) {
  # Filter the DataFrame to get only rows with the current compound
  compound_df <- df_counts_aggregated_by_gene_2[df_counts_aggregated_by_gene_2$Compound == compound, ]
  
  # Create the file name with the desired format and save the DataFrame to a separate file
  file_name <- paste0("4_", gsub("/", "_", compound), ".tsv")
  write.table(compound_df, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
}


# List to store the read DataFrames
df_list <- list()

# Iterate over the list of unique compounds
for (compound in unique_compounds) {
  # Create the file name with the desired format
  file_name <- paste0("4_", gsub("/", "_", compound), ".tsv")
  
  # Read the specific TSV file for the current compound
  compound_df <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t", comment.char = "", fill = TRUE)
  
  # Add the DataFrame to the list
  df_list[[file_name]] <- compound_df
}

# Apply melt to each DataFrame in the list
melted_df_list <- lapply(df_list, function(df) melt(df, id.vars = c("Gene", "Mechanism", "Compound", "Superpathway", "Pathway")))

# Access the DataFrames in the list using the file names as keys
# Example: df_list$"4_Aromatics.tsv"

# You can access the melted DataFrames through melted_df_list
# melted_df_list[[1]] contains the first melted DataFrame, and so on

# Define the new column names
new_column_names <- c("Gene", "Mechanism", "Compound", "Superpathway", "Pathway", "Genome", "Number")

# Iterate over the list of melted DataFrames
for (i in seq_along(melted_df_list)) {
  # Change the column names
  colnames(melted_df_list[[i]]) <- new_column_names
}

# Combine all melted DataFrames into one
combined_df <- bind_rows(melted_df_list, .id = "group")

# Use the sub() function with a single regular expression to remove both the prefix and suffix
combined_df$group <- sub(
  pattern = "^4_(.*)\\.tsv$",
  replacement = "\\1",
  x = combined_df$group
)

# Define the color names and legend titles
colors_heatmap <- c("snow", "#1F988BFF", "#440154FF") # You can add more colors as needed


# Generate the heatmap plot
svg("FAMeDB_plot_Heatmap_Genes.svg", width=15, height=8)
ggplot(combined_df, aes(x = Gene, y = Genome, fill = Number, group = group)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradientn(colors = colors_heatmap, name = "", labels = scales::comma) +
  scale_y_discrete(limits = rev) +
  facet_grid(~ group, scales = "free_x", space = "free") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(hjust = 1, size = 9),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text.x = element_text(size = 8, face = "bold")
  )
dev.off()


# Define the mapping vectors of central compounds names
central_compounds <- c("1,2,3,5‐Tetrahydroxybenzene", "3-Hydroxyanthranilate", "Catechol", "Protocatechuate",
                "Hydroxyquinol", "Gentisate", "Homogentisate", "Methoxyhidroquinone")
cc_abbreviations <- c("1235TeHB", "3HA", "CAT", "PTC", "HQ", "GT", "HGT", "MeH2Q")

# Create vector for renaming
long_2abbreviated <- setNames(cc_abbreviations, central_compounds)

# 3. Create a filtered data frame AND substitute the group names using mutate() and recode()
central_df <- combined_df %>%
  filter(group %in% central_compounds) %>%
  mutate(group = recode(group, !!!long_2abbreviated))

# Generate the heatmap plot using the filtered data
svg("FAMeDB_plot_Heatmap_Genes_Central_Pathways.svg", width=13, height=8)
ggplot(central_df, aes(x = Gene, y = Genome, fill = Number, group = group)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradientn(colors = colors_heatmap, name = "", labels = scales::comma) +
  scale_y_discrete(limits = rev) +
  facet_grid(~ group, scales = "free_x", space = "free") +
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(hjust = 1, size = 9),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text.x = element_text(size = 8, face = "bold")
  )
dev.off()


# Heatmap plotting Function for each unique compound
plot_compound_heatmap <- function(compound_name, data_df, colors_heatmap) {
  
  # Filter the data for the current compound
  filtered_df <- data_df %>%
    filter(group == compound_name)
  
  # Check if the filtered data frame is empty (optional but recommended)
  if (nrow(filtered_df) == 0) {
    message(paste("No data found for compound:", compound_name))
    return(NULL) # Skip plotting
  }
  
  # Construct the output file name dynamically
  output_filename <- paste0("FAMeDB_plot_Heatmap_Genes_", compound_name, ".svg")
  
  # Start the SVG device
  svg(output_filename, width=7, height=8)
  
  # Generate the plot
  p <- ggplot(filtered_df, aes(x = Gene, y = Genome, fill = Number, group = group)) +
    geom_tile(color = "white", lwd = 0.5) +
    scale_fill_gradientn(colors = colors_heatmap, name = "", labels = scales::comma) +
    scale_y_discrete(limits = rev) +
    facet_grid(~ group, scales = "free_x", space = "free") +
    theme_gray() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 9),
      axis.text.y = element_text(hjust = 1, size = 9),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 10, face = "bold")
    )
  
  # Print the plot to the SVG file
  print(p)
  
  # Close the SVG device
  dev.off()
}

# Execute the Function for All Compounds
lapply(unique_compounds, plot_compound_heatmap, 
       data_df = combined_df, 
       colors_heatmap = colors_heatmap)


## BASED ON PATHWAYS

# Group by Pathway and sum the counts across all columns ending with ".faa"
df_counts_aggregated_by_pathway <- df_counts_merged %>%
  group_by(Pathway) %>%
  summarise(across(ends_with(".faa"), \(x) sum(x, na.rm = TRUE))) %>%
  distinct()

# Create a summary data frame for all_pathways
all_pathways_summary <- all_pathways %>%
  select(-Protein_ID) %>%
  distinct()

# Merge df_counts_aggregated with all_pathways_summary based on Pathway
df_counts_aggregated_by_pathway_2 <- merge(df_counts_aggregated_by_pathway, all_pathways_summary, by = "Pathway", all.x = TRUE)

# Remove specific columns
columns_to_remove <- c("Gene")
df_counts_aggregated_by_pathway_2 <- df_counts_aggregated_by_pathway_2[, !(colnames(df_counts_aggregated_by_pathway_2) %in% columns_to_remove)]

# Remove duplicates based on Pathway
df_counts_aggregated_by_pathway_2 <- df_counts_aggregated_by_pathway_2 %>% distinct(Pathway, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_by_pathway_2 <- df_counts_aggregated_by_pathway_2[, c("Mechanism", "Compound", "Superpathway", "Pathway", setdiff(names(df_counts_aggregated_by_pathway_2), c("Mechanism", "Compound", "Superpathway", "Pathway")))]

# Drop columns "central_pathway", "code_compound", and "code_pathway"
df_counts_aggregated_by_pathway_2 <- df_counts_aggregated_by_pathway_2[, !(colnames(df_counts_aggregated_by_pathway_2) %in% c("central_pathway", "code_compound", "code_pathway"))]


# Generate the bubble plot

# Apply melt to obtain a general DataFrame with subpathways
df_counts_aggregated_by_pathway_2_melted <- melt(df_counts_aggregated_by_pathway_2, id.vars = c("Pathway", "Mechanism", "Compound", "Superpathway"))

# Rename the first three columns as ColumnA, ColumnB, and ColumnC
colnames(df_counts_aggregated_by_pathway_2_melted)[1:6] <- c("Pathway", "Mechanism", "Compound", "Superpathway", "Genome", "Hits_number")

# Define breaks for the size scale (using integers)
size_breaks_pathway <- seq(min(df_counts_aggregated_by_pathway_2_melted$Hits_number, na.rm = TRUE),
                           max(df_counts_aggregated_by_pathway_2_melted$Hits_number, na.rm = TRUE),
                           by = 1)
df_counts_aggregated_by_pathway_2_melted[df_counts_aggregated_by_pathway_2_melted == 0] <- NA


# Use the sub() function with a single regular expression to remove the suffix ".faa"
df_counts_aggregated_by_pathway_2_melted$Genome <- sub(
  pattern = "(.*)\\.faa$",
  replacement = "\\1",
  x = df_counts_aggregated_by_pathway_2_melted$Genome
)


svg("FAMeDB_plot_Bubbles_Pathways.svg", width = 18, height = 15) # change width value according to the number of proteomes
ggplot(df_counts_aggregated_by_pathway_2_melted, aes(x = Genome, y = Pathway, size = Hits_number, fill = Compound)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_manual(values = compound_colors_map) +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(breaks = size_breaks_pathway, range = c(2, 10)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 12, face = "bold")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6))) # Adjust the size value as needed for legend of compounds
dev.off()



## BASED ON CENTRAL PATHWAYS

# Filter for Central Pathways (central_pathway == "yes")
all_pathways_central <- all_pathways %>%
  filter(central_pathway == "Yes") %>%
  select(-Protein_ID) %>%
  distinct()


# Merge with central pathway information
# Use the Subathway column as the key for merging
df_counts_aggregated_by_pathway_3 <- merge(df_counts_aggregated_by_pathway, all_pathways_central, by = "Pathway", all.x = TRUE)

# Filter the aggregated data to keep only Central Pathways
df_counts_aggregated_by_pathway_3 <- df_counts_aggregated_by_pathway_3 %>%
  filter(!is.na(central_pathway)) # Keeps only rows that merged with central_pathways_summary (i.e., central_pathway == "Yes")

# Remove specific columns that are no longer needed for grouping
columns_to_remove_pathway <- c("Gene", "Superpathway")
df_counts_aggregated_by_pathway_3 <- df_counts_aggregated_by_pathway_3[, !(colnames(df_counts_aggregated_by_pathway_3) %in% columns_to_remove_pathway)]

# Remove duplicates based on Pathway
df_counts_aggregated_by_pathway_3 <- df_counts_aggregated_by_pathway_3 %>% distinct(Pathway, .keep_all = TRUE)

# Reorder the columns
df_counts_aggregated_by_pathway_3 <- df_counts_aggregated_by_pathway_3[, c("Mechanism", "Compound", "Pathway", setdiff(names(df_counts_aggregated_by_pathway_3), c("Mechanism", "Compound", "Pathway")))]

# Drop code columns
df_counts_aggregated_by_pathway_3 <- df_counts_aggregated_by_pathway_3[, !(colnames(df_counts_aggregated_by_pathway_3) %in% c("central_pathway", "code_compound", "code_pathway"))]


# Generate the bubble plot for Central Pathways

# Apply melt to obtain a general DataFrame with pathways
df_counts_aggregated_by_pathway_3_melted <- melt(df_counts_aggregated_by_pathway_3, id.vars = c("Pathway", "Mechanism", "Compound"))

# Rename the columns
colnames(df_counts_aggregated_by_pathway_3_melted)[1:5] <- c("Pathway", "Mechanism", "Compound", "Genome", "Hits_number")

# Define breaks for the size scale (using integers)
size_breaks_pathway <- seq(min(df_counts_aggregated_by_pathway_3_melted$Hits_number, na.rm = TRUE),
                           max(df_counts_aggregated_by_pathway_3_melted$Hits_number, na.rm = TRUE),
                           by = 1)
df_counts_aggregated_by_pathway_3_melted[df_counts_aggregated_by_pathway_3_melted == 0] <- NA


# Use the sub() function with a single regular expression to remove the suffix ".faa"
df_counts_aggregated_by_pathway_3_melted$Genome <- sub(
  pattern = "(.*)\\.faa$",
  replacement = "\\1",
  x = df_counts_aggregated_by_pathway_3_melted$Genome
)

svg("FAMeDB_plot_Bubbles_Central_Pathways.svg", width = 18, height = 9) # change width value according to the number of proteomes
ggplot(df_counts_aggregated_by_pathway_3_melted, aes(x = Genome, y = Pathway, size = Hits_number, fill = Compound)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_manual(values = compound_colors_map) +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(breaks = size_breaks_pathway, range = c(2, 10)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1, size = 15),
    axis.text.y = element_text(hjust = 1, size = 15),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 12, face = "bold")
  ) +
  guides(fill = guide_legend(override.aes = list(size = 6))) # Adjust the size value as needed for legend of compounds
dev.off()

