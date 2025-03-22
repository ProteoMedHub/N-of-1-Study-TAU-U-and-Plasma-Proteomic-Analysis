library(data.table)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(stringr)

getwd()
setwd("~/Desktop/psy/PatA")
df <- fread("All_identified_proteins_PatA_filtered_norm.csv", sep = ";", header = TRUE)


df[, Protein_Name := str_extract(Description, "^[^O]+")]
df[, Gene := gsub(".*GN=([^ ]+).*", "\\1", Description)]

head(df)
str(df)
summary(df)
colSums(is.na(df))
df[, `Confidence score` := as.numeric(`Confidence score`)]
df[, `Anova (p)` := as.numeric(`Anova (p)`)]
df[, `q Value` := as.numeric(`q Value`)]
df[, `Max fold change` := as.numeric(`Max fold change`)]
df[, `Mass` := as.numeric(`Mass`)]

# Define the treatment and placebo weeks
kali_weeks <- c("W2", "W8", "W16", "W20")
placebo_weeks <- c("W4", "W12", "W24", "W28")

# Assign treatment labels
df[, Treatment := ifelse(`Highest mean condition` %in% kali_weeks, "kali",
                         ifelse(`Highest mean condition` %in% placebo_weeks, "Placebo", "Other"))]
top20 <- df[order(-`Max fold change`)][1:5]

df <- df[(df$Treatment != "Other"), ]

# Define the columns for each replicate for kali and placebo groups.
kali_cols_rep1 <- c("W2_20241203_LSanders_A2_01", "W8_20241203_LSanders_A4_01",
                       "W16_20241203_LSanders_A6_01", "W20_20241203_LSanders_A7_01")
kali_cols_rep2 <- c("W2_20241203_LSanders_A2_02", "W8_20241203_LSanders_A4_02",
                       "W16_20241203_LSanders_A6_02", "W20_20241203_LSanders_A7_02")
kali_cols_rep3 <- c("W2_20241203_LSanders_A2_03", "W8_20241203_LSanders_A4_03",
                       "W16_20241203_LSanders_A6_03", "W20_20241203_LSanders_A7_03")

placebo_cols_rep1 <- c("W4_20241203_LSanders_A3_01", "W12_20241203_LSanders_A5_01",
                       "W24_20241203_LSanders_A8_01", "W28_20241203_LSanders_A9_01")
placebo_cols_rep2 <- c("W4_20241203_LSanders_A3_02", "W12_20241203_LSanders_A5_02",
                       "W24_20241203_LSanders_A8_02", "W28_20241203_LSanders_A9_02")
placebo_cols_rep3 <- c("W4_20241203_LSanders_A3_03", "W12_20241203_LSanders_A5_03",
                       "W24_20241203_LSanders_A8_03", "W28_20241203_LSanders_A9_03")

# Calculate the mean for each replicate across the weeks within each condition.
df[, kali_rep1 := rowMeans(.SD, na.rm = TRUE), .SDcols = kali_cols_rep1]
df[, kali_rep2 := rowMeans(.SD, na.rm = TRUE), .SDcols = kali_cols_rep2]
df[, kali_rep3 := rowMeans(.SD, na.rm = TRUE), .SDcols = kali_cols_rep3]

df[, placebo_rep1 := rowMeans(.SD, na.rm = TRUE), .SDcols = placebo_cols_rep1]
df[, placebo_rep2 := rowMeans(.SD, na.rm = TRUE), .SDcols = placebo_cols_rep2]
df[, placebo_rep3 := rowMeans(.SD, na.rm = TRUE), .SDcols = placebo_cols_rep3]

# For each protein, perform a paired t-test comparing the three paired replicates.
df[, p_value := {
  s_reps <- as.numeric(.SD[, .(kali_rep1, kali_rep2, kali_rep3)])
  p_reps <- as.numeric(.SD[, .(placebo_rep1, placebo_rep2, placebo_rep3)])
  t.test(s_reps, p_reps, paired = TRUE)$p.value
}, by = Accession]

# Calculate overall mean abundances for each condition and compute log2 fold change.
df[, kali_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("kali_rep1", "kali_rep2", "kali_rep3")]
df[, placebo_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = c("placebo_rep1", "placebo_rep2", "placebo_rep3")]
df[, log2_FC := log2(kali_mean / placebo_mean)]

# Define significance based on a fold change threshold of ±1 (i.e. 2-fold change) and p-value cutoff of 0.05.
df[, significance := "NS"]
df[log2_FC > 1  & p_value < 0.05, significance := "Up"]
df[log2_FC < -1 & p_value < 0.05, significance := "Down"]

# Create volcano plot.
volcano_plot <- ggplot(df, aes(x = log2_FC, y = -log10(p_value), color = significance)) +
  geom_point() +
  # Optionally label top proteins by fold change
  geom_text_repel(data = df[order(-abs(log2_FC))][1:10],
                  aes(label = Protein_Name),
                  size = 4) +
  theme_minimal() +
  labs(title = "",
       x = "log2 Fold Change",
       y = "-log10 p-value") +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray"))
print(volcano_plot)

# Subset the top 25 up- and top 25 downregulated proteins
top_up <- df[significance == "Up"][order(p_value)][1:50]
top_down <- df[significance == "Down"][order(p_value)][1:50]
top_proteins <- rbind(top_up, top_down)

# Re-create the per-week mean matrix if not already done.
df[, W2 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W2_", names(df), value = TRUE)]
df[, W8 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W8_", names(df), value = TRUE)]
df[, W24 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W24_", names(df), value = TRUE)]
df[, W28 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W28_", names(df), value = TRUE)]
df[, W4 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W4_", names(df), value = TRUE)]
df[, W12 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W12_", names(df), value = TRUE)]
df[, W16 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W16_", names(df), value = TRUE)]
df[, W20 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^W20_", names(df), value = TRUE)]

# Create the heatmap matrix (proteins in rows, week means in columns)
heatmap_data <- as.matrix(df[, .(W2, W8, W16, W20,
                                 W4, W12, W24, W28)])
rownames(heatmap_data) <- df$Protein_Name

# Scale the data by row to highlight relative differences.
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# Subset the matrix to include only the selected proteins.
selected_heatmap_data <- heatmap_data_scaled[rownames(heatmap_data_scaled) %in% top_proteins$Protein_Name, ]

# Create an annotation for the columns to indicate treatment group.
annotation_col <- data.frame(
  Treatment = rep(c("kali", "Placebo"), each = 4),
  row.names = colnames(selected_heatmap_data)
)

# Generate the heatmap.
pheatmap(selected_heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col,
         main = "")

