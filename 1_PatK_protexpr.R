library(data.table)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# Assign treatment labels
df[, Treatment := ifelse(`Highest mean condition` %in% sulphur_weeks, "Sulphur",
                         ifelse(`Highest mean condition` %in% placebo_weeks, "Placebo", "Other"))]
top20 <- df[order(-`Max fold change`)][1:20]
top20 <- top20[`Max fold change` != "Inf"]

# Ensure the key columns have correct types
df <- df[!is.na(`Max fold change`)]
df <- df[Treatment %in% c("Sulphur", "Placebo")]


# Boxplot for fold change comparison with labeled top proteins
ggplot(df[Treatment %in% c("Sulphur", "Placebo")],
       aes(x = Treatment, y = `Max fold change`, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +  # Boxplot without extreme outliers
  geom_jitter(width = 0.2, alpha = 0.5) +  # Adds jittered points for better visibility
  geom_text_repel(data = top20,
                  aes(label = Protein_Name, y = `Max fold change`),
                  size = 4, max.overlaps = 10) +
  theme_minimal() +
  scale_fill_manual(values = c("Sulphur" = "#1b9e77", "Placebo" = "#d95f02")) +  # Custom colors
  labs(title = "Effect of Sulphur vs Placebo on Protein Expression",
       y = "Max Fold Change",
       x = "Treatment Group")

# Outlier removal for 'Max fold change' using IQR method
Q1 <- quantile(df$`Max fold change`, 0.25, na.rm = TRUE)
Q3 <- quantile(df$`Max fold change`, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

# Assign treatment labels
sulphur_weeks <- c("W2", "W8", "W24", "W28")
placebo_weeks <- c("W4", "W12", "W16", "W20")

ggplot(df, aes(x=p_value)) + geom_histogram(bins=30, fill="blue", alpha=0.7) + theme_minimal()
ggplot(df, aes(x=log2_FC)) + geom_histogram(bins=30, fill="red", alpha=0.7) + theme_minimal()

# Define the treatment and placebo weeks
sulphur_weeks <- c("W2", "W8", "W24", "W28")
placebo_weeks <- c("W4", "W12", "W16", "W20")

######
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "gprofiler2", "enrichR"))
install.packages("data.table")

library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichR)
library(gprofiler2)
library(biomaRt)

# Define the organism for annotation using Ensembl
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Check available attributes (if needed)
# biomaRt::listAttributes(mart)

# Use correct attribute for UniProt or Ensembl protein IDs
gene_mapping <- biomaRt::getBM(
  attributes = c("uniprotswissprot", "external_gene_name"), # Try "uniprotswissprot"
  filters = "uniprotswissprot",
  values = unique(df$Accession),  # Ensure unique accessions
  mart = mart
)

# Merge Gene Symbols back into the dataset
df <- merge(df, gene_mapping, by.x = "Accession", by.y = "uniprotswissprot", all.x = TRUE)
df_match <- df[, c("Accession","Gene","external_gene_name")] # CHECK if is it ok?

# GO Enrichment for each week
go_results <- df[, {
  enriched <- enrichGO(
    gene         = Gene,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    pAdjustMethod = "fdr",
    pvalueCutoff  = 0.05
  )
  as.data.table(enriched@result)
}, by = Treatment]

# Select top 20 GO terms per Treatment group based on adjusted p-value
top_go <- go_results[order(p.adjust), .SD[1:20], by = Treatment]

# For a better y-axis ordering we can order the GO terms within each treatment.
# Here we create a combined label that includes Treatment and Description.
top_go[, term := paste0(Description)]
top_go[, term := factor(term, levels = unique(term[order(p.adjust)]))]

ggplot(top_go, aes(x = Count, y = term)) +
  geom_point(aes(color = p.adjust, size = FoldEnrichment)) +
  facet_wrap(~ Treatment, scales = "free_y") +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  labs(title = "",
       x = "Gene Count",
       y = "GO Term",
       color = "FDR-adjusted p-value") +
  theme(axis.text.y = element_text(size = 8))

library(gprofiler2)
library(ggplot2)

# Suppose you have a vector of gene symbols from your proteomics data:


df[, significance := "NS"]
df[log2_FC > 0  & p_value < 0.05, significance := "Up"]
df[log2_FC < 0 & p_value < 0.05, significance := "Down"]

gene_list <- unique(df$Gene)  # adjust to your gene column
top_up <- df[significance == "Up"][order(p_value)][1:100]
top_down <- df[significance == "Down"][order(p_value)][1:100]

gene_up <- gost(top_up$Gene, organism = "hsapiens", correction_method = 'fdr')
gene_down <- gost(top_down$Gene, organism = "hsapiens", correction_method = 'fdr')

# Run g:Profiler enrichment using selected databases (sources).
gostres <- gost(query = gene_list,
                organism = "hsapiens",
#                sources = c("GO:BP", "GO:MF", "REAC", "KEGG", "TF"),
                correction_method = "fdr")

# Check the results structure
str(gostres)
gostplot(gene_up, interactive = T)
p <- gostplot(gene_up, interactive = F)
p1 <- gostplot(gene_down, interactive = F)

top15_terms_up <- gene_up$result[order(gene_up$result$p_value), "term_id"][1:15]
top15_terms_down <- gene_down$result[order(gene_down$result$p_value), "term_id"][1:15]

publish_gostplot(p, highlight_terms = top15_terms_up, filename = NULL)
publish_gostplot(p1, highlight_terms = top15_terms_down, filename = NULL)

# Load EnrichR databases
dbs <- enrichR::listEnrichrDbs()
dbs <- c("DisGeNET", "DrugBank")
dbs <- c("Reactome_Pathways_2024")
dbs <- c("Proteomics_Drug_Atlas_2023")

top_up <- df[significance == "Up"][order(p_value)][1:100]
top_down <- df[significance == "Down"][order(p_value)][1:100]

# Run enrichment analysis for disease and drug associations
enrichr_results <- top_down[, {
  enriched <- enrichR::enrichr(Gene, dbs)
  as.data.table(enriched$Proteomics_Drug_Atlas_2023)  # Change the DB here
}, by = Treatment]

top15_down <- enrichr_results[order(enrichr_results$Adjusted.P.value)][1:20]

enrichr_results <- top_up[, {
  enriched <- enrichR::enrichr(Gene, dbs)
  as.data.table(enriched$Proteomics_Drug_Atlas_2023)  # Change the DB here
}, by = Treatment]

top15_up <- enrichr_results[order(enrichr_results$Adjusted.P.value)][1:20]

top15 <- rbind(top15_up, top15_down)

ggplot(top15, aes(x = Odds.Ratio, y = Term)) +
  geom_point(aes(color = Adjusted.P.value, size = Combined.Score)) +
  facet_wrap(~ Treatment, scales = "free_y") +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  labs(title = "",
       x = "Odds Ratio",
       y = "Proteomics Drug Term", #Change the DB
       color = "FDR-adjusted p-value") +
  theme(axis.text.y = element_text(size = 8))


########

library(monarchr)
library(tidygraph)
library(dplyr)
library(igraph)

dep_search <- monarch_search("Depression", limit = 40)



diseases <- cbind(V(dep_search)$id, V(dep_search)$name, V(dep_search)$description)
diseases <- as.data.frame(diseases)

dep_diseases <- monarch_engine() |>
  fetch_nodes(query_ids = c("HP:0000716", "MONDO:0012048", "NBO:0000085",
                            "MP:0014414")) |>
  expand(predicates = "biolink:subclass_of", direction = "in", transitive = TRUE) |>
  expand(predicates = c("biolink:causes", "biolink:gene_associated_with_condition")) |>
  expand(categories = "biolink:Gene")

plot(dep_diseases, layout = "sugiyama")

gene_names <- V(dep_diseases)$name
print(head(gene_names))


###########

library(dplyr)
library(tidyr)
library(ggplot2)

# Filter df to only include genes from gene_names and the desired treatments.
df_filtered <- df %>%
  filter(tolower(Gene) %in% tolower(gene_names),
         Treatment %in% c("placebo", "sulphur"))


# Select the week columns that contain the mean values. Adjust week_cols if necessary.
week_cols <- c("W2", "W4", "W8", "W12", "W16", "W20", "W24", "W28")

# Pivot the data from wide to long format.
df_long <- df_filtered %>%
  pivot_longer(cols = all_of(week_cols), names_to = "Week", values_to = "mean_value")

# Reshape the data so that each row represents a unique Gene and Week with separate columns for each Treatment.
df_wide <- df_long %>%
  pivot_wider(names_from = Treatment, values_from = mean_value)

# Calculate the log2 fold change for each gene and week.
df_wide <- df_wide %>%
  mutate(log2_FC = log2(sulphur / placebo))

# Plot the log2 fold change across weeks for each gene.
ggplot(df_wide, aes(x = Week, y = log2_FC, group = Gene, color = Gene)) +
  geom_line() +
  geom_point() +
  labs(title = "Log2 Fold Change (Sulphur vs Placebo) Across Weeks",
       x = "Week",
       y = "Log2 Fold Change") +
  theme_minimal()
