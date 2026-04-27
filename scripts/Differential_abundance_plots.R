### This differential abundance analysis is performed using Maaslin2 and was used to see if there are any differential abundant taxa in Slow and Fast groups. 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Differential abundance with Maaslin2
library(tidyverse)
library(phyloseq)
library(readxl)
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")
library(Maaslin2)
library(writexl)
library(forcats)
library(ggrepel)


ps <- readRDS("datasets/ps_filtered.rds")

# Tax_glom if desired. DON'T rarefy or use a compositional transformation!
ps_glom <- tax_glom(ps, taxrank = "Genus")

# Subset to hindgut samples
ps_hindgut <- subset_samples(ps_glom, sample_type == "hindgut")


sample_data(ps_hindgut)$swim_group <- fct_collapse(
  sample_data(ps_hindgut)$swim_performance,
  fast = c("accelerator", "cruiser sprinter", "manoeuvrer"),
  moderate = c("generalist"),
  slow = c("flow refuging", "burrowing")
)

# Make sure moderate is the reference group
sample_data(ps_hindgut)$swim_group <- factor(
  sample_data(ps_hindgut)$swim_group,
  levels = c("moderate", "slow", "fast")
)

out <- as(out_table(ps_hindgut), "matrix")

# Ensure rows = samples and columns = taxa
if (taxa_are_rows(ps_hindgut)) {
  out <- t(out)
}

otu_df <- as.data.frame(otu_mat)

meta_df <- as(sample_data(ps_hindgut), "data.frame")

# Keep metadata in the same sample order as abundance table
meta_df <- meta_df[rownames(otu_df), , drop = FALSE]


# Run MaAsLin2 differential abundance
set.seed(421)

fit_data <- Maaslin2(
  input_data = otu_df,
  input_metadata = meta_df,
  output = "maaslin2_hindgut_results",
  fixed_effects = c("swim_group"),
  reference = "swim_group,moderate",
  normalization = "TSS",
  transform = "AST",
  min_abundance = 0.001,   # 0.1% relative abundance
  min_prevalence = 0.10,   # present in at least 10% of samples
  max_significance = 0.05,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

maaslin_results <- fit_data$results
write_xlsx(maaslin_results, "Fish_hindgut_MaAsLin2_results.xlsx")

# Taxonomy information
tax_df <- as.data.frame(tax_table(ps_hindgut))
tax_df$feature <- rownames(tax_df)

tax_df <- tax_df %>%
  mutate(
    Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", as.character(Genus))
  )

maaslin_results$feature <- sub("^X", "", maaslin_results$feature)

merged_maaslin <- maaslin_results %>%
  left_join(tax_df, by = "feature")


# Prepare volcano plot data
volcano_data <- merged_maaslin %>%
  filter(
    metadata == "swim_group",
    value %in% c("slow", "fast")
  ) %>%
  mutate(
    # Prevent issues if qval is exactly 0
    plot_qval = if_else(qval == 0, 1e-300, qval),
    
    significance = case_when(
      qval < 0.05 & coef > 0 ~ "Higher than moderate",
      qval < 0.05 & coef < 0 ~ "Lower than moderate",
      TRUE ~ "Not significant"
    ),
    
    taxon_label = ifelse(
      !is.na(Genus) & Genus != "Unclassified",
      Genus,
      feature
    )
  )

# Split into two comparisons
slow_volcano_data <- volcano_data %>%
  filter(value == "slow")

fast_volcano_data <- volcano_data %>%
  filter(value == "fast")


# Volcano plot slow vs moderate
slow_volcano_plot <- ggplot(
  slow_volcano_data,
  aes(x = coef, y = -log10(plot_qval), color = significance)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = subset(slow_volcano_data, qval < 0.05),
    aes(label = taxon_label),
    size = 5,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Higher than moderate" = "#E07A9A",
    "Lower than moderate" = "#fcba65",
    "Not significant" = "#8e7ca6"
  )) +
  theme(text = element_text(family = "Arial")) +
  labs(
    title = "Slow vs. moderate-swimming fish",
    x = "MaAsLin2 coefficient",
    y = "-log10(q-value)",
    color = "Association",
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )

slow_volcano_plot

ggsave(
  filename = "Volcano_slow_vs_moderate.png",
  plot = slow_volcano_plot,
  width = 8,
  height = 7,
  dpi = 300
  )

fast_volcano_plot <- ggplot(
  fast_volcano_data,
  aes(x = coef, y = -log10(plot_qval), color = significance)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = subset(fast_volcano_data, qval < 0.05),
    aes(label = taxon_label),
    size = 3,
    max.overlaps = 25,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Higher than moderate" = "#E07A9A",
    "Lower than moderate" = "#fcba65",
    "Not significant" = "#8e7ca6"
  )) +
  theme(text = element_text(family = "Arial")) +
  labs(
    title = "Fast vs. moderate-swimming fish",
    x = "MaAsLin2 coefficient",
    y = "-log10(q-value)",
    color = "Association"
  )

fast_volcano_plot

ggsave(
  filename = "Volcano_fast_vs_moderate.png",
  plot = fast_volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

slow_sig <- slow_volcano_data %>%
  filter(qval < 0.05) %>%
  arrange(qval)

fast_sig <- fast_volcano_data %>%
  filter(qval < 0.05) %>%
  arrange(qval)

write_xlsx(
  list(
    slow_vs_moderate = slow_sig,
    fast_vs_moderate = fast_sig
  ),
  "Fish_hindgut_significant_taxa.xlsx"
)

