#
# Load libraries
library(tidyverse)

# Create dataset
data <- tribble(
  ~Dataset, ~Species, ~Complete, ~Missing,

  "CompOmics", "Homo sapiens", 1974, 1576,
  "CompOmics", "Saccharomyces cerevisiae", 552, 896,
  "CompOmics", "Escherichia coli", 57, 292,
  "CompOmics", "NAs", 29, 18,

  "MaxQuant", "Homo sapiens", 3227, 303,
  "MaxQuant", "Saccharomyces cerevisiae", 1044, 246,
  "MaxQuant", "Escherichia coli", 181, 140,
  "MaxQuant", "NAs", 48, 6,

  "Proline", "Homo sapiens", 2608, 758,
  "Proline", "Saccharomyces cerevisiae", 797, 472,
  "Proline", "Escherichia coli", 103, 214,
  "Proline", "NAs", 39, 8,

  "TPP", "Homo sapiens", 2913, 746,
  "TPP", "Saccharomyces cerevisiae", 888, 532,
  "TPP", "Escherichia coli", 117, 243,
  "TPP", "NAs", 58, 6
)

# Convert to long format
data_long <- data %>%
  pivot_longer(cols = c(Complete, Missing),
               names_to = "Type",
               values_to = "Count")

# Plot
library(ggplot2)

# Plot
p <- ggplot(data_long, aes(x = Species, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_wrap(~ Dataset, nrow = 2) +
  labs(
    title = "Species Composition Across Proteomics Software Outputs",
    subtitle = "Comparison of Complete vs Missing Submatrices",
    x = "Species",
    y = "Protein Count",
    fill = "Data Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  ) +
  scale_fill_manual(
    values = c(
      "Complete" = "#0072B2",  # blue
      "Missing"  = "#D55E00"   # orange
    )
  )

# Print plot
print(p)
