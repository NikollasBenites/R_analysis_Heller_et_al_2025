library(dplyr)
library(vegan)
library(tidyr)
library(purrr)
library(ggplot2)
library(reshape2)

setwd("~/Documents/R_PCA_files/R_PCA_Heller_et_al_2025/CSV_files")

# =======================
#  Load and Clean Dataset
# =======================
df <- read.csv("Combined_projection_all.csv")
# Replace underscores with spaces in feature names
df <- df %>% rename_all(~ gsub("_", " ", .))
colnames(df) <- c("PC1", "PC2", "PC3", "ID", "Type", "Firing Pattern","Cell ID")



run_pairwise_adonis <- function(data, id_col = "ID", pc_cols = c("PC1", "PC2", "PC3"), permutations = 10000, method = "euclidean") {
  group_combinations <- combn(unique(data[[id_col]]), 2, simplify = FALSE)
  
  results <- map_dfr(group_combinations, function(pair) {
    # Subset data
    pair_data <- data %>% filter(.data[[id_col]] %in% pair)
    
    # Extract response matrix and group vector
    response <- pair_data %>% select(all_of(pc_cols)) %>% as.matrix()
    group <- pair_data[[id_col]]
    
    # Run adonis2
    adonis_result <- adonis2(
      response ~ group,
      permutations = permutations,
      method = method
    )
    
    # Extract results
    tibble(
      group1 = pair[1],
      group2 = pair[2],
      F_value = adonis_result$F[1],
      R2 = adonis_result$R2[1],
      p_value = adonis_result$`Pr(>F)`[1]
    )
  })
  
  return(results)
}


pairwise_results <- run_pairwise_adonis(df)
print(pairwise_results)


print(ggplot(pairwise_results, aes(x = reorder(paste(group1, group2, sep = " vs "), -p_value), y = p_value)) +
  geom_col(fill = "steelblue") +
  theme_minimal() +
  labs(x = "Group Pair", y = "F-value", title = "PERMANOVA p-values by Group Pair") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)))

colnames(pairwise_results) <- c("Group 1","Group 2","F-Value","Rsqrd", "p-Value")
pairwise_results_cleaned <- pairwise_results %>% select(!c("Rsqrd"))

