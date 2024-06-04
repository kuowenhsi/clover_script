# Load necessary libraries for data manipulation and visualization
library(tidyverse)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(paletteer)

# Set working directory to the specified path
setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/Molecular_ecology/clover_script")

# Get a sorted list of data files matching pattern "Q" in the data directory
files_data <- str_sort(list.files("./data")[str_detect(list.files("./data"), "Q")], numeric = TRUE)

# Initialize an empty list to store admix data
admix_data_list <- list()

# Loop through the first 8 files and process them
for (i in 1:8) {
  admix_data_list[[i]] <- read_delim(paste0("./data/", files_data[[i]]), col_names = FALSE, delim = " ") %>%
    mutate(K = i) %>%  # Add column K with the current loop index
    bind_cols(read_tsv("./data/figure_3_bwa_415_20230527_500kb_0.2_passed.psam") %>%
                select(ID = `#IID`, Pop, Latitude, Longitude)) %>%
    pivot_longer(cols = starts_with("X"), names_to = "ancestry", values_to = "ancestry_prop")  # Reshape data to long format
}

# Display the third element in admix_data_list for verification
admix_data_list[[3]]

# Combine all admix data into one dataframe
admix_data <- bind_rows(admix_data_list) %>%
  group_by(K, Pop, ancestry) %>%
  summarise(mean_ancestry_prop = mean(ancestry_prop), Latitude = median(Latitude)) %>%
  group_by(K) %>%
  arrange(Latitude, .by_group = TRUE) %>%
  mutate(Pop = factor(Pop, levels = unique(Pop)))  # Order populations by their latitude

# Create a bar plot for mean ancestry proportion
p <- ggplot(filter(admix_data, K > 1, K < 6), aes(x = Pop, y = mean_ancestry_prop)) +
  geom_bar(aes(fill = ancestry), stat = "identity", position = "stack", show.legend = FALSE, alpha = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_paletteer_d("ggthemes::Classic_10") +
  facet_grid(K ~ .)

# Display the plot
print(p)

# Save the plot to a file
ggsave("./figures/figure_3_A.png", width = 10, height = 5, dpi = 900)

##########################################

# Initialize a list to store PCA data
pca_data <- list()

# Read PCA data with different LD stringency settings
pca_data[["pca_no_LD"]] <- read_tsv("./data/figure_3_bwa_415_20230527_pca.eigenvec")[, 1:11]
pca_data[["pca_100k_0.8"]] <- read_tsv("./data/figure_3_bwa_415_20230527_100kb_0.8_passed_pca.eigenvec")[, 1:11]
pca_data[["pca_200k_0.5"]] <- read_tsv("./data/figure_3_bwa_415_20230527_200kb_0.5_passed_pca.eigenvec")[, 1:11]
pca_data[["pca_500k_0.2"]] <- read_tsv("./data/figure_3_bwa_415_20230527_500kb_0.2_passed_pca.eigenvec")[, 1:11]

# Read additional data for merging
Sara_final <- read_csv("./data/figure_3_WildCloverInfo_Final.csv") %>%
  select(1, 8) %>%
  separate(`Nearest City`, into = c("City", "State"), sep = ", ") %>%
  left_join(read_csv("./data/figure_3_US_states.csv"), by = c("State" = "TWO_LETTERS"))

cyano_climate_data <- read_csv("./data/figure_3_cyano_buf_climate_data_419.csv") %>%
  left_join(Sara_final, by = c("genotype" = "Accession"))

# Combine all PCA data and merge with climate data
pca_data <- bind_rows(pca_data, .id = "LD_stringency") %>%
  left_join(cyano_climate_data, by = c("#IID" = "genotype")) %>%
  arrange(State) %>%
  mutate(Pop = factor(Pop, levels = unique(Pop)))

# Print column names of the combined PCA data
print(colnames(pca_data))

# Create a named vector for population-state mapping
Pop_state <- paste(pca_data$Pop, pca_data$full_name, sep = " - ")
names(Pop_state) <- pca_data$Pop

# Read and process PCA eigenvalues
pca_PVE <- bind_cols(
  read_tsv("./data/figure_3_bwa_415_20230527_pca.eigenval", col_names = "no_LD"),
  read_tsv("./data/figure_3_bwa_415_20230527_100kb_0.8_passed_pca.eigenval", col_names = "100k_0.8"),
  read_tsv("./data/figure_3_bwa_415_20230527_200kb_0.5_passed_pca.eigenval", col_names = "200k_0.5"),
  read_tsv("./data/figure_3_bwa_415_20230527_500kb_0.2_passed_pca.eigenval", col_names = "500k_0.2")
) %>%
  mutate(PVE_no_LD = no_LD / sum(no_LD),
         PVE_100k_0.8 = `100k_0.8` / sum(`100k_0.8`),
         PVE_200k_0.5 = `200k_0.5` / sum(`200k_0.5`),
         PVE_500k_0.2 = `500k_0.2` / sum(`500k_0.2`))

# Summarize PCA data by population and region
pca_data_mean <- pca_data %>%
  group_by(LD_stringency, Pop, full_name, State, Region) %>%
  summarise(Longitude = median(Longitude), Latitude = median(Latitude), PC1 = median(PC1), PC2 = median(PC2), PC3 = median(PC3), PC4 = median(PC4))

# Create PCA plot for PC1 vs PC2 with 100k_0.8 stringency
p2_a <- ggplot(data = filter(pca_data_mean, LD_stringency == "pca_100k_0.8"), aes(x = PC1, y = PC2)) +
  geom_point(data = filter(pca_data, LD_stringency == "pca_100k_0.8"), color = "gray90") +
  geom_point(aes(color = Region, shape = Region)) +
  geom_text_repel(aes(label = State, color = Region), max.overlaps = 25, show.legend = FALSE, size = 3.5) +
  scale_shape_manual(values = rep(c(15, 17, 19), 3)) +
  xlab(paste0("PC1 (", round(pca_PVE$PVE_100k_0.8[[1]] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(pca_PVE$PVE_100k_0.8[[2]] * 100, digits = 2), "%)")) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Create PCA plot for PC1 vs Latitude with 100k_0.8 stringency
p4_a <- ggplot(data = filter(pca_data_mean, LD_stringency == "pca_100k_0.8"), aes(y = PC1, x = Latitude)) +
  geom_point(data = filter(pca_data, LD_stringency == "pca_100k_0.8"), color = "gray90") +
  geom_point(aes(color = Region, shape = Region)) +
  geom_text_repel(aes(label = State, color = Region), max.overlaps = 20, show.legend = FALSE, size = 3.5) +
  scale_shape_manual(values = rep(c(15, 17, 19), 3)) +
  ylab(paste0("PC1 (", round(pca_PVE$PVE_100k_0.8[[1]] * 100, digits = 2), "%)")) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Combine PCA plots and add legend
p <- plot_grid(plot_grid(p2_a + theme(legend.position = "none"), p4_a + theme(legend.position = "none")),
               get_legend(p2_a + theme(legend.position = "bottom")), nrow = 2, align = "h", axis = "l", rel_heights = c(5, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

# Save the combined plot to a file
ggsave("./figures/figure_3_C_D.png", width = 10, height = 6, dpi = 900)

               