# Load necessary libraries
library(tidyverse)
library(broom)
library(vegan)
library(FactoMineR)
library(ggrepel)
library(ggforce)
library(sp)

# Set working directory
setwd("./clover_script")

# Read and process FST data
fst_input <- read_tsv("./data/figure_4_bwa_415_20230527_100kb_0.8_fst.fst.summary") %>%
  mutate(linearized_fst = HUDSON_FST / (1 - HUDSON_FST))

# Rename columns
colnames(fst_input) <- c("POP1", "POP2", "FST", "L_FST")

# Calculate quantiles for FST
quantile(fst_input$FST)

# Get unique population names
pop_names <- unique(c(fst_input$POP1, fst_input$POP2))

# Split FST data by population and process outliers
fst_by_pop_list <- list()
for (i in pop_names) {
  fst_by_pop_list[[i]] <- fst_input[(fst_input$POP1 == i) | (fst_input$POP2 == i), ]
}
fst_by_pop <- bind_rows(fst_by_pop_list, .id = "REF") %>%
  group_by(REF) %>%
  mutate(OUTLIER = FST > (quantile(FST, probs = 0.75) + 1.5 * IQR(FST))) %>%
  mutate(OUTLIER_LABEL = case_when((OUTLIER == TRUE) & (REF == POP1) ~ POP2,
                                   (OUTLIER == TRUE) & (REF != POP1) ~ POP1,
                                   TRUE ~ as.character(NA))) %>%
  arrange(desc(FST), .by_group = TRUE) %>%
  mutate(count_REF = cumsum(OUTLIER)) %>%
  mutate(LABEL_x = case_when((OUTLIER == TRUE) & (count_REF %% 2 == 1) ~ 0.7,
                             (OUTLIER == TRUE) & (count_REF %% 2 == 0) ~ 1.3,
                             TRUE ~ as.numeric(NA)))

# Read and process climate data
data_input <- read_csv("./data/figure_4_cyano_buf_limei_climate_data_415.csv")[, c(2:4, 10:44)]

# Clean column names
colnames(data_input) <- colnames(data_input) %>% str_remove("current_30arcsec_|wc2.1_30s_")
colnames(data_input)[20:38] <- colnames(data_input)[26:44] <- c(
  "Bio1 - Annual Mean Temperature",
  "Bio10 - Mean Temperature of Warmest Quarter",
  "Bio11 - Mean Temperature of Coldest Quarter",
  "Bio12 - Annual Precipitation",
  "Bio13 - Precipitation of Wettest Month",
  "Bio14 - Precipitation of Driest Month",
  "Bio15 - Precipitation Seasonality",
  "Bio16 - Precipitation of Wettest Quarter",
  "Bio17 - Precipitation of Driest Quarter",
  "Bio18 - Precipitation of Warmest Quarter",
  "Bio19 - Precipitation of Coldest Quarter",
  "Bio2 - Mean Diurnal Range",
  "Bio3 - Isothermality",
  "Bio4 - Temperature Seasonality",
  "Bio5 - Max Temperature of Warmest Month",
  "Bio6 - Min Temperature of Coldest Month",
  "Bio7 - Temperature Annual Range",
  "Bio8 - Mean Temperature of Wettest Quarter",
  "Bio9 - Mean Temperature of Driest Quarter"
)

# Inspect unique values for growing degree days
unique(data_input$growingDegDays5)

# Perform PCA on climate data
clim_PCA <- PCA(data_input, quali.sup = 1, quanti.sup = 2:3)
clim_PCA$eig
clim_PCA$var$coord
PCA_cor <- as.data.frame(clim_PCA$var$cor) %>%
  rownames_to_column("env") %>%
  arrange(desc(abs(Dim.1))) %>%
  mutate(ID = seq(n()))

# Extract PCA individual coordinates and combine with data
clim_PCA_indv <- as.data.frame(clim_PCA$ind$coord) %>%
  bind_cols(data_input[, 2:3])

# Aggregate PCA results by population
clim_PCA_Pop <- as.data.frame(clim_PCA$quali.sup$coord) %>%
  rownames_to_column("Pop") %>%
  bind_cols(data_input[, 1:3] %>% group_by(Pop) %>% summarise_all(mean) %>% select(2:3))

# Combine climate data with PCA results
clim_input <- bind_cols(data_input, clim_PCA$ind$coord) %>%
  group_by(Pop) %>%
  summarise_all(mean)

# Convert to spatial format
clim_input_sf <- st_as_sf(clim_input, coords = c("Longitude", "Latitude"), crs = 4326)

# Calculate geographic distances
geo_dist <- st_distance(clim_input_sf, clim_input_sf)
rownames(geo_dist) <- clim_input$Pop
colnames(geo_dist) <- clim_input$Pop
geo_dist_df <- tidy(as.dist(geo_dist, diag = FALSE, upper = FALSE))
colnames(geo_dist_df) <- c("POP1", "POP2", "geo_dist")

# Calculate climate distances
clim_dist <- bind_cols(lapply(clim_input[4:40], function(x) as.numeric(dist(scale(x)))))

# Combine geographic, climate, and FST data
geo_clim_dist <- bind_cols(geo_dist_df, clim_dist) %>%
  left_join(fst_input, by = c("POP1", "POP2"))

# Prepare FST distance matrices
fst_matrix <- matrix(NA, nrow = 43, ncol = 43)
fst_matrix[lower.tri(fst_matrix)] <- geo_clim_dist$FST
fst_dist <- as.dist(fst_matrix, diag = FALSE, upper = FALSE)

L_fst_matrix <- matrix(NA, nrow = 43, ncol = 43)
L_fst_matrix[lower.tri(L_fst_matrix)] <- geo_clim_dist$L_FST
L_fst_dist <- as.dist(L_fst_matrix, diag = FALSE, upper = FALSE)

# Define functions for Mantel tests
# get_mantel_result <- function(x) {
#   result <- mantel(L_fst_dist, dist(x), permutations = 10e3)
#   return(tibble(statistic = result$statistic, signif = result$signif))
# }
# 
# get_mantel_result_geoDist <- function(x) {
#   result <- mantel(L_fst_dist, geo_dist, permutations = 10e3)
#   return(tibble(statistic = result$statistic, signif = result$signif, env_var = "geo_dist"))
# }
# 
# get_mantel.partial_result <- function(x) {
#   result <- mantel.partial(L_fst_dist, dist(x), geo_dist, permutations = 10e3, parallel = 8)
#   return(tibble(statistic = result$statistic, signif = result$signif))
# }

# Perform Mantel tests and save results
# mantel_results <- bind_rows(lapply(clim_input[4:40], get_mantel_result)) %>%
#   mutate(env_var = colnames(clim_input)[4:40]) %>%
#   bind_rows(get_mantel_result_geoDist())
# 
# mantel_results_sig <- mantel_results %>%
#   mutate(sig_color = case_when(signif >= 0.05 ~ "p >= 0.05",
#                                signif < 0.05 & signif >= 0.001 ~ "p < 0.05",
#                                signif < 0.001 ~ "p < 0.001"))
# 
# write_excel_csv(mantel_results_sig, "mantel_results_sig_20230621.csv")
# mantel_results_sig <- read_csv("mantel_results_sig_20230621.csv")

# Perform partial Mantel tests and save results
# mantel.partial_results <- bind_rows(lapply(clim_input[4:40], get_mantel.partial_result)) %>%
#   mutate(env_var = colnames(clim_input)[4:40])
# 
# mantel.partial_results_sig <- mantel.partial_results %>%
#   mutate(sig_color = case_when(signif >= 0.05 ~ "p >= 0.05",
#                                signif < 0.05 & signif >= 0.001 ~ "p < 0.05",
#                                signif < 0.001 ~ "p < 0.001"))
# 
# write_excel_csv(mantel.partial_results_sig, "mantel.partial_results_sig_20230621.csv")
# mantel.partial_results_sig <- read_csv("mantel.partial_results_sig_20230621.csv")


# Plot geographic distance vs. FST
p <- ggplot(data = geo_clim_dist, aes(y = L_FST, x = geo_dist / 1000)) +
  geom_point(alpha = 0.8, size = 0.5) +
  stat_smooth(method = "lm", color = "red") +
  geom_text(x = 0, y = 0.14, label = "Mantel r = 0.21\np = 0.023", hjust = 0, vjust = 1, check_overlap = TRUE) +
  xlab("Geographic distance (km)") +
  ylab(expression(F[ST] / (1 - F[ST]))) +
  theme_bw() +
  theme(panel.grid = element_blank())
p
ggsave("./figures/figure_4_A.png", width = 4, height = 4, dpi = 900)

# Plot growing degree days vs. FST
p <- ggplot(data = geo_clim_dist, aes(y = L_FST, x = growingDegDays5)) +
  geom_point(alpha = 0.8, size = 0.5) +
  stat_smooth(method = "lm", color = "red") +
  geom_text(x = 0, y = 0.14, label = "Mantel r = 0.78\np < 0.001", hjust = 0, vjust = 1, check_overlap = TRUE) +
  xlab(expression("Growing Degree Days (>5"*degree*C*") (standardized)")) +
  ylab(expression(F[ST] / (1 - F[ST]))) +
  theme_bw() +
  theme(panel.grid = element_blank())
p
ggsave("./figures/figure_4_B.png", width = 4, height = 4, dpi = 900)
