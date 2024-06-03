# Load necessary libraries
library(tidyverse)     # For data manipulation and visualization
library(cowplot)       # For combining multiple plots
library(sf)            # For handling spatial data
library(rnaturalearth) # For accessing natural earth data
library(rnaturalearthdata)
library(scatterpie)    # For creating scatter pie plots
library(raster)        # For raster data manipulation
library(ggspatial)     # For spatial data visualization

# Set working directory
setwd("./clover_script")

# Read input data
pca_data_mean <- read_csv("./data/figure_1_cyanotypes.csv")
pca_data_mean_pie_3857 <- read_csv("./data/figure_1_plotting.csv")

# Download this file from https://envirem.github.io/, which is not available on GitHub
env_layer <- raster("./data/current_30arcsec_growingDegDays5.tif")

# Get and crop country boundaries
usa_state <- ne_states(country = "United States of America", returnclass = "sf") %>%
  st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60) %>%
  st_geometry()

canada_state <- ne_states(country = "canada", returnclass = "sf") %>%
  st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60) %>%
  st_geometry()

mexico_state <- ne_states(country = "mexico", returnclass = "sf") %>%
  st_crop(xmin = -135, ymin = -55, xmax = -65, ymax = 60) %>%
  st_geometry()

# Process environmental layer
env_layer_croped <- crop(env_layer, extent(-135, -65, 25, 50)) / 10
env_layer_df <- as.data.frame(env_layer_croped, xy = TRUE) %>%
  drop_na()

# Create PCA plot with scatter pies
p <- ggplot(data = st_as_sf(pca_data_mean, coords = c("Longitude", "Latitude"), agr = "constant", crs = 4326)) +
  geom_sf(data = usa_state, fill = NA, color = "gray75") +
  geom_sf(data = canada_state, fill = NA, color = "gray75") +
  geom_sf(data = mexico_state, fill = NA, color = "gray75") +
  geom_sf(aes(color = Region, shape = Region), size = 2) +
  geom_scatterpie(data = pca_data_mean_pie_3857, aes(x = X, y = Y, r = total_n * 1.2e4), cols = c("AcLi", "acLi", "Acli", "acli"), color = NA) +
  scale_fill_manual("", values = c("cyan4", "pink", "skyblue", "orange")) +
  geom_text(data = pca_data_mean_pie_3857, aes(x = X + 1.5e5, y = Y + 1.4e5, label = Pop), size = 3, show.legend = FALSE) +
  scale_color_discrete("") +
  scale_shape_manual("", values = rep(c(15, 17, 19), 3)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical", panel.grid = element_blank()) +
  coord_sf(ylim = c(2e6, 7e6), expand = FALSE, crs = 3857) +
  annotation_scale(location = "br", width_hint = 0.25)

# Extract legend from the PCA plot
p_legend <- get_legend(p + guides(color = "none", shape = "none") + theme(legend.background = element_blank()))

# Create environmental layer plot
p_env <- ggplot() +
  geom_raster(data = env_layer_df, aes(x = x, y = y, fill = current_30arcsec_growingDegDays5)) +
  geom_sf(data = usa_state, fill = NA) +
  geom_sf(data = canada_state, fill = NA) +
  geom_sf(data = mexico_state, fill = NA) +
  scale_fill_viridis_c(option = "inferno", name = "") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.key.height = unit(0.1, "in"), legend.key.width = unit(0.6, "in"), legend.margin = margin(t = -0.02, b = -0.02, unit = "in"), legend.position = "top", panel.background = element_rect(fill = "white")) +
  xlab(expression("Growing Degree Days (>5" * degree * C * ")")) +
  ylab("") +
  coord_sf(xlim = c(-130, -65), ylim = c(25, 50), expand = FALSE, label_axes = "")

# Combine PCA plot and environmental plot
p_comb <- p + guides(fill = "none") +
  annotation_custom(p_legend, xmin = -20e6, ymax = 3e6) +
  annotation_custom(ggplotGrob(p_env), xmin = -16e6, ymin = 2.75e6, xmax = -11e6, ymax = 4.75e6)

# Save the combined plot
ggsave("./figures/figure_1.png", width = 10, height = 8, dpi = 600)

