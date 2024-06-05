# Load necessary libraries
library(tidyverse)
library(FactoMineR)
library(ungeviz)
library(cowplot)
library(ggpubr)
library(data.table)
library(grid)

# Set working directory
setwd("./clover_script")

# Define function to get p-value from regression
get_regression_p <- function(x, y) {
  if (max(y) == 1) {
    out <- glm(y ~ x, family = "binomial")
    return(-log10(summary(out)$coefficients[2, 4]))
  } else {
    out <- lm(y ~ x)
    return(-log10(summary(out)$coefficients[2, 4]))
  }
}

# Read and process input data
data_input <- read_csv("./data/figure_4_cyano_buf_limei_climate_data_415.csv")

# Clean column names
colnames(data_input) <- colnames(data_input) %>% str_remove("current_30arcsec_|wc2.1_30s_")
colnames(data_input)

# Perform PCA on selected columns
cyano_data_pca <- data_input %>%
  select(10:44)

cyano_data_pca <- as_tibble(PCA(cyano_data_pca)$ind$coord) %>%
  select(PC1 = Dim.1, PC2 = Dim.2)

# Rename columns for clarity
colnames(data_input)[26:44] <- c(
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

# Combine climate data with PCA results and reshape
cyano_data <- data_input %>%
  select(-c(45:65)) %>%
  bind_cols(cyano_data_pca) %>%
  pivot_longer(cols = 10:46, names_to = "predictor", values_to = "values") %>%
  arrange(predictor) %>%
  group_by(predictor) %>%
  summarise(Ac_p = get_regression_p(values, Ac),
            Li_p = get_regression_p(values, Li),
            Ac_CNV_p = get_regression_p(values, Ac_CNV),
            Li_CNV_p = get_regression_p(values, Li_CNV),
            CNH_p = get_regression_p(values, CNH)) %>%
  pivot_longer(cols = 2:6, names_to = "phenotype", values_to = "p_value") %>%
  filter(phenotype %in% c("Ac_CNV_p", "Li_CNV_p", "CNH_p"))

colnames(cyano_data)

# Create contingency table and perform Chi-squared test
chi_result <- chisq.test(table(data_input$Ac, data_input$Li))
chi_result
chi_result$observed
chi_result$expected

# Process data for Chi-squared test visualization
cyano_data_count <- data_input %>%
  select(1:9) %>%
  mutate(AcLi = case_when((Ac == 1) & (Li == 1) ~ 1, TRUE ~ 0),
         acLi = case_when((Ac == 0) & (Li == 1) ~ 1, TRUE ~ 0),
         Acli = case_when((Ac == 1) & (Li == 0) ~ 1, TRUE ~ 0),
         acli = case_when((Ac == 0) & (Li == 0) ~ 1, TRUE ~ 0)) %>%
  select(1, 10:13) %>%
  pivot_longer(cols = 2:5, names_to = "cyanotype", values_to = "value") %>%
  group_by(cyanotype) %>%
  summarise(observed = sum(value)) %>%
  mutate(expected = c(chi_result$expected[2, 2], chi_result$expected[2, 1], chi_result$expected[1, 2], chi_result$expected[1, 1])) %>%
  mutate(cyanotype = factor(cyanotype, levels = c("AcLi", "acLi", "Acli", "acli")))

# Plot Chi-squared test results
p2 <- ggplot(data = cyano_data_count, aes(x = cyanotype)) +
  geom_col(aes(y = observed, fill = cyanotype)) +
  geom_hpline(aes(y = expected), size = 0.1) +
  geom_text(x = 0.5, y = 172, label = "Pearson's Chi-squared test", hjust = 0, size = 3) +
  geom_text(x = 0.5, y = 160, label = "X-squared = 18.29, df = 1, p-value = 1.897e-05", hjust = 0, size = 3) +
  geom_text(aes(y = expected + 6), label = "exp.", size = 3) +
  scale_fill_manual(values = c("cyan4", "pink", "skyblue", "orange")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_y_continuous(name = "Observed / Expected cyanotype counts", expand = c(0, 0, 0, 35)) +
  scale_x_discrete(name = "")

# Summarize data for plotting regression results
data_input_s <- data_input %>%
  select(1:9, growingDegDays5) %>%
  group_by(Pop) %>%
  summarize(Ac_CNV_mean = mean(Ac_CNV), Li_CNV_mean = mean(Li_CNV), growingDegDays5 = mean(growingDegDays5),
            Ac_CNV_se = sd(Ac_CNV) / sqrt(n()), Li_CNV_se = sd(Li_CNV) / sqrt(n()))

# Plot regression results for Ac CNV
p3 <- ggplot(data = data_input_s, aes(x = growingDegDays5, y = Ac_CNV_mean)) +
  geom_errorbar(aes(ymax = Ac_CNV_mean + Ac_CNV_se, ymin = Ac_CNV_mean - Ac_CNV_se), color = "gray85") +
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(size = 3, label.x.npc = 0, label.y.npc = 0.95) +
  xlab(expression("Growing Degree Days (>5" * degree * C * ")")) +
  ylab(expression(italic("Ac") ~ "CNV")) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Plot regression results for Li CNV
p4 <- ggplot(data = data_input_s, aes(x = growingDegDays5, y = Li_CNV_mean)) +
  geom_errorbar(aes(ymax = Li_CNV_mean + Li_CNV_se, ymin = Li_CNV_mean - Li_CNV_se), color = "gray85") +
  geom_point() +
  stat_smooth(method = "lm") +
  stat_cor(size = 3, label.x.npc = 0, label.y.npc = 0.95) +
  xlab(expression("Growing Degree Days (>5" * degree * C * ")")) +
  ylab(expression(italic("Li") ~ "CNV")) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Combine plots into a single figure
combined_p <- plot_grid(p3, p4, p2, rel_widths = c(3.5, 3.5, 3), nrow = 1)

# Save the combined figure
ggsave("./figures/figure_5_A_B_C.png", width = 12, height = 3.5, dpi = 900)

###########################################################

# Load chromosome information

chr_len_temp <- read_tsv("./data/figure_5_HiFi_HiC_LM_combined_v_1.5.fasta.fai", col_names = c("chr", "end")) %>%
  select(1:2)%>%
  mutate(start = 1)%>%
  filter(!str_detect(chr, "drTriRepe4Chr0c"))%>%
  mutate(chr = as.integer(str_remove(chr, "drTriRepe4Chr")))%>%
  arrange(chr) %>%
  mutate(chr = paste0("chr_", str_pad(chr, 2, pad = "0")))%>%
  mutate(lag_pos = lag(end, default = 0))%>%
  mutate(pos_pad = cumsum(lag_pos))%>%
  mutate(padded_start = start + pos_pad - 1, padded_end = end + pos_pad)%>%
  mutate(padded_chr_pos = (padded_start + padded_end)/2)

#################
# plot GWAS results


# Load GWAS results
GWAS_result_list <- list()

GWAS_result_list[[1]] <- read_tsv("./data/figure_5_bwa_415_20230527_GWAS_PC10.Ac_CNV.glm.linear") %>% mutate(trait = "Ac_CNV")

GWAS_result_list[[2]] <- read_tsv("./data/figure_5_bwa_415_20230527_GWAS_PC10.Li_CNV.glm.linear") %>% mutate(trait = "Li_CNV")

GWAS_result_list <- bind_rows(GWAS_result_list)

GWAS_result_joined <- GWAS_result_list %>%
  filter(TEST == "ADD")%>%
  left_join(chr_len_temp, by = c("#CHROM" = "chr"))%>%
  mutate(padded_pos = POS + pos_pad)%>%
  mutate(dot_color = case_when(LOG10_P > -log10(5e-8) ~ "black", TRUE ~ "gray80"))

chr_len_temp_complete <- chr_len_temp%>%
  mutate(trait = factor("Ac", levels = unique(GWAS_result_joined$trait)))%>%
  expand(trait, chr)%>%
  left_join(chr_len_temp, by = "chr")

p_GWAS_Ac <- ggplot(data = filter(GWAS_result_joined, trait == "Ac_CNV"), aes(x = padded_pos, y = LOG10_P))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = chr), inherit.aes = FALSE)+
  geom_point(aes(color = dot_color), size = 0.5)+
  geom_hline(yintercept = -log10(5e-8), color = "red")+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("Chromosome", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16, oob = scales::squish)+
  scale_y_continuous(expression("-"*log[10]*"(p value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(fill = "none", color = "none")


p_QQ_Ac <- ggplot(data = filter(GWAS_result_list, trait == "Ac_CNV"), aes(sample = -log(10^(-LOG10_P))))+
  stat_qq(aes(x = after_stat(theoretical)/log(10), y = after_stat(sample)/log(10)), distribution = qexp, size = 0.5)+
  stat_qq_line(aes(x = after_stat(x)/log(10), y = after_stat(y)/log(10)), distribution = qexp)+
  geom_hline(yintercept = -log10(5e-8), color = "red")+
  xlab(expression("Theoretical"~"-"*log[10]*"(p)"))+
  ylab(expression("Observed"~"-"*log[10]*"(p)"))+
  theme_bw()+
  theme(panel.grid = element_blank())

combined_Ac <- plot_grid(p_GWAS_Ac, p_QQ_Ac, rel_widths = c(8, 2), align = "h", axis = "l")

png(filename = "./figures/figure_5_D.png", width = 10, height = 2, units = "in", res = 900)
combined_Ac
grid.draw(segmentsGrob(x0 = unit(0.101, "npc"), x1 = unit(0.101, "npc"), y0 = unit(0.13, "npc"), y1 = unit(0.23, "npc"), arrow = arrow(type = "closed",length = unit(0.05, "inches")), gp = gpar(col = "red", lwd = 2, fill = "red")))
grid.draw(textGrob(x = unit(0.101, "npc"),y = unit(0.09, "npc"), label = expression(italic("Ac")), hjust = 0.5, vjust = 0.5, gp = gpar(col = "red", size = 2)))
dev.off()


#############


p_GWAS_Li <- ggplot(data = filter(GWAS_result_joined, trait == "Li_CNV"), aes(x = padded_pos, y = LOG10_P))+
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = chr), inherit.aes = FALSE)+
  geom_point(aes(color = dot_color), size = 0.5)+
  geom_hline(yintercept = -log10(5e-8), color = "red")+
  scale_fill_manual(name = "",values = rep(c("white", "gray95"), 8))+
  scale_x_continuous("Chromosome", expand = c(0,0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16)+
  scale_y_continuous(expression("-"*log[10]*"(p value)"), expand = c(0,0,0.1,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(), plot.margin = margin(t = 0.1, r = 0.1, l = 0.1,b = 0.22, unit = "in"))+
  guides(fill = "none", color = "none")

p_QQ_Li <- ggplot(data = filter(GWAS_result_list, trait == "Li_CNV"), aes(sample = -log(10^(-LOG10_P))))+
  stat_qq(aes(x = after_stat(theoretical)/log(10), y = after_stat(sample)/log(10)), distribution = qexp, size = 0.5)+
  stat_qq_line(aes(x = after_stat(x)/log(10), y = after_stat(y)/log(10)), distribution = qexp)+
  geom_hline(yintercept = -log10(5e-8), color = "red")+
  xlab(expression("Theoretical"~"-"*log[10]*"(p)"))+
  ylab(expression("Observed"~"-"*log[10]*"(p)"))+
  theme_bw()+
  theme(panel.grid = element_blank(), plot.margin = margin(t = 0.1, r = 0.1, l = 0.1,b = 0.22, unit = "in"))

combined_Li <- plot_grid(p_GWAS_Li, p_QQ_Li, rel_widths = c(8, 2), align = "h", axis = "l")

png(filename = "./figures/figure_5_E.png", width = 10, height = 2.2, units = "in", res = 900)
combined_Li
grid.draw(segmentsGrob(x0 = unit(0.586, "npc"), x1 = unit(0.586, "npc"), y0 = unit(0.1, "npc"), y1 = unit(0.2, "npc"), arrow = arrow(type = "closed",length = unit(0.05, "inches")), gp = gpar(col = "red", lwd = 2, fill = "red")))
grid.draw(textGrob(x = unit(0.586, "npc"),y = unit(0.06, "npc"), label = expression(italic("Li")), hjust = 0.5, vjust = 0.5, gp = gpar(col = "red", size = 2)))
dev.off()
