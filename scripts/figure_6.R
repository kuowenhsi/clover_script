library(tidyverse)
library(data.table)
library(qvalue)
library(cowplot)

# Set working directory
setwd("./clover_script")

################

# Load chromosome information from a .fasta index file
chr_len_temp <- read_tsv("./data/figure_5_HiFi_HiC_LM_combined_v_1.5.fasta.fai", col_names = c("chr", "end")) %>%
  select(1:2) %>%
  mutate(start = 1) %>%
  filter(!str_detect(chr, "drTriRepe4Chr0c")) %>%
  mutate(chr = as.integer(str_remove(chr, "drTriRepe4Chr"))) %>%
  arrange(chr) %>%
  mutate(lag_pos = lag(end, default = 0)) %>%
  mutate(pos_pad = cumsum(lag_pos)) %>%
  mutate(padded_start = start + pos_pad - 1, padded_end = end + pos_pad) %>%
  mutate(padded_chr_pos = (padded_start + padded_end) / 2)

#################

# Load QTL mapping result
one_drop_data <- list()

# Read and process QTL mapping result files
one_drop_data[[1]] <- read.csv("./data/figure_6_Drop1_DG_pop_at_DMN_site.csv") %>%
  mutate(dataset = "DG_pop_DMN_site")

one_drop_data[[2]] <- read.csv("./data/figure_6_Drop1_DG_pop_at_GFL_site.csv") %>%
  mutate(dataset = "DG_pop_GFL_site")

one_drop_data[[3]] <- read.csv("./data/figure_6_Drop1_GS_pop_at_GFL_site.csv") %>%
  mutate(dataset = "GS_pop_GFL_site")

one_drop_data[[4]] <- read.csv("./data/figure_6_Drop1_GS_pop_at_STL_site.csv") %>%
  mutate(dataset = "GS_pop_STL_site")

# Combine all datasets into one
one_drop_data <- bind_rows(one_drop_data)

# Process combined QTL mapping result
one_drop_data_pRange <- one_drop_data %>%
  mutate(lowmarker_pos = as.integer(str_remove(lowmarker, "chr_..:")), highmarker_pos = as.integer(str_remove(highmarker, "chr_..:"))) %>%
  left_join(chr_len_temp, by = "chr") %>%
  select(dataset, pheno, chr, maxLod, lowmarker_pos, highmarker_pos, pos_pad) %>%
  mutate(padded_lowmarker_pos = lowmarker_pos + pos_pad, padded_highmarker_pos = highmarker_pos + pos_pad, ID = seq(nrow(.))) %>%
  filter(str_detect(pheno, "sqrt")) %>%
  mutate(pheno_group = case_when(
    str_detect(pheno, "Veg") ~ "Vegetative",
    str_detect(pheno, "Flo") ~ "Reproductive",
    str_detect(pheno, "Survival|Life") ~ "Survival",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(chr = factor(chr, levels = 1:16), name = "QTL mapping") %>%
  filter(maxLod > 4)

# Plot QTL mapping results
p_qtl <- ggplot(data = one_drop_data_pRange, aes(xmin = padded_lowmarker_pos, xmax = padded_highmarker_pos, y = pheno_group)) +
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE) +
  geom_linerange(aes(color = pheno_group), linewidth = 3) +
  scale_color_manual(values = c("#D62728", "#FF7F0E", "#2CA02C")) +
  scale_fill_manual(name = "", values = rep(c("white", "gray95"), 8)) +
  scale_x_continuous(expand = c(0, 0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16) +
  annotate("segment", x=124910454, y=4.5, xend=124910454, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  annotate("segment", x=249995270, y=4.5, xend=249995270, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  annotate("segment", x=582870000, y=4.5, xend=582870000, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  annotate("segment", x=713443783, y=4.5, xend=713443783, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  annotate("segment", x=874139073, y=4.5, xend=874139073, yend=3.61, col="red", arrow=arrow(length=unit(0.2, "cm"), type = "closed"), lineend = "round", linewidth = 1) +
  ylab("") +
  coord_cartesian(ylim = c(1, 3), clip="off") +
  theme_bw() +
  theme(plot.margin = unit(c(2, 0.2, 0.2, 0.2), "lines")) +
  guides(fill = "none", color = "none")

p_qtl

##################

# Load GWAS results
Ac_snps <- read_tsv("./data/figure_5_bwa_415_20230527_GWAS_PC10.Ac_CNV.glm.linear") %>%
  filter(LOG10_P > -log10(5e-8), `#CHROM` == "chr_02")

Li_snps <- read_tsv("./data/figure_5_bwa_415_20230527_GWAS_PC10.Li_CNV.glm.linear") %>%
  filter(LOG10_P > -log10(5e-8), `#CHROM` == "chr_12")

AcLi_snp_ID <- c(Ac_snps$ID, Li_snps$ID)

###################

# Download LFMM_output_20230606.txt from Dryad https://doi.org/10.5061/dryad.s7h44j1fd
LFMM_PCA_data <- fread("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/LFMM_output_20230606.txt", header = FALSE)

unique(LFMM_PCA_data$V5)
unique(LFMM_PCA_data$V6)

LFMM_data_sub <- LFMM_PCA_data[(LFMM_PCA_data$V5 %in% unique(LFMM_PCA_data$V5)[7]) & LFMM_PCA_data$V6 == 3,] %>%
  left_join(chr_len_temp, by = c("V1" = "chr")) %>%
  group_by(V5) %>% # group by ENV
  mutate(padded_pos = V2 + pos_pad, qvalue = qvalue(V4)$qvalues) %>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))

LFMM_data_sub_q <- qvalue(LFMM_data_sub$V4)
hist(LFMM_data_sub_q)

simple_LM_GDD5 <- read_tsv("./data/figure_6_simple_lm_GDD5_20240516.tsv") %>%
  left_join(chr_len_temp, by = c("CHR" = "chr")) %>%
  mutate(padded_pos = POSITION + pos_pad, qvalue = qvalue(adjusted_p)$qvalues) %>%
  mutate(dot_color = case_when(qvalue < 0.05 ~ "black", TRUE ~ "gray80"))

simple_LM_GDD5_q <- qvalue(simple_LM_GDD5$adjusted_p)
plot(simple_LM_GDD5_q)
hist(simple_LM_GDD5_q)
hist(simple_LM_GDD5$qvalue)

nrow(filter(LFMM_data_sub, qvalue < 0.05))
LFMM_data_sub_sig <- filter(LFMM_data_sub, qvalue < 0.05)
simple_LM_GDD5_sig <- filter(simple_LM_GDD5, qvalue < 0.05)

# Plot result of simple linear model
ppop_lm <- ggplot(data = simple_LM_GDD5, aes(x = padded_pos, y = -log10(qvalue))) +
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE) +
  
  
  geom_point(color = simple_LM_GDD5$dot_color, size = 0.5) +
  geom_point(data = filter(simple_LM_GDD5, ID %in% AcLi_snp_ID), color = "red", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "", values = rep(c("white", "gray95"), 8)) +
  scale_x_continuous("", expand = c(0, 0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16) +
  scale_y_continuous(expression("-" * log[10] * "(adjusted p value)"), expand = c(0, 0, 0.1, 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")) +
  guides(fill = "none")

ppop_lm

###################

# Load and process PCadapt data
PCadapt_data <- read_tsv("./data/figure_6_PCadapt_output_20230613.txt") %>%
  left_join(chr_len_temp, by = c("CHR" = "chr")) %>%
  mutate(padded_pos = POSITION + pos_pad, qvalue = qvalue(pvalue)$qvalues) %>%
  mutate(dot_color = case_when(pvalue < 5e-8 ~ "black", TRUE ~ "gray80")) %>%
  mutate(name = "PCadapt")

# Q-Q plot for PCadapt result
p_QQ_PCadapt <- ggplot(data = PCadapt_data, aes(sample = -log(pvalue))) +
  stat_qq(aes(x = after_stat(theoretical) / log(10), y = after_stat(sample) / log(10)), distribution = qexp, size = 0.5) +
  stat_qq_line(aes(x = after_stat(x) / log(10), y = after_stat(y) / log(10)), distribution = qexp) +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  xlab(expression("Theoretical" ~ "-" * log[10] * "(p)")) +
  ylab(expression("Observed" ~ "-" * log[10] * "(p)")) +
  labs(title = "Q-Q plot for the PCadapt result") +
  theme_bw() +
  theme(panel.grid = element_blank())

p_QQ_PCadapt

# Load and process candidate genes data for LFMM and PCadapt
cand_genes_pad <- read_csv("./data/figure_6_cand_genes_lfmm_for_plot_20240109.csv") %>%
  left_join(chr_len_temp, by = "chr") %>%
  mutate(padded_pos = gene_pos + pos_pad, nudge_dist = case_when(HJUST == 0 ~ 50e5, HJUST == 1 ~ -50e5))

cand_genes_PCadapt_pad <- read_csv("./data/figure_6_cand_genes_PCadapt_for_plot_20240109.csv") %>%
  left_join(chr_len_temp, by = "chr") %>%
  mutate(padded_pos = gene_pos + pos_pad, nudge_dist = case_when(HJUST == 0 ~ 50e5, HJUST == 1 ~ -50e5))

cand_genes_pad <- cand_genes_pad %>%
  left_join(cand_genes_PCadapt_pad %>% select(chr, gene_name) %>% mutate(also_PCadapt = "YES"), by = c("chr", "gene_name")) %>%
  mutate(font_color = case_when(also_PCadapt == "YES" ~ "red", TRUE ~ "blue3"))

cand_genes_PCadapt_pad <- cand_genes_PCadapt_pad %>%
  left_join(cand_genes_pad %>% select(chr, gene_name) %>% mutate(also_lfmm = "YES"), by = c("chr", "gene_name")) %>%
  mutate(font_color = case_when(also_lfmm == "YES" ~ "red", TRUE ~ "blue3"))

##########

# Plot result of LFMM analysis
ppop <- ggplot(data = LFMM_data_sub, aes(x = padded_pos, y = -log10(qvalue))) +
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE) +
  geom_segment(data = cand_genes_pad, aes(x = padded_pos, xend = padded_pos, yend = Y), y = 0, linewidth = 1, alpha = 0.5, color = "skyblue") +
  geom_text(data = cand_genes_pad, aes(x = padded_pos + nudge_dist, y = Y, label = gene_name, hjust = HJUST), color = cand_genes_pad$font_color, size = 3) +
  geom_point(color = LFMM_data_sub$dot_color, size = 0.5) +
  geom_point(data = filter(LFMM_data_sub, V3 %in% AcLi_snp_ID), color = "red", size = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "", values = rep(c("white", "gray95"), 8)) +
  scale_x_continuous("", expand = c(0, 0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16) +
  scale_y_continuous(expression("-" * log[10] * "(q value)"), expand = c(0, 0, 0.1, 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")) +
  guides(fill = "none")

ppop

# Plot result of PCadapt analysis
pcap <- ggplot(data = PCadapt_data, aes(x = padded_pos, y = -log10(pvalue))) +
  geom_rect(data = chr_len_temp, aes(xmin = padded_start, xmax = padded_end, ymin = -Inf, ymax = Inf, fill = factor(chr, levels = 1:16)), inherit.aes = FALSE) +
  geom_segment(data = cand_genes_PCadapt_pad, aes(x = padded_pos, xend = padded_pos, yend = Y), y = 0, linewidth = 1, alpha = 0.5, color = "skyblue") +
  geom_text(data = cand_genes_PCadapt_pad, aes(x = padded_pos + nudge_dist, y = Y, label = gene_name, hjust = HJUST), color = cand_genes_PCadapt_pad$font_color, size = 3) +
  geom_point(color = PCadapt_data$dot_color, size = 0.5) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "", values = rep(c("white", "gray95"), 8)) +
  scale_x_continuous("", expand = c(0, 0), breaks = chr_len_temp$padded_chr_pos, labels = 1:16) +
  scale_y_continuous(expression("-" * log[10] * "(p value)"), expand = c(0, 0, 0.1, 0.1)) +
  theme_bw() +
  theme(panel.grid = element_blank(), plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines")) +
  guides(fill = "none")

pcap

# Combine all plots into a single figure
p_combined <- plot_grid(p_qtl, ppop, pcap, align = "v", ncol = 1, rel_heights = c(2, 4, 4), axis = 'l', labels = c("(A)", "(B)", "(C)"))

# Save the combined figure
ggsave("./figures/figure_6.png", width = 10, height = 10, dpi = 900)

