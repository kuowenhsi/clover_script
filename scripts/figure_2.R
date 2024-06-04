## LD decay plotting script
library(tidyverse)
library(cowplot)
library(ggpubr)
library(ggsignif)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/Molecular_ecology/clover_script")


# read in data
ld_bins <- read_tsv("./data/figure_2_bwa_415_20230527_100kb_0.8_passed_LD.ld_decay_bins.tsv")

Ac_snps <- read_tsv("./data/figure_2_bwa_415_20230527_GWAS_PC10.Ac_CNV.glm.linear")%>%
  filter(LOG10_P > -log10(5e-8), `#CHROM` == "chr_02")

Li_snps <- read_tsv("./data/figure_2_bwa_415_20230527_GWAS_PC10.Li_CNV.glm.linear")%>%
  filter(LOG10_P > -log10(5e-8), `#CHROM` == "chr_12")

AcLi_snp_ID <- c(Ac_snps$ID, Li_snps$ID)

AcLi_snp_ID_LD <- read_table("./data/figure_2_bwa_415_20230527_AcLi_snp_ID_extraction_bed_LD.ld")%>%
  filter(!(CHR_A == CHR_B))%>%
  mutate(name = "Ac-Li")

AcLi_snp_excluded_LD <- read_table("./data/figure_2_bwa_415_20230527_AcLi_snp_ID_exclusion_thin1000_bed_LD.ld")%>%
  filter(!(CHR_A == CHR_B))%>%
  mutate(name = "inter-chrom.")%>%
  sample_n(3000)



AcLi_snp_combined <- bind_rows(AcLi_snp_ID_LD, AcLi_snp_excluded_LD)%>%
  mutate(name = factor(name, levels = c("inter-chrom.", "Ac-Li")))

ks.test(x = AcLi_snp_ID_LD$R2, y = AcLi_snp_excluded_LD$R2)


p1 <- ggplot(data = ld_bins, aes(x = distance, avg_R2))+
  geom_point(size = 0.5, color = "grey65")+
  geom_vline(xintercept = 100e3, color = "blue")+
  geom_vline(xintercept = 25e3, color = "red")+
  stat_smooth(color = "red", method = "loess")+
  scale_y_continuous(limits = c(0,0.15))+
  xlab("Intra-chromosomal distance (bp)") + ylab(expression(italic(r)^2))+
  theme_bw()


p2 <- ggplot(data = AcLi_snp_combined, aes(x = name, y = R2))+
  geom_point(position = position_jitter(width = 0.2), aes(color = name), size = 0.5)+
  geom_boxplot(fill = "white", width = 0.2, alpha = 0.7, outlier.shape = NA)+
  scale_x_discrete(labels = c("Inter-chrom.", "Ac-Li"))+
  scale_y_continuous(limits = c(0,0.15))+
  scale_color_manual(values = c("gray80", "red"))+
  geom_signif(comparisons = list(c("inter-chrom.", "Ac-Li")), y_position = 0.12, test = "ks.test", map_signif_level = TRUE, tip_length = 0.01)+
  theme_bw()+
  theme(axis.title = element_blank(), legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(t = 30))


p <- plot_grid(p1, p2, align = "h", rel_widths = c(5,2), nrow = 1, labels = c("(A)","(B)"))

ggsave("./figures/figure_2.png", width = 6, height = 4, dpi = 900)


