library(tidyverse)
library(FactoMineR)
library(ungeviz)
library(cowplot)
library(ggpubr)

setwd("/Users/kuowenhsi/OneDrive - Washington University in St. Louis/GBS_wild_population/Figure_5")

get_regression_p <- function(x, y){
  if (max(y) == 1){
    out <- glm(y ~ x, family = "binomial")
    return(-log10(summary(out)$coefficients[2,4]))
  } else{
    out <- lm(y ~ x)
    return(-log10(summary(out)$coefficients[2,4]))
  }
}

data_input <- read_csv("cyano_buf_limei_climate_data_415.csv")

colnames(data_input) <- colnames(data_input)%>% str_remove("current_30arcsec_|wc2.1_30s_")
colnames(data_input)

cyano_data_pca <- data_input%>%
  select(10:44)

cyano_data_pca <- as_tibble(PCA(cyano_data_pca)$ind$coord)%>%
  select(PC1 = Dim.1, PC2 = Dim.2)

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

cyano_data <- data_input%>%
  select(-c(45:65))%>%
  bind_cols(cyano_data_pca)%>%
  pivot_longer(cols = 10:46, names_to = "predictor", values_to = "values")%>%
  arrange(predictor)%>%
  group_by(predictor)%>%
  summarise(Ac_p = get_regression_p(values, Ac),
            Li_p = get_regression_p(values, Li),
            Ac_CNV_p = get_regression_p(values, Ac_CNV),
            Li_CNV_p = get_regression_p(values, Li_CNV),
            CNH_p = get_regression_p(values, CNH))%>%
  pivot_longer(cols = 2:6, names_to = "phenotype", values_to = "p_value")%>%
  filter(phenotype %in% c("Ac_CNV_p", "Li_CNV_p", "CNH_p"))

colnames(cyano_data)

p1 <- ggplot(data = cyano_data, aes(y = reorder(predictor, p_value),x = p_value))+
  geom_point(aes(shape = phenotype, color = phenotype))+
  geom_point(data = filter(cyano_data, p_value < -log10(5e-8)), color = "gray85")+
  geom_vline(xintercept = -log10(5e-8), color = "red")+
  scale_y_discrete("")+
  scale_x_continuous(name = expression("-"*log[10]*"(p value)"), expand = c(0,0,0,10))+
  scale_color_discrete("",limits = c("Ac_CNV_p", "Li_CNV_p", "CNH_p"), labels = c("Ac CNV (numeric)","Li CNV (numeric)", "Ac + Li (binary)"))+
  scale_shape_discrete("",limits = c("Ac_CNV_p", "Li_CNV_p", "CNH_p"),labels = c("Ac CNV (numeric)", "Li CNV (numeric)", "Ac + Li (binary)"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = c(0.8,0.1), legend.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray97"))
p1

ggsave("cyanotype_regression_20230613.png", width = 8, height = 10)

cyano_data_new <- data_input

# https://sites.radford.edu/~rsheehy/Gen_flash/Tutorials/Chi-Square_tutorial/x2-tut.htm

table(cyano_data_new$Ac, cyano_data_new$Li)

chi_result <- chisq.test(table(cyano_data_new$Ac, cyano_data_new$Li))
chi_result$observed
chi_result$expected


cyano_data_count <- data_input%>%
  select(1:9)%>%
  mutate(AcLi = case_when((Ac == 1) & (Li == 1) ~ 1, TRUE ~ 0),
         acLi = case_when((Ac == 0) & (Li == 1) ~ 1, TRUE ~ 0),
         Acli = case_when((Ac == 1) & (Li == 0) ~ 1, TRUE ~ 0),
         acli = case_when((Ac == 0) & (Li == 0) ~ 1, TRUE ~ 0))%>%
  select(1, 10:13)%>%
  pivot_longer(cols = 2:5, names_to = "cyanotype", values_to = "value")%>%
  group_by(cyanotype)%>%
  summarise(observed = sum(value))%>%
  mutate(expected = c(chi_result$expected[2,2], chi_result$expected[2,1], chi_result$expected[1,2], chi_result$expected[1,1]))%>%
  mutate(cyanotype = factor(cyanotype, levels = c("AcLi", "acLi", "Acli", "acli")))



p2 <- ggplot(data = cyano_data_count, aes(x = cyanotype))+
  geom_col(aes(y = observed, fill = cyanotype))+
  geom_hpline(aes(y = expected), size = 0.1)+
  geom_text(x = 0.5, y = 172, label = "Pearson's Chi-squared test", hjust = 0, size = 3)+
  geom_text(x = 0.5, y = 160, label = "X-squared = 18.29, df = 1, p-value = 1.897e-05", hjust = 0, size = 3)+
  geom_text(aes(y = expected + 6), label = "exp.", size = 3)+
  scale_fill_manual(values = c("cyan4","pink","skyblue","orange"))+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_y_continuous(name = "Observed / Expected cyanotype counts", expand = c(0,0,0,35))+
  scale_x_discrete(name = "")
p2

data_input_s <- data_input %>%
  select(1:9, growingDegDays5)%>%
  group_by(Pop)%>%
  summarize(Ac_CNV_mean = mean(Ac_CNV), Li_CNV_mean = mean(Li_CNV), growingDegDays5 = mean(growingDegDays5),
            Ac_CNV_se = sd(Ac_CNV)/sqrt(n()), Li_CNV_se = sd(Li_CNV)/sqrt(n()))

p3 <- ggplot(data = data_input_s, aes(x = growingDegDays5, y = Ac_CNV_mean))+
  geom_errorbar(aes(ymax = Ac_CNV_mean + Ac_CNV_se, ymin = Ac_CNV_mean - Ac_CNV_se), color = "gray85")+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_cor(size = 3, label.x.npc =0, label.y.npc = 0.95)+
  xlab(expression("Growing Degree Days (>5"*degree*C*")"))+
  ylab(expression(italic("Ac") ~ "CNV"))+
  theme_bw()+
  theme(panel.grid = element_blank())


p4 <- ggplot(data = data_input_s, aes(x = growingDegDays5, y = Li_CNV_mean))+
  geom_errorbar(aes(ymax = Li_CNV_mean + Li_CNV_se, ymin = Li_CNV_mean - Li_CNV_se), color = "gray85")+
  geom_point()+
  stat_smooth(method = "lm")+
  stat_cor(size = 3, label.x.npc =0, label.y.npc = 0.95)+
  xlab(expression("Growing Degree Days (>5"*degree*C*")"))+
  ylab(expression(italic("Li") ~ "CNV"))+
  theme_bw()+
  theme(panel.grid = element_blank())

p4
p3

combined_p <- plot_grid(p3, p4, p2, rel_widths = c(3.5, 3.5, 3), nrow = 1)

ggsave("lm_count_top_20240529.png", width = 12, height = 3.5)

##############################################


data_input <- read_csv("cyano_buf_limei_climate_data_415.csv")%>%
  filter(Pop != "BWA" & Pop != "VBC")


cyano_data_pca <- data_input%>%
  select(-c(10:44))%>%
  select(10:28)

cyano_data_pca <- as_tibble(PCA(cyano_data_pca)$ind$coord)%>%
  select(PC1 = Dim.1, PC2 = Dim.2)


cyano_data <- data_input%>%
  select(-c(10:44))%>%
  bind_cols(cyano_data_pca)%>%
  pivot_longer(cols = 10:32, names_to = "predictor", values_to = "values")%>%
  arrange(predictor)%>%
  group_by(predictor)%>%
  summarise(Ac_p = get_regression_p(values, Ac),
            Li_p = get_regression_p(values, Li),
            Ac_CNV_p = get_regression_p(values, Ac_CNV),
            Li_CNV_p = get_regression_p(values, Li_CNV),
            CNH_p = get_regression_p(values, CNH))%>%
  pivot_longer(cols = 2:6, names_to = "phenotype", values_to = "p_value")%>%
  filter(phenotype %in% c("Ac_CNV_p", "Li_CNV_p", "CNH_p"))%>%
  mutate(predictor = factor(predictor, levels = c("PC1", "PC2", paste0("BIO", 1:19), "AI", "Apet")))

colnames(cyano_data)

p1 <- ggplot(data = cyano_data, aes(x = predictor, p_value))+
  geom_point(aes(shape = phenotype, color = phenotype))+
  geom_point(data = filter(cyano_data, p_value < -log10(5e-8)), color = "gray85")+
  geom_hline(yintercept = -log10(5e-8), color = "red")+
  scale_x_discrete("")+
  scale_y_continuous(name = expression("-"*log[10]*"(p value)"), expand = c(0,0,0,10))+
  scale_color_discrete("",limits = c("Ac_CNV_p", "Li_CNV_p", "CNH_p"), labels = c("Ac CNV (numeric)","Li CNV (numeric)", "Ac + Li (binary)"))+
  scale_shape_discrete("",limits = c("Ac_CNV_p", "Li_CNV_p", "CNH_p"),labels = c("Ac CNV (numeric)", "Li CNV (numeric)", "Ac + Li (binary)"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = c(0.9,0.9), legend.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray97"))
p1

cyano_data_new <- data_input

# https://sites.radford.edu/~rsheehy/Gen_flash/Tutorials/Chi-Square_tutorial/x2-tut.htm

table(cyano_data_new$Ac, cyano_data_new$Li)

chi_result <- chisq.test(table(cyano_data_new$Ac, cyano_data_new$Li))
chi_result$observed
chi_result$expected


cyano_data_count <- data_input%>%
  select(1:9)%>%
  mutate(AcLi = case_when((Ac == 1) & (Li == 1) ~ 1, TRUE ~ 0),
         acLi = case_when((Ac == 0) & (Li == 1) ~ 1, TRUE ~ 0),
         Acli = case_when((Ac == 1) & (Li == 0) ~ 1, TRUE ~ 0),
         acli = case_when((Ac == 0) & (Li == 0) ~ 1, TRUE ~ 0))%>%
  select(1, 10:13)%>%
  pivot_longer(cols = 2:5, names_to = "cyanotype", values_to = "value")%>%
  group_by(cyanotype)%>%
  summarise(observed = sum(value))%>%
  mutate(expected = c(chi_result$expected[2,2], chi_result$expected[2,1], chi_result$expected[1,2], chi_result$expected[1,1]))%>%
  mutate(cyanotype = factor(cyanotype, levels = c("AcLi", "acLi", "Acli", "acli")))



p2 <- ggplot(data = cyano_data_count, aes(x = cyanotype))+
  geom_col(aes(y = observed, fill = cyanotype))+
  geom_hpline(aes(y = expected), size = 0.1)+
  geom_text(x = 0.5, y = 172, label = "Pearson's Chi-squared test", hjust = 0, size = 3)+
  geom_text(x = 0.5, y = 160, label = "X-squared = 18.29, df = 1, p-value = 1.897e-05", hjust = 0, size = 3)+
  geom_text(aes(y = expected + 6), label = "exp.", size = 3)+
  scale_fill_manual(values = c("cyan4","pink","skyblue","orange"))+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = "none")+
  scale_y_continuous(name = "Observed / Expected cyanotype counts", expand = c(0,0,0,35))+
  scale_x_discrete(name = "")
p2


combined_p <- plot_grid(p1, p2, rel_widths = c(8, 3))

ggsave("lm_count_top_noBWA_noVBC.png", width = 12, height = 3.5)

