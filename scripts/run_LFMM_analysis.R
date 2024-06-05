# Load necessary libraries
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(qqman)
library(tidyverse)
library(snpStats)
library(FactoMineR)
library(parallel)

# Set working directory
setwd("./clover_script")

# Read genotype data
geno_data <- read.plink("./data/bwa_415_20230606.bed")
dim(geno_data$genotypes@.Data)

# Convert genotype data to a numeric matrix and then to a data frame
geno_data_num <- matrix(data = as.integer(geno_data$genotypes@.Data), nrow = 415, ncol = 345762,
                        dimnames = list(rownames(geno_data$genotypes@.Data), colnames(geno_data$genotypes@.Data))) %>%
  as.data.frame()

# Print first 5 rows and columns of the genotype data for verification
geno_data$genotypes@.Data[1:5, 1:5]
geno_data_num[1:5, 1:5]

# Check row names of the genotype data
rownames(geno_data$fam)
rownames(geno_data_num)

# Read and filter climate data
cyano_climate_data <- read_csv("./data/figure_3_cyano_buf_climate_data_419.csv") %>%
  filter(genotype %in% rownames(geno_data$fam)) %>%
  arrange(genotype)

# Verify if row names of genotype data match with climate data genotypes
identical(rownames(geno_data_num), cyano_climate_data$genotype)

colnames(cyano_climate_data)

#############

dim(geno_data_num)

# Setup parallel backend to use multiple cores
no_cores <- detectCores()
cl <- makeCluster(no_cores)



# Calculate simple linear model
compute_stats <- function(i) {
  model <- lm(geno_data_num[, i] ~ as.matrix(cyano_climate_data[, 16]))
  rss1 <- sum(resid(model)^2)
  
  # Null model for F-statistic
  null_model <- lm(geno_data_num[, i] ~ 1)
  rss0 <- sum(resid(null_model)^2)
  
  # Calculate statistics
  p <- length(coef(model)) - 1
  n <- nrow(geno_data_num)
  f_statistic <- ((rss0 - rss1) / p) / (rss1 / (n - p - 1))
  p_value <- pf(f_statistic, p, n - p - 1, lower.tail = FALSE)
  
  return(data.frame(f_statistic = f_statistic, p_value = p_value, df1 = p, df2 = n - p - 1))
}

# Export necessary objects and functions to the cluster
clusterExport(cl, list("geno_data_num", "cyano_climate_data", "lm", "resid", "sum", "pf", "nrow", "coef"))

# Compute in parallel
lm_result <- parLapply(cl, 1:ncol(geno_data_num), compute_stats)

# Stop the cluster
stopCluster(cl)

# Combine results into a data frame
lm_result_c <- bind_rows(lm_result)
median_f = qf(0.5, df1 = 1, df2 = 413, lower.tail = FALSE)
ob_median_f <- median(lm_result_c$f_statistic)
GIF <- ob_median_f / median_f
lm_result_c <- lm_result_c %>%
  mutate(adjusted_F = f_statistic / GIF) %>%
  mutate(adjusted_p = pf(adjusted_F, df1, df2, lower.tail = FALSE))

# Prepare and save simple linear model output
simpleLM_out <- tibble(CHR = as.integer(str_remove(geno_data$map$chromosome, "chr_")),
                       POSITION = geno_data$map$position,
                       ID = geno_data$map$snp.name,
                       climate.var = colnames(cyano_climate_data[, 16]),
                       k.value = 0) %>%
  bind_cols(lm_result_c)

write_tsv(simpleLM_out, "simple_lm_GDD5_20240516.tsv")

############ LFMM ####################

# Run LFMM for multiple K values and climate variables
for (k in 1:3) {
  for (i in 10:44) {
    print(k)
    print(i)
    
    clover.lfmm <- lfmm_ridge(Y = geno_data_num, X = as.matrix(cyano_climate_data[, i]), K = k)
    clover.pv <- lfmm_test(Y = geno_data_num, X = as.matrix(cyano_climate_data[, i]), lfmm = clover.lfmm, calibrate = "gif")
    
    LFMM_out <- tibble(CHR = as.integer(str_remove(geno_data$map$chromosome, "chr_")),
                       POSITION = geno_data$map$position,
                       ID = geno_data$map$snp.name,
                       pvalue = clover.pv$calibrated.pvalue[, 1],
                       climate.var = colnames(cyano_climate_data[, i]),
                       k.value = k)
    
    write_tsv(LFMM_out, "LFMM_output_20230606.txt", append = TRUE)
  }
}
