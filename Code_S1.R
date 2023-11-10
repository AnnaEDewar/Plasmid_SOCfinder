knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
options(tinytex.compile.min_times = 3) # Allows correct page numbers in contents page

library(tidyverse)
library(data.table)
library(MCMCglmm)
library(ape)
library(phylotools)
library(phytools)
library(phylopath)
library(ggsignif)
library(patchwork)
library(treedataverse)
library(phylolm)
library(kableExtra)
library(FactoMineR)
library(factoextra)


# Function to plot 0s in axes nicely
prettyZero <- function(l){
  max.decimals = max(nchar(str_extract(l, "\\.[0-9]+")), na.rm = T)-1
  lnew = formatC(l, replace.zero = T, zero.print = "0",
                 digits = max.decimals, format = "f", preserve.width=T)
  return(lnew)
}

lm_summary_table <- function(model) {
  coeffs <- coef(summary(model))
  df <- c(summary(model)$df[1], summary(model)$df[2])
  
  # Add stars for significant p-values
  p_values <- coeffs[, 4]
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*", "")))
  
  table_data <- data.frame(
    "Estimate" = format(coeffs[, 1], digits = 4, scientific = FALSE),
    "Std.Error" = format(coeffs[, 2], digits = 4, scientific = FALSE),
    "t value" = format(coeffs[, 3], digits = 4, scientific = FALSE),
    "p-value" = signif(coeffs[, 4], digits = 4),
    "signif." = stars
    
  )
  
  return(table_data)
  
}

summary_mcmc_glmm <- function(mcmc_model) {
  summary <- summary(mcmc_model)
  summary_subset <- as.data.frame(summary$solutions)
  
  p_values <- summary_subset$pMCMC
  stars <- ifelse(p_values < 0.001, "***",
                  ifelse(p_values < 0.01, "**",
                         ifelse(p_values < 0.05, "*",
                                ifelse(p_values < 0.1, ".",""))))
  summary_subset$signif. <- stars
  
  # Extract numeric columns
  num_cols <- sapply(summary_subset, is.numeric)
  numeric_table <- summary_subset[, num_cols]
  
  # Format numeric columns to 4 significant figures
  formatted_table <- format(numeric_table, scientific = FALSE, digits = 4)
  
  # Replace original numeric columns with formatted columns
  summary_subset[, num_cols] <- formatted_table
  
  return(summary_subset)
}


#### Data 
SOC_summary_data_wide_146 <- fread("SOC_summary_data_wide_146.csv")

# Replace NAs with 0
SOC_summary_data_wide_146[is.na(SOC_summary_data_wide_146)] <- 0

# Calculate proportions of social on plasmid(s) and chromosome(s)
genome_SOCs <- SOC_summary_data_wide_146 %>%
  mutate(plas_soc_prop=(social_plasmid_TRUE/(social_plasmid_TRUE+social_plasmid_FALSE)),
         chr_soc_prop=(social_chromosome_TRUE/(social_chromosome_TRUE+social_chromosome_FALSE)))


# Arcsine square root transform the proportions 
# Calculate the difference
genome_SOCs_asin_diff <- genome_SOCs %>%
  mutate(plas_soc_prop_asin=asin(sqrt(plas_soc_prop)),
         chr_soc_prop_asin=asin(sqrt(chr_soc_prop)),
         diff_plas_chr=(plas_soc_prop_asin-chr_soc_prop_asin))

## Calculate mean difference & proportions
genome_SOCs_asin_diff_mean <- genome_SOCs_asin_diff %>%
  group_by(species_name) %>% summarise(mean_diff = mean(diff_plas_chr), 
                                       se_diff = (sd(diff_plas_chr)/sqrt(length(diff_plas_chr))),
                                       mean_plas_soc_prop = mean(plas_soc_prop), 
                                       mean_chr_soc_prop = mean(chr_soc_prop))


genome_B_SOCs <- SOC_summary_data_wide_146 %>%
  mutate(plas_soc_prop=(psortb_plasmid_TRUE/(psortb_plasmid_TRUE+psortb_plasmid_FALSE)),
         chr_soc_prop=(psortb_chromosome_TRUE/(psortb_chromosome_TRUE+psortb_chromosome_FALSE)))

genome_B_SOCs_asin_diff <- genome_B_SOCs %>%
  mutate(plas_soc_prop_asin=asin(sqrt(plas_soc_prop)),
         chr_soc_prop_asin=asin(sqrt(chr_soc_prop)),
         diff_plas_chr=(plas_soc_prop_asin-chr_soc_prop_asin))

genome_B_SOCs_diff <- genome_B_SOCs %>%
  mutate(diff_plas_chr=(plas_soc_prop-chr_soc_prop))

## Mean of psortb only
genome_B_SOCs_asin_diff_mean <- genome_B_SOCs_asin_diff %>%
  group_by(species_name) %>% summarise(mean_diff = mean(diff_plas_chr), 
                                       se_diff = (sd(diff_plas_chr)/sqrt(length(diff_plas_chr))),
                                       mean_plas_soc_prop = mean(plas_soc_prop),
                                       mean_chr_soc_prop = mean(chr_soc_prop))

genome_B_SOCs_asin_diff_mean <- as.data.frame(genome_B_SOCs_asin_diff_mean)


genome_K_SOCs <- SOC_summary_data_wide_146 %>%
  mutate(plas_soc_prop=(kofam_plasmid_TRUE/(kofam_plasmid_TRUE+kofam_plasmid_FALSE)),
         chr_soc_prop=(kofam_chromosome_TRUE/(kofam_chromosome_TRUE+kofam_chromosome_FALSE)))

genome_K_SOCs_asin_diff <- genome_K_SOCs %>%
  mutate(plas_soc_prop_asin=asin(sqrt(plas_soc_prop)),
         chr_soc_prop_asin=asin(sqrt(chr_soc_prop)),
         diff_plas_chr=(plas_soc_prop_asin-chr_soc_prop_asin))

genome_K_SOCs_diff <- genome_K_SOCs %>%
  mutate(diff_plas_chr=(plas_soc_prop-chr_soc_prop))

## Mean of kofam only
genome_K_SOCs_asin_diff_mean <- genome_K_SOCs_asin_diff %>%
  group_by(species_name) %>% summarise(mean_diff = mean(diff_plas_chr), 
                                       se_diff = (sd(diff_plas_chr)/sqrt(length(diff_plas_chr))),
                                       mean_plas_soc_prop = mean(plas_soc_prop),
                                       mean_chr_soc_prop = mean(chr_soc_prop))

genome_K_SOCs_asin_diff_mean <- as.data.frame(genome_K_SOCs_asin_diff_mean)


genome_A_SOCs <- SOC_summary_data_wide_146 %>%
  mutate(plas_soc_prop=(antismash_plasmid_TRUE/(antismash_plasmid_TRUE+antismash_plasmid_FALSE)),
         chr_soc_prop=(antismash_chromosome_TRUE/(antismash_chromosome_TRUE+antismash_chromosome_FALSE)))

genome_A_SOCs_asin_diff <- genome_A_SOCs %>%
  mutate(plas_soc_prop_asin=asin(sqrt(plas_soc_prop)),
         chr_soc_prop_asin=asin(sqrt(chr_soc_prop)),
         diff_plas_chr=(plas_soc_prop_asin-chr_soc_prop_asin))

genome_A_SOCs_diff <- genome_A_SOCs %>%
  mutate(diff_plas_chr=(plas_soc_prop-chr_soc_prop))

genome_A_SOCs_asin_diff_mean <- genome_A_SOCs_asin_diff %>%
  group_by(species_name) %>% summarise(mean_diff = mean(diff_plas_chr), 
                                       se_diff = (sd(diff_plas_chr)/sqrt(length(diff_plas_chr))),
                                       mean_plas_soc_prop = mean(plas_soc_prop),
                                       mean_chr_soc_prop = mean(chr_soc_prop))

genome_A_SOCs_asin_diff_mean <- as.data.frame(genome_A_SOCs_asin_diff_mean)

data_subset_A_SOC <- subset(genome_A_SOCs_asin_diff_mean,mean_diff!=0)

data_subset_A_SOC <- as.data.frame(data_subset_A_SOC)


#### Tree 
dataTree<-read.nexus("SOC_plas_species_tree.nex")
gtdb_tree <- read.tree("species_gtdb.tree")
is.ultrametric(dataTree) #SHould say 'TRUE'
is.ultrametric(gtdb_tree)

## Format tree for bayestraits
is.rooted(dataTree)
is.rooted(gtdb_tree)

## check ALL species in data set 'pangenome_lifestyles' are in tree (should be 126)
genome_B_SOCs_asin_diff_mean$species_name[which((genome_SOCs_asin_diff_mean$species_name %in% dataTree$tip.label)==FALSE)]

## Set priors
# Uninformative prior for when 1 random effect
prior <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002))) 
# Uninformative prior for when 2 random effects
prior2 <- list(R=list(V = 1, nu = 0.002), G=list(G1=list(V=1,nu=0.002), G2=list(V=1,nu=0.002))) 

# 'dataTree' is our phylogeny as an ultrametric tree in nexus format
dataTreeNode <- makeNodeLabel(dataTree, method = "number")

## Make matrix of tree called 'Ainv'
INtree <- inverseA(dataTreeNode, nodes="TIPS") # Converts phylogeny into a covariance matrix
Ainv <- INtree$Ainv #Extracts covariance values

# 'dataTree' is our phylogeny as an ultrametric tree in nexus format
gtdb_treeNode <- makeNodeLabel(gtdb_tree, method = "number")

## Make matrix of tree called 'Ainv'
INtree_gtdb <- inverseA(gtdb_treeNode, nodes="TIPS") # Converts phylogeny into a covariance matrix
Ainv_gtdb <- INtree_gtdb$Ainv #Extracts covariance values

genome_SOCs_asin_diff_mean <- as.data.frame(genome_SOCs_asin_diff_mean)

## GTDB Subset
species_gtdb_SOC_data <- fread("species_in_gtdb_SOC.csv")

species_gtdb_SOC_data[is.na(species_gtdb_SOC_data)] <- 0


species_gtdb_SOC_data <- species_gtdb_SOC_data %>%
  mutate(plas_soc_prop=(social_plasmid_TRUE/(social_plasmid_TRUE+social_plasmid_FALSE)),
         chr_soc_prop=(social_chromosome_TRUE/(social_chromosome_TRUE+social_chromosome_FALSE)))

species_gtdb_SOC_asin_diff <- species_gtdb_SOC_data %>%
  mutate(plas_soc_prop_asin=asin(sqrt(plas_soc_prop)),
         chr_soc_prop_asin=asin(sqrt(chr_soc_prop)),
         diff_plas_chr=(plas_soc_prop_asin-chr_soc_prop_asin))

species_gtdb_SOC_asin_diff_count_10 <- species_gtdb_SOC_asin_diff %>%
  group_by(gtdb_genome_representative) %>%
  count() %>% filter(n >=10)

species_gtdb_SOC_asin_diff_count <- species_gtdb_SOC_asin_diff %>%
  group_by(gtdb_genome_representative) %>%
  count() 

species_gtdb_SOC_asin_diff_subset <- inner_join(species_gtdb_SOC_asin_diff, species_gtdb_SOC_asin_diff_count_10, join_by("gtdb_genome_representative"))

species_gtdb_SOC_asin_diff_subset_mean <- species_gtdb_SOC_asin_diff_subset %>%
  group_by(gtdb_genome_representative)%>% summarise(mean_diff = mean(diff_plas_chr), 
                                                    se_diff = (sd(diff_plas_chr)/sqrt(length(diff_plas_chr))),
                                                    mean_plas_soc_prop = mean(plas_soc_prop),
                                                    mean_chr_soc_prop = mean(chr_soc_prop)) 

species_gtdb_SOC_asin_diff_subset_mean <- as.data.frame(species_gtdb_SOC_asin_diff_subset_mean)

species_gtdb_SOC_asin_diff <- as.data.frame(species_gtdb_SOC_asin_diff)

species_gtdb_SOC_asin_diff$species_n_groups <- species_gtdb_SOC_asin_diff$gtdb_genome_representative


mcmc_model_1 <- MCMCglmm(mean_diff ~ 1, random=~species_name, 
                         data=genome_SOCs_asin_diff_mean,
                         nitt=100000, 
                         prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_mcmc_model_1 <- summary_mcmc_glmm(mcmc_model_1)

#random effect variance/(random effect variance + residual variance)

random_r2_model_1 <- (sum(apply(mcmc_model_1$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_1$VCV,2,mean))))

R2_table_1 <- matrix(c(random_r2_model_1), ncol=1)
R2_table_1 <- format(R2_table_1, digits=4)
rownames(R2_table_1) <- c("Random effect")
colnames(R2_table_1) <- c("R-squared value")

knitr::kable(list(summary_mcmc_model_1,R2_table_1),
             caption="Results from the above MCMCglmm with species' mean difference in plasmid and chromosome proportion of genes coding for cooperative traits as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


mcmc_model_2 <- MCMCglmm(mean_diff ~ 1, random=~gtdb_genome_representative, 
                         data=species_gtdb_SOC_asin_diff_subset_mean,
                         nitt=100000, 
                         prior=prior, ginverse = list(gtdb_genome_representative = Ainv_gtdb),verbose = FALSE)



summary_mcmc_model_2 <- summary_mcmc_glmm(mcmc_model_2)

#random effect variance/(random effect variance + residual variance)

random_r2_model_2 <- (sum(apply(mcmc_model_2$VCV,2,mean)[-2]))/((sum(apply(mcmc_model_2$VCV,2,mean))))

R2_table_2 <- matrix(c(random_r2_model_2), ncol=1)
R2_table_2 <- format(R2_table_2, digits=4)
rownames(R2_table_2) <- c("Random effect")
colnames(R2_table_2) <- c("R-squared value")

knitr::kable(list(summary_mcmc_model_2,R2_table_2),
             caption="Results from the above MCMCglmm with species' mean difference in plasmid and chromosome proportion of genes coding for cooperative traits as the response variable and phylogeny (GTDB tree) as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


mcmc_model_3 <- MCMCglmm(diff_plas_chr ~ 1, random=~gtdb_genome_representative+species_n_groups, 
                         data=species_gtdb_SOC_asin_diff,
                         nitt=100000, 
                         prior=prior2, ginverse = list(gtdb_genome_representative = Ainv_gtdb),verbose = FALSE)



summary_mcmc_model_3 <- summary_mcmc_glmm(mcmc_model_3)

#random effect variance/(random effect variance + residual variance)

random_r2_model_3_phylogeny <- (sum(apply(mcmc_model_3$VCV,2,mean)[1]))/((sum(apply(mcmc_model_3$VCV,2,mean))))

random_r2_model_3_sample_size <- (sum(apply(mcmc_model_3$VCV,2,mean)[2]))/((sum(apply(mcmc_model_3$VCV,2,mean))))

total_r2_model_3 <- (sum(apply(mcmc_model_3$VCV,2,mean)[2])+sum(apply(mcmc_model_3$VCV,2,mean)[1]))/((sum(apply(mcmc_model_3$VCV,2,mean))))

R2_table_3 <- matrix(c(random_r2_model_3_phylogeny, random_r2_model_3_sample_size, total_r2_model_3), ncol=1)
R2_table_3 <- format(R2_table_3, digits=4)
rownames(R2_table_3) <- c("Random effect: phylogeny", "Random effect: species", "Total model")
colnames(R2_table_3) <- c("R-squared value")

knitr::kable(list(summary_mcmc_model_3,R2_table_3),
             caption="Results from the above MCMCglmm with genome-level difference in plasmid and chromosome proportion of genes coding for cooperative traits as the response variable, and phylogeny (GTDB tree) and species name as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


## Average genome calculation
genome_SOCs_asin_diff_average_genome <- genome_SOCs_asin_diff %>%
  group_by(species_name) %>%
  summarise(mean_SOC_plas_T = mean(social_plasmid_TRUE), mean_SOC_chr_T=mean(social_chromosome_TRUE),
            mean_SOC_plas_F = mean(social_plasmid_FALSE), mean_SOC_chr_F=mean(social_chromosome_FALSE))

genome_SOCs_asin_diff_social_plas <- genome_SOCs_asin_diff %>%
  mutate(prop_social_plas = asin(sqrt((social_plasmid_TRUE/(social_chromosome_TRUE+social_plasmid_TRUE)))),
         prop_social = asin(sqrt(((social_plasmid_TRUE+social_chromosome_TRUE)/
                                    (social_chromosome_TRUE+social_plasmid_TRUE+social_chromosome_FALSE+social_plasmid_FALSE)))))

genome_SOCs_asin_diff_social_plas_mean <- genome_SOCs_asin_diff_social_plas %>%
  group_by(species_name) %>%
  summarise(mean_prop_soc_plas = mean(prop_social_plas), mean_prop_soc = mean(prop_social))

## Average social genes on plasmid
soc_plas_genes <- mean(genome_SOCs_asin_diff_average_genome$mean_SOC_plas_T)
non_soc_plas_genes <- mean(genome_SOCs_asin_diff_average_genome$mean_SOC_plas_F)
soc_chr_genes <- mean(genome_SOCs_asin_diff_average_genome$mean_SOC_chr_T)
non_soc_chr_genes <- mean(genome_SOCs_asin_diff_average_genome$mean_SOC_chr_F)


genome_SOCs_asin_diff_average_genome <- as.data.frame(genome_SOCs_asin_diff_average_genome)
genome_SOCs_asin_diff_social_plas_mean <- as.data.frame(genome_SOCs_asin_diff_social_plas_mean)


soc_plas_genes_mcmc <- MCMCglmm(mean_SOC_plas_T ~ 1, random=~species_name, 
                         data=genome_SOCs_asin_diff_average_genome,
                         nitt=100000,
                         prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_soc_plas_genes_mcmc <- summary_mcmc_glmm(soc_plas_genes_mcmc)

#random effect variance/(random effect variance + residual variance)

random_r2_soc_plas_genes_mcmc <- (sum(apply(soc_plas_genes_mcmc$VCV,2,mean)[-2]))/((sum(apply(soc_plas_genes_mcmc$VCV,2,mean))))

R2_table_soc_plas <- matrix(c(random_r2_soc_plas_genes_mcmc), ncol=1)
R2_table_soc_plas <- format(R2_table_soc_plas, digits=4)
rownames(R2_table_soc_plas) <- c("Random effect")
colnames(R2_table_soc_plas) <- c("R-squared value")

knitr::kable(list(summary_soc_plas_genes_mcmc,R2_table_soc_plas),
             caption="Results from the above MCMCglmm with species' mean number of genes coding for cooperative traits on plasmid(s) as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


non_soc_plas_genes_mcmc <- MCMCglmm(mean_SOC_plas_F ~ 1, random=~species_name, 
                                data=genome_SOCs_asin_diff_average_genome,
                                nitt=100000, 
                                prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_non_soc_plas_genes_mcmc <- summary_mcmc_glmm(non_soc_plas_genes_mcmc)

#random effect variance/(random effect variance + residual variance)

random_r2_non_soc_plas_genes_mcmc <- (sum(apply(non_soc_plas_genes_mcmc$VCV,2,mean)[-2]))/((sum(apply(non_soc_plas_genes_mcmc$VCV,2,mean))))

R2_table_non_soc_plas <- matrix(c(random_r2_non_soc_plas_genes_mcmc), ncol=1)
R2_table_non_soc_plas <- format(R2_table_non_soc_plas, digits=4)
rownames(R2_table_non_soc_plas) <- c("Random effect")
colnames(R2_table_non_soc_plas) <- c("R-squared value")

knitr::kable(list(summary_non_soc_plas_genes_mcmc,R2_table_non_soc_plas),
             caption="Results from the above MCMCglmm with species' mean number of genes coding for non-cooperative traits on plasmid(s) as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


soc_chr_genes_mcmc <- MCMCglmm(mean_SOC_chr_T ~ 1, random=~species_name, 
                                    data=genome_SOCs_asin_diff_average_genome,
                                    nitt=100000, 
                                    prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_soc_chr_genes_mcmc <- summary_mcmc_glmm(soc_chr_genes_mcmc)

#random effect variance/(random effect variance + residual variance)

random_r2_soc_chr_genes_mcmc <- (sum(apply(soc_chr_genes_mcmc$VCV,2,mean)[-2]))/((sum(apply(soc_chr_genes_mcmc$VCV,2,mean))))

R2_table_soc_chr <- matrix(c(random_r2_soc_chr_genes_mcmc), ncol=1)
R2_table_soc_chr <- format(R2_table_soc_chr, digits=4)
rownames(R2_table_soc_chr) <- c("Random effect")
colnames(R2_table_soc_chr) <- c("R-squared value")

knitr::kable(list(summary_soc_chr_genes_mcmc, R2_table_soc_chr),
             caption="Results from the above MCMCglmm with species' mean number of genes coding for cooperative traits on chromosome(s) as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


non_soc_chr_genes_mcmc <- MCMCglmm(mean_SOC_chr_F ~ 1, random=~species_name, 
                               data=genome_SOCs_asin_diff_average_genome,
                               nitt=100000,
                               prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE) 


summary_non_soc_chr_genes_mcmc <- summary_mcmc_glmm(non_soc_chr_genes_mcmc)

#random effect variance/(random effect variance + residual variance) 

random_r2_non_soc_chr_genes_mcmc <- (sum(apply(non_soc_chr_genes_mcmc$VCV,2,mean)[-2]))/((sum(apply(non_soc_chr_genes_mcmc$VCV,2,mean))))

R2_table_non_soc_chr <- matrix(c(random_r2_non_soc_chr_genes_mcmc), ncol=1)
R2_table_non_soc_chr <- format(R2_table_non_soc_chr, digits=4)
rownames(R2_table_non_soc_chr) <- c("Random effect")
colnames(R2_table_non_soc_chr) <- c("R-squared value")

knitr::kable(list(summary_non_soc_chr_genes_mcmc,R2_table_non_soc_chr),
             caption="Results from the above MCMCglmm with species' mean number of genes coding for non-cooperative traits on chromosome(s) as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


## Average genes on a plasmid
genome_SOCs_average_plasmid <- genome_SOCs_asin_diff %>%
  mutate(prop_plas = asin(sqrt(((social_plasmid_TRUE+social_plasmid_FALSE)/
                                  (social_chromosome_TRUE+social_chromosome_FALSE+social_plasmid_TRUE+social_plasmid_FALSE)))))

genome_SOCs_average_plasmid_mean <- genome_SOCs_average_plasmid %>%
  group_by(species_name) %>%
  summarise(mean_prop_plas = mean(prop_plas))


genome_SOCs_average_plasmid_mean <- as.data.frame(genome_SOCs_average_plasmid_mean)


prop_plas_mean_mcmc <- MCMCglmm(mean_prop_plas ~ 1, random=~species_name, 
                                    data=genome_SOCs_average_plasmid_mean,
                                    nitt=100000, 
                                    prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_prop_plas_mean_mcmc <- summary_mcmc_glmm(prop_plas_mean_mcmc)

## Backtransform so can be read as a poportion
summary_prop_plas_mean_mcmc$post.mean <- sin(as.numeric(summary_prop_plas_mean_mcmc$post.mean))^2
summary_prop_plas_mean_mcmc[,2] <- sin(as.numeric(summary_prop_plas_mean_mcmc[,2]))^2
summary_prop_plas_mean_mcmc[,3] <- sin(as.numeric(summary_prop_plas_mean_mcmc[,3]))^2

#random effect variance/(random effect variance + residual variance) 

random_r2_prop_plas_mean_mcmc <- (sum(apply(prop_plas_mean_mcmc$VCV,2,mean)[-2]))/((sum(apply(prop_plas_mean_mcmc$VCV,2,mean))))

R2_table_prop_plas_mean_mcmc <- matrix(c(random_r2_prop_plas_mean_mcmc), ncol=1)
R2_table_prop_plas_mean_mcmc <- format(R2_table_prop_plas_mean_mcmc, digits=4)
rownames(R2_table_prop_plas_mean_mcmc) <- c("Random effect")
colnames(R2_table_prop_plas_mean_mcmc) <- c("R-squared value")

knitr::kable(list(summary_prop_plas_mean_mcmc,R2_table_prop_plas_mean_mcmc),
             caption="Results from the above MCMCglmm with species' mean proportion of genes carried on plasmid(s) as the response variable and phylogeny as a random effect. We arcsine square-root transformed the proportions before analysis, and we then backtransformed the output of the MCMCglmm to a proportion.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


prop_soc_mcmc <- MCMCglmm(mean_prop_soc ~ 1, random=~species_name, 
                                    data=genome_SOCs_asin_diff_social_plas_mean,
                                    nitt=100000, 
                                    prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_prop_soc_mcmc <- summary_mcmc_glmm(prop_soc_mcmc)

## Backtransform so can be read as a poportion
summary_prop_soc_mcmc$post.mean <- sin(as.numeric(summary_prop_soc_mcmc$post.mean))^2
summary_prop_soc_mcmc[,2] <- sin(as.numeric(summary_prop_soc_mcmc[,2]))^2
summary_prop_soc_mcmc[,3] <- sin(as.numeric(summary_prop_soc_mcmc[,3]))^2

#random effect variance/(random effect variance + residual variance) 
random_r2_prop_soc_mcmc <- (sum(apply(prop_soc_mcmc$VCV,2,mean)[-2]))/((sum(apply(prop_soc_mcmc$VCV,2,mean))))

R2_table_prop_soc_mcmc <- matrix(c(random_r2_prop_soc_mcmc), ncol=1)
R2_table_prop_soc_mcmc <- format(R2_table_prop_soc_mcmc, digits=4)
rownames(R2_table_prop_soc_mcmc) <- c("Random effect")
colnames(R2_table_prop_soc_mcmc) <- c("R-squared value")

knitr::kable(list(summary_prop_soc_mcmc,R2_table_prop_soc_mcmc),
             caption="Results from the above MCMCglmm with species' mean proportion of genes coding for cooperative traits as the response variable and phylogeny as a random effect. We arcsine square-root transformed the proportions before analysis, and we then backtransformed the output of the MCMCglmm to a proportion.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


prop_soc_plas_mean_mcmc <- MCMCglmm(mean_prop_soc_plas ~ 1, random=~species_name, 
                                   data=genome_SOCs_asin_diff_social_plas_mean,
                                   nitt=100000, 
                                   prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_prop_soc_plas_mean_mcmc <- summary_mcmc_glmm(prop_soc_plas_mean_mcmc)

## Backtransform so can be read as a poportion
summary_prop_soc_plas_mean_mcmc$post.mean <- sin(as.numeric(summary_prop_soc_plas_mean_mcmc$post.mean))^2
summary_prop_soc_plas_mean_mcmc[,2] <- sin(as.numeric(summary_prop_soc_plas_mean_mcmc[,2]))^2
summary_prop_soc_plas_mean_mcmc[,3] <- sin(as.numeric(summary_prop_soc_plas_mean_mcmc[,3]))^2

#random effect variance/(random effect variance + residual variance) 

random_r2_prop_soc_plas_mean_mcmc <- (sum(apply(prop_soc_plas_mean_mcmc$VCV,2,mean)[-2]))/((sum(apply(prop_soc_plas_mean_mcmc$VCV,2,mean))))

R2_table_prop_soc_plas_mean_mcmc <- matrix(c(random_r2_prop_soc_plas_mean_mcmc), ncol=1)
R2_table_prop_soc_plas_mean_mcmc <- format(R2_table_prop_soc_plas_mean_mcmc, digits=4)
rownames(R2_table_prop_soc_plas_mean_mcmc) <- c("Random effect")
colnames(R2_table_prop_soc_plas_mean_mcmc) <- c("R-squared value")

knitr::kable(list(summary_prop_soc_plas_mean_mcmc,R2_table_prop_soc_plas_mean_mcmc),
             caption="Results from the above MCMCglmm with species' mean proportion of genes for cooperative traits carried on plasmid(s) as the response variable and phylogeny as a random effect. We arcsine square-root transformed the proportions before analysis, and we then backtransformed the output of the MCMCglmm to a proportion.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


B_SOCS_diff_mcmc <- MCMCglmm(mean_diff ~ 1, random=~species_name, 
                                data=genome_B_SOCs_asin_diff_mean,
                                nitt=100000, 
                                prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_B_SOCS_diff_mcmc <- summary_mcmc_glmm(B_SOCS_diff_mcmc)

#random effect variance/(random effect variance + residual variance) 

random_r2_B_SOCS_diff_mcmc <- (sum(apply(B_SOCS_diff_mcmc$VCV,2,mean)[-2]))/((sum(apply(B_SOCS_diff_mcmc$VCV,2,mean))))

R2_table_B_SOCS_diff_mcmc <- matrix(c(random_r2_B_SOCS_diff_mcmc), ncol=1)
R2_table_B_SOCS_diff_mcmc <- format(R2_table_B_SOCS_diff_mcmc, digits=4)
rownames(R2_table_B_SOCS_diff_mcmc) <- c("Random effect")
colnames(R2_table_B_SOCS_diff_mcmc) <- c("R-squared value")

knitr::kable(list(summary_B_SOCS_diff_mcmc,R2_table_B_SOCS_diff_mcmc),
             caption="Results from the above MCMCglmm with species' mean difference in plasmid and chromosome proportion of genes coding for extracellular proteins as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


K_SOCS_diff_mcmc <- MCMCglmm(mean_diff ~ 1, random=~species_name, 
                                    data=genome_K_SOCs_asin_diff_mean,
                                    nitt=100000, 
                                    prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_K_SOCS_diff_mcmc <- summary_mcmc_glmm(K_SOCS_diff_mcmc)

#random effect variance/(random effect variance + residual variance) 

random_r2_K_SOCS_diff_mcmc <- (sum(apply(K_SOCS_diff_mcmc$VCV,2,mean)[-2]))/((sum(apply(K_SOCS_diff_mcmc$VCV,2,mean))))

R2_table_K_SOCS_diff_mcmc <- matrix(c(random_r2_K_SOCS_diff_mcmc), ncol=1)
R2_table_K_SOCS_diff_mcmc <- format(R2_table_K_SOCS_diff_mcmc, digits=4)
rownames(R2_table_K_SOCS_diff_mcmc) <- c("Random effect")
colnames(R2_table_K_SOCS_diff_mcmc) <- c("R-squared value")

knitr::kable(list(summary_K_SOCS_diff_mcmc,R2_table_K_SOCS_diff_mcmc),
             caption="Results from the above MCMCglmm with species' mean difference in plasmid and chromosome proportion of genes with a cooperative functional annotation as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


A_SOCS_diff_mcmc_subset <- MCMCglmm(mean_diff ~ 1, random=~species_name, 
                             data=data_subset_A_SOC,
                             nitt=100000, 
                             prior=prior, ginverse = list(species_name = Ainv),verbose = FALSE)


summary_A_SOCS_diff_mcmc_subset <- summary_mcmc_glmm(A_SOCS_diff_mcmc_subset)

#random effect variance/(random effect variance + residual variance) 

random_r2_A_SOCS_diff_mcmc_subset <- (sum(apply(A_SOCS_diff_mcmc_subset$VCV,2,mean)[-2]))/((sum(apply(A_SOCS_diff_mcmc_subset$VCV,2,mean))))

R2_table_A_SOCS_diff_mcmc_subset <- matrix(c(random_r2_A_SOCS_diff_mcmc_subset), ncol=1)
R2_table_A_SOCS_diff_mcmc_subset <- format(R2_table_A_SOCS_diff_mcmc_subset, digits=4)
rownames(R2_table_A_SOCS_diff_mcmc_subset) <- c("Random effect")
colnames(R2_table_A_SOCS_diff_mcmc_subset) <- c("R-squared value")

knitr::kable(list(summary_A_SOCS_diff_mcmc_subset,R2_table_A_SOCS_diff_mcmc_subset),
             caption="Results from the above MCMCglmm with species' mean difference in plasmid and chromosome proportion of genes which are part of a cooperative secondary metabolite cluster as the response variable and phylogeny as a random effect.", digits=4) %>% kable_styling(latex_options = "HOLD_position")


ggplot(genome_SOCs_asin_diff_mean, aes(y=fct_reorder(species_name, mean_diff, .desc=TRUE), x=mean_diff)) +
  geom_col(aes(fill=mean_diff >=0), width=0.7) +
  geom_errorbarh(aes(xmin=(mean_diff-se_diff), xmax=(mean_diff+se_diff)), height=0.5) +
  labs(x="Difference in proportion",y="Species") +
  geom_vline(xintercept=0, linetype="dashed", colour="black") +
  scale_fill_manual(name = "Overrepresented on:", values = setNames(c("#56B4E9","#D55E00"),c(T, F)), labels=c("Chromosome","Plasmid")) +
  theme_classic() +
  theme(plot.title = element_text(size=13, hjust=0.5), axis.title = element_text(size=14),
        axis.text.x = element_text(size=13), #axis.title.y=element_text(angle=360,vjust=0.5), 
        legend.text = element_text(size=13), legend.title = element_text(size=14),
        axis.ticks.y=element_blank(), axis.text.y = element_text(size=5.5, face="italic"), plot.margin = margin(t = 25, unit = "pt"))


ggtree(dataTree,layout='fan',branch.length = "none",open.angle = 25) +
  geom_rootedge(T) +
  geom_tiplab(size=2, fontface=3)


ggtree(gtdb_tree,layout='fan',branch.length = "none",open.angle = 25) +
  geom_rootedge(T) +
  geom_tiplab(size=2)


species_table <- full_join(genome_SOCs_asin_diff_average_genome, genome_SOCs_asin_diff_mean, by="species_name") 

species_table <- species_table %>%
  select(species_name,mean_SOC_plas_T,mean_SOC_chr_T,mean_SOC_plas_F,mean_SOC_chr_F)

species_table$mean_SOC_plas_T <- round(species_table$mean_SOC_plas_T)
species_table$mean_SOC_chr_T <- round(species_table$mean_SOC_chr_T)
species_table$mean_SOC_plas_F <- round(species_table$mean_SOC_plas_F)
species_table$mean_SOC_chr_F <- round(species_table$mean_SOC_chr_F)

species_table <- species_table %>%
    mutate(species_name = str_replace(species_name, "_", " "))

species_table <- as.data.frame(species_table)

knitr::kable(species_table,longtable=TRUE,booktabs=TRUE, linesep="",
             caption="List of all all 146 species in our dataset, and their mean number of genes for cooperative and non-cooperative traits carried on plasmid(s) and chromosome(s) in a genome.",
             col.names = c("Species","Plas. & coop.","Chr. & coop.",
                           "Plas. & non-coop.", "Chr & non-coop."))%>%
  kable_styling(latex_options = c("HOLD_position","repeat_header"), font_size=8) %>%
  column_spec(1,italic=T)

