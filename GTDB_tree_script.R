# GTDB tree 
## Get GTDB tree, based on genomes 

library(tidyverse)
library(data.table)
library(ape)


## Load up data ----
GTDB_metadata <- fread("bac120_metadata_r214.tsv", header=TRUE)

GTDB_metadata_key_info <- GTDB_metadata %>%
  select(accession, gtdb_genome_representative)

GTDB_tree<- read.tree("bac120_r214.tree")

SOC_summary_data_wide_146 <- fread("SOC_summary_data_wide_146.csv")
## Match genome data with GTDB database to get species info

# First create column with genome name
GTDB_metadata_key_info <- GTDB_metadata_key_info %>%
  mutate(genome_accession_no_version = str_sub(accession, 4, -3))

SOC_summary_data_wide_146 <- SOC_summary_data_wide_146 %>%
  mutate(genome_accession_no_version = str_sub(genome_accession, 1, -3))

species_GTDB_info <- left_join(SOC_summary_data_wide_146, GTDB_metadata_key_info, join_by("genome_accession_no_version"))

# filter out genomes without match
species_GTDB_subset <- inner_join(SOC_summary_data_wide_146, GTDB_metadata_key_info, join_by("genome_accession_no_version"))

#write.csv(species_GTDB_subset, "species_in_gtdb_SOC.csv")

# Drop tips that aren't in data set
# Creates what you don't want
dropTip2<-GTDB_tree$tip.label[which(is.na(match(GTDB_tree$tip.label, species_GTDB_subset$gtdb_genome_representative)))]
# Gets rid of it
GTDB_tree_subset<-drop.tip(GTDB_tree, dropTip2, trim.internal=T)

plot(GTDB_tree_subset)

# Finalise
is.ultrametric(GTDB_tree_subset)
tree <- chronoMPL(GTDB_tree_subset) # Make it ultrametric if it isn't

## Remove 0 branch lengths by changing to 1 then remaking the tree ultrametric - repeat until no 0 branch lengths
GTDB_tree_subset$edge.length[GTDB_tree_subset$edge.length<=0]<-1
GTDB_tree_subset<-chronoMPL(GTDB_tree_subset)
which(GTDB_tree_subset$edge.length<=0)

summary(GTDB_tree_subset)
is.rooted(GTDB_tree_subset)

plot(GTDB_tree_subset)

# Write tree
write.tree(GTDB_tree_subset,"species_gtdb.tree")


  






