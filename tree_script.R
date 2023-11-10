#### Packages ----
library(ape)
library(phylotools)
library(doBy)
library(tidyverse)
library(phangorn)

library(MCMCglmm)
library(metafor)
library(dplyr)
library(reshape2)
library(phytools)

## Import tree into R  ----
tree_newick <- read.tree("Original_tree.txt") # In newick format

write.nexus(tree_newick, file="Original_tree.nex")

tree_nexus <- read.nexus("Original_tree.nex")
summary(tree_nexus)

tree <- tree_nexus

is.ultrametric(tree)
tree <- chronoMPL(tree) # Make it ultrametric if it isn't

species_data <- read.csv("species_146.csv", header=T)

## Make data frame of tip labels of tree
#write.csv(data.frame(treeSp = tree$tip.label), file = "treeSp.csv")

#write.csv(data.frame(treeSp = dataTree$tip.label), file = "treeSp_dataTree.csv")

#### Replacing with proper tip name
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodospirillales_Acetobacteraceae_Acetobacter_pasteurianus_386B")] <- "Acetobacter_pasteurianus"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pseudomonadales_Moraxellaceae_Acinetobacter_baumannii_AB30")] <- "Acinetobacter_baumannii"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Aeromonadales_Aeromonadaceae_Aeromonas_hydrophila_AL09_71")] <- "Aeromonas_hydrophila"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Alteromonadales_Alteromonadaceae_Alteromonas_macleodii_Black_Sea_11")] <- "Alteromonas_macleodii"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Bacillaceae_Bacillus_anthracis_52_G")] <- "Bacillus_anthracis"
tree$tip.label[which(tree$tip.label == "Bacteria_Bacteroidetes_Chlorobi_group_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides_fragilis_NCTC_9343")] <- "Bacteroides_fragilis"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Buchnera_Buchnera_aphidicola_str._APS_Acyrthosiphon_pisum")] <- "Buchnera_aphidicola"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Burkholderiales_Burkholderiaceae_Burkholderia_cenocepacia_AU_1054")] <- "Burkholderia_cenocepacia"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_deltaepsilon_subdivisions_Epsilonproteobacteria_Campylobacterales_Campylobacteraceae_Campylobacter_jejuni_jejuni_NCTC_11168")] <- "Campylobacter_jejuni"	
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Micrococcineae_Microbacteriaceae_Clavibacter_michiganensis_subsp._michiganensis_NCPPB_382")] <- "Clavibacter_michiganensis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Clostridium_difficile_CF5")] <- "Clostridioides_difficile"	#same species
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Legionellales_Coxiellaceae_Coxiella_burnetii_CbuK_Q154")] <- "Coxiella_burnetii"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Cronobacter_sakazakii_CMCC_45402")] <- "Cronobacter_sakazakii"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Enterobacter_cloacae_SCF1")] <- "Enterobacter_cloacae"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Enterococcaceae_Enterococcus_faecalis_D32")] <- "Enterococcus_faecalis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Erwinia_amylovora_CFBP1430")] <- "Erwinia_amylovora"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Escherichia_coli_K12__W3110")] <- "Escherichia_coli"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Listeriaceae_Listeria_innocua_sv._6a_Clip11262")] <- "Listeria_innocua"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Morganella_morganii_morganii_KT")] <- "Morganella_morganii"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Neisseriales_Neisseriaceae_Neisseria_gonorrhoeae_FA_1090")] <- "Neisseria_gonorrhoeae"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Pantoea_ananatis_LMG_5342")] <- "Pantoea_ananatis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Bacteroidetes_Chlorobi_group_Bacteroidetes_Bacteroidia_Bacteroidales_Porphyromonadaceae_Parabacteroides_distasonis_ATCC_8503")] <- "Parabacteroides_distasonis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pasteurellales_Pasteurellaceae_Pasteurella_multocida_43137")] <- "Pasteurella_multocida"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Betaproteobacteria_Burkholderiales_Burkholderiaceae_Ralstonia_solanacearum_GMI1000")] <- "Ralstonia_solanacearum"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Raoultella_ornithinolytica_B6")] <- "Raoultella_ornithinolytica"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_RhizobiumAgrobacterium_group_Rhizobium_leguminosarum_bv._trifolii_CB782")] <- "Rhizobium_leguminosarum"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Salmonella_enterica_Enterica_sv._Heidelberg_41578")] <- "Salmonella_enterica"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Shigella_flexneri_2003036")] <- "Shigella_flexneri"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_SinorhizobiumEnsifer_group_Ensifer_meliloti_GR4")] <- "Sinorhizobium_meliloti"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Bacillales_Staphylococcaceae_Staphylococcus_aureus_502A")] <- "Staphylococcus_aureus"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Xanthomonas_oryzae_ATCC_35933")] <- "Xanthomonas_oryzae" #checked genus phylogeny
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Xanthomonadales_Xanthomonadaceae_Xylella_fastidiosa_sandyi_Ann_1")] <- "Xylella_fastidiosa"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Yersinia_pestis_Nepal516")] <- "Yersinia_pestis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Sphingomonadales_Sphingomonadaceae_Zymomonas_mobilis_mobilis_CP4")] <- "Zymomonas_mobilis"	



#### Add species to tree, replace other species in same genus
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhizobiales_Rhizobiaceae_RhizobiumAgrobacterium_group_Agrobacterium_sp._H13_3")] <- "Agrobacterium_fabrum"
tree$tip.label[which(tree$tip.label == "Bacteria_Spirochaetes_Spirochaetia_Spirochaetales_Spirochaetaceae_Borrelia_garinii_PBi_linear")] <- "Borrelia_miyamotoi"
tree$tip.label[which(tree$tip.label == "Bacteria_Chlamydiae_Verrucomicrobia_group_Chlamydiae_Chlamydiia_Chlamydiales_Chlamydiaceae_ChlamydiaChlamydophila_group_Chlamydophila_abortus_S263")] <- "Chlamydia_psittaci"
tree$tip.label[which(tree$tip.label == "Bacteria_Chlamydiae_Verrucomicrobia_group_Chlamydiae_Chlamydiia_Chlamydiales_Chlamydiaceae_ChlamydiaChlamydophila_group_Chlamydophila_pneumoniae_AR39")] <- "Chlamydia_trachomatis"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Citrobacter_rodentium_ICC168")] <- "Citrobacter_freundii"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_Clostridium_acetobutylicum_EA_2018")] <- "Clostridium_botulinum" 
tree$tip.label[which(tree$tip.label == "Bacteria_Deinococcus_Thermus_Deinococci_Deinococcales_Deinococcaceae_Deinococcus_geothermalis_DSM_11300")] <- "Deinococcus_radiodurans" 
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_deltaepsilon_subdivisions_Epsilonproteobacteria_Campylobacterales_Helicobacteraceae_Helicobacter_cetorum_MIT_99_5656")] <- "Helicobacter_pylori"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Klebsiella_oxytoca_HKOPL1")] <- "Klebsiella_oxytoca"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Streptococcaceae_Lactococcus_garvieae_Lg2")] <- "Lactococcus_lactis"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Lactobacillaceae_Lactobacillus_acidophilus_NCFM")] <- "Lacticaseibacillus_paracasei"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Legionellales_Legionellaceae_Legionella_micdadei_ATCC_33218")] <- "Legionella_pneumophila"
tree$tip.label[which(tree$tip.label == "Bacteria_Spirochaetes_Spirochaetia_Spirochaetales_Leptospiraceae_Leptospira_biflexa_serovar_Patoc_strain_Patoc_1_Paris_I")] <- "Leptospira_interrogans"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Leuconostocaceae_Leuconostoc_kimchii_IMSNU11154")] <- "Leuconostoc_mesenteroides"
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Corynebacterineae_Mycobacteriaceae_Mycobacterium_abscessus_bolletii_50594")] <- "Mycobacteroides_abscessus"	
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Lactobacillaceae_Pediococcus_claussenii_P06_ATCC_BAA_344")] <- "Pediococcus_acidilactici"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodobacterales_Rhodobacteraceae_Phaeobacter_gallaeciensis_2.10")] <- "Phaeobacter_inhibens"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Vibrionales_Vibrionaceae_Photobacterium_profundum_SS9")] <- "Photobacterium_damselae"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Thiotrichales_Piscirickettsiaceae_Piscirickettsia_salmonis_LF_89")] <- "Piscirickettsia_salmonis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Corynebacterineae_Nocardiaceae_Rhodococcus_sp._BCP1")] <- "Prescottella_equi"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Proteus_mirabilis_WGLW4")] <- "Proteus_mirabilis"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Providencia_stuartii_MRSN_2154")] <- "Providencia_rettgeri"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Pseudomonadales_Pseudomonadaceae_Pseudomonas_aeruginosa_B136_33")] <- "Pseudomonas_aeruginosa"	
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Alphaproteobacteria_Rhodobacterales_Rhodobacteraceae_Roseobacter_litoralis_Och_149")] <- "Pseudosulfitobacter_pseudonitzschiae"	# See note, new genus
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacteriales_Enterobacteriaceae_Serratia_fonticola_RB_25")] <- "Serratia_marcescens"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Streptococcaceae_Streptococcus_mutans_GS_5")] <- "Streptococcus_thermophilus"
tree$tip.label[which(tree$tip.label == "Bacteria_Actinobacteria_Actinobacteria_Actinobacteridae_Actinomycetales_Streptomycineae_Streptomycetaceae_Streptomyces_coelicolor_A32")] <- "Streptomyces_clavuligerus"
tree$tip.label[which(tree$tip.label == "Bacteria_Deinococcus_Thermus_Deinococci_Thermales_Thermaceae_Thermus_oshimai_JL_2")] <- "Thermus_thermophilus"
tree$tip.label[which(tree$tip.label == "Bacteria_Proteobacteria_Gammaproteobacteria_Vibrionales_Vibrionaceae_Vibrio_anguillarum_775_I")] <- "Vibrio_cholerae"
tree$tip.label[which(tree$tip.label == "Bacteria_Firmicutes_Bacilli_Lactobacillales_Leuconostocaceae_Weissella_ceti_WS74")] <- "Weissella_cibaria"


#### REMOVE EVERYTHING ELSE FROM TREE

# Drop tips that aren't in data set
# Creates what you don't want
dropTip2<-tree$tip.label[which(is.na(match(tree$tip.label, species_data$species)))]
# Gets rid of it
dataTree<-drop.tip(tree, dropTip2, trim.internal=T)
plot(dataTree)

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


#### ADDING SPECIES TO GENUS

## Agrobacterium_tumefaciens
which(dataTree$tip.label=="Agrobacterium_fabrum")		# 33 (edge number of this tip)
which(dataTree$edge[,2] == 33)	  #  74 (row in edge matrix with this end node)
dataTree$edge.length[74]									# 0.06 (edge length)
dataTree <- bind.tip(dataTree, "Agrobacterium_tumefaciens", edge.length=NULL, where=33, position=dataTree$edge.length[74]/4) #4 if adding in first species, /2 if clustering further

## Bacteroides_thetaiotaomicron
which(dataTree$tip.label=="Bacteroides_fragilis")		# 47 (edge number of this tip)
which(dataTree$edge[,2] == 47)	  #  94 (row in edge matrix with this end node)
dataTree$edge.length[94]									# 0.22 (edge length)
dataTree <- bind.tip(dataTree, "Bacteroides_thetaiotaomicron", edge.length=NULL, where=47, position=dataTree$edge.length[94]/4) #4 if adding in first species, /2 if clustering further

## Borreliella_burgdorferi
which(dataTree$tip.label=="Borrelia_miyamotoi")		# 42 (edge number of this tip)
which(dataTree$edge[,2] == 42)	  #  87 (row in edge matrix with this end node)
dataTree$edge.length[87]									# 0.76 (edge length)
dataTree <- bind.tip(dataTree, "Borreliella_burgdorferi", edge.length=NULL, where=42, position=dataTree$edge.length[87]/4) #4 if adding in first species, /2 if clustering further

## Citrobacter_portucalensis
which(dataTree$tip.label=="Citrobacter_freundii")		# 1 (edge number of this tip)
which(dataTree$edge[,2] == 1)	  #  26 (row in edge matrix with this end node)
dataTree$edge.length[26]									# 0.004 (edge length)
dataTree <- bind.tip(dataTree, "Citrobacter_portucalensis", edge.length=NULL, where=1, position=dataTree$edge.length[26]/4) #4 if adding in first species, /2 if clustering further

## Lactococcus_cremoris
which(dataTree$tip.label=="Lactococcus_lactis")		# 62 (edge number of this tip)
which(dataTree$edge[,2] == 62)	  #  127 (row in edge matrix with this end node)
dataTree$edge.length[127]									# 0.117 (edge length)
dataTree <- bind.tip(dataTree, "Lactococcus_cremoris", edge.length=NULL, where=62, position=dataTree$edge.length[127]/4) #4 if adding in first species, /2 if clustering further

## Leptospira_noguchii
which(dataTree$tip.label=="Leptospira_interrogans")		# 45 (edge number of this tip)
which(dataTree$edge[,2] == 45)	  #  92 (row in edge matrix with this end node)
dataTree$edge.length[92]									# 0.76 (edge length)
dataTree <- bind.tip(dataTree, "Leptospira_noguchii", edge.length=NULL, where=45, position=dataTree$edge.length[92]/4) #4 if adding in first species, /2 if clustering further

## Listeria_monocytogenes
which(dataTree$tip.label=="Listeria_innocua")		# 66 (edge number of this tip)
which(dataTree$edge[,2] == 66)	  #  133 (row in edge matrix with this end node)
dataTree$edge.length[133]									# 0.35 (edge length)
dataTree <- bind.tip(dataTree, "Listeria_monocytogenes", edge.length=NULL, where=66, position=dataTree$edge.length[133]/4) #4 if adding in first species, /2 if clustering further

## Pantoea_agglomerans
which(dataTree$tip.label=="Pantoea_ananatis")		# 10 (edge number of this tip)
which(dataTree$edge[,2] == 10)	  #  37 (row in edge matrix with this end node)
dataTree$edge.length[37]									# 0.032 (edge length)
dataTree <- bind.tip(dataTree, "Pantoea_agglomerans", edge.length=NULL, where=10, position=dataTree$edge.length[37]/4) #4 if adding in first species, /2 if clustering further

## Pediococcus_pentosaceus
which(dataTree$tip.label=="Pediococcus_acidilactici")		# 61 (edge number of this tip)
which(dataTree$edge[,2] == 61)	  #  126 (row in edge matrix with this end node)
dataTree$edge.length[126]									# 0.26 (edge length)
dataTree <- bind.tip(dataTree, "Pediococcus_pentosaceus", edge.length=NULL, where=61, position=dataTree$edge.length[126]/4) #4 if adding in first species, /2 if clustering further

## Rhizobium_phaseoli
which(dataTree$tip.label=="Rhizobium_leguminosarum")		# 34 (edge number of this tip)
which(dataTree$edge[,2] == 34)	  #  77 (row in edge matrix with this end node)
dataTree$edge.length[77]									# 0.060 (edge length)
dataTree <- bind.tip(dataTree, "Rhizobium_phaseoli", edge.length=NULL, where=34, position=dataTree$edge.length[77]/4) 

## Streptococcus_suis
which(dataTree$tip.label=="Streptococcus_thermophilus")		# 68 (edge number of this tip)
which(dataTree$edge[,2] == 68)	  #  138 (row in edge matrix with this end node)
dataTree$edge.length[138]									# 0.117 (edge length)
dataTree <- bind.tip(dataTree, "Streptococcus_suis", edge.length=NULL, where=68, position=dataTree$edge.length[138]/4) #4 if adding in first species, /2 if clustering further



# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]

plot(dataTree)
#### Adding multiple species to  genus ----

## Acinetobacter ----
## Acinetobacter_indicus
which(dataTree$tip.label=="Acinetobacter_baumannii")		# 25 (edge number of this tip)
which(dataTree$edge[,2] == 25)	  #  59 (row in edge matrix with this end node)
dataTree$edge.length[59]									# 0.32 (edge length)
dataTree <- bind.tip(dataTree, "Acinetobacter_indicus", edge.length=NULL, where=25, position=dataTree$edge.length[59]/4) #4 if adding in first species, /2 if clustering further

## Acinetobacter_johnsonii
which(dataTree$tip.label=="Acinetobacter_indicus")		# 26 (edge number of this tip)
which(dataTree$edge[,2] == 26)	  #  61 (row in edge matrix with this end node)
dataTree$edge.length[61]									# 0.08 (edge length)
dataTree <- bind.tip(dataTree, "Acinetobacter_johnsonii", edge.length=NULL, where=26, position=dataTree$edge.length[61]/2) #4 if adding in first species, /2 if clustering further

## Acinetobacter_haemolyticus
which(dataTree$tip.label=="Acinetobacter_baumannii")		# 25 (edge number of this tip)
which(dataTree$edge[,2] == 25)	  #  60 (row in edge matrix with this end node)
dataTree$edge.length[60]									# 0.08 (edge length)
dataTree <- bind.tip(dataTree, "Acinetobacter_haemolyticus", edge.length=NULL, where=25, position=dataTree$edge.length[60]/2) #4 if adding in first species, /2 if clustering further

## Acinetobacter_pittii
which(dataTree$tip.label=="Acinetobacter_baumannii")		# 25 (edge number of this tip)
which(dataTree$edge[,2] == 25)	  #  61 (row in edge matrix with this end node)
dataTree$edge.length[61]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Acinetobacter_pittii", edge.length=NULL, where=25, position=dataTree$edge.length[61]/2) #4 if adding in first species, /2 if clustering further

## Acinetobacter_seifertii
which(dataTree$tip.label=="Acinetobacter_baumannii")		# 25 (edge number of this tip)
which(dataTree$edge[,2] == 25)	  #  62 (row in edge matrix with this end node)
dataTree$edge.length[62]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Acinetobacter_seifertii", edge.length=NULL, where=25, position=dataTree$edge.length[62]/2) #4 if adding in first species, /2 if clustering further

## Aeromonas 
## Aeromonas_caviae
which(dataTree$tip.label=="Aeromonas_hydrophila")		# 22 (edge number of this tip)
which(dataTree$edge[,2] == 22)	  #  55 (row in edge matrix with this end node)
dataTree$edge.length[55]									# 0.26 (edge length)
dataTree <- bind.tip(dataTree, "Aeromonas_caviae", edge.length=NULL, where=22, position=dataTree$edge.length[55]/4) #4 if adding in first species, /2 if clustering further

## Aeromonas_veronii
which(dataTree$tip.label=="Aeromonas_hydrophila")		# 22 (edge number of this tip)
which(dataTree$edge[,2] == 22)	  #  56 (row in edge matrix with this end node)
dataTree$edge.length[56]									# 0.066 (edge length)
dataTree <- bind.tip(dataTree, "Aeromonas_veronii", edge.length=NULL, where=22, position=dataTree$edge.length[56]/2) #4 if adding in first species, /2 if clustering further

## Aeromonas_salmonicida
which(dataTree$tip.label=="Aeromonas_veronii")		# 23 (edge number of this tip)
which(dataTree$edge[,2] == 23)	  #  58 (row in edge matrix with this end node)
dataTree$edge.length[58]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Aeromonas_salmonicida", edge.length=NULL, where=23, position=dataTree$edge.length[58]/2) #4 if adding in first species, /2 if clustering further


## Bacillus ----
## Add Priestia first, actualy nested within Bacillus genera (subtilis and cereus clades form separate monophyletic groups)
## Bacillus_subtilis
which(dataTree$tip.label=="Bacillus_anthracis")		# 81 (edge number of this tip)
which(dataTree$edge[,2] == 81)	  #  161 (row in edge matrix with this end node)
dataTree$edge.length[161]									# 0.37 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_subtilis", edge.length=NULL, where=81, position=dataTree$edge.length[161]/4) #4 if adding in first species, /2 if clustering further

## Priestia_megaterium
which(dataTree$tip.label=="Bacillus_anthracis")		# 81 (edge number of this tip)
which(dataTree$edge[,2] == 81)	  #  162 (row in edge matrix with this end node)
dataTree$edge.length[162]									# 0.09 (edge length)
dataTree <- bind.tip(dataTree, "Priestia_megaterium", edge.length=NULL, where=81, position=dataTree$edge.length[162]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_cytotoxicus
which(dataTree$tip.label=="Bacillus_anthracis")		# 81 (edge number of this tip)
which(dataTree$edge[,2] == 81)	  #  163 (row in edge matrix with this end node)
dataTree$edge.length[163]									# 0.046 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_cytotoxicus", edge.length=NULL, where=81, position=dataTree$edge.length[163]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_mycoides
which(dataTree$tip.label=="Bacillus_anthracis")		# 81 (edge number of this tip)
which(dataTree$edge[,2] == 81)	  #  164 (row in edge matrix with this end node)
dataTree$edge.length[164]									# 0.023 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_mycoides", edge.length=NULL, where=81, position=dataTree$edge.length[164]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_toyonensis
which(dataTree$tip.label=="Bacillus_anthracis")		# 81 (edge number of this tip)
which(dataTree$edge[,2] == 81)	  #  165 (row in edge matrix with this end node)
dataTree$edge.length[165]									# 0.011 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_toyonensis", edge.length=NULL, where=81, position=dataTree$edge.length[165]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_cereus
which(dataTree$tip.label=="Bacillus_anthracis")		# 81 (edge number of this tip)
which(dataTree$edge[,2] == 81)	  #  166 (row in edge matrix with this end node)
dataTree$edge.length[166]									# 0.006 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_cereus", edge.length=NULL, where=81, position=dataTree$edge.length[166]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_thuringiensis
which(dataTree$tip.label=="Bacillus_cereus")		# 82 (edge number of this tip)
which(dataTree$edge[,2] == 82)	  #  168 (row in edge matrix with this end node)
dataTree$edge.length[168]									# 0.0003 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_thuringiensis", edge.length=NULL, where=82, position=dataTree$edge.length[168]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_licheniformis
which(dataTree$tip.label=="Bacillus_subtilis")		# 88 (edge number of this tip)
which(dataTree$edge[,2] == 88)	  #  175 (row in edge matrix with this end node)
dataTree$edge.length[175]									# 0.093 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_licheniformis", edge.length=NULL, where=88, position=dataTree$edge.length[175]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_velezensis
which(dataTree$tip.label=="Bacillus_licheniformis")		# 89 (edge number of this tip)
which(dataTree$edge[,2] == 89)	  #  177 (row in edge matrix with this end node)
dataTree$edge.length[177]									# 0.045 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_velezensis", edge.length=NULL, where=89, position=dataTree$edge.length[177]/2) #4 if adding in first species, /2 if clustering further

## Bacillus_amyloliquefaciens
which(dataTree$tip.label=="Bacillus_velezensis")		# 90 (edge number of this tip)
which(dataTree$edge[,2] == 90)	  #  179 (row in edge matrix with this end node)
dataTree$edge.length[179]									# 0.023 (edge length)
dataTree <- bind.tip(dataTree, "Bacillus_amyloliquefaciens", edge.length=NULL, where=90, position=dataTree$edge.length[179]/2) #4 if adding in first species, /2 if clustering further

plot(dataTree)

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


## Burkholderia ----
## Burkholderia_glumae
which(dataTree$tip.label=="Burkholderia_cenocepacia")		# 40 (edge number of this tip)
which(dataTree$edge[,2] == 40)	  #  86 (row in edge matrix with this end node)
dataTree$edge.length[86]									# 0.15 (edge length)
dataTree <- bind.tip(dataTree, "Burkholderia_glumae", edge.length=NULL, where=40, position=dataTree$edge.length[86]/4) 

## Burkholderia_pseudomallei
which(dataTree$tip.label=="Burkholderia_cenocepacia")		# 40 (edge number of this tip)
which(dataTree$edge[,2] == 40)	  #  87 (row in edge matrix with this end node)
dataTree$edge.length[87]									# 0.037 (edge length)
dataTree <- bind.tip(dataTree, "Burkholderia_pseudomallei", edge.length=NULL, where=40, position=dataTree$edge.length[87]/2) 

## Burkholderia_multivorans
which(dataTree$tip.label=="Burkholderia_cenocepacia")		# 40 (edge number of this tip)
which(dataTree$edge[,2] == 40)	  #  88 (row in edge matrix with this end node)
dataTree$edge.length[88]									# 0.018 (edge length)
dataTree <- bind.tip(dataTree, "Burkholderia_multivorans", edge.length=NULL, where=40, position=dataTree$edge.length[88]/2) 

plot(dataTree)


## Campylobacter ----
## Campylobacter_fetus
which(dataTree$tip.label=="Campylobacter_jejuni")		# 55 (edge number of this tip)
which(dataTree$edge[,2] == 55)	  #  113 (row in edge matrix with this end node)
dataTree$edge.length[113]									# 0.27 (edge length)
dataTree <- bind.tip(dataTree, "Campylobacter_fetus", edge.length=NULL, where=55, position=dataTree$edge.length[113]/4) 

## Campylobacter_coli
which(dataTree$tip.label=="Campylobacter_jejuni")		# 55 (edge number of this tip)
which(dataTree$edge[,2] == 55)	  #  114 (row in edge matrix with this end node)
dataTree$edge.length[114]									# 0.06 (edge length)
dataTree <- bind.tip(dataTree, "Campylobacter_coli", edge.length=NULL, where=55, position=dataTree$edge.length[114]/2) 

plot(dataTree)


## Clostridium ----
## Clostridium_butyricum
which(dataTree$tip.label=="Clostridium_botulinum")		# 98 (edge number of this tip)
which(dataTree$edge[,2] == 98)	  #  194 (row in edge matrix with this end node)
dataTree$edge.length[194]									# 0.51 (edge length)
dataTree <- bind.tip(dataTree, "Clostridium_butyricum", edge.length=NULL, where=98, position=dataTree$edge.length[194]/4) 

## Clostridium_perfringens
which(dataTree$tip.label=="Clostridium_butyricum")		# 99 (edge number of this tip)
which(dataTree$edge[,2] == 99)	  #  196 (row in edge matrix with this end node)
dataTree$edge.length[196]									# 0.127 (edge length)
dataTree <- bind.tip(dataTree, "Clostridium_perfringens", edge.length=NULL, where=99, position=dataTree$edge.length[196]/2) 


# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


## Enterobacter
## Enterobacter_hormaechei
which(dataTree$tip.label=="Enterobacter_cloacae")		# 7 (edge number of this tip)
which(dataTree$edge[,2] == 7)	  #  33 (row in edge matrix with this end node)
dataTree$edge.length[33]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Enterobacter_hormaechei", edge.length=NULL, where=7, position=dataTree$edge.length[33]/4) 

## Enterobacter_kobei
which(dataTree$tip.label=="Enterobacter_cloacae")		# 7 (edge number of this tip)
which(dataTree$edge[,2] == 7)	  #  34 (row in edge matrix with this end node)
dataTree$edge.length[34]									# 0.0056 (edge length)
dataTree <- bind.tip(dataTree, "Enterobacter_kobei", edge.length=NULL, where=7, position=dataTree$edge.length[34]/2) 

## Enterobacter_roggenkampii
which(dataTree$tip.label=="Enterobacter_kobei")		# 8 (edge number of this tip)
which(dataTree$edge[,2] == 8)	  #  36 (row in edge matrix with this end node)
dataTree$edge.length[36]									# 0.002 (edge length)
dataTree <- bind.tip(dataTree, "Enterobacter_roggenkampii", edge.length=NULL, where=8, position=dataTree$edge.length[36]/2) 

## Enterobacter_asburiae
which(dataTree$tip.label=="Enterobacter_roggenkampii")		# 9 (edge number of this tip)
which(dataTree$edge[,2] == 9)	  #  38 (row in edge matrix with this end node)
dataTree$edge.length[38]									# 0.001 (edge length)
dataTree <- bind.tip(dataTree, "Enterobacter_asburiae", edge.length=NULL, where=9, position=dataTree$edge.length[38]/2) 


## Enterococcus ----
## Enterococcus_faecium
which(dataTree$tip.label=="Enterococcus_faecalis")		# 82 (edge number of this tip)
which(dataTree$edge[,2] == 82)	  #  167 (row in edge matrix with this end node)
dataTree$edge.length[167]									# 0.15 (edge length)
dataTree <- bind.tip(dataTree, "Enterococcus_faecium", edge.length=NULL, where=82, position=dataTree$edge.length[167]/4) 

## Enterococcus_hirae
which(dataTree$tip.label=="Enterococcus_faecium")		# 83 (edge number of this tip)
which(dataTree$edge[,2] == 83)	  #  169 (row in edge matrix with this end node)
dataTree$edge.length[169]									# 0.075 (edge length)
dataTree <- bind.tip(dataTree, "Enterococcus_hirae", edge.length=NULL, where=83, position=dataTree$edge.length[169]/2) 



## Escherichia ----
## Escherichia_albertii
which(dataTree$tip.label=="Escherichia_coli")		# 4 (edge number of this tip)
which(dataTree$edge[,2] == 4)	  #  30 (row in edge matrix with this end node)
dataTree$edge.length[30]									# 0.018 (edge length)
dataTree <- bind.tip(dataTree, "Escherichia_albertii", edge.length=NULL, where=4, position=dataTree$edge.length[30]/4) 

## Escherichia_fergusonii
which(dataTree$tip.label=="Escherichia_coli")		# 4 (edge number of this tip)
which(dataTree$edge[,2] == 4)	  #  31 (row in edge matrix with this end node)
dataTree$edge.length[31]									# 0.004 (edge length)
dataTree <- bind.tip(dataTree, "Escherichia_fergusonii", edge.length=NULL, where=4, position=dataTree$edge.length[31]/2) 

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


## Klebsiella ----
## Klebsiella_aerogenes
which(dataTree$tip.label=="Klebsiella_oxytoca")		# 15 (edge number of this tip)
which(dataTree$edge[,2] == 15)	  #  47 (row in edge matrix with this end node)
dataTree$edge.length[47]									# 0.026 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_aerogenes", edge.length=NULL, where=15, position=dataTree$edge.length[47]/4) 

## Klebsiella_michiganensis
which(dataTree$tip.label=="Klebsiella_oxytoca")		# 15 (edge number of this tip)
which(dataTree$edge[,2] == 15)	  #  48 (row in edge matrix with this end node)
dataTree$edge.length[48]									# 0.006 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_michiganensis", edge.length=NULL, where=15, position=dataTree$edge.length[48]/2) 

## Klebsiella_grimontii
which(dataTree$tip.label=="Klebsiella_michiganensis")		# 16 (edge number of this tip)
which(dataTree$edge[,2] == 16)	  #  50 (row in edge matrix with this end node)
dataTree$edge.length[50]									# 0.003 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_grimontii", edge.length=NULL, where=16, position=dataTree$edge.length[50]/2) 

## Klebsiella_variicola
which(dataTree$tip.label=="Klebsiella_aerogenes")		# 18 (edge number of this tip)
which(dataTree$edge[,2] == 18)	  #  53 (row in edge matrix with this end node)
dataTree$edge.length[53]									# 0.006 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_variicola", edge.length=NULL, where=18, position=dataTree$edge.length[53]/2) 

## Klebsiella_pneumoniae
which(dataTree$tip.label=="Klebsiella_variicola")		# 19 (edge number of this tip)
which(dataTree$edge[,2] == 19)	  #  55 (row in edge matrix with this end node)
dataTree$edge.length[55]									# 0.003 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_pneumoniae", edge.length=NULL, where=19, position=dataTree$edge.length[55]/2) 

## Klebsiella_quasipneumoniae
which(dataTree$tip.label=="Klebsiella_pneumoniae")		# 20 (edge number of this tip)
which(dataTree$edge[,2] == 20)	  #  57 (row in edge matrix with this end node)
dataTree$edge.length[57]									# 0.0016 (edge length)
dataTree <- bind.tip(dataTree, "Klebsiella_quasipneumoniae", edge.length=NULL, where=20, position=dataTree$edge.length[57]/2) 

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


## Lactobacillus ----
## Latilactobacillus_sakei
which(dataTree$tip.label=="Lacticaseibacillus_paracasei")		# 89 (edge number of this tip)
which(dataTree$edge[,2] == 89)	  #  181 (row in edge matrix with this end node)
dataTree$edge.length[181]									# 0.3 (edge length)
dataTree <- bind.tip(dataTree, "Latilactobacillus_sakei", edge.length=NULL, where=89, position=dataTree$edge.length[181]/4) 

## Lactobacillus_helveticus
which(dataTree$tip.label=="Lacticaseibacillus_paracasei")		# 89 (edge number of this tip)
which(dataTree$edge[,2] == 89)	  #  182 (row in edge matrix with this end node)
dataTree$edge.length[182]									# 0.078 (edge length)
dataTree <- bind.tip(dataTree, "Lactobacillus_helveticus", edge.length=NULL, where=89, position=dataTree$edge.length[182]/2) 

## Lacticaseibacillus_rhamnosus
which(dataTree$tip.label=="Lacticaseibacillus_paracasei")		# 89 (edge number of this tip)
which(dataTree$edge[,2] == 89)	  #  183 (row in edge matrix with this end node)
dataTree$edge.length[183]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Lacticaseibacillus_rhamnosus", edge.length=NULL, where=89, position=dataTree$edge.length[183]/2) 

## Ligilactobacillus_salivarius
which(dataTree$tip.label=="Latilactobacillus_sakei")		# 92 (edge number of this tip)
which(dataTree$edge[,2] == 92)	  #  187 (row in edge matrix with this end node)
dataTree$edge.length[187]									# 0.078 (edge length)
dataTree <- bind.tip(dataTree, "Ligilactobacillus_salivarius", edge.length=NULL, where=92, position=dataTree$edge.length[187]/2) 

## Lactiplantibacillus_plantarum
which(dataTree$tip.label=="Ligilactobacillus_salivarius")		# 93 (edge number of this tip)
which(dataTree$edge[,2] == 93)	  #  189 (row in edge matrix with this end node)
dataTree$edge.length[189]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Lactiplantibacillus_plantarum", edge.length=NULL, where=93, position=dataTree$edge.length[189]/2) 

## Limosilactobacillus_reuteri
which(dataTree$tip.label=="Lactiplantibacillus_plantarum")		# 94 (edge number of this tip)
which(dataTree$edge[,2] == 94)	  #  191 (row in edge matrix with this end node)
dataTree$edge.length[191]									# 0.019 (edge length)
dataTree <- bind.tip(dataTree, "Limosilactobacillus_reuteri", edge.length=NULL, where=94, position=dataTree$edge.length[191]/2) 

## Levilactobacillus_brevis
which(dataTree$tip.label=="Limosilactobacillus_reuteri")		# 95 (edge number of this tip)
which(dataTree$edge[,2] == 95)	  #  193 (row in edge matrix with this end node)
dataTree$edge.length[193]									# 0.001 (edge length)
dataTree <- bind.tip(dataTree, "Levilactobacillus_brevis", edge.length=NULL, where=95, position=dataTree$edge.length[193]/2) 

## Latilactobacillus_curvatus
which(dataTree$tip.label=="Latilactobacillus_sakei")		# 92 (edge number of this tip)
which(dataTree$edge[,2] == 92)	  #  188 (row in edge matrix with this end node)
dataTree$edge.length[188]									# 0.039 (edge length)
dataTree <- bind.tip(dataTree, "Latilactobacillus_curvatus", edge.length=NULL, where=92, position=dataTree$edge.length[188]/2) 

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


## Mycobacterium ----
## Mycobacterium_avium
which(dataTree$tip.label=="Mycobacteroides_abscessus")		# 84 (edge number of this tip)
which(dataTree$edge[,2] == 84)	  #  167 (row in edge matrix with this end node)
dataTree$edge.length[167]									# 0.3 (edge length)
dataTree <- bind.tip(dataTree, "Mycobacterium_avium", edge.length=NULL, where=84, position=dataTree$edge.length[167]/4) 

## Mycobacterium_intracellulare
which(dataTree$tip.label=="Mycobacterium_avium")		# 85 (edge number of this tip)
which(dataTree$edge[,2] == 85)	  #  169 (row in edge matrix with this end node)
dataTree$edge.length[169]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Mycobacterium_intracellulare", edge.length=NULL, where=85, position=dataTree$edge.length[169]/2) 


## Pseudomonas ----
## Pseudomonas_putida
which(dataTree$tip.label=="Pseudomonas_aeruginosa")		# 39 (edge number of this tip)
which(dataTree$edge[,2] == 39)	  #  88 (row in edge matrix with this end node)
dataTree$edge.length[88]									# 0.3 (edge length)
dataTree <- bind.tip(dataTree, "Pseudomonas_putida", edge.length=NULL, where=39, position=dataTree$edge.length[88]/4) 

## Pseudomonas_syringae
which(dataTree$tip.label=="Pseudomonas_putida")		# 40 (edge number of this tip)
which(dataTree$edge[,2] == 40)	  #  90 (row in edge matrix with this end node)
dataTree$edge.length[90]									# 0.08 (edge length)
dataTree <- bind.tip(dataTree, "Pseudomonas_syringae", edge.length=NULL, where=40, position=dataTree$edge.length[90]/2) 


## Rhodococcus ----
## Rhodococcus_pyridinivorans
which(dataTree$tip.label=="Prescottella_equi")		# 85 (edge number of this tip)
which(dataTree$edge[,2] == 85)	  #  170 (row in edge matrix with this end node)
dataTree$edge.length[170]									# 0.12 (edge length)
dataTree <- bind.tip(dataTree, "Rhodococcus_pyridinivorans", edge.length=NULL, where=85, position=dataTree$edge.length[170]/4) 

## Rhodococcus_erythropolis
which(dataTree$tip.label=="Rhodococcus_pyridinivorans")		# 86 (edge number of this tip)
which(dataTree$edge[,2] == 86)	  #  172 (row in edge matrix with this end node)
dataTree$edge.length[172]									# 0.03 (edge length)
dataTree <- bind.tip(dataTree, "Rhodococcus_erythropolis", edge.length=NULL, where=86, position=dataTree$edge.length[172]/2) 

## Rhodococcus_qingshengii
which(dataTree$tip.label=="Rhodococcus_erythropolis")		# 87 (edge number of this tip)
which(dataTree$edge[,2] == 87)	  #  174 (row in edge matrix with this end node)
dataTree$edge.length[174]									# 0.015 (edge length)
dataTree <- bind.tip(dataTree, "Rhodococcus_qingshengii", edge.length=NULL, where=87, position=dataTree$edge.length[174]/2) 

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]

## Shigella ----
## Shigella_sonnei
which(dataTree$tip.label=="Shigella_flexneri")		# 7 (edge number of this tip)
which(dataTree$edge[,2] == 7)	  #  35 (row in edge matrix with this end node)
dataTree$edge.length[35]									# 0.019 (edge length)
dataTree <- bind.tip(dataTree, "Shigella_sonnei", edge.length=NULL, where=7, position=dataTree$edge.length[35]/4) 

## Shigella_dysenteriae
which(dataTree$tip.label=="Shigella_sonnei")		# 8 (edge number of this tip)
which(dataTree$edge[,2] == 8)	  #  37 (row in edge matrix with this end node)
dataTree$edge.length[37]									# 0.004 (edge length)
dataTree <- bind.tip(dataTree, "Shigella_dysenteriae", edge.length=NULL, where=8, position=dataTree$edge.length[37]/2) 

## Shigella_boydii
which(dataTree$tip.label=="Shigella_dysenteriae")		# 9 (edge number of this tip)
which(dataTree$edge[,2] == 9)	  #  39 (row in edge matrix with this end node)
dataTree$edge.length[39]									# 0.002 (edge length)
dataTree <- bind.tip(dataTree, "Shigella_boydii", edge.length=NULL, where=9, position=dataTree$edge.length[39]/2) 


## Staphylococcus ----
## Staphylococcus_simulans
which(dataTree$tip.label=="Staphylococcus_aureus")		# 117 (edge number of this tip)
which(dataTree$edge[,2] == 117)	  #  234 (row in edge matrix with this end node)
dataTree$edge.length[234]									# 0.37 (edge length)
dataTree <- bind.tip(dataTree, "Staphylococcus_simulans", edge.length=NULL, where=117, position=dataTree$edge.length[234]/4) 

## Staphylococcus_saprophyticus
which(dataTree$tip.label=="Staphylococcus_aureus")		# 117 (edge number of this tip)
which(dataTree$edge[,2] == 117)	  #  235 (row in edge matrix with this end node)
dataTree$edge.length[235]									# 0.093 (edge length)
dataTree <- bind.tip(dataTree, "Staphylococcus_saprophyticus", edge.length=NULL, where=117, position=dataTree$edge.length[235]/2) 

## Staphylococcus_haemolyticus
which(dataTree$tip.label=="Staphylococcus_aureus")		# 117 (edge number of this tip)
which(dataTree$edge[,2] == 117)	  #  236 (row in edge matrix with this end node)
dataTree$edge.length[236]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Staphylococcus_haemolyticus", edge.length=NULL, where=117, position=dataTree$edge.length[236]/2) 

## Staphylococcus_epidermidis
which(dataTree$tip.label=="Staphylococcus_aureus")		# 117 (edge number of this tip)
which(dataTree$edge[,2] == 117)	  #  237 (row in edge matrix with this end node)
dataTree$edge.length[237]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Staphylococcus_epidermidis", edge.length=NULL, where=117, position=dataTree$edge.length[237]/2) 

## Staphylococcus_argenteus
which(dataTree$tip.label=="Staphylococcus_aureus")		# 117 (edge number of this tip)
which(dataTree$edge[,2] == 117)	  #  238 (row in edge matrix with this end node)
dataTree$edge.length[238]									# 0.04 (edge length)
dataTree <- bind.tip(dataTree, "Staphylococcus_argenteus", edge.length=NULL, where=117, position=dataTree$edge.length[238]/2) 

# Check which species aren't in tree
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]

## Vibrio ----
## Vibrio_campbellii
which(dataTree$tip.label=="Vibrio_cholerae")		# 36 (edge number of this tip)
which(dataTree$edge[,2] == 36)	  #  84 (row in edge matrix with this end node)
dataTree$edge.length[84]									# 0.11 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_campbellii", edge.length=NULL, where=36, position=dataTree$edge.length[84]/4) 

## Vibrio_parahaemolyticus
which(dataTree$tip.label=="Vibrio_campbellii")		# 37 (edge number of this tip)
which(dataTree$edge[,2] == 37)	  #  86 (row in edge matrix with this end node)
dataTree$edge.length[86]									# 0.02 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_parahaemolyticus", edge.length=NULL, where=37, position=dataTree$edge.length[86]/2) 

## Vibrio_alginolyticus
which(dataTree$tip.label=="Vibrio_parahaemolyticus")		# 38 (edge number of this tip)
which(dataTree$edge[,2] == 38)	  #  88 (row in edge matrix with this end node)
dataTree$edge.length[88]									# 0.01 (edge length)
dataTree <- bind.tip(dataTree, "Vibrio_alginolyticus", edge.length=NULL, where=38, position=dataTree$edge.length[88]/2) 


## Xanthamonas ----
## Xanthomonas_campestris
which(dataTree$tip.label=="Xanthomonas_oryzae")		# 57 (edge number of this tip)
which(dataTree$edge[,2] == 57)	  #  121 (row in edge matrix with this end node)
dataTree$edge.length[121]									# 0.12 (edge length)
dataTree <- bind.tip(dataTree, "Xanthomonas_campestris", edge.length=NULL, where=57, position=dataTree$edge.length[121]/4) 

## Xanthomonas_citri
which(dataTree$tip.label=="Xanthomonas_oryzae")		# 57 (edge number of this tip)
which(dataTree$edge[,2] == 57)	  #  122 (row in edge matrix with this end node)
dataTree$edge.length[122]									# 0.029 (edge length)
dataTree <- bind.tip(dataTree, "Xanthomonas_citri", edge.length=NULL, where=57, position=dataTree$edge.length[122]/2) 

## Xanthomonas_phaseoli
which(dataTree$tip.label=="Xanthomonas_citri")		# 58 (edge number of this tip)
which(dataTree$edge[,2] == 58)	  #  124 (row in edge matrix with this end node)
dataTree$edge.length[124]									# 0.014 (edge length)
dataTree <- bind.tip(dataTree, "Xanthomonas_phaseoli", edge.length=NULL, where=58, position=dataTree$edge.length[124]/2) 


## Yersinia ----
## Yersinia_ruckeri
which(dataTree$tip.label=="Yersinia_pestis")		# 29 (edge number of this tip)
which(dataTree$edge[,2] == 29)	  #  74 (row in edge matrix with this end node)
dataTree$edge.length[74]									# 0.07 (edge length)
dataTree <- bind.tip(dataTree, "Yersinia_ruckeri", edge.length=NULL, where=29, position=dataTree$edge.length[74]/4) 

## Yersinia_enterocolitica
which(dataTree$tip.label=="Yersinia_pestis")		# 29 (edge number of this tip)
which(dataTree$edge[,2] == 29)	  #  75 (row in edge matrix with this end node)
dataTree$edge.length[75]									# 0.018 (edge length)
dataTree <- bind.tip(dataTree, "Yersinia_enterocolitica", edge.length=NULL, where=29, position=dataTree$edge.length[75]/2) 

## Yersinia_pseudotuberculosis
which(dataTree$tip.label=="Yersinia_pestis")		# 29 (edge number of this tip)
which(dataTree$edge[,2] == 29)	  #  76 (row in edge matrix with this end node)
dataTree$edge.length[76]									# 0.01 (edge length)
dataTree <- bind.tip(dataTree, "Yersinia_pseudotuberculosis", edge.length=NULL, where=29, position=dataTree$edge.length[76]/2) 

# Check which species aren't in tree
# Should now say 'character(0)'
species_data$species[which((species_data$species %in% dataTree$tip.label)==FALSE)]


plot(dataTree)


## Remove 0 branch lengths by changing to 1 then remaking the tree ultrametric - repeat until no 0 branch lengths
dataTree$edge.length[dataTree$edge.length<=0]<-1
dataTree<-chronoMPL(dataTree)
which(dataTree$edge.length<=0)

summary(dataTree)
is.rooted(dataTree)

write.nexus(dataTree, file="SOC_plas_species_tree.nex")


