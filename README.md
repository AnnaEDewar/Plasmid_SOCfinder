# Genes for cooperation are not more likely to be carried by plasmids
Authors: Anna Dewar, Laurence Belcher, Thomas Scott, Stuart West

Affiliation: Department of Biology, University of Oxford, United Kingdom

## Overview
This repository contains code and data for the research article 'Genes for cooperation are not more likely to be carried by plasmids'. For any queries please contact Anna Dewar at anna.dewar@biology.ox.ac.uk.

## Supplementary Material 1

S1 contains full details of all results. The document was compiled from an Rmarkdown file. The original .Rmd file, including code for all models and results in the S1 document, is available in this repository as: 'Code_S1.Rmd'. The code is also available as a regular R file: 'Code_S1.R'.

### Data & tree
The .Rmd file requires the data file 'SOC_summary_data_wide_146.csv', which contains all data used in our analyses, and the phylogenies 'SOC_plas_species_tree.nex' and 'species_gtdb.tree'. For ease of use we recommend saving all files, including the data, tree and .Rmd file, into a folder alongside an Rstudio project. 

### Compiling 'Code_S1.Rmd' as a pdf
If you wish to locally compile the file 'Code_S1.pdf' into a PDF identical to the S1 document, please also download the file 'figure_order_header.tex' and include this within the same folder as the .Rmd file.

## Phylogenies
### Supplementary Material 2: supertree phylogeny
To build our phylogeny, we used a recently published maximum likelihood tree generated with 16S ribosomal protein data as the basis for our phylogeny (Hug et al. (2016), 'A new view of the tree of life', Nat. Microbiol.). We used the R package ‘ape’ to identify all branches that matched either a species or a genus in our dataset. In cases where we had multiple species within a single genus, we used the R package ‘phytools’ to add these species as additional branches in the tree. We used published phylogenies from the literature to add any within-genus clustering of species’ branches (details and references of these phylogenies are available in 'Supplementary_material_2.pdf' and 'species_tree_refs.xlsx').

The code for how did this is available to download as 'tree_script.R'. It requires the files 'Original_tree.txt' and 'species_146.csv'. Each line of the script edits the tree to produce the final tree, so the lines of code should be run in the order they are in the script.

### GTDB tree
We also used the GTDB bacterial reference tree as an alternative phylogeny (version 214.1). The code we used to make a subset of the tree with only species clusters corresponding to our genomes is in: 'GTDB_tree_script.R'. It requires: 'SOC_summary_data_wide_146.csv', 'bac120_r214.tree' and 'bac120_metadata_r214.tsv' [note: this file is too large to upload to GitHub, but is available to download from GTDB at: https://data.gtdb.ecogenomic.org/releases/release214/214.1/.

# Scripts to use SOCfinder on large genomic datasets
Finally, we have included a html guide with all of the code we used to conduct our analysis with SOCfinder, including selecting and downloading genomes, and extracting data. We hope this will be a helpful example of how SOCfinder could be used in future studies. This is available as: 'plas_chr_SOCfinder_analysis.html'.
