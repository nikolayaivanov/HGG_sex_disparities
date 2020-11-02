#
# Prepubescent -----------------
rm(list=ls())
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_cases.rda")
DMPs_prepubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_prepubescent_controls.rda")
DMPs_prepubescent_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_prepubescent_cases, 
	controls=DMPs_prepubescent_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename= '/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_prepubescent.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_prepubescent_ORA.csv')

#Postpubescent -----------------
rm(list=ls())
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_cases.rda")
DMPs_postpubescent_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_postpubescent_controls.rda")
DMPs_postpubescent_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_postpubescent_cases, 
	controls=DMPs_postpubescent_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename= '/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename= '/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_postpubescent.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_postpubescent_ORA.csv')

#0-21 yo -----------------
rm(list=ls())
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_cases.rda")
DMPs_allPeds_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_allPeds_controls.rda")
DMPs_allPeds_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_allPeds_cases, 
	controls=DMPs_allPeds_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_allPeds.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_allPeds_ORA.csv')

#>21 yo -----------------
rm(list=ls())
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

#load genes
load('/athena/masonlab/scratch/users/nai2008/gencode_v35_hg19_genes_and_genomicState.rda')

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_cases.rda")
DMPs_adult_cases=all_dmps

load("/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/rdas/annotated_DMPs_adult_controls.rda")
DMPs_adult_controls=all_dmps

lost_sex_DMPs_interrogation(cases=DMPs_adult_cases, 
	controls=DMPs_adult_controls, 
	lost_DMPs_overlapping_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults.csv', 
	lost_DMPs_overlapping_genes_enrichR_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults_for_enrichR.csv',
	lost_DMPs_overlapping_imprinted_genes_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_imprinted_genes_adults.csv',
	genes=genes,
	ORA_table_filename='/athena/masonlab/scratch/users/nai2008/DNAm/third_round_of_scripts/data_tables/lost_DMPs_overlapping_genes_adults_ORA.csv')



