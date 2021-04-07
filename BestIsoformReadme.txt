Readme Doc for bestIsoformFunction_v1.R
2021-03-11

Directory Structure:
My code expects the project to have its own directory with multiple subdirectories

Project Directory (Parent Directory):
/ProbeDesign

Subdirectories:
/ProbeDesign/Data
	* This subdirectory will house (1) the gene list and the (2) RNA seq dataset, more info on these below

/ProbeDesign/R
	* The R code lives in this folder


Necessary Files:

(1) A list of Genes, including ENSEMBL IDs saved as a .CSV file. An example file is provided in the repository and is called "geneListForBestIsoformFunction.csv"

(2) An RNA seq data set. I have used the OV and PAAD TCGA RNA seq data sets which were analyzed to give isoform level transcript counts.
These data can be downloaded from a pilot study using google cloud for RNA-seq analysis using data from TCGA and CCLE 

https://osf.io/gqrz9/files

Information from the repository regarding who did this analysis and how:
Google Cloud Pilot RNA-Sequencing for CCLE and TCGA
Contributors: PJ Tatlow
Date created: 2016-07-06 12:59 PM | Last Updated: 2016-07-19 11:32 AM
Category:  Project
Description: We processed over 10,000 RNA-Seqeuencing samples from the Cancer Cell Line Encyclopedia and The Cancer Genome Atlas using kallisto
License: CC-By Attribution 4.0 International

I downloaded the files at this location:
for HGSOC
matrices/TCGA/TCGA_OV_tpm.tsv.gz

for pancreatic cancer
matrices/TCGA/TCGA_PAAD_tpm.tsv.gz

You can download these or other datasets in this location (or from the analogous matrices/CCLE/ directory)


** Important**
The bestIsoformFunction_v1.R code is written to expect RNA seq data sets in this format. If you use data sets from other sources, either wrangle it into this format OR adjust the bestIsoformFunction_v1.R code to recognize the new formatting

