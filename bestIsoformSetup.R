# If you don't have Bioconductor yet, you will need to run this block
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install("Biostrings")
BiocManager::install("biomaRt")


# Need these packages
library("Biostrings")
library("biomaRt")
library("tidyverse")

# To make the files for IDT orders, the files need to be in .xlsx format
install.packages("writexl")
library("writexl")
