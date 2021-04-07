# Load the source code for the function bestIsoform
source("R/bestIsoformFunction_v1.R")

# The arguments of this function take: 
# 1) a csv with your gene list, 
# 2) the tcga dataset you want to use, 
# 3) and a string that you want to name the project

# The function saves graphs of each ENSEMBL Id and an output .csv file
# The the data frame that is output as a .csv is also returned to the variable
# "genes" if you want to look at the data in the R console
genes <- bestIsoform("Data/geneListForBestIsoformFunction.csv", "Data/TCGA_OV_tpm.tsv", "geneList_test")


# Load the source code for the function makeCodebook
source("R/makeCodebook_v1.R")

# Format the list of ENSEMBL Gene IDs and ENSEMBL Transcript IDs for the Matlab MerFISH code
makeCodebook("Output/geneList_test/geneList_test_best_isoforms.csv","geneList_test")



