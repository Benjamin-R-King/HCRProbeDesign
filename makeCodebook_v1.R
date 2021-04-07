# A function to format the output of the function "bestIsoform" into a .csv
# that will be recognized by the Matlab MerFISH probe design code

# input_table is the output of the function bestIsoform
makeCodebook <- function(input_table,codebook_name){
  # open the input_table .csv file
  input_table <- read_csv(input_table)
  
  # Check to see if there is a codebook output directory. If not, make one
  if (!dir.exists("Output/Codebook")){
    dir.create("Output/Codebook") 
  }
  
  # Open a blank text file
  codebook <- file(description=paste0("Output/Codebook/",codebook_name,"_codebook.txt"), open="w")
  # Write the header lines
  writeLines(paste0("version,","1"),codebook)
  writeLines(paste0("codebook_name,",codebook_name),codebook)
  writeLines(paste("bit_names", "RS0015", "RS0083","RS0095","RS0109","RS0175", "RS0237","RS0247", "RS0255","RS0307","RS0332","RS0343","RS0384","RS0406","RS0451","RS0468","RS0548",sep=","),codebook)
  writeLines("name,id,barcode",codebook)
  
  # Write the ENSEMBL Id, ENSEMBL Transcript, barcode (irrelevant for our purposes, but required by MerFISH code),
  # and Gene Symbol
  for(i in 1:nrow(input_table)){
    print(i)
    print(Genes[i,])
    writeLines(paste(input_table[i,2],input_table[i,3],"110000000001010", input_table[i,4], sep=","), codebook)
  }
  # Close the new text file
  close(codebook)
}
