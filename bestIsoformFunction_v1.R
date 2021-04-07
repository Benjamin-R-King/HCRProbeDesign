# To run this function, you need:
# 1) a directory called "Data" with: an RNA seq data set, a list of genes in a .csv
# more info on the appropriate formatting of these is available in the read me

bestIsoform = function(GenePath, DataPath, Name){
  # Check to make sure there is an "Output" and "Graphics" directory
  if (!dir.exists("Output")){
    dir.create("Output") 
  }
  if (!dir.exists("Graphics")){
    dir.create("Graphics") 
  }
  
  # These are the genes that I am interested in
  GeneList <- read_csv(GenePath, col_names=TRUE)
  ids <- GeneList$ensembl_gene_id
  graphics_path <- paste0(getwd(),"/Graphics/", Name)
  print(GeneList)
  
  # This is all the RNA seq data for all the patients (n ~ 430)
  SeqData <- read_tsv(DataPath, col_names=TRUE)
  
  # The first column of the RNA seq data set is a composite of a lot of different accession ids (ensembl, havana, etc.)
  # These are separated by |
  # Separate this string at the | character and save each in its own column of a tibble
  # Some list manipulations are necessary to wrangle the tibble into the correct form
  Info <- str_split(SeqData$X1, fixed("|"))
  Info <- data.table::transpose(Info)
  Info <- as_tibble(Info, .name_repair = "unique")
  Info <- dplyr::select(Info,c(1,2,6,7,8))
  
  # Naming the columns on this newly created tibble
  names(Info)[1:5] <- c("ensembl_isoform","ensembl_gene","symbol", "length", "type")
  print(Info)
  
  # I didn't comment this line when I wrote it, I can't remember what it does, if anything
  SeqData <- select_if(SeqData, is.numeric)
  
  # Make output vectors
  geneVector <- c()
  isoformVector <- c()
  symbolVector <- c()
  lengthVector <- c()
  rankVector <- c()
  sumTPMvector <- c()
  
  # Data wrangling to make a data frame that has the data that I want, in  the format I want
  
  # iterate over each ensembl id
  for (g in seq_along(ids)){
    gene <- ids[g]
    # Which rows match ensembl_id[g]?
    indices <- which(substr(Info$ensembl_gene, 1, 15) == gene)
    # Creating an empty variable that will be updated below
    bestIndex <- NULL
    # "rank" is a column from the input csv file, it's really just an arbitrary
    # number. Could consider removing this in a future version
    rank <- GeneList$Rank[g]
    
    # an if test. As long as there were matching ensembl ids found, proceed
    if (!is_empty(indices)){
      # Make empty vectors. They have a slot for each isoform (row) that matched the ensembl id
      sumVectorPlot <- vector("double", length=length(indices))
      isoformVectorPlot <- vector("character", length=length(indices))
      symbolVectorPlot <- vector("character", length=length(indices))
      lengthVectorPlot <- vector("double", length=length(indices))
      typeVectorPlot <- vector("character", length=length(indices))
      
      # Now start to fill these vectors that I just created with relevant data
      # iterate over each isoform in the indices vector
      for (i in seq_along(indices)){
        # using the slice function to just grab the row from that isoform
        # sum the TPM from all patients for that isoform
        sumVectorPlot[i] <- sum(dplyr::slice(SeqData,indices[i]))
        isoformVectorPlot[i] <- substr(Info$ensembl_isoform[indices[i]],1,15)
        symbolVectorPlot[i] <- Info$symbol[indices[i]]
        lengthVectorPlot[i] <- Info$length[indices[i]]
        typeVectorPlot[i] <- Info$type[indices[i]]
        
        # Need to know what the best isoform is
        # Our criteria for "best" is the one with highest sum(TPM)
        # You could change this criterium
        # bestIndex is empty when you encounter the first isoform, so initially
        # assign bestIndex as the first isoform
        if (i == 1){
          bestIndex <- indices[i]
          # creating a temp variable previousSum to keep track of the current
          # "best isoform"
          previousSum <- sumVectorPlot[i]
        }
        
        # After the first isoform has been processed, do this block of code
        else {
          currentSum <- sumVectorPlot[i]
          
          # Check: does the current isoform have a higher sum(TPM) than the previous
          # best isoform? if it does, overwrite the variables bestIndex and previousSum
          if (currentSum > previousSum){
            bestIndex <- indices[i]
            previousSum <- currentSum
          }
        }
      }
      # At this point we have iterated over each isoform, have filled all of the 
      # data vectors, and have decided which of the isoforms is the best
      
      
      # Adding info to the vectors that will be used to construct the .csv output below
      geneVector <- c(geneVector,substr(Info$ensembl_gene[bestIndex], 1, 15))
      isoformVector <- c(isoformVector,substr(Info$ensembl_isoform[bestIndex],1,15))
      symbolVector <- c(symbolVector, Info$symbol[bestIndex])
      lengthVector <- c(lengthVector, Info$length[bestIndex])
      rankVector <- c(rankVector, rank)
      
      bestName <- substr(Info$ensembl_isoform[bestIndex],1,15)
      bestName <- paste0("Best Isoform: ",bestName)
      geneSymbol <- symbolVectorPlot[1]
      
      # Creating a tibble that will be used to construct the summary plots for 
      # each gene
      tempTibble <- tibble(isoform=isoformVectorPlot, type = typeVectorPlot, symbol = symbolVectorPlot, sum_transcripts=sumVectorPlot)
      
      # this first plot is sumTPM on y axis and each isoform on x axis, bar graph
      # color of the bar is the type of transcript "protein coding", etc...
      # plotX <- ggplot(tempTibble, aes(x= isoform, y=sum_transcripts, fill=type)) + geom_bar(stat="identity") +
      #   ggtitle(paste0(geneSymbol,"\n",gene,"\n",bestName)) + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) + theme_classic()

      plotX <- ggplot(tempTibble, aes(x= isoform, y=sum_transcripts, fill=type)) + geom_bar(stat="identity") +
        ggtitle(paste0(geneSymbol,"\n",gene,"\n",bestName)) + theme_classic() + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
        scale_x_discrete(expand=c(0,0)) +
        
        # ylim(0,max(tempTibble$sum_transcripts)) +
        theme(panel.grid = element_blank(),panel.border = element_blank())
      
      
      # Managing output directories
      # Checking to see if the graphics path exists; if it doesn't, create it
      if (!dir.exists(graphics_path)){
        dir.create(graphics_path) 
      }
      # Save the first plot in that directory
      ggsave(paste0(rank,"_",geneSymbol,"_Best.png"), width = 8, height = 4, path=graphics_path)
      
      # Make another plot 
      plotX <- ggplot(tempTibble, aes(x= isoform, y=sum_transcripts, fill=type)) + geom_bar(stat="identity") +
        ggtitle(paste0(geneSymbol,"\n",gene,"\n",bestName)) + theme_classic() + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
        scale_x_discrete(expand=c(0,0)) +
        
        # ylim(0,max(tempTibble$sum_transcripts)) +
        theme(panel.grid = element_blank(),panel.border = element_blank())
      
      # print(plotX)
      ggsave(paste0(rank,"_",geneSymbol,"_Best.pdf"), width = 8, height = 4, path=graphics_path)
      
    }
    
    # carrying out this block of code only if we found isoforms that matched
    # the ENSEMBL ID that we are iterating on
    if (!is_empty(indices)){
      # Subset the part of the RNAseq tibble that is for 
      tempTibble <- slice(SeqData,indices)
      tempIsoformTibble <- slice(SeqData,bestIndex)
      tempTibble <- summarize_all(tempTibble, .funs=sum)
      sumTPM <- sum(tempTibble)
      sumTPMvector <- c(sumTPMvector,sumTPM)
      # print(ncol(tempTibble))
      # print(ncol(tempIsoformTibble))
      # print(tempTibble)
      # print(tempIsoformTibble)
      normalizeTibble <- as_tibble(tempIsoformTibble/tempTibble)
      normalizeTibble <- t(normalizeTibble)
      normalizeTibble <- as_tibble(normalizeTibble)
      
      # print(normalizeTibble)
      
      # Creating a histogram showing the the best isoform only for each gene
      # a value of 1 on the x axis in this histogram represents an observation
      # of a patient where all reads were assigned to this isoform
      # a value of 0 is an observation where none of the reads for that patient were
      # assigned to this isoform
      # This gives an idea of how "dominant" this isoform is
      plotX <- ggplot(normalizeTibble, aes(x= V1)) + geom_histogram(bins=30) +
        ggtitle(paste0(geneSymbol,"\n",gene,"\n",bestName)) + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) + theme_classic()
     
      # print(plotX)
      ggsave(paste0(rank,"_",geneSymbol,"_Histogram.png"), width = 6, height = 4, path=graphics_path)
      
      plotX <- ggplot(normalizeTibble, aes(x= V1)) + geom_histogram(bins=30) +
        ggtitle(paste0(geneSymbol,"\n",gene,"\n",bestName)) + theme_classic() + theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
        coord_cartesian(xlim=c(0,1))
      
      # print(plotX)
      ggsave(paste0(rank,"_",geneSymbol,"_Histogram.pdf"), width = 6, height = 4, path=graphics_path)
      
    }
    
  }
  
  # Creating a tibble with all the information about what the best isoform is
  # Using the vectors that were built above
  Output <- tibble(gene_rank=rankVector,ensembl_gene=geneVector,ensembl_isoform=isoformVector,symbol=symbolVector,length=lengthVector,sum_TPM=sumTPMvector)
  
  output_path <- paste0(getwd(),"/Output/", Name)
  if (!dir.exists(output_path)){
    dir.create(output_path) 
  }
  output_path <- paste0(output_path,"/",Name, "_best_isoforms.csv")
  write.csv(Output, file = output_path, row.names = FALSE)
  
  return(Output)
}