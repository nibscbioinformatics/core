
writeLines("

######################################
# HELP message #
######################################

this function performs a hypergeometric test of sets enrichment
based on list of names only
the function returns (and prints) a matrix object with the results
of the tests.

syntax:
results <- hypergeoEnrichment(geneList, setsDataFrame, genomeGenes)
where:

A - geneList
the list of elements you want to test the enrichment for, i.e. your list of 
significant genes, or variants, or group of things you want to test if enriched for 
some set

B - setsDataFrame
this is a dataframe with two columns
one with the name of the set (this will be repeated for every element in the set)
one with the element member of the set (gene, variant, anything)
in this way the data frame will contain all sets you want to test the enrichment
of your group for.

C - genomeGenes
this object is a vector, and represents the Universe, i.e. all genes, or all variants,
or all elements existing. it is important to create the null distribution of your test.

")




hypergeoEnrichment <- function(geneList, setsDataFrame, genomeGenes){
  
  sets <- unique(setsDataFrame$set)
  writeLines("Hypergeometric test of Enrichment")
  header <- paste("Set Name","Count in Set", "Set total", "unadjusted p-value")
  writeLines(header)
  results <- data.frame()
  
  for (set in sets){
        
    setData <- setsDataFrame[setsDataFrame$set == set,]
    setGenes <- setData$gene
    setLength <- length(setGenes)
    
    #writeLines(set)
    #head(setData)
    
    # calculate values
    listInSet <- length(which(geneList %in% setGenes))
    genomeInSet <- length(which(genomeGenes %in% setGenes))
    genomeNotInSet <- length(genomeGenes) - genomeInSet
    listSize <- length(geneList)
    
    pval<-phyper(
      listInSet-1,
      genomeInSet,
      genomeNotInSet,
      listSize,
      lower.tail=F
    )
    
    results <- rbind(results, cbind(set, listInSet, setLength, pval))
  }
  names(results)<-c("SetName", "AnnotatedInSet", "SetSize", "pvalue")
  print(results)
  return(results)
}

