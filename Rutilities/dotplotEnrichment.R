dotplotEnrichment <- function(enrichResultsTable, totTestedGenes, maxToDisplay = 10, pvalTreshold = 0.05){
  # geneRatio is the ratio of significant genes belonging to the set (k) compared to total significant genes (n)
  # symbols k and n refers to the standard hypergeometric distribution where N is background and M is setSize
  # i.e. how many of N belong to the Set
  enrichResultsTable$AnnotatedInSet = as.numeric(as.character(enrichResultsTable$AnnotatedInSet))
  enrichResultsTable$GeneRatio = enrichResultsTable$AnnotatedInSet / totTestedGenes
  if(pvalTreshold) {
    enrichResultsTable = enrichResultsTable[enrichResultsTable$padjust <= pvalTreshold,]
  }
  # limit results to established threshold if any
  enrichResultsTable = enrichResultsTable[1:maxToDisplay,]
  enrichResultsTable = enrichResultsTable[!is.na(enrichResultsTable$SetName),]
  idx <- order(enrichResultsTable$GeneRatio, decreasing=FALSE)
  enrichResultsTable$SetName <- factor(enrichResultsTable$SetName, levels = unique(enrichResultsTable$SetName[idx]))
  plot <- ggplot(enrichResultsTable, aes(x=GeneRatio, y=SetName, size=AnnotatedInSet, color=padjust))+
    geom_point()+
    scale_color_gradient(low="red", high="blue") +
    labs(x = "Gene Ratio", title = "Enrichment results")+ 
    theme(plot.title=element_text(hjust = 0.5, 
                                  size =14, 
                                  face = "bold"))
  return(plot)
}