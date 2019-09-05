## ----setup, include=FALSE------------------------------------------------
setwd("~/PROJECTS_WORK/132")
library(knitr)
opts_chunk$set(message=FALSE, error=FALSE, warning=FALSE, tidy=TRUE, tidy.opts=list(blank=FALSE), cache=TRUE, cache.lazy = FALSE, echo = TRUE,  results = 'asis', fig.pos='H', fig.wide=TRUE)
library(tidyverse)
library(reshape2)
library(pander)
library(grid)
library(VennDiagram)
library(DESeq2)
library(ImpulseDE2)
library(splines)
library(pheatmap)
library(TCseq)
library(biomaRt)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(pander)
library(eulerr)
library(qvalue)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
# ddsHTSeq <- readRDS("XXX_ddsHTSeq_timeCourse_spline.RDS")
# sampleTable <- readRDS("XXX_sampleTable_metadata.RDS")
# objectImpulseDE2 <- readRDS("XXX_objectImpulseDE2.RDS")
load("XXX_RNAseq_time-course_analysis_20190830.RData")


## ----dataLoadCluster, eval=FALSE-----------------------------------------
## ### this is run on the cluster while the rest is run locally
## directory <- "/usr/share/sequencing/projects/132/deseq"
## sampleFiles <- grep("Bio",list.files(directory),value=TRUE)
## 
## ### build metadata table
## ### starting from name - 132_Bio-1-D_wk-1_S15.htseq.a5.txt
## 
## sampleTable <- data.frame(sampleName = sub("^132_(Bio-.+)-.{1}_wk-.{1}_S.+.txt$",
##                                            "\\1", sampleFiles),
##                           fileName = sampleFiles,
##                           condition = ifelse(grepl("-D_", sampleFiles), "DMEM", "XVIVO"),
##                           time = sub("^132_Bio-.+-.{1}_wk-(.{1})_S.+.txt$",
##                                            "\\1", sampleFiles),
##                           stringsAsFactors = FALSE)
## ## sample names need to be unique
## sampleTable$sampleName <- paste0(sampleTable$sampleName, "_",
##                                  sampleTable$condition, "_",
##                                  sampleTable$time)
## sampleTable <- numerize(sampleTable, c("time"))


## ----splineCluster, eval=FALSE-------------------------------------------
## matTimeSplineBasis <- ns(
##     sampleTable$time, df=4)
## colnames(matTimeSplineBasis) <-
##     paste0("spline", seq(1, dim(matTimeSplineBasis)[2]))


## ----splineSamples, eval=FALSE-------------------------------------------
## sampleTable <- cbind(sampleTable, matTimeSplineBasis)


## ----countLoadCluster, eval=FALSE----------------------------------------
## ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
##                                        directory = directory,
##                                        design = ~condition + condition:spline1 + condition:spline2 +
##                                          condition:spline3 + condition:spline4)


## ----preFiltering, eval=FALSE--------------------------------------------
## keep <- rowSums(counts(ddsHTSeq)) >= 10
## ddsHTSeq <- ddsHTSeq[keep,]
## ddsHTSeq$condition <- relevel(ddsHTSeq$condition, ref = "DMEM")


## ----dispersion, eval=FALSE----------------------------------------------
## ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
## ddsHTSeq <- estimateDispersions(ddsHTSeq)


## ----DELRT, eval=FALSE---------------------------------------------------
## ddsHTSeq <- nbinomLRT(ddsHTSeq,
##                       full = ~condition + condition:spline1 + condition:spline2 +
##                                          condition:spline3 + condition:spline4,
##                  reduced = ~spline1 + spline2 + spline3 + spline4)


## ----resSummary, results='asis'------------------------------------------
res <- results(ddsHTSeq)
res@elementMetadata@listData$description[2]
res@elementMetadata@listData$description[5]


## ----resOrdered----------------------------------------------------------
resOrdered <- as.data.frame(res[order(res$pvalue),])
resOrdered$gene <- row.names(resOrdered)
pandoc.table(head(resOrdered, 10L)[c(2,3,6)])


## ----MAplot, fig.cap="MA Plot. A significant amount of dispersion can be noticed"----
plotMA(res, ylim=c(-4,4))


## ----contrastNames-------------------------------------------------------
resultsNames(ddsHTSeq)


## ----shrinkage, eval=FALSE-----------------------------------------------
## resLFC <- lfcShrink(ddsHTSeq, coef="conditionXVIVO.spline4", type="apeglm")


## ----MAplotShrinked, fig.cap="Shrinked Fold Change. MA plot after running an apeglm fold change shirnkage"----
plotMA(resLFC, ylim=c(-4,4))


## ----extractNormalised---------------------------------------------------
vsd <- vst(ddsHTSeq, blind=FALSE)
normalisedCounts <- assay(vsd)


## ----PCAvsd--------------------------------------------------------------
plotPCA(vsd, intgroup=c("condition", "time"))


## ----heatMap-------------------------------------------------------------
select <- order(rowMeans(counts(ddsHTSeq,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsHTSeq)[,c("condition","time")])
pheatmap(normalisedCounts[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


## ----geneCoordinates, eval=FALSE-----------------------------------------
## genes <- row.names(normalisedCounts)
## genesUniverseCoordinates <- readRDS("~/CODE/core/Rutilities/geneSymbols_coordinates_biomart.RDS")
## 
## ## removing non standard chromosomes
## 
## genesUniverseCoordinates <- genesUniverseCoordinates %>%
##   filter(!grepl("CHR", chr))
## 
## ### there will be several symbols in this list which are not in the database
## ### for some reasons
## 
## exprGenesCoordinates <- data.frame(symbol = genes)
## exprGenesCoordinates <- exprGenesCoordinates %>%
##   left_join(genesUniverseCoordinates, by = "symbol")
## 
## ## removing duplicated entries
## exprGenesCoordinates = exprGenesCoordinates[!duplicated(exprGenesCoordinates$symbol),]
## 
## ## fixing those not present in database with a dummy coordinate
## exprGenesCoordinates[is.na(exprGenesCoordinates$chr),"chr"]<-1
## exprGenesCoordinates[is.na(exprGenesCoordinates$start),"start"]<-0
## exprGenesCoordinates[is.na(exprGenesCoordinates$end),"end"]<-1
## 
## ## tcseq needs specific identifiers
## names(exprGenesCoordinates)<-c("id", "chr", "start", "end")
## 


## ----tcSamples-----------------------------------------------------------
samplesAnnotationTC <- sampleTable[,c(1,3,4)]
names(samplesAnnotationTC) <- c("sampleid", "group", "timepoint")



## ----tcaSetup------------------------------------------------------------
tca <- TCA(
  design = samplesAnnotationTC,
  counts = counts(ddsHTSeq),
  genomicFeature = exprGenesCoordinates
)


## ----tcacomputation------------------------------------------------------
tca <- DBanalysis(tca)


## ----timecoursetable-----------------------------------------------------
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE)


## ----timeclustering------------------------------------------------------
tca <- timeclust(tca, algo = "cm", k = 10, standardize = TRUE)


## ----clustplot, include=FALSE--------------------------------------------
clustplot <- timeclustplot(tca, value = "z-score(PRKM)", cols = 5)


## ----clustprint, fig.cap="Gene Clusters. Clusters of genes identified by fuzzy clustering of their time course trend.", fig.height=15, fig.width=10----
grid.newpage()
pushViewport(viewport(layout=grid.layout(5,2)))
vplayout<-function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
print(clustplot[[1]], vp=vplayout(1,1))
print(clustplot[[2]], vp=vplayout(1,2))
print(clustplot[[3]], vp=vplayout(2,1))
print(clustplot[[4]], vp=vplayout(2,2))
print(clustplot[[5]], vp=vplayout(3,1))
print(clustplot[[6]], vp=vplayout(3,2))
print(clustplot[[7]], vp=vplayout(4,1))
print(clustplot[[8]], vp=vplayout(4,2))
print(clustplot[[9]], vp=vplayout(5,1))
print(clustplot[[10]], vp=vplayout(5,2))


## ----membership----------------------------------------------------------
clusters <- data.frame(
  gene = names(clustResults(tca)@cluster),
  cluster = unname(clustResults(tca)@cluster)
)
clusteredDEresults <- as.data.frame(resOrdered) %>%
  left_join(clusters, by = "gene")
write_tsv(clusteredDEresults, "XXX_results_DESeq2_spline_clusters.txt")


## ----DEclusters----------------------------------------------------------
pandoc.table(
  clusteredDEresults %>%
  filter(padj < 1e-3 & abs(log2FoldChange)>2) %>%
  group_by(cluster) %>%
  summarise(sigGenes = n())
)


## ----clusterThreePlot----------------------------------------------------
print(clustplot[[3]])


## ----impulseCalc, eval=FALSE---------------------------------------------
## ## also impulse needs specific column names for the annotations
## sampleAnnoationImpulse <- samplesAnnotationTC
## names(sampleAnnoationImpulse) <- c("Sample", "Condition", "Time")
## ## Time needs to be numeric
## sampleAnnoationImpulse$Time <- as.numeric(as.character(sampleAnnoationImpulse$Time))
## ## Then case/control needs to be named as such (very picky...)
## ## We set DMEM as control and XVIVO as cases
## sampleAnnoationImpulse$Condition <- ifelse(sampleAnnoationImpulse$Condition == "DMEM", "control", "case")
## 
## objectImpulseDE2 <- runImpulseDE2(
##   matCountData           = counts(ddsHTSeq),
##   dfAnnotation           = sampleAnnoationImpulse,
##   boolCaseCtrl           = TRUE,
##   boolIdentifyTransients = TRUE,
##   scaNProc               = 1 )


## ----impulseResults------------------------------------------------------
impulseResults <- objectImpulseDE2$dfImpulseDE2Results[order(objectImpulseDE2$dfImpulseDE2Results$padj),]
write_tsv(impulseResults, "XXX_results_ImpulseDE2.txt")
pandoc.table(head(impulseResults, 10L)[c(4,3,17)])


## ----ImpDeseqOverlap-----------------------------------------------------
methodsOverlap <- euler(
  list(
    A=unique(impulseResults$Gene[impulseResults$padj<1e-04]),
    B=unique(resOrdered$gene[resOrdered$padj<1e-04])
    ))
plot(methodsOverlap,
     quantities = TRUE,
     legend = list(labels = c("ImpulseDE2", "DESeq2")))


## ----impulseGenes--------------------------------------------------------
lsgplotsGenes <- plotGenes(
  vecGeneIDs       = NULL,
  scaNTopIDs       = 10,
  objectImpulseDE2 = objectImpulseDE2,
  boolCaseCtrl     = TRUE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL)


## ----impulseGenesPlots, fig.cap="Top10 most significant. Top significant genes resulting from the ImpulseDE2 analysis, showing monotonous trend over time compared to controls.", fig.height=15, fig.width=10----
grid.newpage()
pushViewport(viewport(layout=grid.layout(5,2)))
vplayout<-function(x,y) viewport(layout.pos.row=x, layout.pos.col=y)
print(lsgplotsGenes[[1]], vp=vplayout(1,1))
print(lsgplotsGenes[[2]], vp=vplayout(1,2))
print(lsgplotsGenes[[3]], vp=vplayout(2,1))
print(lsgplotsGenes[[4]], vp=vplayout(2,2))
print(lsgplotsGenes[[5]], vp=vplayout(3,1))
print(lsgplotsGenes[[6]], vp=vplayout(3,2))
print(lsgplotsGenes[[7]], vp=vplayout(4,1))
print(lsgplotsGenes[[8]], vp=vplayout(4,2))
print(lsgplotsGenes[[9]], vp=vplayout(5,1))
print(lsgplotsGenes[[10]], vp=vplayout(5,2))


## ----impulseResultsTransient---------------------------------------------
impulseTransient <- impulseResults %>%
  filter(isTransient=="TRUE") %>%
  arrange(padj, impulseTOsigmoid_padj)
pandoc.table(head(impulseTransient))


## ----impulseTransientPlots, fig.cap="Top10 most significant transient genes. Top significant genes with a transient expression time course, resulting from the ImpulseDE2 analysis.", fig.height=15, fig.width=10----
lsgplotsTransientGenes <- plotGenes(
  vecGeneIDs       = head(impulseTransient$Gene, 10L),
  objectImpulseDE2 = objectImpulseDE2,
  boolCaseCtrl     = TRUE,
  dirOut           = NULL,
  strFileName      = NULL,
  vecRefPval       = NULL, 
  strNameRefMethod = NULL)
grid.newpage()
pushViewport(viewport(layout=grid.layout(5,2)))
print(lsgplotsTransientGenes[[1]], vp=vplayout(1,1))
print(lsgplotsTransientGenes[[2]], vp=vplayout(1,2))
print(lsgplotsTransientGenes[[3]], vp=vplayout(2,1))
print(lsgplotsTransientGenes[[4]], vp=vplayout(2,2))
print(lsgplotsTransientGenes[[5]], vp=vplayout(3,1))
print(lsgplotsTransientGenes[[6]], vp=vplayout(3,2))
print(lsgplotsTransientGenes[[7]], vp=vplayout(4,1))
print(lsgplotsTransientGenes[[8]], vp=vplayout(4,2))
print(lsgplotsTransientGenes[[9]], vp=vplayout(5,1))
print(lsgplotsTransientGenes[[10]], vp=vplayout(5,2))


## ----universe------------------------------------------------------------
#### creating the Universe from ENSEMBL
universe <- keys(org.Hs.egSYMBOL)
# all symbols
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
universeSymbols <- unlist(as.list(x[mapped_genes]))

### universe can also be created more appropriately from the analysis
universeFromAnalysis = unlist(mget(impulseResults$Gene[!is.na(impulseResults$Gene)], org.Hs.egSYMBOL2EG, ifnotfound=NA))

### however, one must get rid of those symbols without easy correspondence to entrezIDs
universeFromAnalysisCleaned <- universeFromAnalysis[which(universeFromAnalysis!="NA")]



## ----geneLists-----------------------------------------------------------

sigImpulseAll <- bitr(impulseResults %>%
                        filter(padj<1e-04) %>%
                        pull(Gene), 
                      fromType = "SYMBOL", 
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = org.Hs.eg.db)


sigImpulseTransient <- bitr(impulseResults %>%
                        filter(padj<1e-03) %>%
                        filter(isTransient=="TRUE") %>%
                        pull(Gene), 
                      fromType = "SYMBOL", 
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = org.Hs.eg.db)



## ----geneListFoldchange--------------------------------------------------
sigImpulseAll <- sigImpulseAll %>%
  left_join(resOrdered %>%
              dplyr::select(gene,log2FoldChange),
            by = c("SYMBOL" = "gene")) %>%
  filter(abs(log2FoldChange) > 2)

sigImpulseAllList <- sigImpulseAll$log2FoldChange
names(sigImpulseAllList) <- sigImpulseAll$ENTREZID
sigImpulseAllList <- sort(sigImpulseAllList, decreasing = TRUE)

sigImpulseTransient <- sigImpulseTransient %>%
  left_join(resOrdered %>%
              dplyr::select(gene,log2FoldChange),
            by = c("SYMBOL" = "gene"))

sigImpulseTransientList <- sigImpulseTransient$log2FoldChange
names(sigImpulseTransientList) <- sigImpulseTransient$ENTREZID
sigImpulseTransientList <- sort(sigImpulseTransientList, decreasing = TRUE)


## ------------------------------------------------------------------------
sigImpulseAllBP <- enrichGO(gene=unique(sigImpulseAll$ENTREZID),
                              OrgDb = 'org.Hs.eg.db',
                              ont="BP",
                              universe = universeFromAnalysisCleaned,
                              qvalueCutoff = 0.3,
                              readable=T)
# head(sigImpulseAllBP@result[c(1:7)], 50L)
# emapplot(sigImpulseAllBP)

sigImpulseAllBPX <- setReadable(sigImpulseAllBP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(sigImpulseAllBPX)
enrichplot::upsetplot(sigImpulseAllBPX, 8)



## ------------------------------------------------------------------------
sigImpulseTransientegoBP <- enrichGO(gene=unique(sigImpulseTransient$ENTREZID),
                              OrgDb = 'org.Hs.eg.db',
                              ont="BP",
                              universe = universeFromAnalysisCleaned,
                              qvalueCutoff = 0.3,
                              readable=T)
head(sigImpulseTransientegoBP@result[c(1:7)], 50L)
emapplot(sigImpulseTransientegoBP)

sigImpulseTransientegoBPX <- setReadable(sigImpulseTransientegoBP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(sigImpulseTransientegoBPX)
enrichplot::upsetplot(sigImpulseTransientegoBP, 8)



sigImpulseTransientegoKEGG <- enrichKEGG(gene=unique(sigImpulseTransient$ENTREZID),
                              organism = 'hsa',
                              pvalueCutoff = 0.05)


sigImpulseTransientGSEA <- gseGO(geneList = sigImpulseTransientList,
              OrgDb = org.Hs.eg.db,
              ont = "BP",
              nPerm = 10000,
              minGSSize = 10,
              maxGSSize = 1000,
              pvalueCutoff = 0.05)


## ----save, eval=FALSE----------------------------------------------------
## save.image("XXX_RNAseq_time-course_analysis_20190830.RData")

