## ----setup, include=FALSE------------------------------------------------
setwd("C:/Users/tbleazar/OneDrive - MHRA/Documents/R/rnaseq/272")
library(knitr)
opts_chunk$set(message=FALSE, error=FALSE, warning=FALSE, tidy=TRUE, tidy.opts=list(blank=FALSE), cache=TRUE, cache.lazy = FALSE, echo = TRUE,  results = 'asis', fig.pos='H', fig.wide=TRUE)
library(tidyverse)
library(reshape2)
library(pander)
library(grid)
library(VennDiagram)
library(DESeq2)
#BiocManager::install("ImpulseDE2")
library(ImpulseDE2)
library(splines)
library(pheatmap)
#BiocManager::install("TCseq")
library(TCseq)
library(biomaRt)
library(RSQLite)
library(backports)
library(digest)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(eulerr)
#BiocManager::install("qvalue")
library(qvalue)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(DOSE)
library(tidyr)
library(ggplot2)
library(scales)
library(dplyr)
library(Glimma)
library(limma)
library(apeglm)
library(gplots)


directory <- "X:/272/htseq"
sampleFiles <- grep(".htseq.out", list.files(directory), value=TRUE)

sampleFilesFrame <- data.frame(sampleFiles)
sampleFilesBits <- separate(data=sampleFilesFrame, col=sampleFiles, sep="_", into=c("272","t","sample","rest"))
sampleFilesBits[8,"sample"] <- "24cont_a"
sampleFilesBits[9,"sample"] <- "24cont_b"
sampleFilesBits[10,"sample"] <- "24cont_c"
sampleFilesBits[11,"sample"] <- "24cont_d"
sampleNamesFiltered <- sampleFilesBits$sample[c(1:7,12:32)]
sampleFilesFiltered <- sampleFiles[c(1:7,12:32)]
sampleTimesFiltered <- c(0,0,0,0,24,24,24,24,2,2,2,2,0.5,0.5,0.5,0.5,3,3,3,3,1,1,1,1,6,6,6,6)

sampleTable <- data.frame(sampleName=sampleNamesFiltered, fileName=sampleFilesFiltered, condition=rep("case", times=length(sampleNamesFiltered)), time=sampleTimesFiltered, stringsAsFactors=FALSE)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~1)
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
ddsHTSeq <- estimateDispersions(ddsHTSeq)

#Interactive html of MDS plot
glMDSPlot(ddsHTSeq)

samplesAnnotationTC <- sampleTable[,c(1,3,4)]
names(samplesAnnotationTC) <- c("Sample", "Condition", "Time")

#ImpulseDE2
objectImpulseDE2 <- runImpulseDE2(matCountData=counts(ddsHTSeq), dfAnnotation=samplesAnnotationTC, boolCaseCtrl=FALSE, vecConfounders=NULL, boolIdentifyTransients=TRUE, scaNProc=3)
impulseResults <- objectImpulseDE2$dfImpulseDE2Results[order(objectImpulseDE2$dfImpulseDE2Results$padj), ]
write.csv(impulseResults, "272-ImpulseDE2.csv")

#PCA Plot
vsd <- vst(ddsHTSeq, blind=TRUE)
normalisedCounts <- assay(vsd)
plotPCA(vsd, intgroup=c("time")) + 
  scale_colour_gradientn(colors=c("blue", "red", "white"), values=rescale(c(0,6,24)), limits=c(0,24)) +
  labs(color="Time (h)")

#Heatmap
distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- sampleNamesFiltered
pdf("272-heatmap.pdf", width=14, height=14)
heatmap.2(mat, trace = "none", margin = c(22,22))
dev.off()

#Gene co-ordinates for TCseq
genes <- row.names(normalisedCounts)
genesUniverseCoordinates <- readRDS("geneSymbols_coordinates_biomart.RDS")
## removing non standard chromosomes
genesUniverseCoordinates <- genesUniverseCoordinates %>% filter(!grepl("CHR",chr))
### there will be several symbols in this list which are not in the database
### for some reasons
exprGenesCoordinates <- data.frame(symbol = genes)
exprGenesCoordinates <- exprGenesCoordinates %>% left_join(genesUniverseCoordinates, by = "symbol")
## removing duplicated entries
exprGenesCoordinates = exprGenesCoordinates[!duplicated(exprGenesCoordinates$symbol), ]
## fixing those not present in database with a dummy coordinate
exprGenesCoordinates[is.na(exprGenesCoordinates$chr), "chr"] <- 1
exprGenesCoordinates[is.na(exprGenesCoordinates$start), "start"] <- 0
exprGenesCoordinates[is.na(exprGenesCoordinates$end), "end"] <- 1
## tcseq needs specific identifiers
names(exprGenesCoordinates) <- c("id", "chr", "start", "end")

#TCseq for gene clustering
names(samplesAnnotationTC) <- c("sampleid", "group", "timepoint")
samplesAnnotationTC$timepoint <- as.character(samplesAnnotationTC$timepoint)
tca <- TCA(design = samplesAnnotationTC, counts = counts(ddsHTSeq), genomicFeature = exprGenesCoordinates)
tca <- DBanalysis(tca, categories=c("timepoint"))
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE)
tca <- timeclust(tca, algo = "cm", k = 8, standardize = TRUE)
clustplot <- timeclustplot(tca, value = "z-score(PRKM)", cols = 5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(5, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(clustplot[[1]], vp = vplayout(1, 1))
print(clustplot[[2]], vp = vplayout(1, 2))
print(clustplot[[3]], vp = vplayout(2, 1))
print(clustplot[[4]], vp = vplayout(2, 2))
print(clustplot[[5]], vp = vplayout(3, 1))
print(clustplot[[6]], vp = vplayout(3, 2))
print(clustplot[[7]], vp=vplayout(4,1))
print(clustplot[[8]], vp=vplayout(4,2))

## ----membership----------------------------------------------------------
clusters <- data.frame(
  Gene = names(clustResults(tca)@cluster),
  cluster = unname(clustResults(tca)@cluster)
)
clusters <- left_join(clusters, impulseResults, by="Gene")
clusters <- clusters[clusters$padj<0.05, c("Gene", "cluster", "padj")]
write.csv(clusters, file="272-cluster-8-members.csv")

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

sigImpulseAll <- bitr(impulseResults %>%
                        filter(padj<1e-04) %>%
                        pull(Gene), 
                      fromType = "SYMBOL", 
                      toType = c("ENSEMBL", "ENTREZID"),
                      OrgDb = org.Hs.eg.db)

clusterfour <- bitr(clusters %>%
                      filter(cluster==4) %>%
                      pull(Gene), 
                    fromType = "SYMBOL", 
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

clustersix <- bitr(clusters %>%
                      filter(cluster==6) %>%
                      pull(Gene), 
                    fromType = "SYMBOL", 
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

clusterseven <- bitr(clusters %>%
                      filter(cluster==7) %>%
                      pull(Gene), 
                    fromType = "SYMBOL", 
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

## ------------------------------------------------------------------------
sigImpulseAllBP <- enrichGO(gene=unique(sigImpulseAll$ENTREZID),
                            OrgDb = 'org.Hs.eg.db',
                            ont="BP",
                            universe = universeFromAnalysisCleaned,
                            qvalueCutoff = 0.3,
                            readable=T)
sigImpulseAllBPX <- setReadable(sigImpulseAllBP, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(sigImpulseAllBPX)
enrichplot::upsetplot(sigImpulseAllBPX, 8)

enrichfour <- enrichGO(gene=unique(clusterfour$ENTREZID),
                            OrgDb = 'org.Hs.eg.db',
                            ont="BP",
                            universe = universeFromAnalysisCleaned,
                            qvalueCutoff = 0.3,
                            readable=T)
enrichfourX <- setReadable(enrichfour, 'org.Hs.eg.db', 'ENTREZID')
enrichplot::upsetplot(enrichfourX, 8)

enrichsix <- enrichGO(gene=unique(clustersix$ENTREZID),
                       OrgDb = 'org.Hs.eg.db',
                       ont="BP",
                       universe = universeFromAnalysisCleaned,
                       qvalueCutoff = 0.3,
                       readable=T)
enrichsixX <- setReadable(enrichsix, 'org.Hs.eg.db', 'ENTREZID')
enrichplot::upsetplot(enrichsixX, 8)

enrichseven <- enrichGO(gene=unique(clusterseven$ENTREZID),
                       OrgDb = 'org.Hs.eg.db',
                       ont="BP",
                       universe = universeFromAnalysisCleaned,
                       qvalueCutoff = 0.3,
                       readable=T)
enrichsevenX <- setReadable(enrichseven, 'org.Hs.eg.db', 'ENTREZID')
enrichplot::upsetplot(enrichsevenX, 8)
