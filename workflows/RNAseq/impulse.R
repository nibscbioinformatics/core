#! /bin/Rscript

writeLines("########## loading packages")

setwd("/usr/share/sequencing/projects/XXX/analysis/time-course")
library(DESeq2)
library(ImpulseDE2)
library(tidyverse)

writeLines("############### loading data")

load("PXXX_RNAseq_time-course_analysis_20190830.RData")

## also impulse needs specific column names for the annotations
sampleAnnoationImpulse <- samplesAnnotationTC
names(sampleAnnoationImpulse) <- c("Sample", "Condition", "Time")
## Time needs to be numeric
sampleAnnoationImpulse$Time <- as.numeric(as.character(sampleAnnoationImpulse$Time))
## Then case/control needs to be named as such (very picky...)
## We set DMEM as control and XVIVO as cases
sampleAnnoationImpulse$Condition <- ifelse(sampleAnnoationImpulse$Condition == "DMEM", "control", "case")


writeLines("############## running impulse2")

objectImpulseDE2 <- runImpulseDE2(
  matCountData           = counts(ddsHTSeq),
  dfAnnotation           = sampleAnnoationImpulse,
  boolCaseCtrl           = TRUE,
  boolIdentifyTransients = TRUE,
  scaNProc               = 1 )


writeLines("############### impulse done")
writeLines("############### saving impulse only as RDS")

saveRDS(objectImpulseDE2, "PXXX_objectImpulseDE2.RDS")

writeLines("############### script done")
