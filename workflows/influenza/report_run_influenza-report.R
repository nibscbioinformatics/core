#!/bin/Rscript

args <- commandArgs(TRUE)

library(rmarkdown)

workdir <- getwd()

markdown<-args[1]
report <-args[2]
results<-args[-c(1,2)]

render(input = markdown, output_file = report, params = list(data = results), knit_root_dir=workdir, output_dir=workdir)