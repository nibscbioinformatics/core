#!/bin/Rscript

args <- commandArgs(TRUE)

library(rmarkdown)

markdown<-args[1]
report <-args[2]
results<-args[-c(1,2)]

render(input = markdown, output_file = report, params = list(data = results))