#!/usr/bin/env Rscript

#Load libraries
library(Seurat)
library(tidyverse)
library(SeuratWrappers)
library("Matrix")
library("readr")

#### Reading in output fresh for 50% mito cutoff filtering
#Rscript 00readoutput.R /data/paezbaenal2/projectOutputs/mligDavies/STARsolo/adultUnsorted/adultUnsortedSolo.out/Gene/filtered/ /data/paezbaenal2/projectOutputs/mligDavies/STARsolo/adultUnsorted/seuratoutputs/pctest/adultUns.rds adultUns 
#Rscript 00readoutput.R /data/paezbaenal2/projectOutputs/mligDavies/STARsolo/nucspotFACS/nucspotFACSSolo.out/Gene/filtered/ /data/paezbaenal2/projectOutputs/mligDavies/STARsolo/nucspotFACS/seuratoutputs/pctest/nucspotFACS.rds nucspotFACS

#Establish parameters: directory, path to generated seurat object
args_c = commandArgs(trailingOnly = TRUE)
dir = args_c[1]
seu_path = args_c[2]
proj = args_c[3]

#Read in manually
mtx <- readMM(paste0(dir, "matrix.mtx"))
rownames(mtx) <- read_tsv(paste0(dir, "features.tsv"), col_names=FALSE)[, 1, drop=TRUE]
colnames(mtx) <- read_tsv(paste0(dir, "barcodes.tsv"), col_names=FALSE)[, 1, drop=TRUE]

#Create seurat object from manual read in method
seu <- CreateSeuratObject(counts=mtx, min.cells=3, min.feature=200, project = proj)
print(Idents(seu))
saveRDS(seu, file = seu_path)  # seu.rds is saved

print(seu)
head(seu)