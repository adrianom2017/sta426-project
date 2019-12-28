#Adriano Martinelli, adrianom@student.ethz.ch
#Script which installs packages
#Inspired by renku lab script

#Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Update to given version
if(BiocManager::version() != '3.10'){
  BiocManager::install(version='3.10', update=TRUE, ask=FALSE, quiet=TRUE)}

#BiocManager::install('remotes', quiet=TRUE)
#remotes::install_github('HelenaLC/muscat', quiet=FALSE)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(gridExtra)
library(grid)
library(edgeR)
library(irlba)
library(muscat)
pkgs <- c("SingleCellExperiment",
          "scater",
          "ggplot2","gridExtra","grid",
          "edgeR",
          "irlba",
          "muscat")

BiocManager::install(pkgs, update=FALSE, ask=FALSE, quiet=TRUE)
