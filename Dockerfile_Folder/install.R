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
pkgs <- c("SingleCellExperiment",
          "scater",
          "ggplot2","gridExtra","grid",
          "edgeR",
          "irlba",
          "countsimQC",
          "tidyverse",
          "TMB",
          "cowplot",
          "scDblFinder",
          "dplyr",
          "UpSetR",
          "uwot",
          "Rtsne"
          )

BiocManager::install(pkgs, update=FALSE, ask=FALSE, quiet=TRUE)
remotes::install_github("adrianom2017/sta426-project/sta426")
#remotes::install_github("HelenaLC/muscat", ref="dream2", quiet = TRUE)
