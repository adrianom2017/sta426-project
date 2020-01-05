set.seed(146)

plotData<-function(sce){
  plotColData(sce, x = "sum", y="detected",colour_by="cluster_id")
}
plotdist<-function(sce){
  par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
  hist(sce$detected, xlab="Detected Genes", 
       ylab="Number of cells", breaks=20, main="", col="blue")
}
plothighexpression<-function(sce){
  plotHighestExprs(sce, exprs_values = "counts")
}
plotVar<-function(vars){
  plotExplanatoryVariables(vars)
}
plotExpressionProfile<-function(sce){
  plotExpression(sce, rownames(sce)[1:6], x = "cluster_id",colour_by="cluster_id")
}

#Preprocessing
prep_steps <- function(sce){
  # Doublet removal
  sce <- scDblFinder(sce, samples = sce$sample_id)
  sce <- sce[, sce$scDblFinder.class != "doublet"]
  # 1. remove undetected genes
  sce <- sce[rowSums(counts(sce) > 0) > 0, ]
  # `scater` QC & filtering
  per.cell <- perCellQCMetrics(sce, subsets=list(Mito=grep("mt-", rownames(sce))))
  #print(plotData(sce))
  colData(sce) <- cbind(colData(sce), per.cell)
  #remove putative low-quality cells that have very few or many detected genes.
  #print(plotdist(sce))
  total.drop <- isOutlier(sce$total, nmads = 2, type = "both", log = TRUE, batch = sce$sample_id)
  detected.drop <- isOutlier(sce$detected, nmads=2, type="lower", log=TRUE)
  subsets_Mt_percent.drop <- isOutlier(sce$subsets_Mt_percent, nmads = 2, type = "higher") &   sce$subsets_Mt_percent > 0.08
  sce <- sce[, !(total.drop | detected.drop |   subsets_Mt_percent.drop)]
  data.frame(ByDetectedGenes=sum(detected.drop),ByTotal=sum(total.drop),Bysubsets_Mt_percent=sum(subsets_Mt_percent.drop),Remaining=ncol(sce))
  # Remove low-abundance/lowly expressed genes
  #print(plothighexpression(sce))
  sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
  # Variable-level QC metrics
  counts <- assay(sce, "counts")
  
  return(sce)
}