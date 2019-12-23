
#data, n_treated, n_untreated, samples_treated, samples_untreated, n_comp / or auto, logFC magnitude + portion / list

create_dataset = function(sce, n_comp, n_cells, kNN, kNN_subsample, n_samples, logFC){
  
  samples = unique(colData(sce)$sample_id)
  
  #Sanity checks
  if(length(samples) < n_samples){
    stop("Number of samples in original data set smaller than n_samples")
  }
  
  if(length(n_cells) == 1){
    n_cells = rep(n_cells, length(samples))
  }else{
    if(!(length(n_cells) == n_samples)){
      stop("Lenght of n_cells vector not the same as n_samples")
    }
  }
  
  if(!typeof(logFC) == "list"){
    stop("logFC argument no list")
  }else if(length(logFC) == 2){
    if(logFC[[2]] > 1){
      stop("Specified proportion in logFC is lager than 1. Invalid.")
    }
    
    if(logFC[[1]] == 0){
      stop("Specified magnitude of logFC is 0")
    }
  }
  else if(length(logFC) != 1 | length(logFC) != 2){
    stop("length of logFC list argument has to be 2 (with magnitue and proportion) or 1 (vector of logFCs)")
  }
  
  #TODO
  sce_sim = SingleCellExperiment(list(logcounts_sim = matrix(nrow = nrow(sce), ncol = 0)))
  colD = data.frame(matrix(nrow = 0, ncol = 3))
  names(colD) = c("cluster_id","sample_id","cell_id")
  
  for(i in 1:length(samples)){
    sample_sce = create_sample(sce[, colData(sce)$sample_id == samples[i]], n_comp, n_cells[i], kNN, kNN_subsample)
    sce_sim = SingleCellExperiment(list(logcounts_sim = cbind(assay(sce_sim, 'logcounts_sim'), assay(sample_sce, 'logcounts_sim'))))
    colD = rbind(colD, colData(sample_sce))
  }
  
  colData(sce_sim) = cbind(colData(sce_sim),colD)
  
  #Compute dispersion of genes
  dispresion = compute_dispresion(sce)
  
  #Compute logFC for each gene
  lFC = compute_logFC(sce, logFC)
  
  #Create SingleCellExperiment object
  rowD = data.frame(gene_id = rownames(sce), logFC = lFC, dispersion = dispersion)
  rowData(sce_sim) = cbind(rowData(sce_sim), rowD)
  
  #Add counts
  #TODO how to trainsfom
  tmp = 2^logcounts(sce_sim) - 1
  tmp[tmp < 0 ] = 0
  counts(sce_sim) = tmp
  
  #Compute the counts from NB
  sce_sim = nb_counts(sce_sim)
  
  return(sce_sim)
}

nb_counts = function(sce_sim){
  
  #TODO incoorperate library.size
  ds = rep(rowData(sce_sim)$dispersion, each = ncol(sce_sim)) 
  
  mu = c(counts(sce_sim)) * 2^rowData(sce_sim)$logFC
  
  counts =  rnbinom(n=nrow(sce_sim)*ncol(sce_sim), size=1/ds, mu=mu)
  counts = matrix(counts, nrow = nrow(sce_sim),ncol = ncol(sce_sim))
  assay(sce_sim, 'counts_sim') = counts
  return(sce_sim)
}

compute_logFC = function(sce, logFC){

  if(length(logFC) == 1){
    return(logFC[[1]])
  }else{
    l = rep(0, nrow(sce))
    n_DE_genes = floor(nrow(sce)*logFC[[2]])
    l[sample(1:nrow(sce), size = n_DE_genes)] = sample(c(logFC[[1]], 1/logFC[[1]]), size = n_DE_genes, replace = TRUE)
    return(l)
  }
  
}

compute_dispresion = function(sce){
  dds <- DGEList(counts(sce))
  mm <- model.matrix(~1, data=as.data.frame(colData(sce)))
  dds <- estimateDisp(dds, mm)
  dispersion <- dds$tagwise.dispersion
  return(dispersion)
}

create_sample = function(sce, n_comp, n_cells, kNN, kNN_subsample){
  
  n.cells.per.types = n_cells_per_types(sce, n_cells)
  cell_types = unique(colData(sce)$cluster_id)
  
  #Initialise matrix
  sample_expression_matrix = matrix(nrow = nrow(sce), ncol = 0)
  for(i in 1:length(n.cells.per.types)){
    show(cell_types[i])
    if(!n.cells.per.types[i] == 0){
      tmp = expression_matrix(sce[,colData(sce)$cluster_id == cell_types[i]], n_comp, n.cells.per.types[i], kNN, kNN_subsample)
      sample_expression_matrix = cbind(sample_expression_matrix, tmp)}
  }
  
  colD = data.frame(cluster_id = rep(cell_types, n.cells.per.types),
                    sample_id = colData(sce)$sample_id[1],
                    cell_id = colnames(sample_expression_matrix))
  
  #Create a SingleCellExperiment object with the mean logcounts for this sample
  sample_sce = SingleCellExperiment(list(logcounts_sim = sample_expression_matrix))
  colData(sample_sce) =  cbind(colData(sample_sce), colD)
  return(sample_sce)
}

#Computes the number of cells we need to sample for each cell type
n_cells_per_types = function(sce, n_cells_per_sample){
  #Vector of weights per cell type
  w = c(table(colData(sce)$cluster_id ))
  
  #number of differenet cell types in sample
  n_cell_types = length(unique(colData(sce)$cluster_id))
  
  cells = sample(c(1:n_cell_types), size = n_cells_per_sample,  prob = w, replace = TRUE)
  return(tabulate(cells))
}

expression_matrix = function(sce, n_comp, n_cells_per_celltype, kNN, kNN_subsample){

  if(!c("logcounts") %in% assayNames(sce)){
    stop("logcounts not defined. Call logNormCounts(sce) to add them to your sce")
  }
  
  if(n_cells_per_celltype <= dim(sce)[2]){
    idx_cells = sample(c(1:dim(sce)[2]), n_cells_per_celltype, replace = FALSE)
  }else{
    idx_cells = sample(c(1:dim(sce)[2]), n_cells_per_celltype, replace = TRUE)
  }
  
  #PCA
  pca_data <- prcomp_irlba(t(logcounts(sce)),n = n_comp)
  
  #Create distance matrix
  distance_matrix<-as.matrix(dist(pca_data$x,method="euclidean"))
  diag(distance_matrix) = Inf
  
  #Sort distance matric by row
  sorted_distances <- apply(distance_matrix,1,function(x) sort(x, index.return = T)$ix)
  
  #Sample cells and average PC's loading of the N nearest neighbours
  w_clustered <- sapply(idx_cells, function (idx) {
    d = sorted_distances[idx,1:kNN]
    cell_neighbors <- sample(d, kNN_subsample, replace = TRUE)
    colMeans(pca_data$x[cell_neighbors,])
  })
  
  w_clustered = t(w_clustered)
  
  #TODO check transpose https://stackoverflow.com/questions/29783790/how-to-reverse-pca-in-prcomp-to-get-original-data
  expression_matrix <- t(w_clustered  %*% t(pca_data$rotation)) + pca_data$center
  colnames(expression_matrix) = rownames(colData(sce))[idx_cells]
  return(expression_matrix)
}