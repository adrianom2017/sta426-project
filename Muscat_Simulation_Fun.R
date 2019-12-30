#Aggregation of single-cell to pseudobulk data
compute.pd <- function(sim,assay,fun){
  pb <- aggregateData(sim,assay = assay, fun = fun,
                      by = c("cluster_id", "sample_id"))
  
  return(pb)
}


Myfunc <- function(data, true_DE_genes) {
  DE_genes_res <- (true_DE_genes %in% data["Ground_Truth"]) & (DE_genes %in% data[method_name])
}



DS.analysis.pd <- function(sim,assay,fun,ds_method,topnumber,gene_info_pos){
  
  #Aggregation of single-cell to pseudobulk data for sum counts
  pb <- compute.pd(sim,assay,fun) 
  # Pseudobulk-level multidimensional scaling (MDS) plot. Each point represents a cluster-sample instance; points are colored by cluster ID and shaped by group ID
  pb_mds <- pbMDS(pb)
  
  # run edgeR/limma-voom on pseudobulks
  res <- pbDS(pb, method = ds_method, verbose = FALSE)
  # reformat results
  tbl <- resDS(sim, res, bind = "col")  
  
  # filter & sort
  #tbl <- dplyr::filter(tbl, p_adj.loc < 0.05,logFC > 1)
  
  tbl <- arrange(tbl, tbl$p_adj.loc)

  
  sim_df<- metadata(sim)[gene_info_pos]

  if (gene_info_pos == 3){
    sim_map <-select(sim_df$gene_info, "gene","cluster_id","category")
    #print("metadata(sim)$gene_info$logFC")
    #print(length(metadata(sim)$gene_info$logFC)) #19122
  }
  else{
    sim_map <-select(sim_df$gene_info2, "gene","cluster_id","category")
    #print("metadata(sim)$gene_info2$logFC")
    #print(length(metadata(sim)$gene_info2$logFC)) #19122
  }
  
  
  #Incomporate category type into tbl object
  tbl <- merge(tbl,sim_map,by=c("gene","cluster_id"))
  
  # no. of DS genes per cluster
  res_by_k <- split(tbl, tbl$cluster_id)
  vapply(res_by_k, nrow, numeric(1))
  # top hits in each cluster
  top <- do.call("rbind", lapply(res_by_k, head, 3))
  top <- select(top, -c("contrast", "p_adj.glb"))
  top$gene <- gsub("^.*\\.", "", top$gene)
  format(data.frame(top, row.names = NULL), digits = 3)
  
  tbl$category <- lapply(tbl$category, function(x) {gsub("de", "1", x) })
  tbl$category <- lapply(tbl$category, function(x) {gsub("ee", "0", x) })
  
  tbl$category <- as.integer(unlist(tbl$category))
  

  ## Empty COBRAData object:
  COBRAData()
  
  method_name = paste0(assay,"-",fun,"-",ds_method) 
  pval <- as.data.frame(tbl$p_val)
  colnames(pval) = method_name
  padj <- as.data.frame(tbl$p_adj.loc)
  colnames(padj) = method_name
  truth = as.data.frame(as.integer(unlist(tbl$category)))
  colnames(truth) = method_name
  
  print("TBL")
  print(tbl)
  
  tbl$category_method[tbl$p_val < 0.05] <- 1
  tbl$category_method[tbl$p_val >= 0.05] <- 0
  
  DE_genes <- tbl$gene[tbl$category_method == 1]
  print(length(DE_genes))
  #print(DE_genes)
  EE_genes <- tbl$gene[tbl$category_method == 0]
  print(length(EE_genes))
  #print(EE_genes)
  
  true_DE_genes <- tbl$gene[tbl$category == 1]
  print(length(true_DE_genes))
  #print(true_DE_genes)
  true_EE_genes <- tbl$gene[tbl$category == 0]
  print(length(true_EE_genes))
  #print(true_EE_genes)
  
  # example of list input (list of named vectors)
  #listInput <- list(Method = c(DE_genes,EE_genes), Truth  =c(true_DE_genes,true_EE_genes) )
  
  # Create an UpsetR Plot
  data <-tbl %>% select(contains("category"))
  print(head(data,5))
  print("END test")
  colnames(data) <- c("Ground_Truth",method_name)
  

  #print(upset(data,main.bar.color = "black", mainbar.y.label = "DE/EE gene Intersections", 
  #      sets.x.label = "genes per Category",query.legend = "bottom",queries =  list(list(query = Myfunc, 
  #                    params = list(0), color = "orange", active = T),
  #                    list(query = intersects, params = list(method_name), color = "red", active = F), 
  #                    list(query = intersects,params = list("Ground_Truth"), color = "yellow", active = F),
  #                    list(query = intersects,params = list(method_name, "Ground_Truth"), color = "green",active = T))))
  
  
 # print(upset(data,main.bar.color = "black", mainbar.y.label = "DE/EE gene Intersections", 
#    sets.x.label = "genes per Category",query.legend = "bottom",
#    queries =  list( list(query = intersects, params = list(method_name), color = "red", active = F), 
#              list(query = intersects,params = list("Ground_Truth"), color = "yellow", active = F),
 #           list(query = intersects,params = list(method_name, "Ground_Truth"), color = "green",active = T))))


  
  #print(upset(tbl$category),sets = levels(sim$gene))
  #print(upset(data.frame(tbl$category_method)),sets = levels(sim$gene))
  
  return(list(sim = sim,res = res ,tbl = tbl,res_by_k = res_by_k, pval  = pval,padj = padj, truth = truth ))
  
}



compute.mm <- function(sim,ds_method, vst){
  res <- muscat::mmDS(sim, method = ds_method, vst = vst)
  return(res)
}

DS.analysis.mm <- function(sim,ds_method,type,vst,topnumber,gene_info_pos){
  
  res <- compute.mm(sim,ds_method,vst)
  
  tbl <- rbind(res$Neuronal_excit,res$Neuronal_inhib)
  tbl <- arrange(tbl, tbl$p_adj.loc)
  
  sim_df<- metadata(sim)[gene_info_pos]
  if (gene_info_pos == 3){
    sim_map <-select(sim_df$gene_info, "gene","cluster_id","category")
  }
  else{
    sim_map <-select(sim_df$gene_info2, "gene","cluster_id","category")
  }
  
  #Incomporate category type into tbl object
  tbl <- merge(tbl,sim_map,by=c("gene","cluster_id"))
  
  
  # no. of DS genes per cluster
  res_by_k <- split(tbl, tbl$cluster_id)
  vapply(res_by_k, nrow, numeric(1))
  # top hits in each cluster
  top <- do.call("rbind", lapply(res_by_k, head, 3))
  top <- select(top, -c("p_adj.glb"))
  top$gene <- gsub("^.*\\.", "", top$gene)
  format(data.frame(top, row.names = NULL), digits = 3)
  
  #Preparing TPR vs FDR object
  tbl$category <- lapply(tbl$category, function(x) {gsub("de", 1, x) })
  tbl$category <- lapply(tbl$category, function(x) {gsub("ee", 0, x) })
  
  ## Empty COBRAData object:
  COBRAData()
  
  method_name = paste0(type,"-",ds_method,"-",vst)
  pval <- as.data.frame(tbl$p_val)
  colnames(pval) = method_name
  padj <- as.data.frame(tbl$p_adj.loc)
  colnames(padj) = method_name
  truth = as.data.frame(as.integer(unlist(tbl$category)))
  colnames(truth) = method_name

  
  return(list(sim = sim,res = res ,tbl = tbl,res_by_k = res_by_k, pval  = pval,padj = padj, truth = truth))
  
}

# TPR vs FDR Plot
# TPR vs FDR
library(iCOBRA)
#library(purrr)
Plot.TPR.vs.FDR <- function(pval,padj,truth,method_name,num){
  
  #Initialize a COBRA object
  cobradata <- COBRAData(pval = pval, padj = padj,truth = truth)
  #methods <- colnames(cobradata@truth)
  #print(method_name)
  cobraperf <- calculate_performance(cobradata,binary_truth = method_name, aspects =c("fdrtpr", "fdrtprcurve"))
  #cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = "Dark2", incltruth = TRUE)
  
  # p <- plot_fdrtprcurve (prepare_data_for_plot(cobraperf))
  
  color <-  c("red","blue", "green","orange","purple","yellow","brown")
  p <- plot_fdrtprcurve(
    prepare_data_for_plot(cobraperf), 
    linewidth = 0.6, pointsize = 1) +
    guides(color = guide_legend(override.aes = list(size = 1.2, linetype = 0))) +
    scale_color_manual(values = color[num]) 
  #  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0.05, 0)) +
  #  scale_x_sqrt(breaks = c(0.01, 0.05, 0.1, seq(0.2, 1, 0.2)), expand = c(0.05, 0))
  
  p$layers[[1]] <- NULL # remove vertical dashed lightgrey lines 
  # p$layers[[1]]$aes_params$size <- 0.4 # dashed lines at thresholds
  
  #print(p)
  
  #ggsave(p, device = "pdf")
  return(p)
}

# Violin plots
Violin_plots<-function (sim,res_by_k,cs_by_k){
  v1 <- plotExpression(sim[, cs_by_k$Neuronal_excit], 
                       features = res_by_k$Neuronal_excit$gene[seq_len(8)],
                       x = "sample_id", colour_by = "group_id") + theme_classic() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  v2 <- plotExpression(sim[, cs_by_k$Neuronal_inhib], 
                       features = res_by_k$Neuronal_inhib$gene[seq_len(8)],
                       x = "sample_id", colour_by = "group_id") + theme_classic() + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #print(v1)
  #print(v2)
  return(list(v1 = v1, v2 = v2))
}

# wrapper to prettify reduced dimension plots
.plot_dr <- function(sim, dr, col)
  plotReducedDim(sim, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

DR.colored.by.expression<-function(sim,tbl,topnumber,cs_by_k){
  
  cs100 <- unlist(sapply(cs_by_k, function(u) 
    sample(u, min(length(u), 100))))
  
  top <- bind_rows(tbl) %>% 
    top_n(topnumber, dplyr::desc(p_adj.loc)) %>% 
    pull("gene")
  # for ea. gene in 'top', plot t-SNE colored by its expression 
  ps <- lapply(top, function(g).plot_dr(sim[, cs100], "TSNE", g) + 
                 ggtitle(g) + theme(legend.position = "none"))
  # arrange plots
  p <- plot_grid(plotlist = ps, ncol = 4, align = "vh")
  #print(p)
  return(p)
}

# Heatmaps
pbHeatmap_plots<- function(sim,res){
  
  # single gene across all clusters
  # top-20 DS genes for single clusters
  print(typeof(sim))
  print(typeof(res))
  print(pbHeatmap(sim, res, k = "Neuronal_inhib"))
  print(pbHeatmap(sim, res, k = "Neuronal_excit"))
  
  # top-5 DS genes per cluster
  print(pbHeatmap(sim, res, top_n = 5))
  
}


DS.analysis.Visualization.mm <-function(ds,ds_method,type,vst,topnumber,num,sim_type){
  
  sim <- ds$sim
  res <-ds$res
  tbl <- ds$tbl
  res_by_k <- ds$res_by_k
  
  
  ## Violin Plots
  cs_by_k <- split(colnames(sim), sim$cluster_id)
  v <- Violin_plots(sim,res_by_k,cs_by_k)
  v1 <-v$v1
  v2 <- v$v2
  
  ### Between-cluster concordance
  ds_gs <- lapply(res_by_k, pull, "gene")
  upset <- upset(fromList(ds_gs), sets = levels(sim$cluster_id))
  #print(upset)
  
  #decideTests(tbl$)
  
  ### DR colored by expression
  if ((ds_method == "dream")){
    #print("DR colored by expression Plot")
    ds_expression <- DR.colored.by.expression(sim,tbl,topnumber,cs_by_k)
  } 
  else{
    ds_expression <- NULL
  }
  
  method_name = paste0(type,"-",ds_method,"-",vst)
  #print(method_name)
  TPR.vs.FDR <- Plot.TPR.vs.FDR(ds$pval,ds$padj,ds$truth,method_name,num)
  
  ### Write results to .rds
  mainDir <-  paste0(getwd(),"/")
  subDir <- "output/"
  ifelse(!dir.exists(file.path(mainDir, subDir,sim_type)), dir.create(file.path(mainDir, subDir,sim_type), recursive = TRUE), FALSE)
  
  print("DS_METHOD")
  saveRDS(res, file.path(paste0(mainDir,subDir,"/",sim_type), paste0("DS_results_mm","_",ds_method,".rds")))
  
  return(list(v1 = v1, v2 = v2, upset = upset, ds_expression = ds_expression, TPR.vs.FDR = TPR.vs.FDR ))
  
}

DS.analysis.Visualization.pb <-function(ds,assay,fun,ds_method,topnumber,num,sim_type){
  
  sim <- ds$sim
  res <-ds$res
  tbl <- ds$tbl
  res_by_k <- ds$res_by_k
  
  
  if (!(assay == "counts" && fun=="sum" && ds_method == "limma-voom")
      && !(assay=="logcounts" && fun == "mean" && ds_method =="edgeR")){
    ### Pseudobulk-level heatmaps
    #pbHeatmap_plots(sim,res)
    
  }
  
  ### Violin Plots
  cs_by_k <- split(colnames(sim), sim$cluster_id)
  v <- Violin_plots(sim,res_by_k,cs_by_k)
  v1 <- v$v1
  v2 <- v$v2
  
  ### Between-cluster concordance
  
  ds_gs <- lapply(res_by_k, pull, "gene")
  
  upset <- upset(fromList(ds_gs), sets = levels(sim$cluster_id))
  #print(upset)
  
  ### DR colored by expression
  if (!(assay=="logcounts" && fun == "mean" && ds_method =="edgeR")){
    ds_expression <- DR.colored.by.expression(sim,tbl,topnumber,cs_by_k)
  }
  else{
    ds_expression <- NULL
  }
  
  method_name = paste0(assay,"-",fun,"-",ds_method) 
  TPR.vs.FDR <- Plot.TPR.vs.FDR(ds$pval,ds$padj,ds$truth,method_name,num)
  
  ### Write results to .rds
  mainDir <-  paste0(getwd(),"/")
  subDir <- "output/"
  ifelse(!dir.exists(file.path(mainDir, subDir,sim_type)), dir.create(file.path(mainDir, subDir,sim_type), recursive = TRUE), FALSE)
  
  saveRDS(res, file.path(paste0(mainDir,subDir,"/",sim_type), paste0("DS_results_",assay,"_",fun,"_",ds_method,".rds")))
  
  return(list(v1 = v1, v2 = v2, upset = upset, ds_expression = ds_expression, TPR.vs.FDR = TPR.vs.FDR ))
  
}

