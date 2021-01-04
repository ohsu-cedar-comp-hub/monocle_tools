#****************************************************************
# Single cell functions
#****************************************************************

so_to_cds <- function(so, use_pca = "pca", use_umap = "umap", use_cluster = "orig.ident") {
  # Project PC dimensions to whole data set
  so <- ProjectDim(so, reduction = use_pca)
  
  # Create an expression matrix
  expression_matrix <- so@assays$RNA@counts
  
  # Get cell metadata
  cell_metadata <- so@meta.data
  if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
    print(sprintf("Cell identifiers match"))
  } else {
    print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                  ncol(expression_matrix), nrow(cell_metadata)))
    print("If the counts are equal, sort differences will throw this error")
  }
  
  # get gene annotations
  gene_annotation <- data.frame(gene_short_name = rownames(so@assays$RNA), row.names = rownames(so@assays$RNA))
  if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
    print(sprintf("Gene identifiers all match"))
  } else {
    print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                  nrow(expression_matrix), nrow(gene_annotation)))
    print("If the counts are equal, sort differences will throw this error")
  }
  
  # Seurat-derived CDS
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  # Transfer Seurat embeddings
  # Note that these may be calculated on the Integrated object, not the counts
  #   and thus will involve fewer genes
#  reducedDim(cds, type = "PCA") <- so@reductions$pca@cell.embeddings 
  reducedDim(cds, type = "PCA") <- so@reductions[[use_pca]]@cell.embeddings 
#  cds@preprocess_aux$prop_var_expl <- so@reductions$pca@stdev
  cds@preprocess_aux$prop_var_expl <- so@reductions[[use_pca]]@stdev
#  plot_pc_variance_explained(cds)
  
  # Transfer Seurat UMAP embeddings
  cds@int_colData@listData$reducedDims$UMAP <- so@reductions[[use_umap]]@cell.embeddings
  #    plot_cells(cds)
  
  # Copy cluster info from Seurat
  cds@clusters$UMAP_so$clusters <- so@meta.data[[use_cluster]]
  
  cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)
#  DimPlot(so, reduction = "umap")
#  plot_cells(cds, color_cells_by = "partition", group_label_size = 3.5)
#  plot_cells(cds, color_cells_by = cellids, show_trajectory_graph = FALSE, group_label_size = 3.5)
  
  # Fix from https://gitmemory.com/cole-trapnell-lab
  rownames(cds@principal_graph_aux$UMAP$dp_mst) <- NULL
  colnames(cds@int_colData@listData$reducedDims$UMAP) <- NULL
  
  return(cds)
}