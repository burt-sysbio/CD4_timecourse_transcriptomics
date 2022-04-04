
proc_kinetic_pathways <- function(df){
  df_wide <- df %>% pivot_wider(names_from = celltype, values_from = annot)
  
  # genes that are assigned to cluster x in all four celltypes (rowsums=4 criterion)
  df_clust1 <- df_wide[rowSums(df_wide[,2:5]==1)==4,]
  df_clust2 <- df_wide[rowSums(df_wide[,2:5]==2)==4,]
  df_clust3 <- df_wide[rowSums(df_wide[,2:5]==3)==4,]
  gene_vec <- list(df_clust1, df_clust2, df_clust3)
  
  sdir <- paste0("genesets_output/pathway_analysis/")
  
  # df for individual clusters
  
  for(i in seq(1:3)){
    df3 <- gene_vec[[i]]
    sname <- paste0(sdir, "kinetic_genes_clust", i, ".txt")
    
    write.table(df3, sname, row.names = F, col.names = F, quote = F)
  }
}