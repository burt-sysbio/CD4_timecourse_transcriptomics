require(ontologyIndex)
require(clusterProfiler)
require(msigdbr)
require(readr)
require(dplyr)


filenames_full <- list.files(path = "genesets_output/pathway_analysis/", pattern = "*.txt", full.names = T)
filenames_save <- list.files(path = "genesets_output/pathway_analysis/", pattern = "*.txt", full.names = F)
files <- lapply(filenames_full, read.table, col.names = "gene")


background <- read.csv("genesets_output/background_genes_cleaned.txt", header = F, col.names = "gene")
universe <- as.character(background$gene)

# load ontologies
cat1 <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:WIKIPATHWAYS")
cat2 <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "REACTOME")
cat3 <- msigdbr(species = "Mus musculus", category = "H")
cat4 <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
cat5 <- msigdbr(species = "Mus musculus", category = "C3", subcategory = "TFT:TFT_Legacy")

# add slim ontology from PANTHER
go_slim <- get_ontology("genesets_output/PANTHERGOslim.obo", propagate_relationships = c("is_a", "part_of"), extract_tags = "everything")
go_slim <- as.data.frame(go_slim)
go_slim <- go_slim[go_slim$namespace == "biological_process",]

# get all IDs, there are also some alternative (secondary annotations in the SLIM annotations)
# split anternative ids
#go_slim <- go_slim %>% separate(col = alt_id, into = c("alt_id1", "alt_id2", "alt_id3", "alt_id4", "alt_id5"), sep = ";", remove = F)
#go_slim_ids <- go_slim %>% select(starts_with("alt_id")) %>% select(-alt_id)
#go_slim_ids$id <- rownames(go_slim_ids)
#go_slim_vec <- as.character(unlist(go_slim_ids, use.names=FALSE))

# turns out (I checked in console) that alternative ids are not in msigdbr df
cat4_red <- cat4[cat4$gs_exact_source %in% go_slim$id,]


categories <- list(cat1, cat2, cat3, cat4, cat5)
cat_names <- c("WIKIPATHWAYS", "REACTOME", "HALLMARK", "GOBP", "TFT")

all_sets <- bind_rows(categories)
all_sets_nogo <- all_sets[all_sets$gs_subcat != "GO:BP",]
# pathway analysis for each file and each data base

sdir <- "data/ORA_output/ORA_burt/"

genesets1 <- as.data.frame(all_sets) %>% select(gs_name, gene_symbol)
genesets2 <- as.data.frame(all_sets_nogo) %>% select(gs_name, gene_symbol)

for(i in seq_along(files)){
  savename <- filenames_save[[i]]
  degs <- files[[i]]
  
  # run ora for all categories combined
  test <- enricher(degs$gene, TERM2GENE = genesets1, universe = universe, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1)
  res <- test@result
  out <- paste0(sdir, "ORA_combined_", savename)
  write.table(res, out)
  
  test <- enricher(degs$gene, TERM2GENE = genesets2, universe = universe, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1)
  res <- test@result
  out <- paste0(sdir, "ORA_nogo_", savename)
  write.table(res, out)
  # run ORA for each category separately
  for(j in seq_along(categories)){
    genesets <- categories[[j]]

    catname <- cat_names[[j]]
    genesets <- as.data.frame(genesets) %>% select(gs_name, gene_symbol)
    
    test <- enricher(degs$gene, TERM2GENE = genesets, universe = universe, pAdjustMethod = "fdr", maxGSSize = 1000, pvalueCutoff = 1, qvalueCutoff = 1)
    res <- test@result
    
    out <- paste0(sdir, "ORA_", catname, "_", savename)
    write.table(res, out)
  }
}


