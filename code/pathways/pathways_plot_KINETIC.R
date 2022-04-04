# take multiple pathways and plot them together in dotplot
require(readr)
require(tidyr)
require(ggplot2)
require(stringr)
require(dplyr)
require(forcats)
require(RColorBrewer)
require(pheatmap)

source("code/utils.R")

get_files <- function(nogo){
  # list files
  mypath <- "data/ORA_output/ORA_burt/"
  filenames <- list.files(mypath, pattern = ".txt", full.names = T)
  filenames2 <- list.files(mypath, pattern = ".txt", full.names = F)
  
  # get only the deg files (but not for Th0)
  pattern <- "clust_new"
  filenames <- filenames[grepl(pattern, filenames)]
  filenames2 <- filenames2[grepl(pattern, filenames2)]
  
  # bool: use GO annoation or not?
  if(nogo){
    pattern_comb <- "nogo" # all databases combined
  } else {
    pattern_comb <- "combined" # all DB except GO
  }
  
  filenames <- filenames[grepl(pattern_comb, filenames)]
  filenames2 <- filenames2[grepl(pattern_comb, filenames2)]
  
  # read files and combine
  files <- lapply(filenames, read.table, header = T)
  out <- mapply(myfun1, files, filenames2, SIMPLIFY = F)
  out <- bind_rows(out)
  return(out)
}

myfun1 <- function(df, fname){
  ####
  #### add the cluster number to the file. we agreed to switch cluster numbers based on which is upregulated first
  ####
  # add the filename as a column, transform filename first
  #name <- str_split(fname, "_")
  #name <- paste0(name[[1]][5], "vs", name[[1]][6])
  print("using reordered cluster IDs!")
  if(grepl("clust_new1", fname)){
    fname <- "C2"
  } else if (grepl("clust_new2", fname)){
    fname <- "C3"
  } else if (grepl("clust_new3", fname)){
    # cluster reordered
    fname <- "C1"
  }
  df$name <- fname
  return(df)
}

# remove nonsignif categories
apply_fdr_filter <- function(df, alpha, filter = 1){
  df <- df %>% filter(qvalue <= alpha)
  counts <- df %>% group_by(ID) %>% count()
  myselec <- counts$ID[counts$n>=filter]
  df <- df %>% filter(ID %in% myselec)  
  return(df)
}

proc_files <- function(out, user_curated, savedir){
  # add database as column and add -log10 column
  out <- out %>% separate(ID, into = c("DB"), remove = F)
  out <- out %>% mutate(fdrlog10 = -log10(qvalue))
  

  # rename transcription factor targets because the annotation doesnt work in combined feature
  mycategories <- c("GOBP", "WP", "REACTOME", "HALLMARK")
  out$DB[!(out$DB %in% mycategories)] <- "TFT"
  
  # use manual annotation 
  fname <- paste0(savedir, "pathways_curated_new.csv")
  out <- proc_pathways_curated(out, user_curated, fname)

  return(out)  
}

proc_pathways_curated <- function(out, user_curated, fname){
  
  # check if manual curation exists 
  # if it exists and user_curated True, reduce out to manual curation subset and change pathway names and ordering according to file
  # otherwise if it exists and user_curated False, add IDs to full table
  # otherwise return original
  if(file.exists(fname)){
    sep = ","
    print("current separator is")
    print(sep)
    pathways_curated <- read.csv(fname, strip.white = TRUE)
    
    print(pathways_curated)
    stopifnot("New_ID" %in% colnames(pathways_curated))

    # show only pathways in manual curation file or alternatively add index of manual curation to full table    
    if(user_curated){
      # reorder according to provided manual order in csv file
      pathways_curated$myorder <- rownames(pathways_curated)
      print(pathways_curated$myorder)
      out_merge <- inner_join(out, pathways_curated)
      
      print("not using manual IDs currently, uncomment me otherwise")
      #out_merge$ID[!is.na(out_merge$New_ID)] <- out_merge$New_ID
      #out <- out %>% mutate(New_ID = fct_reorder(New_ID, as.numeric(myorder), .desc = T), ID = New_ID)
    } 
    
  } else {
    print("no manual curation found, returning original")
  }  
  return(out_merge)
}
# remove the category stuff from plot
#out <- out %>% mutate(ID = str_sub(ID, (nchar(category)+2),-1))
plot_dotplot <- function(out, alpha, filter, width, height, nogo, user_curated, sname){

  S1 <- ggplot(out, aes(x= name, y=ID, size=Count, color=fdrlog10, group=name)) + 
    geom_point(alpha = 0.8) + 
    theme_bw(base_size = 15) +
    scale_color_gradient(low = "lightgrey",  high = "black", space = "Lab", limits = c(1, 5))+
    scale_size(range = c(2, 8)) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 9))
  
  ggsave(sname, bg = "white", width = width, height = height) 
  
  # also save as table for manual curation that can be read in after processing
  if(!user_curated){
    sname2 <- str_sub(sname, 1, -5)
    sname3 <- paste0(sname2, ".csv")
    sname4 <- paste0(sname2, "_pathway_summary.csv")
    
    ## reorder 
    out_table <- out$ID[order(out$ID, decreasing = T)]
    out_table <- unique(out_table)
    write.table(out_table, sname3, row.names = F, col.names = c("ID"))  
    write.csv(out, sname4)
  }
}


plot_heatmap <- function(out, sname){
  test <- out %>% select(name, fdrlog10, ID) %>% pivot_wider(names_from = name, values_from = fdrlog10)
  
  # all NAs get 0 instead (log10(1) = 0)
  test[is.na(test)]= 0
  
  test <- as.data.frame(test)
  
  test <- test %>% mutate(ID = str_replace(ID, "REACTOME", "R"))
  test <- test %>% mutate(ID = str_replace(ID, "HALLMARK", "H"))
  mypattern <- "IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL"
  test <- test %>% mutate(ID = str_replace(ID, mypattern, "INTERACTION_LYMPHOID_AND_NON_LYMPHOID_CELL"))
  # make coarse grained heatmap only showing signif *,** and *** differences
  rownames(test) <- test$ID
  test <- test[,2:4]
  test <- test[c("C1", "C2", "C3")]
  test2 <- test
  
  thres1 <- -log10(0.05)
  thres2 <- -log10(0.01)
  thres3 <- -log10(0.001)
  
  test2[test<thres1] = 0
  test2[test>thres1] = 1
  test2[test>thres2] = 2
  test2[test>thres3] = 3
  
  cellwidth = 10
  cellheight = 10
  width = 8
  height = 6
  # get three distinct types of red
  mycolors = colorRampPalette(brewer.pal(n = 4, name ="Reds"))(4)  
  mycolors[1] <- "white"

  sname <- str_sub(sname, 1, -5)
  sname2 <- paste0(sname, ".pdf")
  pheatmap(test2,
           color = mycolors,
           cluster_rows = F,
           cluster_cols = F,
           fontsize_row = 10,
           fontsize_col = 11,
           cellwidth = cellwidth,
           cellheight = cellheight,
           legend = F,
           width = width,
           height = height,
           filename = sname2)
  
  p <- pheatmap(test2,
           color = mycolors,
           cluster_rows = F,
           cluster_cols = F,
           fontsize_row = 10,
           fontsize_col = 11,
           cellwidth = cellwidth,
           cellheight = cellheight,
           legend = F,
           width = width,
           height = height,
           filename = sname2)
  
  
  sname3 <- paste0(sname, ".svg")
  
  save_pheatmap(p, sname3, width = width, height = height)
}


pipeline <- function(alpha, filter, width, height, nogo, user_curated){
  
  filedir <- "data/ORA_output/ORA_burt/"
  if(nogo){
    savedir <- "figures/pathway_results/kinetic_clusters/ALLDB_noGO/"
  } else {
    savedir <- "figures/pathway_results/kinetic_clusters/ALLDB/"
  }
  
  # list files in directory, read them and combine into data frame
  out <- get_files(nogo)
  # process df
  out <- proc_files(out, user_curated, savedir)
  # keep only signif categories
  out <- apply_fdr_filter(out, alpha, filter)
  # plot and save results
  
  # manually remove Wikipathways and Transcription factors
  print("removing Wikipath and TFT")
  out <- out %>% filter(!(DB %in% c("TFT", "WP")))
  
 
  if(user_curated){
    sname0 <- "curated"
  } else {
    sname0 <- "uncurated"
  }
  

  sname <- paste0(savedir, "ORA_fdr", alpha, "filter", filter,"_", sname0, ".png")
  sname_heatmap <- paste0(savedir, "ORA_fdr_heatmap_", alpha, "filter", filter,"_", sname0, ".png") 
  # for full table its too large so split into GO and nonGO
  if((!user_curated) & (!nogo)){
    out1 <- out[grepl("GOBP", out$ID),]
    out2 <- out[!grepl("GOBP", out$ID),]
    sname1 <- paste0(savedir, "ORA_fdr", alpha, "filter", filter,"_GO.png")
    sname2 <- paste0(savedir, "ORA_fdr", alpha, "filter", filter,"_other.png")
    plot_dotplot(out1, alpha, filter, 10, 24, nogo, user_curated, sname1)
    plot_dotplot(out2, alpha, filter, 12, 12, nogo, user_curated, sname2)
  }
  
  #plot_dotplot(out, alpha, filter, width, height, nogo, user_curated, sname)
  plot_heatmap(out, sname_heatmap)
  
  return(out)
}

###################################################################################################
###################################################################################################
###################################################################################################
# analysis starts here

#out1 <- pipeline(alpha = 0.05, filter = 1, width = 10, height = 18, nogo = T, user_curated = F)
out1 <- pipeline(alpha = 0.05, filter = 1, width = 10, height = 18, nogo = T, user_curated = T)
