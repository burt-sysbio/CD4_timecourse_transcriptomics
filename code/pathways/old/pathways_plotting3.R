# take multiple pathways and plot them together in dotplot
require(readr)
require(tidyr)
require(ggplot2)
require(stringr)
require(dplyr)
require(forcats)

get_files <- function(nogo){
  # list files
  mypath <- "data/ORA_output/ORA_burt/"
  filenames <- list.files(mypath, pattern = ".txt", full.names = T)
  filenames2 <- list.files(mypath, pattern = ".txt", full.names = F)
  
  
  # get only the deg files (but not for Th0)
  pattern_degs <- "cor0.5"
  filenames <- filenames[grepl(pattern_degs, filenames)]
  filenames2 <- filenames2[grepl(pattern_degs, filenames2)]
  
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
  # add the filename as a column, transform filename first
  name <- str_split(fname, "_")
  name <- paste0(name[[1]][5], "vs", name[[1]][6])
  df$name <- name
  return(df)
}

# remove nonsignif categories
apply_fdr_filter <- function(df, alpha, filter){
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
  
  # manually change some names
  out$ID[grepl("SOMATIC_RECOMBINATION", out$ID)] <- "GOBP_ADAPTIVE_IMMUNE_RESPONE_BASED_ON_SOMATIC_RECOMBINATION"
  
  # filter for user curated pathways in same directory
  if(user_curated){
    fname <- paste0(savedir, "pathways_curated.csv")
    pathways_curated <- read.csv(fname)
    pathways_curated$myorder <- rownames(pathways_curated)

    out <- inner_join(out, pathways_curated)
    # reorder according to provided manual order in csv file
    print(pathways_curated)
    out <- out %>% mutate(New_ID = fct_reorder(New_ID,as.numeric(myorder), .desc = T), ID = New_ID)

  }
  
  return(out)  
}
# remove the category stuff from plot
#out <- out %>% mutate(ID = str_sub(ID, (nchar(category)+2),-1))
plot_dotplot <- function(out, alpha, filter, width, height, nogo, user_curated, savedir){

  S1 <- ggplot(out, aes(x= name, y=ID, size=Count, color=fdrlog10, group=name)) + 
    geom_point(alpha = 0.8) + 
    theme_bw(base_size = 15) +
    scale_color_gradient(low = "lightgrey",  high = "black", space = "Lab")+
    scale_size(range = c(2, 8)) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 9))
  
  if(user_curated){
    sname0 <- "curated"
  } else {
    sname0 <- "uncurated"
  }
  
  sname <- paste0(savedir, "ORA_fdr", alpha, "filter", filter,"_", sname0, ".png")
  ggsave(sname, bg = "white", width = width, height = height) 
  
  # also save as table for manual curation that can be read in after processing
  if(!user_curated){
    sname2 <- str_sub(sname, 1, -5)
    sname2 <- paste0(sname2, ".csv")
    
    ## reorder 
    out_table <- out$ID[order(out$ID, decreasing = T)]
    out_table <- unique(out_table)
    write.table(out_table, sname2, row.names = F, col.names = c("ID"))    
  }
}

pipeline <- function(alpha, filter, width, height, nogo, user_curated){
  
  filedir <- "data/ORA_output/ORA_burt/"
  if(nogo){
    savedir <- "figures/pathway_results/ALLDB_noGO/"
  } else {
    savedir <- "figures/pathway_results/ALLDB/"
  }
  
  # list files in directory, read them and combine into data frame
  out <- get_files(nogo)
  # process df
  out <- proc_files(out, user_curated, savedir)
  # keep only signif categories
  out <- apply_fdr_filter(out, alpha, filter)
  # plot and save results
  plot_dotplot(out, alpha, filter, width, height, nogo, user_curated, savedir)
}


###################################################################################################
###################################################################################################
###################################################################################################
# analysis starts here


# filter data if user curation criterion is true
filters <- c("INTERFERON", "IMMUN", "INFLAMM", "INTERLEUKIN",
             "DEFENSE", "VIRUS", "TCELL", "T_CELL", "CYTOKINE",
             "EFFECTOR", "SIGNAL", "BACTERIUM", "STIMULUS", "IRF", "VIRAL", "IFN",
             "LYMPHO", "TRANSCRIPTION", "INNATE", "ILC", "LEUCO")

filters <- paste(filters, collapse = "|")
# kick out covid stuff
#out <- out[grepl(filters, out$ID),]
#out <- out[!grepl("COVID|CORONA|SARS", out$ID),]

out1 <- pipeline(alpha = 0.1, filter = 2, width = 15, height = 18, nogo = T, user_curated = F)
out2 <- pipeline(alpha = 0.1, filter = 2, width = 15, height = 30, nogo = F, user_curated = F)
#out3 <- pipeline(alpha = 0.1, filter = 2, width = 12, height = 6, nogo = T, user_curated = T)
out4 <- pipeline(alpha = 0.1, filter = 2, width = 8, height = 6, nogo = F, user_curated = T)


#out2 <- pipeline(out, alpha = 0.1, filter = 2, width = 15, height = 30, nogo = use_GO, curated = user_curated)
#out3 <- pipeline(out, alpha = 0.05, filter = 1, width = 15, height = 30, nogo = use_GO, curated = user_curated)
#out4 <- pipeline(out, alpha = 0.05, filter = 2, width = 15, height = 30, nogo = use_GO, curated = user_curated)
