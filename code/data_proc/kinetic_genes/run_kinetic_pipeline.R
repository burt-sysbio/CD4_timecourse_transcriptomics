require(maSigPro)
require(readr)
require(dplyr)
require(tidyr)

# get and save all the kinetic genes as a table
# set wd explicitly in case script is used in server call
#setwd("~/Documents/projects/2020/kinetic_paper/")

source("code/data_proc/utils_data_proc.R")
source("code/data_proc/kinetic_genes/masigpro_run_kinetic1.R")
source("code/data_proc/kinetic_genes/masigpro_run_kinetic2.R")
source("code/data_proc/kinetic_genes/kinetic_2_timepoints.R")
source("code/data_proc/kinetic_genes/kinetic_combined.R")
source("code/data_proc/kinetic_genes/make_correlation_df.R")
source("code/data_proc/kinetic_genes/classswitch_dataprep.R")
source("code/data_proc/kinetic_genes/cluster_pathways.R")

# load data
df <- read_csv("data/peine_data_log2fc_annot.csv")

# criterion for kinetic fold change at two consecutive time points
fold_change <- 1

# uncomment to run masigpro again
#run_masig(df)
proc_masig(df)
run_fc_crit(df, fold_change)
out <- summarize_kinetic(fold_change)

# also compute correlation
#make_correlation_df(df)


print("make sure that cluster ordering is correct based on pyhton script figure 2 cluster quantification")


clust_annot <- read.csv("data/clustered_genes/kinetic_genes_reordered_3clust.csv")

proc_classswitch_data(clust_annot)
proc_kinetic_pathways(clust_annot)
