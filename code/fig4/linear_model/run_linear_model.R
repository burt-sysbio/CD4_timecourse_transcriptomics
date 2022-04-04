# run the whole linear model pipeline
require(dplyr)
require(tidyr)
require(pheatmap)
require(readr)
#require(ontologyIndex)
require(clusterProfiler)

#load output from linear model and perform pathway analysis
#needs to process df and derive clean gene sets, then match to stat and TF targets

require(msigdbr)
require(stringr)
require(tidyr)

#use as background all genes between Th1 and Th2 DEGs
source("code/paper_theme.R")
source("code/fig4/linear_model/fit_linear_model.R")
source("code/fig4/linear_model/proc_linear_model.R")
source("code/fig4/chipseq/run_chipseq_pathways.R")
source("code/fig4/linear_model/heatmap_linear_model_categories.R")
source("code/fig4/chipseq/barplot_chipseq_enrichment.R")
source("code/fig4/linear_model/timecourse_linear_model_categories.R")


#baselines <- c("kinetic_th12", "degs_Th1_Th2_cor0.5",  "degs_Th1_Th2_cor0","degs_Th1_Th2_cor0.3")
#baselines <- c("degs_Th1_Th2_cor0.3")

threshold_array <- seq(0.24,0.5,0.01)
baselines <- paste0("degs_Th1_Th2_cor", threshold_array)

for(baseline in baselines){
  fit_model(baseline)
  proc_model(baseline)
  #plot_heatmap_categories(baseline)
  out <- run_pathways(baseline)
  #plot_pathways(baseline)
  #plot_timecourse(baseline)
}
