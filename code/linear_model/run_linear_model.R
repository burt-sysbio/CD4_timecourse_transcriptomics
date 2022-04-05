# run the whole linear model pipeline
#load output from linear model and perform pathway analysis
#needs to process df and derive clean gene sets, then match to stat and TF targets

source("code/paper_theme.R")
source("code/linear_model/fit_linear_model.R")
source("code/linear_model/proc_linear_model.R")
source("code/linear_model/run_chipseq_pathways.R")

baselines <- c("degs_Th1_Th2_cor0.3")

# uncomment if analysis should be run for different correlation thresholds
# threshold_array <- seq(0.24,0.5,0.01) # provide input array for correlation thresholds
# baselines <- paste0("degs_Th1_Th2_cor", threshold_array)

for(baseline in baselines){
  fit_model(baseline)
  proc_model(baseline)
  run_pathways(baseline)
}
