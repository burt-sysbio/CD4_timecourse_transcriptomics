# Data analysis code
Code underlying the manuscript "Dissecting the dynamic transcriptional landscape of early T helper cell differentiation into Th1, Th2, and Th1/2 hybrid cells", Burt et al., Frontiers Immunology (2022).

# Usage
Download/clone repository to local drive. Data to reproduce figures are contained in /data dir. Corresponding code is provided in /code dir. genesets_input contains sets for candidate genes and genesets for enrichment analysis. Utility functions called by scripts are provided in utils.R. To reproduce analysis for chipseq enrichment first run code/linear_model/run_linear_model.R.

# Requirements
All analysis were run in R version 4.1.2
packages:
- dplyr 1.0.8
- tidyr 1.2.0
- stringr
- readr 2.1.2
- clusterprofiler 4.2.2
- RColorBrewer 1.1-2
- ggplot2 3.3.5
- masigpro 1.66.0
- pheatmap 1.0.12

PCA timecourse simulations were run in python v. 3.10.4
packages:
- numpy 1.21.5
- pandas 1.4.1
- seaborn 0.11.2
- matplotlib 3.5.1
