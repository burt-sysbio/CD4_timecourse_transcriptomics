# run order

#/kinetic genes
masigpro_run_kinetic1 --> run masigpro on each celltype separately
masigpro_run_kinetic2 --> combine output from masigpro_run_kinetic1 into one table
kinetic_2_timepoints --> identify genes with FC>1 to day 0 at two consecutive timepoints
kinetic_combined --> get intersection of kinetic genes from both approaches

# 
masigpro_run_deg_kinetic_baseline --> run masigpro only on kinetic genes with full output

# clustering
fig2/fig2_1_heatmap_clusters --> cluster kinetic genes
cluster_pathways --> save genes assigned to clusters

classswitch_dataprep --> based on cluster assignment of kinetic genes identify class switches
quantDEG_volcano --> add pval / corr idx filters
quantDEG_intersect --> venn diagram approach for DEGs
quantDEG_summary --> combine output from quantDEG_volcano into one table
classswitch_pathways --> check classswitches for quantDEGs