#Source code files
source("riddler_fit.R")
source("copy_call.R")
source("cnv_heatmap.R")
#Required R packages in source files
#library(robustbase)
#library(MASS)
#library(stringr)
#library(fitdistrplus)
#library(reshape)
#library(ggplot2)
#library(ComplexHeatmap)

#First time run: create window by cell matrix from fragment file
source("util/fast_window_matrix.R")
window_matrix(frag_file,w_file="window_1MB.bed",out_file="test_data.1MB.bed")

#Run algorithm using: 
#	matricized reads (mat_file) 
#	GC and mappability feature file (feat_file)
#	peak feature file (peak_file)
#  See util scripts and examples to create these inputs for different window sizes.
#  Custom covariates can be used, so long as they are at the same resolution.

#Fit robust GLM and compute expected values.  Coverage minimum per cell set to 10,000 reads
#  Dropout values in observed and expected matrices are assigned values of -1, for downstream filtering.
mat_file = "test_data.1MB.bed"
rid_obj = riddler_fit(mat_file,feat_file="features.1MB.csv",peak_file="multi_tissue_peaks.1MB.bed",cov_min=10000)

#Call CNVs from fitted model.
#  'w_size' is the size of the sliding window for p-value calculation
#  'prob_thresh' is the 2-sided FDR threshold value.  
#  'mult' is a list of multipliers of base expected values to predict as copy numbers.
rid_obj = copy_call(rid_obj,w_size=5,prob_thresh=.025,mult=c(0,.5,1,1.5,2,3))

#Create heatmap plot of gains and losses
#  'ref_bar': An optional vector of values, at the same resolution, that will be displayed above heatmap.
#		Intended for use as reference distribution of CNVs.
#  'ref_title': A title that will be displayed next the the reference vector, if used.
heat_plot = cnv_heatmap(rid_obj,heat_title="RIDDLER CNVs",ref_bar=NULL,ref_title=NULL)
#plot can be saved as a png, jpg, pdf, etc, using ggasve.
ggsave(plot=heat_plot,"example_plot.png",device="png",height=5,width=8)

