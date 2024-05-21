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

#Run algorithm using: 
#	matricized reads (mat_file) 
#	GC and mappability feature file (feat_file)
#	peak feature file (peak_file)
#  See util scripts and examples to create these inputs for different window sizes.
#  Custom covariates can be used, so long as they are at the same resolution and format.

#Fit robust GLM and compute expected values.  
#	cov_min: Coverage minimum per cell, default set to 10,000 reads
#	peak_corr_filter: If using a multi-columned peak file, like tissue peak set, filter columns who's correlation with bulk signal is greater than this value.
#		Set to -1 to include all peak types.
#  Dropout values in observed and expected matrices are assigned values of -1, for downstream filtering.
#  Parameters to adjust robust glm fitting:
#	'maxit': maximum iterations, default 100
#	'huber': parameter for huber loss, default 1.345 (allows for 95% efficiency in absense of outliers)
#	'acc': threshold for convergence of coeficients, default 1e-04
mat_file = "data/SNU601_1MB.bed"
rid_obj = riddler_fit(mat_file,feat_file="features/features.1MB.csv",peak_file="features/tissue_peaks.1MB.bed",cov_min=10000,peak_corr_filter=.5)

#Call CNVs from fitted model.
#  'w_size' is the size of the sliding window for p-value calculation
#  'prob_thresh' is the 2-sided FDR threshold value.  
#  'mult' is a list of multipliers of base expected values to predict as copy numbers.
rid_obj = copy_call(rid_obj,w_size=5,prob_thresh=.025,mult=c(0,.5,1,1.5,2,3))

#Create heatmap plot of gains and losses
#  'ref_bar': An optional vector of values, at the same resolution, that will be displayed above heatmap.
#		Intended for use as reference distribution of CNVs.
#  'ref_title': A title that will be displayed next to the reference vector, if used.
heat_plot = cnv_heatmap(rid_obj,heat_title="RIDDLER CNVs",ref_bar=NULL,ref_title=NULL,cluster=TRUE)
#plot can be saved as a png, jpg, pdf, etc, using ggasve.
ggsave(plot=heat_plot,"example_plot.png",device="png",height=5,width=8)

