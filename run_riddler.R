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

#### Creation of input files for desired resolution ####
## Input files for 1MB resolution are provided.  
## If you want to run RIDDLER at different resolutions, see util/create_features.sh to created required files.

#### Creation of Peak and Data matrices for window resolution ####
#Create window by cell matrix from sorted fragment file. Example shown for SNU601 data example
#	Uses window location file, minimum fragments for cell filtering 
source("util/fast_window_matrix.R")
window_matrix("SNU601_subset_frags_sorted.tsv",w_file="window_1MB.bed",out_file="SNU601_1MB.bed",frag_min=10000)

#Create peak feature file from source with multiple peak typs, such as DNase multi-tissue file.
#       Input file should have a bed like format, with 'chr start stop' as the first three columns.
#       'tid' specifies the column number that labels the peak types in the file.
#               Each peak type label will create a different covariate in the output.
source("util/multi_peak_sum.R")
multi_peak_sum(file="util/Multisample_peaks.tsv",wfile="window_1MB.bed",output="tissue_peaks.1MB.bed",chr_sort=FALSE,tid=10)

#Create peak feature file, from a standard peak caller output like MACS2
#	Input file should have a bed like format, with 'chr start stop' as the first three columns.
#	Set 'chr_sort' to TRUE is file is not already sorted.
source("util/peak_sum.R")
peak_sum(file,wfile="window_1MB.bed",output,chr_sort=FALSE)
####

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
mat_file = "SNU601_1MB.bed"
rid_obj = riddler_fit(mat_file,feat_file="features.1MB.csv",peak_file="tissue_peaks.1MB.bed",cov_min=10000,peak_corr_filter=.5)

#Call CNVs from fitted model.
#  'w_size' is the size of the sliding window for p-value calculation
#  'prob_thresh' is the 2-sided FDR threshold value.  
#  'mult' is a list of multipliers of base expected values to predict as copy numbers.
rid_obj = copy_call(rid_obj,w_size=5,knn=2,prob_thresh=.025,mult=c(0,.5,1,1.5,2,3))

#Perform clustering to define subgroups. Choose your favorite clustering algorithm.
library(Rphenograph)
library(Rcpp)
pheno = Rphenograph(as.matrix(t(rid_obj$CNV)),k=20)[[2]]
clusters = membership(pheno)

#Create and plot phylogenetic tree from cluster CNV profiles
source("phylo_tree.R")
tree_obj = phylo_tree(rid_obj$CNV,clusters,norm=1,uniq_thresh=0.05,cons_thresh=.35)
node_plot = tree_plot(tree_obj)
#Including a list of cell labels will create a pie-chart tree instead
fake_labs = sample(c("A","B","C","D"),ncol(rid_obj$CNV),replace=TRUE)
pie_plot = tree_plot(tree_obj,labels=fake_labs,lab_name="Fake Labels")
#plot can be saved as a png, jpg, pdf, etc, using ggasve.
ggsave(plot=node_plot,"example_phylotree.png",device="png",height=5,width=8)
ggsave(plot=pie_plot,"example_pie_phylotree.png",device="png",height=5,width=8)

#Create heatmap plot of gains and losses
#  'ref_bar': An optional vector of values, for each cell, that will be displayed to the left of heatmap.
#		Can be used for cluster groups.
#  'ref_title': A title that will be displayed next to the reference vector, if used.
#		Defaults to 'cluster'
heat_plot = cnv_heatmap(rid_obj,heat_title="RIDDLER CNVs",ref_bar=clusters,ref_title="Cluster")
#plot can be saved as a png, jpg, pdf, etc, using ggasve.
ggsave(plot=heat_plot,"example_heatmap.png",device="png",height=5,width=8)

