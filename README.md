# RIDDLER

# Running Code
Algorithm is run in R, with a bedfile of cell barcoded de-duplicated reads as the primary input.  See run_riddler.R for an example use case and function inputs.

## Feature files and peak covariates
RIDDLER requires a feature and window location file to run, and optionally a peak covariate file as well.  We have included examples feature and window files for window sizes of 1MB, 500KB, 250KB, and 100KB, located in features/features.[size].csv and features/window_[size].bed.  These feature files contain GC and mappability scores for hg38, aggegated at the specified resolutions, which are useful for a wide variety of assays.  See later sections for steps on how to generate these for different window sizes.

An additional covariate file can also be passed to the algorithm, such as ATAC-seq peak regions for various tissues.  Examples for different window sizes are also included in features/tissue_peaks.[size].bed. The 5th and onward columns of this file will be used as separate covariates, with the option to set a parameter for minimum correlation with data in order to be included (most correlated column is always included).

### Creating new feature files for different window resolutions
To run the algorithm at different resolutions, the feature files will need to re-generated at the desired resolution.  This will require bedtools (https://bedtools.readthedocs.io/en/latest/) for GC content generation, a mappability kmer big-wig file (we use the 50kmer Umap option from https://bismap.hoffmanlab.org/), and the script bigWigAverageOverBed (https://www.encodeproject.org/software/bigwigaverageoverbed/).

We have included a bash script util/create_features.sh to creatue new feature files once the above dependencies have been included.  The script will need to be run with the files in that folder.

Script usage:
./create_features.sh [size of windows] [suffix] [path to assembly] [assembly chromosome sizes] [mappability big wig]

Chromosome sizes for hg38 are included as util/hg38.fa.sizes.  The suffix parameter will name the outputs as features.suffix.csv and window_suffix.bed.

Example usage:
./create_features.sh 10000000 10MB ~/hg38.fa hg38.fa.sizes k50.Umap.MultiTrackMappability.bw

### Creating new peak covariates
The peak covariate input for the algorithm is the relative coverage of each window by the input peak set.  We have included tissue_peaks.1MB.bed (and other resolution versions), as a default covariate for 1MB windows derived from a multi-tissue atlas for ATAC-seq data.  Each tissue in the atlas is included as a different peak covariate, making this a reliable default for many cell types.

To generate the multi-tissue peak file for a different resolution, you can run util/multi_peak_sum.R in R, using the window_[suffix].bed output from the feature generation step above.  You will need to download the DNase site data from https://zenodo.org/records/3838751 (DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz)

Example usage in R:

multi_peak_sum("DHS_Index_and_Vocabulary_hg38_WM20190703.txt","window_10MB.bed","tissue_peaks.10MB.bed")

If you wish to use your own called peaks as covariates, you can use util/peak_sum.R to create a peak feature file from any bed format input, such as the output from MACS2

## Workflow
Once the necessary feature and window files are created, the algorithm is run in three steps, corresponding to three function calls: 1. Bin the reads per-cell at the desired window resolution. 2. Fit the robust poisson GLM model to the data and features. 3. Call CNVs from model residuals.

### Bin reads
The function window_matrix in util/fast_window_matrix.R will bin reads in a bed-like file format, with columns corresponding to chr, start, end, barcode. It matches read locations to a window file, such as those generated in the above steps, creating a window by cell matrix of reads.  This matrix is saved to a specified file.

source("util/fast_window_matrix.R")

window_matrix("data/SNU601_subset_frags_sorted.tsv",w_file="features/window_1MB.bed",out_file="data/SNU601_1MB.bed",frag_min=10000)

### Fit GLM model
After the read matrix file is created, the GLM model can be fit. The function riddler_fit takes the matrix file location, feature file location, peak file location (or NULL if not used), and additional parameters and returns a riddler object containing the normalized version of the matrix (obs_norm), expected values under the model (expected), the fitted model (model), and window locations (window)

source("riddler_fit.R")

rid_obj = riddler_fit(mat_file,feat_file="features/features.1MB.csv",peak_file="features/tissue_peaks.1MB.bed",cov_min=10000,peak_corr_filter=.5)

### Call copy numbers
The function copy_call takes a riddler object, and identifies significant deviations from the expected model values via a sliding window test across each cell. These significant deviations are called as CNVs, and assigned a copy number based on the deviation from normal (RIDDLER uses the copy number of 1 for a normal cell). The returned riddler object contains a matrix of copy number values for each cell and window (CNV) and the p-value calculated for each window's deviation from the expected value (prob).
The list of copy number values to be considered is passed as an argument to the function, expressed as multipliers of the base copy number (1).

source("copy_call.R")

rid_obj = copy_call(rid_obj,w_size=5,prob_thresh=.025,mult=c(0,.5,1,1.5,2,3))

## Plotting CNVs
We include the function cnv_heatmap, a wrapper for ComplexHeatmap, to create plots of the reported CNV matrix from RIDDLER. 

source("cnv_heatmap.R")

heat_plot = cnv_heatmap(rid_obj,heat_title="RIDDLER CNVs",ref_bar=NULL,ref_title=NULL,cluster=TRUE)

ggsave(plot=heat_plot,"example_plot.png",device="png",height=5,width=8)
