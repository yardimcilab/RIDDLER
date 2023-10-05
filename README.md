# RIDDLER

# Running Code
RIDDLER requires feature and window files to run.  We have included examples for 1MB windows, see later sections for steps on how to generate these for different window sizes.
Algorithm is run in R, with the fragment file for scATAC-seq data as the primary input.  See run_riddler.R for an example use case and function inputs.  Fragment file will need to be converted into a window by cell matrix for the main functions, which is also covered in the example.

# Creating new feature files for different window resolutions
To run the algorithm at different resolutions, the feature files will need to re-generated.  We have included a bash script create_features.sh to do this in the util/ folder.  The script will need to be run with the files in that folder.

Script usage:
./create_features.sh [size of windows] [suffix] [path to assembly] [assembly chromosome sizes]

Chromosome sizes for hg38 are included as util/hg38.fa.sizes.  The suffix parameter will name the outputs as features.suffix.csv and window_suffix.bed.

Example usage:
./create_features.sh 10000000 10MB ~/hg38.fa hg38.fa.sizes

# Creating new peak covariates
The peak covariate input the algorithm are the relative coverage of each window by the input peak set.  We have included multi_tissue_peaks.1MB.bed, as a default covariate for 1MB windows derived from a multi-tissue atlas.  Peak covariates for custom peak sets and window sizes can be created using util/peak_sum.R


