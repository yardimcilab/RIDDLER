# RIDDLER

# Running Code
RIDDLER requires feature and window files to run.  We have included examples for window sizes of 1MB, 500KB, 250KB, and 100KB.  See later sections for steps on how to generate these for different window sizes.
Algorithm is run in R, with the fragment file for scATAC-seq data as the primary input.  See run_riddler.R for an example use case and function inputs.  Fragment file will need to be converted into a window by cell matrix for the main functions, which is also covered in the example.

# Creating new feature files for different window resolutions
To run the algorithm at different resolutions, the feature files will need to re-generated at the desired resolution.  This will require bedtools (https://bedtools.readthedocs.io/en/latest/) for GC content generation, a mappability kmer big-wig file (we use the 50kmer Umap option from https://bismap.hoffmanlab.org/), and the encode script bigWigAverageOverBed (https://www.encodeproject.org/software/bigwigaverageoverbed/).

We have included a bash script util/create_features.sh to creatue new feature files once the above dependencies have been included.  The script will need to be run with the files in that folder.

Script usage:
./create_features.sh [size of windows] [suffix] [path to assembly] [assembly chromosome sizes] [mappability big wig]

Chromosome sizes for hg38 are included as util/hg38.fa.sizes.  The suffix parameter will name the outputs as features.suffix.csv and window_suffix.bed.

Example usage:
./create_features.sh 10000000 10MB ~/hg38.fa hg38.fa.sizes k50.Umap.MultiTrackMappability.bw

# Creating new peak covariates
The peak covariate input for the algorithm is the relative coverage of each window by the input peak set.  We have included multi_tissue_peaks.1MB.bed (and other resolution versions), as a default covariate for 1MB windows derived from a multi-tissue atlas.  Each tissue in the atlas is included as a different peak covariate, making this a reliable default for many cell types.

To generate the multi-tissue peak file for a different resolution, you can run util/multi_peak_sum.R in R, using the window_[suffix].bed output from the feature generation step above.  You will need to download the DNase site data from https://zenodo.org/records/3838751 (DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz)

Example usage in R:
multi_peak_sum("DHS_Index_and_Vocabulary_hg38_WM20190703.txt","window_10MB.bed","tissue_peaks.10MB.bed")

If you wish to use your own called peaks as covariates, you can use util/peak_sum.R to create a peak feature file from any bed format input.
