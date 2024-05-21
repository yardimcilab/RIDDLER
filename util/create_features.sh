#/bin/bash
size=$1
name=$2
genome=$3
sizes=$4
map=$5
#Create window file
Rscript window_file.R $size $name $sizes
#Get gc content
./gc_content.sh $name $genome window_$name.bed
#Get mappability
./bigWigAverageOverBed $map window_$name.bed mappability_50mer_$name.tab
#create feature file
Rscript combine_gc_map.R mappability_50mer_$name.tab GCcontent.$name.nuc.bed features.$name.csv
#remove intermediate file
rm GCcontent.$name.nuc.bed
rm mappability_50mer_$name.tab
