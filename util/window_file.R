window_bed = function(window,outfile,size_file){
  options(scipen = 100)
  #create list of chromosomes
  chr_list = c(paste("chr",1:22,sep=""),"chrX","chrY")
  chr_sz = chrom_lengths(chr_list,size_file)
  #Create window psuedo-bed table
  window_peaks = data.frame(chr=NULL,start=NULL,end=NULL)
  id = 1
  for(i in 1:length(chr_list)){
    end_loc = chr_sz[i]
    chr = chr_list[i]
    num = ceiling(end_loc/window)
    starts = ((1:num)-1)*window
    ends = (1:num)*window-1
    ends[num]=end_loc
    w_name = paste("window",id:(id+num-1),sep="_")
    id = id+num
    window_peaks = rbind(window_peaks,data.frame(chr=rep(chr,num),start=starts,end=ends,name=w_name,stringsAsFactors = FALSE))
  }
  write.table(window_peaks,outfile,sep="\t",row.names = FALSE,quote=FALSE,col.names=FALSE)
}

chrom_lengths = function(chr_list,file){
  genome_vals = read.csv(file,sep="\t",header=FALSE)
  ids = match(chr_list,genome_vals[,1])
  val = genome_vals[ids,2]
  return(val)
}

