#Sum total coverage of peaks over different windows
#Works for DNase multi-tissue file, separating each tissue type as a different column in output file.
# file: DNase multi-tissue file, or similar.
#	File must start with peak location as 'chr start stop', have a column denoting types of peaks
# wfile: bed file of windows
# output: output file for peak coverage
# tid: column index of peak type information. Each unique value becomes a column in the output file.
multi_peak_sum = function(file,wfile,output,chr_sort=FALSE,tid=10){
  library(Matrix)
  library(data.table)
  print(file)
  tid=10
  #Read in bed file
  atac = fread(file,data.table=FALSE)
  types = unique(atac[,tid])
  chr_list = c(paste("chr",1:22,sep=""),"chrX","chrY")
  #atac = atac[atac[,1] %in% chr_list,]
  if(chr_sort){
    atac_n = NULL
    for(chr in chr_list){
      atac_n = rbind(atac_n,atac[atac[,1]==chr,])
    }
    atac = atac_n
  }
  an = dim(atac)[1]

  #Read in window file
  window = read.table(wfile,header=FALSE,stringsAsFactors = FALSE)
  n = dim(window)[1]
  #Create matrix of total reads per window
  total = matrix(0,nrow=n,ncol=length(types))
  colnames(total) = types

  #Sum intersections with window
  #Leverages both bed files being sorted
  x=1
  j=1
  for(i in 1:an){
    at = atac[i,]
    if(at[1] %in% chr_list){
      id = which(types == unlist(at[tid]))
      for(j in x:n){
        wd = window[j,]
        if(at[1]==wd[1] && at[2] < wd[3]){
          total[j,id] = unlist(total[j,id]+at[3]-at[2])
          break
        }
      }
    }
    x=j
  }

  #save sums
  write.table(cbind(window,total),file=output,row.names=FALSE,col.names = FALSE,quote=FALSE)
}
