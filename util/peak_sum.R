#Sum total coverage of peaks over different windows
# file: bed file of peak locations
# wfile: bed file of windows
# output: output file for peak coverage
# chr_sort: should peak location file be sorted (set as TRUE unless file is already sorted)
peak_sum = function(file,wfile,output,chr_sort=TRUE){
  library(Matrix)
  print(file)
  #Read in bed file
  atac = read.table(file,header=TRUE,stringsAsFactors = FALSE,sep="\t")
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
  #Create vector of total reads per window
  total = rep(0,n)

  #Sum intersections with window
  #Leverages both bed files being sorted
  x=1
  j=1
  for(i in 1:an){
    at = atac[i,]
    if(at[1] %in% chr_list){
      for(j in x:n){
        wd = window[j,]
        if(at[1]==wd[1] && at[2] < wd[3]){
          total[j] = unlist(total[j]+at[3]-at[2])
          break
        }
      }
    }
    x=j
  }

  #save sums
  write.table(cbind(window,total),file=output,row.names=FALSE,col.names = FALSE,quote=FALSE)
}
