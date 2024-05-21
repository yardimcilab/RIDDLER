#Turn bed reads from all cells into window x cell matrix
window_matrix = function(bed_file,w_file,out_file,frag_min=5000){
  library(data.table)
  print(paste0("Reading ",bed_file))
  print(paste0("Write output to ",out_file))
  chr_list = c(paste("chr",1:22,sep=""),"chrX","chrY")

  #Load fragments, get barcodes with at least frag_min reads
  frags = fread(bed_file,data.table=FALSE)
  print("Fragments loaded")
  b_tab = table(frags[,4])
  barcodes = names(b_tab[b_tab >= frag_min])
  frags = frags[frags[,4] %in% barcodes,]
  print(paste0("Barcodes filtered: ",length(barcodes)))
  print(paste0("Fragments to read: ",nrow(frags)))
  #Read in window file
  window = read.table(w_file,header=FALSE,stringsAsFactors = FALSE,)
  n = nrow(window)
  window = cbind(window,1:n)
  
  hash_window = list()
  for(i in 1:length(chr_list)){
    hash_window[[i]] = window[window[,1]==chr_list[i],]
  }
  names(hash_window) = chr_list

  
  #Create matrix of total reads per window
  counts = matrix(0,nrow=n,ncol=length(barcodes))
  colnames(counts) = barcodes
  

  t = proc.time()[3]
  for(chr in chr_list){
    frags_c = frags[frags[1] == chr,]
    chr_wnd = hash_window[[chr]]
    #quickly find intersecting windows
    s_ids = findInterval(frags_c[,2],chr_wnd[,2]) 
    e_ids = findInterval(frags_c[,3],chr_wnd[,3])+1
    c_ids = which(s_ids != e_ids)
    if(length(c_ids) > 0){
      #find best overlapping window
      offset = apply(abs(frags_c[c_ids,c(2,3)]-chr_wnd[s_ids[c_ids],3]),1,which.max)-1
      s_ids[c_ids] = s_ids[c_ids] + offset
    }
    #Adjust chr ids to global ids
    w_ids = chr_wnd[s_ids,5]
    #Add to counts
    tab_frag = table(w_ids,frags_c[,4])
    counts[as.numeric(rownames(tab_frag)),colnames(tab_frag)] = tab_frag
    
    print(paste(chr,proc.time()[3]-t))
    t = proc.time()[3]
  }


  #save sums
  write.table(cbind(window[,1:4],counts),file=out_file,row.names=FALSE,col.names = TRUE,quote=FALSE)
}
