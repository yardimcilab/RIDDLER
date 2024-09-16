cnv_heatmap = function(rid_obj,heat_title,states=TRUE,ref_bar=NULL,ref_title=NULL,cluster=TRUE){
  library(reshape)
  library(stringr)
  library(ggplot2)
  library(ComplexHeatmap)

  #load data, meta-data
  rv = rid_obj$CNV
  wnd = rid_obj$window

  #wnd = read.csv(infile,header=TRUE,sep=" ")
  #wnd = wnd[,1:3]
  Yids = wnd[,1] == "chrY" | wnd[,1] == "chrX" | wnd[,1] == "chrM"
  wnd = wnd[!Yids,]

  #dendrogram clustering
  if(cluster){
    hc = hclust(dist(t(rv)))
    rv = rv[,hc$order]
  }

  if(!states){
    #Switch to 3 state values
    cnv_names=c("Loss","Normal","Gain")
    temp = rv
    temp[rv < 1] = cnv_names[1]
    temp[rv == 1] = cnv_names[2]
    temp[rv > 1] = cnv_names[3]
    rv = temp
    if(!is.null(ref_bar)){
      temp = ref_bar
      temp[ref_bar < 1] = cnv_names[1]
      temp[ref_bar == 1] = cnv_names[2]
      temp[ref_bar > 1] = cnv_names[3]
      ref_bar = temp
    }
    heat_color = c("dodgerblue","gray90","firebrick1")
    names(heat_color) = c(cnv_names)
  } else{
    cnv_names = c(0,0.5,1,1.5,2,3)
    heat_color = c("blue","dodgerblue","gray90","palevioletred1","firebrick1","red4")
    names(heat_color) = c(cnv_names)
  }
  
  if(!is.null(ref_bar)){
    anno_df = data.frame(x = ref_bar)
    colnames(anno_df) = ref_title
    col_list = list(heat_color)
    names(col_list) = ref_title
    ha = HeatmapAnnotation(df = anno_df,col=col_list,annotation_name_side = "left",
        show_legend=FALSE,border=TRUE)
    ht = Heatmap(t(rv),cluster_rows=TRUE,cluster_columns=FALSE,name="CNV",col=heat_color,show_row_dend = FALSE,
        show_column_names=FALSE,show_row_names=FALSE,column_split=factor(wnd[,1],levels=unique(wnd[,1])),
        column_title_rot = 45,column_title_side = "bottom",row_title = "scATAC",
        top_annotation = ha,border=TRUE,column_gap=unit(2,"mm"))
  } else{
    ht = Heatmap(t(rv),cluster_rows=TRUE,cluster_columns=FALSE,name="CNV",col=heat_color,show_row_dend = FALSE,
	show_column_names=FALSE,show_row_names=FALSE,column_split=factor(wnd[,1],levels=unique(wnd[,1])),
	column_title_rot = 45,column_title_side = "bottom",row_title = "scATAC",
	border=TRUE,column_gap=unit(2,"mm"))
  }
  p = grid.grabExpr(draw(ht,column_title=heat_title))


  return(p)
}
