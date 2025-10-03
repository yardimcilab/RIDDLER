## Create plots from tree object
# tree_obj: output of phylo_tree()
# labels: For pie chart version, labels of each cell. Cells must align with clusters in tree_obj
# lab_name: Title of labels. i.e. 'Cell Types'
# slice_colors: Color vector for pie chart labels
# node_colors: Color vector for cluster names
# pie_plot: return pie-chart version of tree
tree_plot = function(tree_obj,labels=NULL,lab_name=NULL,slice_colors=NULL,node_colors=NULL){
  library(igraph)
  library(ggnetwork)
  library(PieGlyph)
  library(ggplot2)

  Adj = tree_obj$Adj
  clusters = tree_obj$clusters
  if(is.null(node_colors)){
    node_colors = tree_obj$node_colors
  }

  #Make tree network
  rid = which(row.names(Adj) == "Normal")
  g = graph_from_adjacency_matrix(Adj)
  layout = layout_as_tree(g,root=rid)
  row.names(layout) = row.names(Adj)

  gf = ggnetwork(Adj,layout=layout)
  node_ids = (nrow(gf)-nrow(Adj)+1):nrow(gf)
  nf = gf[node_ids,]
  ef = gf[-node_ids,]

  if(!is.null(labels)){
    labs = unique(labels)
    D = matrix(0,nrow=nrow(nf),ncol=length(labs))
    colnames(D) = labs
    nf = cbind(nf,D)
    for(i in 1:nrow(nf)){
      v = nf$vertex.names[i]
      for(ll in labs){
        nf[i,ll] = sum(labels[clusters==i] == ll)
      }
    }
    #node_names = names(node_colors)[nf$vertex.names]
    node_names = nf$vertex.names
    lf = cbind(nf,node_names=node_names)
    rng = max(lf$y)-min(lf$y)
    lf$y = lf$y+rng*.1

    if(is.null(slice_colors)){
      slice_colors = colorRampPalette(c("royalblue","cyan","paleturquoise","plum2","purple","deeppink","red3","lightcoral","orange","gold"))(length(labs))
      names(slice_colors) = labs
    }

    p=ggplot()+geom_segment(data=ef,aes(x=x,y=y,xend=xend,yend=yend))+theme_blank()+
	geom_pie_glyph(slices=labs,data=nf,aes(x=x,y=y),radius=.4)+scale_fill_manual(values=slice_colors)+labs(fill=lab_name)+
	geom_label(data=lf,aes(label=node_names,x=x,y=y))+xlim(-.05,1.05)
  } else{  
    node_names = nf$vertex.names
    lf = cbind(nf,node_names=node_names)
    rng = max(lf$y)-min(lf$y)
    lf$y = lf$y+rng*.1

    p=ggplot()+geom_segment(data=ef,aes(x=x,y=y,xend=xend,yend=yend))+theme_blank()+
        geom_point(data=nf,aes(x=x,y=y,color=node_names),size=4,show.legend=FALSE)+
	scale_color_manual(values=node_colors)+
        geom_label(data=lf,aes(label=node_names,x=x,y=y))+xlim(-.05,1.05)
  }
  return(p)
}

## Create phylegenetic tree of CNVs based on clusters
# norm: What copy number value is considered normal
# uniq_thresh: percentage of non-inherited CNV differences allowed between parent and child nodes
# cons_thresh: percent of cells CNV must be expressed in for cluster consensus profile
phylo_tree = function(cnvs,clusters,norm=1,uniq_thresh=0.05,cons_thresh=.35){
  library(cluster)
  library(ggplot2)
  library(igraph)
  library(ggplotify)
  if(length(unique(clusters))<2){
    cons = consensus_cnv(cnvs,cons_thresh,norm)
    nnum = sum(cons != norm)
    if(nnum > length(cons)*uniq_thresh){
      return(c("red"))
    } else{
      return(c("blue"))
    }
  }
  ### Identify consensus CNVs
  cons = matrix(0,nrow=nrow(cnvs),ncol=max(clusters))
  for(i in 1:ncol(cons)){
    cons[,i] = consensus_cnv(cnvs[,clusters==i],cons_thresh,norm)
  } 

  ### Construct tree
  nodes = paste0("Clone ",1:ncol(cons))
  Adj = matrix(0,nrow=ncol(cons),ncol=ncol(cons))
  colnames(Adj) = nodes
  row.names(Adj) = nodes
  used = rep(FALSE,length(nodes))
  #If a cluster is normal, remove it from merging
  nnum = min(colSums(cons != norm))
  if(nnum <= nrow(cons)*uniq_thresh){ 
    norm_id = which.min(colSums(cons != norm))
    used[norm_id] = TRUE 
  }
  #Connect valid parent nodes first
  D = parent_dist(cons[,!used],norm)
  colnames(D) = nodes[!used]
  row.names(D) = colnames(D)
  for(i in 1:ncol(D)){
    j = which.min(D[,i])
    if(D[j,i] < uniq_thresh*nrow(cons) && D[j,i] <= D[i,j]){
      id1 = which(nodes==row.names(D)[i])
      id2 = which(nodes==row.names(D)[j])
      used[id1] = TRUE
      Adj[id1,id2] = 1
      Adj[id2,id1] = 1
    }
  }
  #Connect closest matching remaining nodes with intermediate states
  num=1
  while(sum(!used)>1){
    #find closest matching nodes
    D = parent_dist(cons[,!used],norm)
    colnames(D) = nodes[!used]
    row.names(D) = colnames(D)
    ids = which(D == min(D), arr.ind = TRUE)[1,]
    id1 = which(nodes==row.names(D)[ids[1]])
    id2 = which(nodes==row.names(D)[ids[2]])
    parent = D[ids[1],ids[2]] < uniq_thresh*nrow(cons)
    #replace nodes with ancestor
    used[id1] = TRUE
    used[id2] = TRUE
    #Add intermediate node if necessary
    if(!parent){
      nodes = c(nodes,paste0("Inter_",num))
      num = num+1
      Adj = rbind(Adj,0)
      Adj = cbind(Adj,0)
      con2 = cons[,c(id1,id2)]
      cids = apply(abs(con2-norm),1,which.min)
      cnv = rep(norm,nrow(cons))
      for(j in 1:length(cnv)){cnv[j] = con2[j,cids[j]]}
      cons = cbind(cons,cnv)
      used = c(used,FALSE)
      n = ncol(cons)
      Adj[id1,n] = 1
      Adj[n,id1] = 1
      Adj[id2,n] = 1
      Adj[n,id2] = 1
    } else{
      used[id1] = FALSE
      Adj[id1,id2] = 1
      Adj[id2,id1] = 1
    }
  }

  ### Add normal clone if there isn't one already
  if(nnum > nrow(cons)*uniq_thresh){
    #Check if a normal clone has been added
    nnum = sum(cons[,!used] != norm)
    if(nnum <= nrow(cons)*uniq_thresh){
      norm_id = which(!used)
      nodes[norm_id] = "Normal"
    } else{
      id = which.min(colSums(cons != norm))
      cons = cbind(cons,norm)
      nodes = c(nodes,"Normal")
      Adj = rbind(Adj,0)
      Adj = cbind(Adj,0)
      n= ncol(cons)
      Adj[id,n] = 1
      Adj[n,id] = 1
      norm_id = n
    }
  } else{
    nodes[norm_id] = "Normal"
    id = order(colSums(abs(cons - cons[,norm_id])))[2]
    Adj[norm_id,id] = 1
    Adj[id,norm_id] = 1
  }

  row.names(Adj) = nodes
  colnames(Adj) = nodes

  ### Prune unncessary intermediate nodes
  x = which(rowSums(Adj)==2 & grepl("Inter",nodes))
  while(length(x) > 0){
    for(i in x){
      ids = which(Adj[i,]==1)
      Adj[ids[1],ids[2]] = 1
      Adj[ids[2],ids[1]] = 1
    }
    Adj = Adj[-x,-x]
    nodes = nodes[-x]
    cons = cons[,-x]
    x = which(rowSums(Adj)==2 & grepl("Inter",nodes))
  }

  ### Color scale for plotting
  colnames(cons) = nodes
  node_colors = rep("",length(nodes))
  names(node_colors) = nodes
  cids = grepl("Clone",nodes)
  tids = grepl("Inter",nodes)
  nid = nodes == "Normal"
  node_colors[nid] = "forestgreen"
  node_colors[cids] = colorRampPalette(c("royalblue","cyan","paleturquoise","plum2","purple","deeppink","red3","lightcoral","orange","gold"))(sum(cids))
  node_colors[tids] = colorRampPalette(c("gray0","gray90"))(sum(tids))


  tree_obj = list(Adj=Adj,node_colors=node_colors,node_names=nodes,cons=cons,clusters=clusters) 

  return(tree_obj)
}

tree_order = function(A,id,prev=0){
  nids = which(A[id,]==1)
  nids = nids[nids != prev]
  if(length(nids) == 0){ return(NULL) }
  res = c()
  for(n in nids){
    res = c(res,n,tree_order(A,n,id))
  }
  return(res)
}

parent_dist = function(cons,norm){
  D = matrix(Inf,nrow=ncol(cons),ncol=ncol(cons))
  for(i in 1:nrow(D)){
    for(j in 1:ncol(D)){
      if(i == j){next}
      D[i,j] = cnv_deviations(cons[,i],cons[,j],norm)
    }
  }
  colnames(D) = colnames(cons)
  row.names(D) = colnames(D)
  return(D)
}

cnv_deviations = function(cnv1,cnv2,norm){
  lids = cnv1 < norm
  gids = cnv1 > norm
  dev = sum(cnv1[lids] < cnv2[lids]) + sum(cnv1[gids] > cnv2[gids])
  return(dev)
}

consensus_cnv = function(M,cutoff,norm){
  cons = rep(norm,nrow(M))
  for(i in 1:nrow(M)){
    nids = M[i,] == norm
    if(sum(!nids)==0){next}
    t = table(unlist(M[i,!nids]))
    t = t/(sum(t)+sum(nids))
    if(max(t) > cutoff){
      cons[i] = as.numeric(names(t)[which.max(t)])
    }
  }
  return(cons)
}

cons_heatmap = function(rv,wnd,colors,norm,show_leg=TRUE){
  library(ComplexHeatmap)
  
  cnv_names = c(0,0.5,1,1.5,2,3)
  heat_color = c("blue","dodgerblue","gray90","palevioletred1","firebrick1","red4")
  names(heat_color) = c(cnv_names)

  ha = rowAnnotation(clone=colnames(rv),col=list(clone=colors),show_legend=FALSE)
  ht = Heatmap(t(rv),cluster_rows=FALSE,cluster_columns=FALSE,name="CNV",col=heat_color,show_row_dend = FALSE,
        show_column_names=FALSE,show_row_names=FALSE,column_split=factor(wnd[,1],levels=unique(wnd[,1])),
        column_title_rot = 45,column_title_side = "bottom",row_title = "",
        #row_split=colnames(rv),row_title_side="left",row_title_rot=0,
        left_annotation=ha,border=TRUE,column_gap=unit(2,"mm"),show_heatmap_legend=show_leg)

  p = grid.grabExpr(draw(ht,column_title="Average CNVs"))

  return(p)
}

combine_clusters = function(cnvs,clusters,thresh=0.025,norm=1,cons_thresh=.25){
  cons = matrix(0,nrow=nrow(cnvs),ncol=max(clusters))
  for(i in 1:ncol(cons)){
    cons[,i] = consensus_cnv(cnvs[,clusters==i],cons_thresh,norm)
  }
  for(i in 1:(ncol(cons)-1)){
    for(j in (i+1):ncol(cons)){
      diff = sum(cons[,i] != cons[,j])
      if(diff < thresh*nrow(cons)){
        #Merge labels, re-calculate consensus CNVs
        clusters[clusters==j] = i
        cons[,i] = consensus_cnv(cnvs[,clusters==i],cons_thresh,norm)
        cons[,j] = -1
      }
    }
  }
  nums = unique(clusters)
  clusters_o = clusters
  for(i in 1:length(nums)){
    clusters[clusters_o==nums[i]] = i
  }

  return(clusters)
}


