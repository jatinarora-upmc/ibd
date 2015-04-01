#-------------------------------------------------------------------------------#
#Author: Jatin Arora, CeMM, Vienna
#-------------------------------------------------------------------------------#

rm(list=ls())

##-----------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------
##--------------------------------------- Functions ---------------------------------------

##----------------------------------- calculate distance ------------------------------------------
prioritize_predict <- function(trans, seeds, cands){

  # calculate probabilities at convergence after random walk iterations
  lx=limitProb(P=trans,nodes=n.full.G,seeds=seeds)
  probs=lx[order(lx,decreasing=TRUE)] # return the probabilities in decreasing order
  orig=data.frame(genes=names(probs),values=as.numeric(probs))

  ## get the prob values for candidate genes
  indexs=as.numeric(sapply(cands, function(x) which(orig$genes %in% x)))
  cands.prob=orig[indexs[!is.na(indexs)],]
  cands.prob=cands.prob[order(cands.prob$values, decreasing = T),]
  return(list(cands.prob=cands.prob,lx=lx))

}

##----------------------------------- Plot network ------------------------------------------

plot_subnetwork <- function(pp, g, lx, seeds, cands, n, lay) {
#plot_subnetwork <- function(pp, n.nodes, seeds, cands, lay) {

  ## get neighbors n.add in addition to seeds after random walk
  neighb=getNeighborhood(g,lx,seeds,n.add=n)
  n.nodes=sort(nodes(neighb$neighb))
  ## get pairs from database whose both elements are in neighbor nodes
  n.edges=c()
  for (i in 1:nrow(pp)){
    if(is.element(as.character(pp[i,"from"]), n.nodes) & is.element(as.character(pp[i,"to"]), n.nodes)){
      n.edges=rbind(n.edges,pp[i,])
    }
  }
  ## get the edges from pp to be plotted
  edges.toplot=data.frame(from=n.edges$from,to=n.edges$to)
  ind.nodes=data.frame(nodes=unique(c(as.character(edges.toplot$from), as.character(edges.toplot$to))))
  ## determine which nodes are seeds, cands and others
  grp.seeds=c()
  for(s in seeds){
    grp.seeds=c(grp.seeds,which(ind.nodes$nodes == s))
  }
  grp.cands=c()
  for(s in cands){
    grp.cands=c(grp.cands,which(ind.nodes$nodes == s))
  }
  grp.others=as.numeric(setdiff(rownames(ind.nodes), c(grp.seeds,grp.cands)))
  grps=list(seeds=grp.seeds, candidates=grp.cands, others=grp.others)
  ## change node labels to hgnc symbols
  dat=read.delim("genes/all")
  all.genes=dat[!is.na(dat$ent),]
  nodes.genes=all.genes[which(all.genes$uni %in% ind.nodes$nodes),]
  edges.toplot.genes=matrix(,nrow(edges.toplot),2)
  colnames(edges.toplot.genes)=colnames(edges.toplot)
  for(i in 1:nrow(edges.toplot)){

      fr=as.character(nodes.genes[nodes.genes$uni == as.character(edges.toplot[i,"from"]),"hgnc"])
      if(length(fr) > 0){
         edges.toplot.genes[i,"from"]=fr[1]
      } else {
         edges.toplot.genes[i,"from"] = as.character(edges.toplot[i,"from"])
      }
      to=as.character(nodes.genes[nodes.genes$uni == as.character(edges.toplot[i,"to"]),"hgnc"])
      if(length(to) > 0){
         edges.toplot.genes[i,"to"]=to[1]
      } else {
         edges.toplot.genes[i,"to"] = as.character(edges.toplot[i,"to"])
      }

  }
  library(qgraph)
  qgraph(edges.toplot.genes[complete.cases(edges.toplot.genes),], groups=grps, directed=F, pastel=T, details=T, layout=lay)

}


plot_heatmap <- function(mat,source){

   library(ggplot2)
   library(reshape)
   library(stats)

   mat.y=melt(mat)
   mat.y$X1 <- factor(as.character(mat.y$X1), levels = unique(as.character(mat.y$X1)), ordered = TRUE)
   mat.y$X2 <- factor(as.character(mat.y$X2), levels = unique(as.character(mat.y$X2)), ordered = TRUE)
   pdf(paste("plots/",source,".pdf", sep=""))
   p = ggplot(mat.y, aes(y=X1,x=X2))
   p = p + geom_tile(aes(fill=value)) +  scale_fill_gradient2(low = "red", mid = "white", high = "blue")
   p = p + ggtitle(source) + xlab("") + ylab("")
   p = p + theme(axis.text.x=element_text(angle=-90, size=10, colour="black", vjust=0.5), axis.text.y=element_text(size=10, colour="black"))
   print(p)
   dev.off()
   print("plot done")
}



##--------------------------------------- PIDs --------------------------------------------

pid <- function(trans,seeds,n.full.G){

   print(length(seeds))
   pid.genes=toupper(read.delim("pid/pid.known", header=F)$V1)
   pid=all.genes[which(all.genes$hgnc %in% pid.genes),"uni"]
   cands=intersect(as.character(pid),n.full.G)
   cands=setdiff(cands,ibd.genes$uni)
   ranks.frame=rep(0,length(cands))
   i.full.G=igraph.from.graphNEL(graphNEL = full.G,name = T,weight = T)

   for(s in seeds){
      deg=degree(i.full.G,s)
      pred=prioritize_predict(trans = trans, seeds = s, cands = cands)
      data=data.frame(rank=rownames(pred$cands.prob), gene=pred$cands.prob$genes)
      colnames(data)=c(paste(deg,s,sep="."),paste(deg,s,sep="."))
      ranks.frame=cbind(ranks.frame,data)
   }

   return(ranks.frame[,-1])
}

##--------------------------------------- Permutation Tests--------------------------------
permutation_test <- function(n.G, trans, seeds, cands, pred){

  no.permut=10000
  rank.mat=matrix(NA,length(cands),no.permut)
  rownames(rank.mat)=pred$cands.prob$genes
  for(i in 1:no.permut){
    rand.seeds=sample(n.G,length(seeds))
    rand.pred=prioritize_predict(trans, rand.seeds, cands)
    for(j in 1:nrow(pred$cands.prob)){
      s=pred$cands.prob$genes[j]
      rank.mat[j,i]=rownames(rand.pred$cands.prob[rand.pred$cands.prob$genes==s,])
    }
  }

  df=data.frame(orig=rownames(pred$cands.prob),rank.mat)
  pval=c()
  for(i in 1:nrow(df)){
      val=(length(which(as.numeric(as.character(unlist(df[i,-1]))) <= as.numeric(as.character(unlist(df[i,"orig"]))))))/no.permut
      pval=c(val,pval)
  }
  write.table(data.frame(df,pval),"pval.txt",row.names=T,col.names=T,append=F,quote=F,sep="\t")
  print("permutation test done")
}

##----------------------------------- Heat Map ------------------------------------------
format_mat <- function(x.mat,source){
   ind=c()
   for(i in 1:nrow(x.mat)){
      if(all(is.na(x.mat[i,]))){
         ind=c(ind,i)
      }
   }
   x.mat=x.mat[-ind,]
   ind=c()
   for(i in 1:ncol(x.mat)){
      if(all(is.na(x.mat[,i]))){
         ind=c(ind,i)
      }
   }
   x.mat=x.mat[,-ind]
   for(i in 1:nrow(x.mat)){
      rownames(x.mat)[i] = as.character(all.genes[all.genes$ent == (rownames(x.mat)[i]),"hgnc"])[1]
   }
   for(i in 1:ncol(x.mat)){
      colnames(x.mat)[i] = as.character(all.genes[all.genes$ent == (colnames(x.mat)[i]),"hgnc"])[1]
   }
   plot_heatmap(x.mat,source)

}

##----------------------------------- Novarino et al. ------------------------------------------
hsp.test <- function(trans,n.full.G){

   path=read.delim("pathways/cpdb.entez.tab copy")
   dat=read.delim("genes/all")
   all.genes=dat[!is.na(dat$ent),]
   hsp.genes=as.character(read.delim("hsp/seeds", header=F)$V1)
   hsp.seeds.uni=as.character(all.genes[which(all.genes$hgnc %in% hsp.genes),"uni"])
   hsp.seeds.ent=as.character(all.genes[which(all.genes$hgnc %in% hsp.genes),"ent"])
   lx=limitProb(P=trans,nodes=n.full.G,seeds=hsp.seeds.uni)
   probs=lx[order(lx,decreasing=TRUE)]
   hsp.cands=as.character(read.delim("hsp/cands", header=F)$V1)
   hsp.cands.ent=as.character(all.genes[which(all.genes$hgnc %in% hsp.cands),"ent"])

   go.mat=matrix(,length(hsp.seeds.ent),length(hsp.cands.ent))
   path.mat=matrix(,length(hsp.seeds.ent),length(hsp.cands.ent))
   rownames(go.mat)=hsp.seeds.ent
   colnames(go.mat)=hsp.cands.ent
   rownames(path.mat)=hsp.seeds.ent
   colnames(path.mat)=hsp.cands.ent

   for(ontologies in c("MF","BP","CC")){

      for(g1 in hsp.seeds.ent){
         path1=as.character(path[which(regexpr(pattern = paste("\\b",g1,"\\b",sep=""), text = path$entrez_gene_ids) > 0), "external_id"])
         for(g2 in hsp.cands.ent){
            sim=geneSim(gene1 = g1, gene2 = g2, ont = ontologies, organism="human", measure = "Wang", drop = NULL, combine="BMA")
            go.mat[g1,g2]=sim[[1]]
            path2=as.character(path[which(regexpr(pattern = paste("\\b",g2,"\\b",sep=""), text = path$entrez_gene_ids) > 0), "external_id"])
            if(length(path1) > 0 & length(path2) > 0){
               path.mat[g1,g2]=length(intersect(path1,path2))/length(union(path1,path2))
            }
         }
      }

      format_mat(go.mat,paste("hsp.go.",ontologies,sep=""))
   }

   #format_mat(path.mat,"hsp.kegg")

}

##-----------------------------------------------------------------------------------------
##--------------------------------------- Main --------------------------------------------
##-----------------------------------------------------------------------------------------

## load libraries
library(SparseM)
library(RBGL)
library(Rgraphviz)
library(igraph)
library(GOSemSim)
library(org.Hs.eg.db)
source("expandGraph.R")

## load interactions and thier weights
#pp=read.delim("interactions/global/bg.go.sc.enriched.0.01.resnik.bwa",header=F, sep="\t")
pp=read.delim("funconet/funconet1sept.tab", sep="\t", header=T)
#colnames(pp)=c("from","to","bp","mf","cc")
#pp[is.na(pp)]=0
#pp$go=pp$bp + pp$mf + pp$cc
#pp$kegg.sc=(pp$kegg.sc-min(pp$kegg.sc))/(max(pp$kegg.sc)-min(pp$kegg.sc))
#pp$ex.sc = abs(pp$ex.sc)
#pp=subset(pp, pp[,"go.sc"] > 0)

#coeff.go=1
#coeff.kegg=1
#coeff.ex=1

## make the graph
#weight= coeff.go*pp$go + 1
#pp$weight=1
ft=data.frame(from=pp$source, to=pp$target) #, weight=pp$weight)
#ft=data.frame(from=pp$from, to=pp$to, weight=weight)
#full.g <- ftM2graphNEL(as.matrix(ft), edgemode="undirected", W = weight)
full.g <- ftM2graphNEL(as.matrix(ft), edgemode="undirected")
#g=graphBAM(ft,edgemode="undirected",ignore_dup_edges=T)
#full.g=as(g, "graphNEL")

## first get all connected components of the graph, take the sub-graph which is most connected
full.cc <- connectedComp(full.g)
full.G <- subGraph(full.cc[[1]],full.g)

## get nodes from most connected sub-graph
n.full.G <- nodes(full.G) # n = no. of nodes
full.am <- graph2SparseM(full.G,useweights=F)
full.P <- adjmat2transprob(full.am)
#for (k in 1:(length(full.am@ia)-1)){
#    s = sum(full.am@ra[full.am@ia[k]:(full.am@ia[k+1]-1)])
#    if(s > 0){
#      full.am@ra[full.am@ia[k]:(full.am@ia[k+1]-1)] = full.am@ra[full.am@ia[k]:(full.am@ia[k+1]-1)]/s
#    }
#}

#trans=t(full.am)
dat=read.delim("genes/all")
all.genes=dat[!is.na(dat$ent),]

## take only those which exist in database
ibd.genes=read.delim("genes/seeds")
seeds=intersect(ibd.genes$uni,n.full.G)
## also for candiate genes
cands.genes=read.delim("genes/cands")
cands=intersect(cands.genes$uni,n.full.G)

## prioritze candiadates
#orig.pred=prioritize_predict(full.P, seeds, cands)

## permutation test
#permutation_test(n.G = n.full.G, trans = full.P, seeds = seeds, cands = cands, pred = orig.pred)

## measure proximity of seeds from pid genes
dist.pid=pid(trans = full.P, seeds = cands, n.full.G = n.full.G)
write.table(dist.pid, "dist.pid.tsv", row.names=F, col.names=T, append=F, quote=F, sep="\t")

## after selecting the best candidates, select the subnetwork
seeds="P05156"
pred=prioritize_predict(trans = full.P, seeds = seeds, cands = cands)
lx=pred$lx
probs=lx[order(lx,decreasing=TRUE)]

i.full.G=igraph.from.graphNEL(graphNEL = full.G,name = T,weight = T)
degree(i.full.G,seeds)
sg=induced.subgraph(graph = i.full.G, v = names(probs[1:50]))
label.comm=label.propagation.community(sg)
clust.res=data.frame(nodes=label.comm$names, clust=label.comm$membership)
nrow(clust.res[clust.res$clust==clust.res[clust.res$nodes==seeds,"clust"],])
plot(label.comm,sg)


## plot the subnetwork
plot_subnetwork(pp = ft, g = full.G, lx = pred$lx, seeds = seeds, cands = cands, n = 20, lay = "spring")
## extract sub network, include top ranking candidates in seeds
#lx=orig.pred$lx
#probs=lx[order(lx,decreasing=TRUE)]
#distribution(probs,length(probs))
#plot_subnetwork(ft,names(probs)[1:25],seeds,cands,lay = "spring")
#subnetwork_go_enrichment(pp = ft,g = full.G,probs = probs,cands = cands)


##-------------------------------------- sub-network GO enrichment --------------------
subnetwork_go_enrichment <- function(pp,g,probs,cands){

   ## calculate how much is average similarity between seed and randomly choosen nodes. this would help to decide if seed subnetwork is more enriched for seed go terms or just by chance
   ##------------ Function
   permutation <- function(g1.ent,n.add){
      mean.en.rand=c()
      ontologies="MF"
      for(i in 1:100){
         rands=sample(intersect(all.genes$uni,n.full.G),n.add)
         go.rand=c()
         for(nod in rands){
            g2.ent=all.genes[all.genes$uni == nod,"ent"][1]
            sim=geneSim(gene1 = g1.ent, gene2 = g2.ent, ont = ontologies, organism="human", measure = "Wang", drop = NULL, combine="BMA")
            if(!is.na(sim)){
               go.rand=c(go.rand,sim$geneSim)
            }
         }
         mean.en.rand=rbind(mean.en.rand,mean(go.rand))
      }
      frame.rand=data.frame(mean.en.rand)
      colnames(frame.rand)=c("go.ran")
      return(mean(frame.rand$go.ran))
   } # permutation

   distribution <- function(probs,limit){

      ## kegg pahtways
      path=read.delim("pathways/cpdb.entez.tab copy")
      ## go ontologies
      dat=read.delim("genes/all")
      all.genes=dat[!is.na(dat$ent),]
      seed=names(probs[1])
      print(seed)
      #other.nodes=intersect(all.genes$uni,names(probs)[-1])
      other.nodes=names(probs)[-1]
      g1.ent=all.genes[all.genes$uni == seed,"ent"]
      g1.ens=as.character(all.genes[all.genes$uni == seed,"ens"])
      path1=as.character(path[which(regexpr(pattern = paste("\\b",g1.ent,"\\b",sep=""), text = path$entrez_gene_ids) > 0), "external_id"])
      go.frame=c()
      ontologies="MF"
      for(n in other.nodes[1:limit]){
         g2.ent=all.genes[all.genes$uni == n, "ent"][1]
         if(!is.na(g2.ent)){
            path.intr=NA
            if(length(path1) > 0){
               path2=as.character(path[which(regexpr(pattern = paste("\\b",g2.ent,"\\b",sep=""), text = path$entrez_gene_ids) > 0), "external_id"])
               path.intr=length(intersect(path1,path2))/length(path1)
            }
            sim=geneSim(gene1 = g1.ent, gene2 = g2.ent, ont = ontologies, organism="human", measure = "Wang", drop = NULL, combine="BMA")
            go.frame=rbind(go.frame, c(n,sim[[1]],as.numeric(probs[n]), path.intr))
         }
      }

      frame <- data.frame(node=go.frame[,1],go=as.numeric(go.frame[,2]),prob=as.numeric(go.frame[,3]),path=as.numeric(go.frame[,4]))
      print(dim(frame))
      frame[is.na(frame)]<-0
      library(ggplot2)
      frame$node=factor(as.character(frame$node), levels=unique(as.character(frame$node)), ordered = TRUE)
      pdf(paste("plots/",seed,".pdf",sep=""))
      p <- ggplot(frame)
      p <- p + geom_point(aes(node,go,colour="green"))# + geom_smooth(aes(node,go), method="loess",colour="green",size=2)
      p <- p + geom_point(aes(node,path,colour="blue"))# + geom_smooth(aes(node,path,group=1),method="loess",colour="red",size=2)
      p <- p + theme_bw() + xlab("rank in increasing order") + ylab("similarity") + ggtitle(all.genes[all.genes$uni == seed,"hgnc"])
      p <- p + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=10), axis.ticks.x=element_blank())
      p <- p + scale_colour_discrete(name  ="Legend",labels=c("Pathway", "GO"))
      print(p)
      dev.off()
   }
   ##------------ MAIN

   seed=names(probs[1])
   other.nodes=names(probs)[-1]
   ## convert graphNEL to igraph
   i.full.G=igraph.from.graphNEL(graphNEL = g,name = T,weight = T)
   deg=degree(i.full.G,seed)
   print(deg)
   dat=read.delim("genes/all")
   all.genes=dat[!is.na(dat$ent),]

   ## expression data
   #exp=read.delim("expression/atlas/e-mtab-1733.tab")
   #colnames(exp)=sub("[.].*","",sub("FPKM.","",colnames(exp)))

   ## get neighbours and thier nodes
   neigh.nodes=names(probs)[-(1:(deg+1))]
   ## convert ids
   g1.ent=all.genes[all.genes$uni == seed,"ent"]
   g1.ens=as.character(all.genes[all.genes$uni == seed,"ens"])
   ontologies="MF"
   ## set arrays
   go=c()
   ex=c()
   sig.neighs=names(probs)[2:(deg+1)]
   for(n in sig.neighs){
      g2.ent=all.genes[all.genes$uni == n, "ent"][1]
      sim=geneSim(gene1 = g1.ent, gene2 = g2.ent, ont = ontologies, organism="human", measure = "Wang", drop = NULL, combine="BMA")
      if(!is.na(sim)){
         go=c(go,sim$geneSim)
      }
   }
   ## loop over all neighbours
   add=1
   for(n in neigh.nodes){

      go.old=go
      print(mean(go.old))
      g2.ent=all.genes[all.genes$uni == n,"ent"][1]
      g2.ens=as.character(all.genes[all.genes$uni == n,"ens"][1])
      sim=geneSim(gene1 = g1.ent, gene2 = g2.ent, ont = ontologies, organism="human", measure = "Wang", drop = NULL, combine="BMA")
      if(!is.na(sim)){
         go=c(go,sim$geneSim)
      }
      #ex=c(ex,cor(as.numeric(exp[exp$ensembl == g1.ens,-c(1)]), as.numeric(exp[exp$ensembl == g2.ens,-c(1)]), method="spearman"))

      go.rand=permutation(g1.ent,(deg+add))
      print(mean(go))
      print(go.rand)
      print("---")
      if((mean(go) > go.rand) & ((mean(go) > mean(go.old)) | (mean(go.old) - mean(go) < 0.05))){
         #& abs(mean(go)-mean(go.old)) < 0.01
         sig.neighs=c(sig.neighs,n)
         add=add+1
      }
      else{
         print("convergence")
         break
      }
   }
   plot_subnetwork(pp,c(seed,sig.neighs),seed,cands,"spring")

}

##--------------------------------------- call functions--------------------------------
##--------------------------------------------------------------------------------------

#permutation_test(n.full.G, trans, seeds, cands, orig.pred)


