#-------------------------------------------------------------------------------#
#Author: Jatin Arora, CeMM, Vienna
#-------------------------------------------------------------------------------#

rm(list=ls())

library(SparseM)
library(RBGL)
library(Rgraphviz)
library(igraph)
#library(GOSemSim)
#library(org.Hs.eg.db)

source("expandGraph.R")
pp=read.delim("funconet/funconet1sept.tab", sep="\t", header=T)
ft=data.frame(from=pp$source, to=pp$target)

## make the graph
g=graph.data.frame(d = ft, directed = F)
#trans.g=transitivity(graph = g, type = "global")
#print(trans.g)

## read mapping files and known pid genes
dat=read.delim("genes/all")
all.genes=dat[!is.na(dat$ent),]
pid.genes=toupper(read.delim("pid/pid.known", header=F)$V1)
pid=as.character(all.genes[which(all.genes$hgnc %in% pid.genes),"uni"])

cands.all=read.delim("genes/cands")
cands=as.character(cands.all$uni)
## --------------------------- perform clustering on line graph ------------------------------------

## convert to line graph
line.g=line.graph(graph = g)
#trans.line.g=transitivity(graph = line.g, type = "global")
#print(trans.line.g)

## ---------------------------- label propagation clustering -----------------------------------

label.clust=label.propagation.community(graph = line.g)
clust.info=data.frame(get.edgelist(g),as.numeric(V(line.g)),label.clust$membership)
colnames(clust.info)=c("node1","node2","node.id","cluster")

## iterate for each node of interest
for(v in cands){

   print(v)
   clust.ids=unique(clust.info[clust.info$node2 == v | clust.info$node1 == v,"cluster"])
   ## iterate over each cluster, however check for the one which has the maximum instances of node of interest
   for(i in 1:length(clust.ids)){
      edges=clust.info[clust.info$cluster == clust.ids[i], c("node1","node2")]
      print(c(nrow(edges),length(intersect(union(edges$node1, edges$node2), pid))))
   }

}



clust.g=graph.data.frame(d = edges, directed = F)
plot(clust.g)

## ---------------------------- MCL clustering ------------------------------------

#write.table(get.edgelist(graph = line.g, names = T),"mcl/line.mcl.input",col.names=F,row.names=F,quote=F,append=F,sep=" ")
#system("/home/jarora/Setups/mcl-14-137/local/bin/mcl mcl/line.mcl.input --abc -I 1.8 -overlap keep -o mcl/line.1.8.clust")

## load cluster file
clust=read.delim("mcl/line.1.8.clust", header = F)

## convert back from line graph to original graph
map = data.frame(get.edgelist(g),as.numeric(V(line.g)))
colnames(map)=c("node1","node2","node.id")

for(ind in 1: nrow(clust)){
   cl = clust[ind,]
   cl = cl[!is.na(cl)]
   clust.ft = map[which(map$node.id %in% cl),c("node1","node2")]
   clust.n = union(clust.ft$node1, clust.ft$node2)
   len.c = length(intersect(clust.n,cands)) # if cluster has any candidate
   if(len.c > 0){
      len.p = length(intersect(clust.n, pid)) # if cluster has both candidate and pid gene
      if(len.p > 0){
         print(paste(len.c, len.p, ind))
      }
   }
}

clust.g=graph.data.frame(d = clust.ft, directed = F)
plot(clust.g)

