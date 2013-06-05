tabla_de_pares <- "MAMMOTH_T0630-D1_noorigconf.out"  ## Input file name

library(apcluster)
library(mefa)
library(compiler)

enableJIT(3)

tabla <- read.table(tabla_de_pares,sep=";",dec=".",stringsAsFactors=FALSE,na.strings="no_data",header=TRUE)
tabla <- subset(tabla,select=c("PDB_1","PDB_2","Target","RMSD"))
niveles <- unique(tabla$Target)

textclust <- function(clust,x,minimo){
  cut  <- cutree(clust,h=minimo)
  nombre <- paste(c(x,".hclusterRMSD",".txt"),collapse="")
  tabla <- data.frame(PDB=character(length(cut)))
  tabla$PDB <- names(cut)
  for(i in seq_along(cut)){
    tabla$cluster[i] <- cut[[i]]
  }
  write.table(x=tabla[order(tabla$cluster),],file=nombre,row.names=FALSE,quote=FALSE,sep="\t")
}

for(x in niveles){tryCatch({
  dos <- subset(tabla, Target==x)
  if(exists("dos") & nrow(dos)<2) next
  tryCatch({dif <- vec2dist(dos$RMSD,size=length(unique(dos$PDB_2))+1,labels=unique(c(dos$PDB_1[1],dos$PDB_2)))
            #sim <- as.matrix(vec2dist(dos$TMscore_TMvalue,size=length(unique(dos$PDB_2))+1,labels=unique(c(dos$PDB_1[1],dos$PDB_2))))
  },error=function(e){NULL})

  if(exists("sim")){ tryCatch({
    ap <- apcluster(as.matrix(sim),q=0.5)
    png(paste(c(x,".apclusterTM.png"),collapse=""), width=600, height=600)
    plot(ap,sim, main=paste("Affinity Propagation Clustering for ",x,". Similarity: TM score. q:0.5",sep=""),mar=c(7,7))
    dev.off()
  },error=function(e){NULL}) }

  if(exists("dif") & sum(is.na(dif))==0){
    tryCatch({
      clust <- hclust(dif,method="complete")
      altura <- clust$height
      cut <- 1
      if(min(altura)<1){
        if(max(altura)<1){ cut <- 0 }
        minimo <- 1
      }else{
	minimo <- min(altura) + 0.1 }
      textclust(clust,x,minimo)
      png(paste(c(x,".hclusterRMSD.png"),collapse=""),width=8000, height=8000)
      if(cut){  plot(clust,hang=-1,main=x,xlab=paste("Distance: RMSD from Mammoth, clusters at ",as.character(minimo)," Ångströms"))
                clust <- rect.hclust(clust,h=minimo) 
              }else{
                plot(clust,hang=-1,main=x,xlab=paste("Distance: RMSD from Mammoth, 1 cluster at ",as.character(minimo)," Ångström")) 
              }
      dev.off()
      rm(list=c("altura","clust","minimo")); gc()
    },error=function(e){NULL})
  }
  rm(list=c("dos","dif","ap")); gc()
},error=function(e){NULL})}



