require(tidyverse)

### arabidopsis orthologous genes to medicago
orthologs <- read.table("mart_export_ara2med.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)
colnames(orthologs)<-c("TAIR_ID","Mt_ID","Mt_ID_percent","TAIR_ID_percent")

### arabidopsis ppi interaction
ppi <-read.table("1TAIR_inter.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)
ppi <- as.data.frame(t(apply(ppi,1,sort)),stringsAsFactors = FALSE)
ppi <- ppi[ppi$V1!=ppi$V2,]



## mutiple to multiple orthologs
Ortho_PPI <-function(at1,at2){
  mt1 <- orthologs$Mt_ID[orthologs$TAIR_ID==at1]
  mt2 <- orthologs$Mt_ID[orthologs$TAIR_ID==at2]
  return(expand.grid(mt1,mt2,stringsAsFactors = FALSE))  
}
ppi_med<-apply(ppi,1,function(x)Ortho_PPI(x[1],x[2]))
ppi_med<-do.call(rbind,ppi_med)

ppi_med<-na.omit(ppi_med)
ppi_med$type<-"PPI"
ppi_med <-unique(ppi_med)
dim(ppi_med)
ara_ppi<-ppi_med

#### rice data
orthologs <- read.table("mart_export_rice2med.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)
colnames(orthologs)<-c("Rice_ID","Mt_ID","Mt_ID_percent","Rice_ID_percent")
orthologs<-na.omit(orthologs)


ppi <-read.table("1rice_inter.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)
ppi <- as.data.frame(t(apply(ppi,1,sort)),stringsAsFactors = FALSE)
ppi <- ppi[ppi$V1!=ppi$V2,]

Ortho_PPI_rice <-function(rc1,rc2){
  mt1 <- orthologs$Mt_ID[orthologs$Rice_ID==rc1]
  mt2 <- orthologs$Mt_ID[orthologs$Rice_ID==rc2]
  return(expand.grid(mt1,mt2,stringsAsFactors = FALSE))  
}
ppi_rice<-apply(ppi,1,function(x)Ortho_PPI_rice(x[1],x[2]))
ppi_rice<-do.call(rbind,ppi_rice)

ppi_rice<-na.omit(ppi_rice)
ppi_rice$type<-"PPI"
ppi_rice <-unique(ppi_rice)
dim(ppi_rice)
rice_ppi<-ppi_rice


ppi_all <-unique(rbind(ara_ppi,rice_ppi))
ppi_all<-ppi_all[ppi_all$Var1!=""&ppi_all$Var2!="",]
write.table(ppi_all,"PPI2.txt",col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE)
