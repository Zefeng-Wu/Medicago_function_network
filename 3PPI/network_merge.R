library(igraph)
library(stringr)

setwd("/data/wzf/Medicago_individule_network/")
### gene info 
gene_info <- read.csv("gene_name_mapping.csv",stringsAsFactors = FALSE) 
gene_info$ID.Mt4.0 <- str_replace_all(gene_info$ID.Mt4.0,pattern = "Medtr",replacement = "MTR_")
gene_info <-subset(gene_info,gene_info$ID.Mt4.0!="#N/A")


## 1.co-expression network

co_exp <-read.table("gene_gene_cor_0.8.txt",header = TRUE,stringsAsFactors=FALSE)
co_exp$row<-str_replace_all(co_exp$row,pattern = "Medtr",replacement = "MTR_")
co_exp$col<-str_replace_all(co_exp$col,pattern = "Medtr",replacement = "MTR_")
co_exp$Co_express <- "Co-express"
co_exp_g<-graph_from_data_frame(d = co_exp,directed = FALSE)
co_exp_g<-simplify(co_exp_g,remove.multiple = TRUE,remove.loops = TRUE,edge.attr.comb = "random")

## 2. ppi
ppi <-read.table("PPI2.txt",header = TRUE,stringsAsFactors = FALSE)
ppi <-ppi[-3]
ppi$PPI<-"PPI"
ppi_g <- graph_from_data_frame(d = ppi,directed = FALSE)
ppi_g<-simplify(ppi_g,remove.multiple = TRUE,remove.loops = TRUE,edge.attr.comb = "random")


## 3. ppi_by_domain
ppi_domain <-read.table("5PPI_by_domain.txt",header = TRUE,stringsAsFactors = FALSE)
ppi_domain$Domain_interact<-"Domain_interact" 
ppi_domain_g <- graph_from_data_frame(d = ppi_domain,directed = FALSE)
ppi_domain_g<-simplify(ppi_domain_g,remove.multiple = TRUE,remove.loops = TRUE,edge.attr.comb = "random")


## 4. co_evo
co_evo <- read.table("2co_evo0.8.txt",header = TRUE,stringsAsFactors = FALSE)
co_evo$Co_evo <- "Co-evo"
co_evo_g <-graph_from_data_frame(d = co_evo,directed = FALSE)

### merge all network
u_g<-union(co_exp_g,ppi_g,ppi_domain_g, co_evo_g)
el<-get.edge.attribute(u_g)


## visulization indvidule gene with their neibour network
Sub_Net_Single_Gene_Vis<-function(focus_gene,highlights){ # della protein
  require(networkD3)
  neib <- neighborhood(u_g,nodes  = focus_gene)
  sub_g <- induced_subgraph(u_g,vids = names(neib[[1]]))

  ### Whether need have ppi evidence (optional)
  sub_g <- subgraph.edges(sub_g, E(sub_g)[!is.na(E(sub_g)$PPI)])  # limit the edge by PPI type
  
  ## replce gene names by synomous names if available
  temp_df <- data.frame(gn=names(V(sub_g)),stringsAsFactors = FALSE)
  temp_df$syno<- ifelse(temp_df$gn%in%gene_info$ID.Mt4.0,
                        gene_info$ACRONYM[match(temp_df$gn,gene_info$ID.Mt4.0)],
                        temp_df$gn)
  V(sub_g)$label <- temp_df$syno
  
  
  #plot.igraph(sub_g,vertex.size=2,vertex.label.cex=0.2,layout=layout.circle)
  ## plot netowrk by d3
  karate_d3 <- igraph_to_networkD3(sub_g,group = ifelse(names(V(sub_g))%in%highlights,"a","b"))
  karate_d3$nodes$label <- temp_df$syno
  karate_d3$nodes$size <- 10
  
  ColourScale <- 'd3.scaleOrdinal()
            .domain(["a", "b"])
           .range(["#FF6900", "#694489"]);'
  
  forceNetwork(Links = karate_d3$links,
               Nodes = karate_d3$nodes, 
               Source = 'source', 
               Target = 'target', 
               NodeID = 'label', 
               Group = 'group',
               zoom = TRUE,
               height = 1000,
               width = 1000,
               #linkDistance = networkD3::JS("function(d) { return 5*d.value; }"),
               colourScale = ColourScale,
               opacity = 0.9,
               opacityNoHover = 1,
               Nodesize = "size",
               fontSize = 10,
               legend = T)
} 
focus_gene ="MTR_3g065980"
focus_gene ="MTR_7g113680" #nf-y
focus_gene = "MTR_8g043970" #ccamk
focus_gene = "MTR_5g099060" #nin
focus_gene = "MTR_5g026850"
Sub_Net_Single_Gene_Vis(focus_gene,focus_gene)

### visulization gene list in a network
Sub_Net_Gene_List_Vis<-function(gene_list,highlights){ # della protein
  
  sub_g <- induced_subgraph(u_g,vids = gene_list[gene_list%in%names(V(u_g))])
  
  ### Whether need have ppi evidence (optional)
  #sub_g <- subgraph.edges(sub_g, E(sub_g)[!is.na(E(sub_g)$PPI)])  # limit the edge by PPI type
  
  ## replce gene names by synomous names if available
  temp_df <- data.frame(gn=names(V(sub_g)),stringsAsFactors = FALSE)
  temp_df$syno<- ifelse(temp_df$gn%in%gene_info$ID.Mt4.0,
                        gene_info$ACRONYM[match(temp_df$gn,gene_info$ID.Mt4.0)],
                        temp_df$gn)
  V(sub_g)$label <- temp_df$syno
  
  
  #plot.igraph(sub_g,vertex.size=2,vertex.label.cex=0.2,layout=layout.circle)
  ## plot netowrk by d3
  karate_d3 <- igraph_to_networkD3(sub_g,group = ifelse(names(V(sub_g))%in%highlights,"a","b"))
  karate_d3$nodes$label <- temp_df$syno
  karate_d3$nodes$size <- 10
  
  ColourScale <- 'd3.scaleOrdinal()
            .domain(["a", "b"])
           .range(["#FF6900", "#694489"]);'
  
  forceNetwork(Links = karate_d3$links,
               Nodes = karate_d3$nodes, 
               Source = 'source', 
               Target = 'target', 
               NodeID = 'label', 
               Group = 'group',
               zoom = TRUE,
               height = 1000,
               width = 1000,
               linkDistance = 100,
               colourScale = ColourScale,
               opacity = 0.9,
               opacityNoHover = 1,
               Nodesize = "size",
               fontSize = 10,
               legend=T,
               charge = -1000)
} 

##########nodule genes#######
#nodule_genes
##################
nodule_genes <-read.table("1other_data/1.5nodule_genes.txt",header = FALSE,sep="\t",stringsAsFactors = FALSE)
nodule_genes$V2<-str_split_fixed(nodule_genes$V2,pattern = "/",n = 2)[,1]
nodule_genes <- subset(nodule_genes,nodule_genes$V3!="NCR")
nodule_genes <- subset(nodule_genes,nodule_genes$V1%in%names(V(u_g)))
nodule_genes$color <- factor(x=nodule_genes$V3,
                             levels = unique(nodule_genes$V3),
                             labels = c("#00008B", "#FF4040", "#53868B", "#E066FF", "#EE9572", "#FFFF00", "#8968CD", "#8B7500", "#458B00", "#8B8378"))
nodule_genes$color<-as.character(nodule_genes$color)

sub_g <- induced_subgraph(u_g,vids = nodule_genes$V1)
#sub_g <- subgraph.edges(sub_g, E(sub_g)[!is.na(E(sub_g)$PPI)]) 

temp_df <- data.frame(gn=names(V(sub_g)),stringsAsFactors = FALSE)
temp_df$syno<- ifelse(temp_df$gn%in%nodule_genes$V1,
                      nodule_genes$V2[match(temp_df$gn,nodule_genes$V1)],
                      temp_df$gn)
V(sub_g)$label <- temp_df$syno
V(sub_g)$color <- nodule_genes$color[match(names(V(sub_g)),nodule_genes$V1)]

E(sub_g)[!is.na(E(sub_g)$Co_express)]$color<-"purple" 
E(sub_g)[!is.na(E(sub_g)$Co_evo)]$color<-"darkgreen" 
E(sub_g)[!is.na(E(sub_g)$Domain_interact)]$color<- "steelblue"
E(sub_g)[!is.na(E(sub_g)$PPI)]$color<-"pink"  
E(sub_g)$width=2

par(mar=c(5,5,5,5))
plot.igraph(sub_g,
            layout=layout_on_sphere,
            vertex.size=8,
            vertex.label.cex=1,
            vertex.label.color="black",
            #vertex.color=
            vertex.label.dist=1,
            vertex.frame.color=NA)

legend("bottom",
       legend=unique(nodule_genes$V3),
       pch=19, #shape
       box.lty=2, # 
       pt.cex= 3, #lines size 
       cex=0.8, #box size
       col=c("#00008B", "#FF4040", "#53868B", "#E066FF", "#EE9572", "#FFFF00", "#8968CD", "#8B7500", "#458B00", "#8B8378"),
       y.intersp=1.8,
       ncol=4,
       yjust=-2,
       inset=-0.1)
## legend for line
legend("top",
       legend=c("Co-express","Co_evo","Domain_interact","PPI"),
       #pch=19, #shape
       lty="solid",
       box.lty=2, # 
       pt.lwd= 10, #lines size 
       cex=0.8, #box size
       col=c("purple", "darkgreen", "steelblue", "pink"),
       y.intersp=1.8,
       ncol=4,
       yjust=-2,
       inset=-0.1)


##################################
### test the tfs whoses DNA binding site are enriched in nodule genes promoters
##################################
focus_tf_list <-c("Medtr3g102980","Medtr1g086250","Medtr2g014300","Medtr8g089895","Medtr5g054300","Medtr2g026725","Medtr3g435480",
                   "Medtr1g017350","Medtr7g118330","Medtr7g118330","Medtr4g094908","Medtr1g086250")

focus_tf_list <-str_replace_all(focus_tf_list,pattern = "Medtr","MTR_")

for (tf in focus_tf_list){
  message(tf)
  neibs <-neighborhood(u_g,nodes = tf)
  
  ## neibs must be in nodule gene lsit
  sub_g <- induced_subgraph(u_g,vids = names(neibs[[1]])[names(neibs[[1]])%in%nodule_genes$V1])
  
  ### color in different interact type (optional)
  E(sub_g)[!is.na(E(sub_g)$Co_express)]$color<-"purple" 
  E(sub_g)[!is.na(E(sub_g)$Co_evo)]$color<-"darkgreen" 
  E(sub_g)[!is.na(E(sub_g)$Domain_interact)]$color<- "steelblue"
  E(sub_g)[!is.na(E(sub_g)$PPI)]$color<-"pink"  
  E(sub_g)$width=2
  message("ok")
  ## replce gene names by synomous names if available
  temp_df <- data.frame(gn=names(V(sub_g)),stringsAsFactors = FALSE)
  temp_df$syno<- ifelse(temp_df$gn%in%gene_info$ID.Mt4.0,
                        gene_info$ACRONYM[match(temp_df$gn,gene_info$ID.Mt4.0)],
                        temp_df$gn)
  V(sub_g)$label <- temp_df$syno
  V(sub_g)$color <- nodule_genes$color[match(names(V(sub_g)),nodule_genes$V1)]
  
  plot.igraph(sub_g,layout=layout.fruchterman.reingold,
              vertex.size=8,
              vertex.label.cex=1,
              vertex.label.color="black",
              #vertex.color=
              vertex.label.dist=1,
              vertex.frame.color=NA)
 
  
  ### legend for nodes
  legend("bottom",
         legend=unique(nodule_genes$V3),
         pch=19, #shape
         box.lty=2, # 
         pt.cex= 3, #lines size 
         cex=0.8, #box size
         col=c("#00008B", "#FF4040", "#53868B", "#E066FF", "#EE9572", "#FFFF00", "#8968CD", "#8B7500", "#458B00", "#8B8378"),
         y.intersp=1.8,
         ncol=4,
         yjust=-2,
         inset=-0.1)
  
  #legend for lines
  legend("top",
         legend=c("Co-express","Co_evo","Domain_interact","PPI"),
         #pch=19, #shape
         lty="solid",
         box.lty=2, # 
         pt.lwd= 10, #lines size 
         cex=0.8, #box size
         col=c("purple", "darkgreen", "steelblue", "pink"),
         y.intersp=1.8,
         ncol=4,
         yjust=-2,
         inset=-0.1)
}

#######################
###### gene go and kegg
#######################





#######################
###### gene go and kegg
#######################



