library(RpsiXML)
library(stringr)
file_list <- dir("1other_data/intact/arath/")

df <- c()
for (file_xml in file_list){
  message(file_xml)
  intactxml <- file.path("1other_data/intact/arath/", file_xml)
  x <- psimi25XML2Graph(intactxml, psimi25source = INTACT.PSIMI25,verbose=FALSE)
  y<-names(x@edgeData)
  temp_int_df <- as.data.frame(str_split_fixed(y,pattern = "\\|",n = 2),stringsAsFactors = FALSE)
  temp_int_df <-as.data.frame(t(apply(temp_int_df,1,sort)),stringsAsFactors = FALSE)
  df <-rbind(df,temp_int_df)  
  }

df<-unique(df)
df<-df[df$V1!=df$V2,]


##
df2<-df
uniprot2gene <-read.table("1other_data/ARATH_3702_idmapping.dat",header = FALSE,sep="\t",quote = "",stringsAsFactors = FALSE) # downlaod from uniprot ftp
uniprot2gene$V1<-str_split_fixed(uniprot2gene$V1,pattern = "-",n = 2)[,1]

######## keep ATXG

uniprot2gene<-subset(uniprot2gene,grepl(pattern = "AT[0-9]G",x = uniprot2gene$V3))
uniprot2gene<-subset(uniprot2gene,startsWith(prefix = "AT",uniprot2gene$V3))

uniprot2gene$V1<-toupper(uniprot2gene$V1)
uniprot2gene$V3<-toupper(uniprot2gene$V3)

##### remove isoform transcripts
uniprot2gene$V3<-str_split_fixed(uniprot2gene$V3,pattern = "\\.",n = 2)[,1]
uniprot2gene<-unique(uniprot2gene[,c(1,3)])

# convert uniprot into tair gene ids 
df2$V1<-uniprot2gene$V3[match(df2$V1,uniprot2gene$V1)]
df2$V2<-uniprot2gene$V3[match(df2$V2,uniprot2gene$V1)]

## small letter to big letter
df2$V1<-toupper(df2$V1)
df2$V2<-toupper(df2$V2)

### remove na 
df2<-df2[!is.na(df2$V1)&!is.na(df$V2),]

## remove "NA"
df2<-na.omit(df2)

write.table(df2,file = "1other_data/ara_ppi_from_intact.txt",sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
