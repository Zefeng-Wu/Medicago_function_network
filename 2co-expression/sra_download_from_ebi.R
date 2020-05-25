#/usr/bin/R

require(stringr)
sra_list_file<-read.table("ena.txt",stringsAsFactors = FALSE,header = TRUE,sep="\t")
for (m in 1:nrow(sra_list_file)){
  system("date")
  sra_download_site<- sra_list_file$fastq_ftp[m]
  message(sra_list_file$run_accession[m])
  
  if (sra_list_file$library_layout[m]=="PAIRED"){
    
    r1_reads <-str_split_fixed(sra_download_site,pattern = ";",n =2 )[,1] #ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949258/SRR949258_1.fastq.gz
    r2_reads <-str_split_fixed(sra_download_site,pattern = ";",n =2 )[,2] #ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949258/SRR949258_1.fastq.gz

    r1_download<-str_replace_all(r1_reads,pattern = "ftp.sra.ebi.ac.uk",replacement = "era-fasp@fasp.sra.ebi.ac.uk:")
    r2_download<-str_replace_all(r2_reads,pattern = "ftp.sra.ebi.ac.uk",replacement = "era-fasp@fasp.sra.ebi.ac.uk:")
    #ascp -QT -l 300m -P 33001  -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR949/SRR949627/SRR949627_1.fastq.gz
    
    download_command1<-paste(c("ascp -QT -l 300m -P 33001  -i /public/home/wuzefeng/softwares/yes/etc/asperaweb_id_dsa.openssh",r1_download,"0paired/"),collapse  = " ")
    download_command2<-paste(c("ascp -QT -l 300m -P 33001  -i /public/home/wuzefeng/softwares/yes/etc/asperaweb_id_dsa.openssh",r2_download,"0paired/"),collapse  = " ")
    #message(download_command1,"---",download_command2)
    system(download_command1)
    system(download_command2)
  }
  else{
    r_download<-str_replace_all(sra_download_site,pattern = "ftp.sra.ebi.ac.uk",replacement = "era-fasp@fasp.sra.ebi.ac.uk:")
    download_command<-paste(c("ascp -QT -l 300m -P 33001  -i /public/home/wuzefeng/softwares/yes/etc/asperaweb_id_dsa.openssh",r_download,"0single"),collapse =  " ")
    
    #message("single...",download_command)
    system(download_command)
    #message(download_command)
  }
}


