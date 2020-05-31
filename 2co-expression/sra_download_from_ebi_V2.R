#/usr/bin/R

require(stringr)
sra_list_file<-read.table("ena.txt",stringsAsFactors = FALSE,header = TRUE,sep="\t")
sra_list_file_md5<-read.table("ena_md5.txt",header = TRUE,sep = "\t")

for (m in 1:nrow(sra_list_file)){ #
  system("date")
  sra_download_site<- sra_list_file$fastq_ftp[m]
  message("Run id is : ",m," ---- ", sra_list_file$run_accession[m])
  if (sra_download_site==""){next;} # empty files
  if (sra_list_file$library_layout[m]=="PAIRED"){
    
    r1_reads <-str_split_fixed(sra_download_site,pattern = ";",n =3 )[,1] #ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949258/SRR949258_1.fastq.gz
    message("r1 reads is: ",basename(r1_reads))
    r1_download<-str_replace_all(r1_reads,pattern = "ftp.sra.ebi.ac.uk",replacement = "era-fasp@fasp.sra.ebi.ac.uk:")
    #ascp -QT -l 300m -P 33001  -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR949/SRR949627/SRR949627_1.fastq.gz
    
    download_command1<-paste(c("ascp -QT -l 50m -P 33001  -i /public/home/wuzefeng/softwares/yes/etc/asperaweb_id_dsa.openssh",r1_download,"0paired/"),collapse  = " ")
    #message(download_command1,"---",download_command2)
    
	  if(! file.exists(paste("0paired/",basename(r1_reads),sep=""))){
     system(download_command1)}
    
    true_md5_paired<- sra_list_file_md5$fastq_md5[sra_list_file_md5$run_accession==sra_list_file$run_accession[m]]
    r1_reads_true_md5 <-str_split_fixed(true_md5_paired,pattern = ";",n =3 )[,1]
    r2_reads_true_md5 <-str_split_fixed(true_md5_paired,pattern = ";",n =3 )[,2]
    
    r1_reads_file <- paste("0paired/",basename(r1_reads),sep="")
    r1_reads_download_md5_command<-paste(c("md5sum",r1_reads_file),collapse = " ")
    r1_reads_download_md5_result <- system(r1_reads_download_md5_command,intern = TRUE)
    r1_reads_download_md5<-str_split_fixed(r1_reads_download_md5_result,pattern = " ",n = 2)[1]
    
    message("r1_reads true md5 is: ",r1_reads_true_md5)
    message("r1_reads download md5 is: ",r1_reads_download_md5)
    
    if (r1_reads_true_md5!= r1_reads_download_md5){
      message("r1 download failed!! Must REMOVE and redownload!")
      system(paste("rm 0paired/",basename(r1_reads),sep=""))
      system(download_command1)
    }
    
####### r2 reads handle #############################################################################
####### r2 reads handle #############################################################################
    r2_reads <-str_split_fixed(sra_download_site,pattern = ";",n =3 )[,2] #ftp.sra.ebi.ac.uk/vol1/fastq/SRR949/SRR949258/SRR949258_2.fastq.gz
    if (r2_reads!=""){
      message("r2 reads is: ", basename(r2_reads))
      r2_download<-str_replace_all(r2_reads,pattern = "ftp.sra.ebi.ac.uk",replacement = "era-fasp@fasp.sra.ebi.ac.uk:")
      download_command2<-paste(c("ascp -QT -l 50m -P 33001  -i /public/home/wuzefeng/softwares/yes/etc/asperaweb_id_dsa.openssh",r2_download,"0paired/"),collapse  = " ")
	    if(! file.exists(paste("0paired/",basename(r2_reads),sep=""))){
	      message("Paired")
        system(download_command2)}
  
      r2_reads_file <- paste("0paired/",basename(r2_reads),sep="")
      r2_reads_download_md5_command<-paste(c("md5sum",r2_reads_file),collapse = " ")
      r2_reads_download_md5_result <- system(r2_reads_download_md5_command,intern = TRUE)
      r2_reads_download_md5<-str_split_fixed(r2_reads_download_md5_result,pattern = " ",n = 2)[1]

      message("r2_reads true md5 is: ",r2_reads_true_md5)
      message("r2_reads download md5 is: ",r2_reads_download_md5)
  
      if (r2_reads_true_md5!= r2_reads_download_md5){
        message("r2 download failed!! Must REMOVE and redownload!")
        system(paste("rm 0paired/",basename(r2_reads),sep=""))
        system(download_command2)
  }
    }
    else(message("No paired 2 reads!"))
    }
  else{ # single reads
    r_download<-str_replace_all(sra_download_site,pattern = "ftp.sra.ebi.ac.uk",replacement = "era-fasp@fasp.sra.ebi.ac.uk:")
    download_command<-paste(c("ascp -QT -l 50m -P 33001  -i /public/home/wuzefeng/softwares/yes/etc/asperaweb_id_dsa.openssh",r_download,"0single/"),collapse =  " ")
    
    #message("single...",download_command)
	   if(! file.exists(paste("0single/",basename(r_download),sep=""))){
	       message("single")
         system(download_command)}
    
    true_md5<- sra_list_file_md5$fastq_md5[sra_list_file_md5$run_accession==sra_list_file$run_accession[m]]
    
    reads_file <- paste("0single/",basename(r_download),sep="")
    reads_download_md5_command<-paste(c("md5sum",reads_file),collapse = " ")
    reads_download_md5_result <- system(reads_download_md5_command,intern = TRUE)
    reads_download_md5<-str_split_fixed(reads_download_md5_result,pattern = " ",n = 2)[1]
    
    message("reads true md5 is: ",true_md5)
    message("reads download md5 is: ",reads_download_md5)
    
    if(reads_download_md5!=true_md5){
      message("download failed!! Must REMOVE and redownload!")
      system(paste("rm 0single/",basename(r_download),sep=""))
      system(download_command)
    }
    #message(download_command)
  }
}
