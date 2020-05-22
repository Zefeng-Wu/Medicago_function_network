#!/bin/bash

## 1 download ensemble plants data
mkdir pep
cd pep
for m in $(ls *.fa); do echo $m; cat $m| cut -d" " -f1,4 | sed 's/gene://' | sed '/^>/ s/$/ /' | tr -s "\s" |  tr -d "\n" | sed 's/>/\n>/g' | sed '1d' | awk -v OFS="\t" '{print $1, $2, $3, length($3)}' | sort -k2,2 -k4nr | awk 'L!=$2 {print  ">"$2,"\n" $3} {L=$2}' > ../1longest_pep/$(basename $m); done

makeblastdb  -dbtype prot -in Medicago_truncatula.MedtrA17_4.0.pep.all.fa
nohup sh -c 'for m in $(ls 1longest_pep/*.fa);do echo $m; blastp -db 2Mt2others/Medicago_truncatula.MedtrA17_4.0.pep.all.fa -out 2Mt2others/$(basename $m).blout.txt -evalue 1e-5 -query $m  -outfmt 6 -num_threads 20 ;done' &

mkdir cds
mkdir gtf

rsync -av rsync://ftp.ensemblgenomes.org/all/pub/plants/current/fasta/*/cds/*.cds.all.fa.gz cds/
## 
## 2 downlaod ensembl fungi pep
rsync -av rsync://ftp.ensemblgenomes.org/all/pub/fungi/current/fasta/*/pep/*.pep.all.fa.gz ./
rsync -av rsync://ftp.ensemblgenomes.org/all/pub/fungi/current/fasta/*/*/pep/*.pep.all.fa.gz ./

## 3 downlaod ensembl bacterial pep
    
    for m in $(seq 1 183); do echo $m; rsync -av rsync://ftp.ensemblgenomes.org/all/pub/bacteria/current/fasta/bacteria_${m}_collection/*/pep/*.pep.all.fa.gz ./ ;done
    python python species_uniq.py
    
    nohup sh -c 'for m in $(ls *.fa);do echo $m; blastp -db ../ensembl_plants/2Mt2others/Medicago_truncatula.MedtrA17_4.0.pep.all.fa -out 2Mt2others/$(basename $m).blout.txt -evalue 1e-5 -query $m  -outfmt 6 -num_threads 25 ;done' &

## 4 metazoa (no alternative splicing)
rsync -av rsync://ftp.ensemblgenomes.org/all/pub/metazoa/current/fasta/*/*/pep/*.pep.all.fa.gz ./
mkdir 1longest_pep
mkdir 2Mt_blout
for m in $(ls *.fa); do echo $m; cat $m| cut -d" " -f1,4 | sed 's/gene://' | sed '/^>/ s/$/ /' | tr -s "\s" |  tr -d "\n" | sed 's/>/\n>/g' | sed '1d' | awk -v OFS="\t" '{print $1, $2, $3, length($3)}' | sort -k2,2 -k4nr | awk 'L!=$2 {print  ">"$2,"\n" $3} {L=$2}' >1longest_pep/$(basename $m); done
nohup sh -c 'for m in $(ls 1longest_pep/*.fa);do echo $m; blastp -db ../ensembl_plants/2Mt2others/Medicago_truncatula.MedtrA17_4.0.pep.all.fa -out 2Mt_blout/$(basename $m).blout.txt -evalue 1e-5 -query $m  -outfmt 6 -num_threads 20 ;done' &

## 5 genome
rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_fasta/*/pep/*.pep.all.fa.gz ./
nohup sh -c 'for m in $(ls 1longest_pep/*.fa);do echo $m; blastp -db ../ensembl_plants/2Mt2others/Medicago_truncatula.MedtrA17_4.0.pep.all.fa -out 2Mt_blout/$(basename $m).blout.txt -evalue 1e-5 -query $m  -outfmt 6 -num_threads 20 ;done' &

## domain data
#用hmmscan软件扫描所有蛋白序列的结构域
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz # 下载pfam结构域数据
hmmpress Pfam-A.hmm
nohup  hmmscan --tblout 1Mt.hmmer.out.txt --noali -E 1e-2 --domE 1e-2 --cpu 28 Pfam-A.hmm Medicago_truncatula.MedtrA17_4.0.pep.all.fa &
grep -v "#" 1Mt.hmmer.out.txt | awk 'OFS="\t",u=split($2,m,"."){$2=m[1];print $0}' > 2Mt.hmmer.out.modify

wget http://3diana.cnb.csic.es/DIMERO/chimera_plugin/DIANA.tar.gz #下载domain结构与互作数据
grep "PF" 3did_flat | cut -f4,5 | sed 's/\@Pfam//g' | sed 's/(//g' | sed 's/)//g'| sed 's/\.[0-9]*//g' > 3pfam_domain_intaction



