#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=1573077420@qq.com

#########################################################################
# File Name: T084V2_bed_add_genesymbol.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 11 Apr 2019 02:49:25 PM CST
# Usage:
#########################################################################

#### Prepare the required files
# ucscgene ucscgene_to_genesymbol: generated from ucsc website
#rm !(ex*)
#mart_export.txt.gz # http://grch37.ensembl.org/biomart/martview/5c1b7aa415ae11ee9c6aff7712d7b7a7 choose uniq
#[ ! -e mart_export.txt ] && gunzip mart_export.txt.gz # 63676 rows, 56639 unique gene names, there are many strange chromosomes
[ ! -e Homo_sapiens.gene_info.gz ] &&  wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
[ ! -e Homo_sapiens.gene_info ] && gunzip Homo_sapiens.gene_info.gz
#wget https://www.omim.org/static/omim/data/mim2gene.txt # omim database
#wget ftp://ftp.ensembl.org/pub/grch37/release-96/mysql/homo_sapiens_cdna_96_37/gene.txt.gz 
#wget https://decipher.sanger.ac.uk/files/downloads/HI_Predictions_Version3.bed.gz # Decipher database
#wget https://decipher.sanger.ac.uk/files/downloads/population_cnv.txt.gz # Decipher database
#wget http://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz # Decipher database

[ ! -e T084V2.ig.bed ] && awk '{OFS="\t"; print $0,"."}' ../doc/T084V2.ig.bed > T084V2.ig.bed # bedmap adds fiels in mapping files to the end of reference file and uses '|' to seperate them, so add a column of '.' to T084V2.ig.bed to avoid directly append '|' and other information at the end column of T084V2.ig.bed
[ ! -e T084V2.ig.sort.bed ] && sort-bed T084V2.ig.bed > T084V2.ig.sort.bed # 72904 rows
[ ! -e T084V2_CNV.igt.bed ] && awk '{OFS="\t"; print $0,"."}' ../doc/T084V2_CNV.igt.bed > T084V2_CNV.igt.bed
[ ! -e T084V2_CNV.igt.sort.bed ] && sort-bed T084V2_CNV.igt.bed > T084V2_CNV.igt.sort.bed

[ ! -e Homo_sapiens_protein-coding.gene ] && grep protein-coding Homo_sapiens.gene_info |cut -f 3 > Homo_sapiens_protein-coding.gene # 20116 genes

#### Use ucsc
cp ../doc/{ucscgene,ucscgene_to_genesymbol} .
join -1 4 -2 1 <(sort -k4,4 ucscgene) <(sort -k1,1 ucscgene_to_genesymbol) |cut -d' ' -f 2-4,13 |tr " " "\t" > ucscgene_genesymbol.bed
sort-bed ucscgene_genesymbol.bed > ucscgene_genesymbol.sort.bed
python ../src/extract_protein_coding.py
sort-bed ucscgene_genesymbol_protein_coding.bed > ucscgene_genesymbol_protein_coding.sort.bed # 18064 unique genes

#### Use ensembl
#join -1 1 -2 1 <(sort -k1,1 mart_export.txt) <(sort -k1,1 Homo_sapiens_protein-coding.gene) |tr " " "\t" |awk '{OFS="\t"; print "chr"$2,$3,$4,$1}' > ensembl_protein-coding.bed
#sort-bed ensembl_protein-coding.bed > ensembl_protein-coding.sort.bed

#### Using bedmap to mapping these two files to add genesymbol to the bed
bedmap --echo --echo-map-id-uniq T084V2.ig.sort.bed ucscgene_genesymbol_protein_coding.sort.bed > T084V2.ig.cnv.bed # T084V2.ig.sort.bed:reference file, ucscgene_genesymbol_protein_coding.sort.bed:mapping file, --echo:T084V2.ig.sort.bed all rows, --echo-map-id-uniq:some genes in the mapping file have several intervals that overlap with each other
bedmap --echo --echo-map-id-uniq T084V2_CNV.igt.sort.bed ucscgene_genesymbol_protein_coding.sort.bed > T084V2_CNV.igt.cnv.bed
# bedmap --echo --echo-map ucscgene_genesymbol_protein_coding.sort.bed T084V2.ig.sort.bed

#### Tongji each gene of Homo_sapiens_protein-coding.gene has how many rows in other.refseq.ucsc.refseq.hg19
#grep SCNN1D ucscgene_genesymbol_protein_coding.sort.bed
#for i in `cat Homo_sapiens_protein-coding.gene`; do grep -cw $i other.refseq.ucsc.refseq.hg19 >> gene_count; done


