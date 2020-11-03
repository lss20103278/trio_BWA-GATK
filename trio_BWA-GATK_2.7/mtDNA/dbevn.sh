
soft_path="/SSD750/PB1/soft1"
DB="/SSD750/PB1/db1/Homo"
ref_genome="$DB/refseq/hg19.fa"
bwa_index="$DB/bwa_index/hg19.fa"
ref_1000G_indel="$DB/GATK/1000G_phase1.indels.hg19.sites.vcf"
ref_1000M_indel="$DB/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
#ref_annotation="/SSD750/PB1/db1/Homo/annovar"
ref_annotation="/home/ana005/data/annovar/new/humandb"
ref_snp="$DB/GATK/dbsnp_138.hg19.vcf"
ref_vqsr_hapmap="/SSD750/PB3/db3/Homo/GATK/hapmap_3.3.hg19.sites.vcf"
ref_vqsr_omni="/SSD750/PB3/db3/Homo/GATK/1000G_omni2.5.hg19.sites.vcf"
ref_vqsr_1000G="/SSD750/PB3/db3/Homo/GATK/1000G_phase1.snps.high_confidence.hg19.sites.vcf"
adapter="$DB/refseq/Illumina_Adapter_Contaminants.fasta"

panel="$DB/refseq/hg19_SeqCap_EZ_Exome_v3_capture_plus20bp.bed"               #WXS
#panel="$DB/design/SCH.capture.hg19.20bp.bed"                                  #PANEL
#panel="$DB/refseq/Agilentinherited.bed"                                   #Agilent
#panel="/home/ana005/data/annovar/bed/IDT/IDT.bed"            #IDT
#panel="/home/ana005/data/annovar/bed/Agilent/S06588914/S06588914_Regions.bed"    # Agilent_wes
