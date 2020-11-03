wget ftp://ftp.kobic.kr/Data/refseq/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_assembly_report.txt
# mart_export.txt.gz select Dataset; select "Ensembl Genes 96"; select "Human genes (GRCh37.p13); select Atrributes; select Structures; select "Transcript stable ID" "Chromosome/scaffold name" "Transcript length (including UTRs and CDS)" "Gene name" "Exon region start (bp)" "Exon region end (bp)" "Exon stable ID"

grep ENST00000374690 ensemble_GRCh37 |awk '{OFS=t; print chr,,,AR-001,}' > AR.bed
awk '{OFS=t; for(i=int();i<=int();i++){print ,i,,}}' AR.bed > AR.loci
grep ENST00000254958 ensemble_GRCh37 |awk '{OFS=t; print chr,,,JAG1-201,}' > JAG1.bed
awk '{OFS=t; for(i=int();i<=int();i++){print ,i,,}}' JAG1.bed > JAG1.loci
grep ENST00000405650 ensemble_GRCh37 |awk '{OFS=t; print chr,,,SRD5A2-001,}' > SRD5A2.bed 
awk '{OFS=t; for(i=int();i<=int();i++){print ,i,,}}' SRD5A2.bed > SRD5A2.loci
grep ENST00000423058 ensemble_GRCh37 |awk '{OFS=t; print chr,,,SCN1A-005,}' > SCN1A.bed 
awk '{OFS=t; for(i=int();i<=int();i++){print ,i,,}}' SCN1A.bed > SCN1A.loci
grep ENST00000416240 ensemble_GRCh37 |awk '{OFS=t; print chr,,,SLC25A13-201,}' > SLC25A13.bed
awk '{OFS=t; for(i=int();i<=int();i++){print ,i,,}}' SLC25A13.bed > SLC25A13.loci
grep ENST00000371100 ensemble_GRCh37 |awk '{OFS=t; print chr,,,GNAS-214,}' > GNAS.bed 
grep ENST00000265715 ensemble_GRCh37 |awk '{OFS=t; print chr,,,SLC26A4-201,}' > SLC26A4.bed
grep ENST00000373366 ensemble_GRCh37 |awk '{OFS=t; print chr,,,GJB3-202,}' > GJB3.bed
grep ENST00000357033 ensemble_GRCh37 |awk '{OFS=t; print chr,,,DMD-203,}' > DMD.bed
grep ENST00000242839 ensemble_GRCh37 |awk '{OFS=t; print chr,,,ATP7B-201,}' > ATP7B.bed
grep ENST00000242592 ensemble_GRCh37 |awk '{OFS=t; print chr,,,ACADS-201,}' > ACADS.bed
grep ENST00000358776 ensemble_GRCh37 |awk '{OFS=t; print chr,,,ACADSB-201,}' > ACADSB.bed
grep ENST00000393562 ensemble_GRCh37 |awk '{OFS="\t"; print "chr"$4,$2,$3,"G6PD-001",$6}' > G6PD.bed
