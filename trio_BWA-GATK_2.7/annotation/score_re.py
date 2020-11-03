import re
import os

p1=re.compile("Pathogenic")
p2=re.compile("Likely pathogenic")
#p3=re.compile("risk factor")
#p4=re.compile("association")
p3=re.compile("Likely benign")
p4=re.compile("Benign")
p4_2 = re.compile("Uncertain significance")
p5=re.compile("CLIN_pathogenic")
p6=re.compile("CLIN_likely_pathogenic")
p7=re.compile("CLIN_risk_factor")
p8=re.compile("CLIN_association")
s1=[4,3,2,1]
if os.path.exists('list_all'):
	infile=open("list_all","r")
else:
	infile=open("list", "r")
#inhgmd=open("hgmd.vcf","r")
#hgmd={}
#for line in inhgmd.readlines():
#  hgmd["chr"+line.split('\t')[0]+','+line.split('\t')[1]]=1
diseaseid={}
disease={}
typedisease={}
inheri={}
indisease=open("/DATA/sslyu/trio_BWA-GATK_2.7/annotation/db/gene_disease.txt","r")
for line in indisease.readlines():
  if line.split('\t')[3] in disease.keys():
    disease[line.split('\t')[3]]=disease[line.split('\t')[3]]+','+line.split('\t')[1]
  else:
    disease[line.split('\t')[3]]=line.split('\t')[1]
  if line.split('\t')[3] in diseaseid.keys():
    diseaseid[line.split('\t')[3]]=diseaseid[line.split('\t')[3]]+','+line.split('\t')[0]
  else:
    diseaseid[line.split('\t')[3]]=line.split('\t')[0]
  if line.split('\t')[3] in typedisease.keys():
    typedisease[line.split('\t')[3]]=typedisease[line.split('\t')[3]]+','+line.split('\t')[2]
  else:
    typedisease[line.split('\t')[3]]=line.split('\t')[2]
  if line.split('\t')[3] in inheri.keys():
    inheri[line.split('\t')[3]]=inheri[line.split('\t')[3]]+','+line.strip().split('\t')[4]
  else:
    inheri[line.split('\t')[3]]=line.strip().split('\t')[4]
for file in infile.readlines():
  print file.strip()
  inline=open(file.strip()+".ann.hg19_multianno.txt","r")
  outfile=open(file.strip()+".score","w")
  outfile.write(open(file.strip()+".ann.hg19_multianno.txt","r").readline().strip()+'\t'+"score_clinvar"+'\t'+"score_ensembl"+'\t'+"score_pre"+'\t'+"score_HGMD"+'\t'+"score_total"+'\t'+"DiseaseID"+'\t'+"Disease"+'\t'+"Disease type"+'\t'+"InheritanceStatus"+'\n')
  for line in inline.readlines()[1:]:
    s2=[0,0,0,0]
    s3=[0,0,0,0]
    scorepre=0
    scorehgmd=0
    if p1.search(line.strip().split('\t')[119]): #
      s2[0]=1
    if p2.search(line.strip().split('\t')[119]): #
      s2[1]=1
    if p3.search(line.strip().split('\t')[119]): #
      s2[2]=1
    if p4.search(line.strip().split('\t')[119]): #
      s2[3]=1
    if p4_2.search(line.strip().split('\t')[119]): #
      s2[3]=1
    scoreclin=max(map(lambda x, y: x*y,  s1, s2))
    if p5.search(line.strip().split('\t')[124]): #
      s3[0]=1
    if p6.search(line.strip().split('\t')[124]): #
      s3[1]=1
    if p7.search(line.strip().split('\t')[124]):
      s3[2]=1
    if p8.search(line.strip().split('\t')[124]):
      s3[3]=1
    scoreensembl=max(map(lambda x, y: x*y,  s1, s3))
    if line.strip().split('\t')[25]=="D": #
      scorepre=scorepre+1
    if line.strip().split('\t')[31]=="D": #
      scorepre=scorepre+1
    if line.strip().split('\t')[37]=="D": #
      scorepre=scorepre+1
    if line.strip().split('\t')[43]=="D": #
      scorepre=scorepre+1
    if line.strip().split('\t')[57]=="D":
    	scorepre=scorepre+1
    if line.strip().split('\t')[89] != 'NA':
      if float(line.strip().split('\t')[89])>0.5:
        scorepre=scorepre+1
#    if line.split('\t')[0]+','+line.split('\t')[1] in hgmd.keys():
    confhgmd=['NA','','.','-']
    if (line.split('\t')[125] not in confhgmd): #
      scorehgmd=5
    if line.split('\t')[6].split(',')[0] in diseaseid.keys():
      outfile.write(line.strip()+'\t'+str(scoreclin)+'\t'+str(scoreensembl)+'\t'+str(scorepre)+'\t'+str(scorehgmd)+'\t'+str(scoreclin+scoreensembl+scorepre+scorehgmd)+'\t'+diseaseid[line.split('\t')[6].split(',')[0]]+'\t'+disease[line.split('\t')[6].split(',')[0]]+'\t'+typedisease[line.split('\t')[6].split(',')[0]]+'\t'+inheri[line.split('\t')[6].split(',')[0]]+'\n')
    else:
      outfile.write(line.strip()+'\t'+str(scoreclin)+'\t'+str(scoreensembl)+'\t'+str(scorepre)+'\t'+str(scorehgmd)+'\t'+str(scoreclin+scoreensembl+scorepre+scorehgmd)+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\n')
      
