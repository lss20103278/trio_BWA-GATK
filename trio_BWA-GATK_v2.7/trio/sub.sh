#!/bin/sh

if [ ! -e ped/peddy.ped ]
then
echo "ERROR:peddy.ped must exist,please provide.Exiting....."
exit
fi

gname=`awk '{print $2}' ped/peddy.ped`
for i in $gname
do
[ ! -e gvcf/$i.raw.g.vcf ] && mv gvcf/$i.g.vcf gvcf/$i.raw.g.vcf
./gvcf.sh $i
done

for i in `cat list`
do
[ ! -e $i ] && mkdir -p $i
cd $i 
cp ../*.sh ./
sed -i "s/^sample.*$/sample=$i/;s/^#SBATCH -J .*$/#SBATCH -J $i/" cmd.sh
sbatch cmd.sh
cd ..
done
