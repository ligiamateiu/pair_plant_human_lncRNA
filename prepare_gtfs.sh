# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt

# download refseq gtf

cat assembly_summary.txt | grep -v ^#|cut -f 8,20|while IFS='\t' read line; do
link=$(echo $line|awk '{print $NF}')
# echo $link
organism=$(echo $line|awk '{$(NF--)=""; print}')
organism=${organism// /_}
#remove last element "_"
organism=${organism::-1}
echo $organism

sampleid=$(basename $link)
wget $link"/"$sampleid"_genomic.gtf.gz"
mv $sampleid"_genomic.gtf.gz" $organism"_genomic.gtf.gz"
done


# extract the lncRNAs coordinates from all plant gtfs in a single bed file

fileout="all_plant.lncRNA.gtf.bed"
test -f $fileout || touch $fileout

for i in $(ls gtfs|grep -e ".gz$");do
sample=$i
sample=${sample/_genomic.gtf.gz/}
echo $sample
zcat gtfs/$i|grep -i lncRNA|awk -F "\t" -v var=$sample '{print $0"\t"var}' | cut -f 1,4,5,7,9- >> $fileout
done
