#wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt

mkdir refseq_fasta

## Download refseq fasta for 175 plant species

cat assembly_summary.txt | grep -v ^#|cut -f 8,20|while IFS='\t' read line; do
link=$(echo $line|awk '{print $NF}')
organism=$(echo $line|awk '{$(NF--)=""; print}')
organism=${organism// /_}
organism=${organism::-1}
echo $organism

sampleid=$(basename $link)
wget $link"/"$sampleid"_genomic.fna.gz"
mv $sampleid"_genomic.fna.gz" $organism"_genomic.fna.gz"
done



## Align plant 10k-mers to hg38

fileout="beds10k/all_plant_vs_humangenes.bed"
test -f $fileout || touch $fileout


for i in $(ls refseq_fasta|grep gz);do

# or via gnu parallel 
## ls refseq_fasta|grep gz|parallel -j 6 './run_alignment_10k.sh {}'
## i=$1
cp refseq_fasta/$i .

sample=$i
sample=${sample/_genomic.fna.gz/}
gunzip $i


mv $sample"_genomic.fna"  $sample".fasta" 

## split fasta in 10k-mers
pyfasta split -n 1 -k 10000 -o 200 $sample".fasta"

## keep track of the chrid and kmer seq coordinates
awk '/>/{print $1"_"$NF} !/>/{print $0}' $sample".split.10Kmer.200overlap.fasta" > $sample".tmp.fasta"
mv $sample".tmp.fasta" $sample".split.10Kmer.fasta"


## first generate the human genome index
## minimap2 -d ../v110/hg38v110.mmi Homo_sapiens.GRCh38.dna.primary_assembly.fa

## run alignment and select "reads" mapped wit a minimum quality over 20 
echo "$sample"
minimap2  -ax map-ont ../v110/Homo_sapiens.GRCh38.dna.primary_assembly.fa $sample".split.10Kmer.fasta"  |samtools view -b -q 20 | samtools sort -o bams10k/aligned_$sample"_to_human.10k.sorted.bam" 

## index bams
samtools index bams10k/aligned_$sample"_to_human.10k.sorted.bam"

## convert bams to bed
bamToBed -i bams10k/aligned_$sample"_to_human.10k.sorted.bam" > beds10k/$sample".10k.bed"

## intersect with human gtf to find the coordinates of the plant 10k-mers overlapping human gene features
intersectBed -a beds10k/$sample".10k.bed" -b ../v110/Homo_sapiens.GRCh38.110.genes.bed -wo |sort|uniq|awk -v var=$sample '{print $0"\t"var}'  >> beds10k/all_plant_vs_humangenes.bed

rm "$sample".*
done



##  find the coordinates of all mapped plant 10k-mers overlapping human lncRNAs
cat beds10k/all_plant_vs_humangenes.bed |grep lncRNA|grep "_sequence"|awk -F"[\t, ]" '{print $4,$NF,$12,$7,$8,$9}'|sed 's/"//g'|sed 's/;//g'|sed 's/_sequence_/\t/g'|awk '{OFS="\t";print $1,$2,$2+10000,$3,$4,$5,$6,$7,$8}' |sed 's/\t$//' >beds10k/all_plant_vs_human_lncRNA.bed


# ##  find the coordinates of all mapped plant 10k-mers of a minimum 1000bp overlap with human lncRNAs and are also lncRNAs in plants 
intersectBed -a beds10k/all_plant_vs_human_lncRNA.bed -b all_plant.lncRNA.gtf.bed -wo |awk -F"\t" '{if($NF>1000)print}' >plant_lncrna_human_lncrna_1000overlap.tsv

