
# work only with lnRNA with min 1000 bp overlap in human lncRNA



###### extract human lncrna fa sequences and make blast database for each of them
for i in $(cat $input |cut -f 5|sort|uniq);do
echo "processing human "$i
hs=$i."fa"
test -f $hs || touch $hs

linehs=$(cat $input|grep $i|head -1 )

hschr=$(echo $linehs|cut -d " " -f 6)
hsstart=$(echo $linehs|cut -d " " -f 7)
hsend=$(echo $linehs|cut -d " " -f 8)
 
samtools faidx ../../v110/Homo_sapiens.GRCh38.dna.primary_assembly.fa $hschr":"$hsstart"-"$hsend > $i.tmp.txt

awk '/>/{print $0"_human"} !/>/{print $0}' $i.tmp.txt > $hs
rm  $i.tmp.txt

#make blastdb
makeblastdb -in $hs -dbtype nucl -out $i
done



###### extract org lncrna seq
input="plant_lncrna_human_lncrna_1000overlap.tsv"

while IFS= read -r line
do 
## ENSG id
hs=$(echo $line|cut -d" " -f 5)
#organism
org=$(echo $line|cut -d" " -f 4)
echo $org
seqfile=$org."fa"

## copy fasta, unpack and modify fasta headers
fastafile=$org"_genomic.fna.gz"

 if [ ! -f "$fastafile" ];then
cp "../refseq_fasta/$fastafile" .
gunzip $org"_genomic.fna.gz"
wait

## append organism name to the header
awk '/>/{print $1"_"$NF} !/>/{print $0}' $org"_genomic.fna" >$org"_genomic.tmp.fna"
wait

## index fasta file
samtools faidx $org"_genomic.tmp.fna"
fi


## coordinates of plant lncRNAs
chr=$(echo $line|cut -d" " -f 9)
start=$(echo $line|cut -d" " -f 10)
end=$(echo $line|cut -d" " -f 11)

echo -e "processing $org\t$hs\t$chr\t$start\t$end"

### extract lncRNA seq from corresponding fasta and append the organism name to fasta header and convert nucleotides to uppercase
samtools faidx  $org"_genomic.tmp.fna" $chr"_sequence:"$start"-"$end > $org"_"$hs"_"$chr"_"$start"-"$end".tmp.fa"
cat $org"_"$hs"_"$chr"_"$start"-"$end".tmp.fa" | awk '/>/{print $0} !/>/{print toupper($0)}' |awk -v var=$org '/>/{print $0"_"var} !/>/{print $0}' |sed 's/_sequence//g' >seqs/$org"_"$hs"_"$chr"_"$start"-"$end".fa" 

## run blastn (megablast) query: plant lncRNAseq vs human lncRNA
blastn -query seqs/$org"_"$hs"_"$chr"_"$start"-"$end".fa"  -db $hs -task megablast -out outblast_10k/$org"_"$hs"_"$chr"_"$start"-"$end".tsv" -outfmt "6 qseqid sseqid pident nident slen length evalue bitscore gapopen gaps frames" 

rm $org"_"$hs"_"$chr"_"$start"-"$end".tmp.fa"
done < "$input"
