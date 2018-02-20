WD=`pwd`
LUMPY=/gs/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/lumpy-sv
YAHA=/gs/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/yaha/bin/yaha
BAMADDRG=/gs/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/bamaddrg/bamaddrg

## Index with YAHA
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -o $WD/yahaIndex.out
#PBS -e $WD/yahaIndex.err
#PBS -V
#PBS -N yahaIndex
module load mugqic/samtools mugqic/bedtools
cd $SCRATCH
ln -s /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa .
$YAHA -g Homo_sapiens.GRCh37.fa
" > $WD/yahaIndex.qsub
# qsub yahaIndex.qsub -A bws-221-ae

## Step 1 : Get discordant/splitters bam files
sed 1d ../../bam-samples.tsv | while read line
do
    SAMPLE=`echo $line | awk '{print $1}'`
    echo $SAMPLE
    BAMFILE=`echo $line | awk '{print $2}'`
    echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=1:00:00:00
#PBS -o $WD/$SAMPLE-1.out
#PBS -e $WD/$SAMPLE-1.err
#PBS -V
#PBS -N $SAMPLE
module load mugqic_dev/samtools/0.1.19 mugqic/bedtools gcc/4.9.1

mkdir $SCRATCH/$SAMPLE

echo \`date\` - Extract unmapped reads
samtools view -b -F 1294 $BAMFILE > $SCRATCH/$SAMPLE/$SAMPLE.discordants.unsorted.bam &
samtools view $BAMFILE | $LUMPY/scripts/split_unmapped_to_fasta.pl -b 20 > $SCRATCH/$SAMPLE/$SAMPLE.unmapped.fastq

wait

echo \`date\` - YAHA alignment
$YAHA -x $SCRATCH/Homo_sapiens.GRCh37.X15_01_65525S -q $SCRATCH/$SAMPLE/$SAMPLE.unmapped.fastq -osh stdout -t 12 | samtools view -Sb - | $BAMADDRG -s $SAMPLE -b - > $SCRATCH/$SAMPLE/$SAMPLE.splitters.unsorted.bam
rm $SCRATCH/$SAMPLE/$SAMPLE.unmapped.fastq

echo DONE
" > $WD/$SAMPLE-1.qsub
    # qsub $SAMPLE-1.qsub -A bws-221-ae
done    


## Step 2 : Sort bam files and run Lumpy
sed 1d ../../bam-samples.tsv | while read line
do
    SAMPLE=`echo $line | awk '{print $1}'`
    echo $SAMPLE
    BAMFILE=`echo $line | awk '{print $2}'`
    echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=1:00:00:00
#PBS -o $WD/$SAMPLE-2.out
#PBS -e $WD/$SAMPLE-2.err
#PBS -V
#PBS -N $SAMPLE
module load mugqic_dev/samtools/1.3 mugqic/bedtools gcc/4.9.1

echo \`date\` - Sort both alignments
samtools sort -@ 6 -o $SCRATCH/$SAMPLE/$SAMPLE.splitters.bam -T $SCRATCH/$SAMPLE/$SAMPLE.splitters.sort $SCRATCH/$SAMPLE/$SAMPLE.splitters.unsorted.bam
samtools sort -@ 6 -o $SCRATCH/$SAMPLE/$SAMPLE.discordants.bam -T $SCRATCH/$SAMPLE/$SAMPLE.discordants.sort $SCRATCH/$SAMPLE/$SAMPLE.discordants.unsorted.bam
rm $SCRATCH/$SAMPLE.discordants.unsorted.bam
rm $SCRATCH/$SAMPLE.splitters.unsorted.bam

echo \`date\` - Lumpy
$LUMPY/bin/lumpyexpress -B $BAMFILE -T $SCRATCH/$SAMPLE/lumpytemp -o $WD/$SAMPLE.vcf \
	-S $SCRATCH/$SAMPLE/$SAMPLE.splitters.bam \
	-D $SCRATCH/$SAMPLE/$SAMPLE.discordants.bam

echo DONE
" > $WD/$SAMPLE-2.qsub
    # qsub $SAMPLE-2.qsub -A bws-221-ae
done    


## Merge results
echo -e "sample\tCHROM\tPOS\tINFO" > svs-Lumpy-twin.tsv
sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $1}'`
echo $SAMP
grep -v "#" $SAMP.vcf | cut -f1,2,8 | awk -v samp=$SAMP '{print samp"\t"$0}' >> svs-Lumpy-twin.tsv
done

