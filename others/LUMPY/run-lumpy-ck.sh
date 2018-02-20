WD=`pwd`
SCRATCH=/lb/scratch
LUMPY=/lb/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/lumpy-sv
YAHA=/lb/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/yaha/bin/yaha
BAMADDRG=/lb/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/bamaddrg/bamaddrg

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
ln -s /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa .
$YAHA -g Homo_sapiens.GRCh37.fa
" > $WD/yahaIndex.qsub
# qsub yahaIndex.qsub

## Run pipeline
sed 1d ../../bam-samples.tsv | grep -v D000G2G | sed 1,10d |  while read line
do
    SAMPLE=`echo $line | awk '{print $2}'`
    echo $SAMPLE
    BAMFILE=`echo $line | awk '{print $3}'`
    echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=3:00:00:00
#PBS -o $WD/$SAMPLE.out
#PBS -e $WD/$SAMPLE.err
#PBS -V
#PBS -N $SAMPLE
source /home/jmonlong/.bashrc
module load mugqic/samtools/0.1.19 mugqic/bedtools

echo \`date\` - Extract unmapped reads
samtools view -b -F 1294 $BAMFILE > $SCRATCH/$SAMPLE.discordants.unsorted.bam &
samtools view $BAMFILE | $LUMPY/scripts/split_unmapped_to_fasta.pl -b 20 > $SCRATCH/$SAMPLE.unmapped.fastq

echo \`date\` - YAHA alignment
$YAHA -x $SCRATCH/Homo_sapiens.GRCh37.X15_01_65525S -q $SCRATCH/$SAMPLE.unmapped.fastq -osh stdout -t 12 | samtools view -Sb - | $BAMADDRG -s $SAMPLE -b - > $SCRATCH/$SAMPLE.splitters.unsorted.bam
# rm $SCRATCH/$SAMPLE.unmapped.fastq

wait

echo \`date\` - Sort both alignments
samtools sort $SCRATCH/$SAMPLE.discordants.unsorted.bam $SCRATCH/$SAMPLE.discordants &
samtools sort $SCRATCH/$SAMPLE.splitters.unsorted.bam $SCRATCH/$SAMPLE.splitters &
wait
# rm $SCRATCH/$SAMPLE.discordants.unsorted.bam
# rm $SCRATCH/$SAMPLE.splitters.unsorted.bam

echo \`date\` - Lumpy
$LUMPY/bin/lumpyexpress -B $BAMFILE -T $SCRATCH/$SAMPLE -o $WD/$SAMPLE.vcf \
        -S $SCRATCH/$SAMPLE.splitters.bam \
        -D $SCRATCH/$SAMPLE.discordants.bam

echo DONE
" > $WD/$SAMPLE.qsub
    qsub $SAMPLE.qsub
done

## Merge results
echo -e "sample\tCHROM\tPOS\tINFO" > svs-Lumpy-cagekid.tsv
sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $2}'`
echo $SAMP
grep -v "#" $SAMP.vcf | cut -f1,2,8 | awk -v samp=$SAMP '{print samp"\t"$0}' >> svs-Lumpy-cagekid.tsv
done
