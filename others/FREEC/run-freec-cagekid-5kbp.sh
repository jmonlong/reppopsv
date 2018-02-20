WD=`pwd`

sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $1}'`
echo $SAMP
BAMF=`echo $line | awk '{print $2}'`
mkdir -p $SAMP
echo "[general]
chrFiles = /lb/project/mugqic/analyste_dev/genomes/Homo_sapiens/hg19/fasta/byChro
chrLenFile = $WD/../../../installArchives/hg19.len
samtools = /sb/programs/analyste/software/samtools/samtools-0.1.19/samtools
contaminationAdjustment = true
window = 5000
gemMappabilityFile = $WD/../../../installArchives/out100m2_hg19.gem
maxThreads = 6
telocentromeric = 0
intercept = 1
outputDir = $WD/$SAMP/
[sample]
mateFile = $BAMF
inputFormat = BAM
mateOrientation = FR
[control]
" > $SAMP/$SAMP.freec.cfg
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -o $WD/$SAMP/$SAMP.out
#PBS -e $WD/$SAMP/$SAMP.err
#PBS -V
#PBS -N freec$SAMP
$WD/../../../installArchives/FREEC/freec -conf $WD/$SAMP/$SAMP.freec.cfg
" > $SAMP/$SAMP.qsub
qsub $SAMP/$SAMP.qsub
done


## Merge results into one file
echo -e "sample\tchr\tstart\tend\tcn\ttype" > cnvs-FREEC-cagekid-5kbp.tsv
sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $1}'`
echo $SAMP
awk -v samp=$SAMP '{print samp"\t"$0}' $SAMP/*_CNVs >> cnvs-FREEC-cagekid-5kbp.tsv
done




sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $1}'`
echo $SAMP
BAMF=`echo $line | awk '{print $2}'`
mkdir -p $SAMP/loose
echo "[general]
chrFiles = /lb/project/mugqic/analyste_dev/genomes/Homo_sapiens/hg19/fasta/byChro
chrLenFile = $WD/../../../installArchives/hg19.len
samtools = /sb/programs/analyste/software/samtools/samtools-0.1.19/samtools
contaminationAdjustment = true
window = 5000
gemMappabilityFile = $WD/../../../installArchives/out100m2_hg19.gem
maxThreads = 6
telocentromeric = 0
intercept = 1
outputDir = $WD/$SAMP/loose/
breakPointThreshold=0.6
[sample]
mateFile = $BAMF
inputFormat = BAM
mateOrientation = FR
[control]
" > $SAMP/loose/$SAMP.freec.cfg
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -o $WD/$SAMP/$SAMP.out
#PBS -e $WD/$SAMP/$SAMP.err
#PBS -V
#PBS -N freec$SAMP
$WD/../../../installArchives/FREEC/freec -conf $WD/$SAMP/loose/$SAMP.freec.cfg
" > $SAMP/loose/$SAMP.qsub
qsub $SAMP/loose/$SAMP.qsub
done

## Merge results into one file
echo -e "sample\tchr\tstart\tend\tcn\ttype" > cnvs-FREEC-cagekid-5kbp-loose.tsv
sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $1}'`
echo $SAMP
awk -v samp=$SAMP '{print samp"\t"$0}' $SAMP/loose/*_CNVs >> cnvs-FREEC-cagekid-5kbp-loose.tsv
done
