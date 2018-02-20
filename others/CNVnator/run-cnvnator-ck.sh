WD=`pwd`
CNVNATOR=/lb/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/CNVnator_v0.3.2/src/cnvnator
SCRATCH=/lb/scratch/
ROOT=/lb/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/root-6.04.14

# Split chr fasta
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -o $WD/splitFasta.out
#PBS -e $WD/splitFasta.err
#PBS -V
#PBS -N splitFasta
module load mugqic/samtools 
for chr in \`seq 1 22\`
do
samtools faidx /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \$chr > $SCRATCH/\$chr.fa
done
" > $WD/splitFasta.qsub
# qsub splitFasta.qsub


# Run CNVnator
sed 1d ../../bam-samples.tsv | grep -v -f samples.done | grep -f samples.counted | while read line
do
    SAMPLE=`echo $line | awk '{print $2}'`
    echo $SAMPLE
    BAMFILE=`echo $line | awk '{print $3}'`
    rm $WD/$SAMPLE.out $WD/$SAMPLE.err
    echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=20:00:00
#PBS -o $WD/$SAMPLE.out
#PBS -e $WD/$SAMPLE.err
#PBS -V
#PBS -N ctor$SAMPLE

# Setup ROOT etc
module load itgenome/gcc/4.8.5
source $ROOT/bin/thisroot.sh

# rm $SCRATCH/$SAMPLE.root

# CNVnator pipeline
# for chr in \`seq 1 22\`
# do
# $CNVNATOR -root $SCRATCH/$SAMPLE.root -chrom \$chr  -tree $BAMFILE -unique
# done

# 5kbp
$CNVNATOR -root $SCRATCH/$SAMPLE.root -his 5000 -d $SCRATCH
$CNVNATOR -root $SCRATCH/$SAMPLE.root -stat 5000
$CNVNATOR -root $SCRATCH/$SAMPLE.root -partition 5000
$CNVNATOR -root $SCRATCH/$SAMPLE.root -call 5000 > $WD/$SAMPLE-5kbp.tsv
" > $WD/$SAMPLE.qsub
    qsub $SAMPLE.qsub
done



## Merge 5kbp results
echo -e "sample\tCNV_type\tcoordinates\tCNV_size\tnormalized_RD\teval1\teval2\teval3\teval4\tq0" > cnvs-CNVnator-cagekid-5kbp.tsv
sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $2}'`
echo $SAMP
awk -v samp=$SAMP '{print samp"\t"$0}' $SAMP-5kbp.tsv >> cnvs-CNVnator-cagekid-5kbp.tsv
done
