WD=`pwd`
CNVNATOR=/gs/project/mugqic/projects/jmonlong/PopSV-forPaper/installArchives/cnvnator/CNVnator_v0.3.2/src/cnvnator

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
samtools faidx /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \$chr > $SCRATCH/\$chr.fa
done
" > $WD/splitFasta.qsub
#qsub splitFasta.qsub -A bws-221-ae

# Run CNVnator
sed 1d ../../bam-samples.tsv | while read line
do
    SAMPLE=`echo $line | awk '{print $1}'`
    echo $SAMPLE
    BAMFILE=`echo $line | awk '{print $2}'`
    echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=36:00:00
#PBS -o $WD/$SAMPLE.out
#PBS -e $WD/$SAMPLE.err
#PBS -V
#PBS -N ctor$SAMPLE

# Setup ROOT
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source \$ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
lsetup root

# CNVnator pipeline
for chr in \`seq 1 22\`
do
$CNVNATOR -root $SCRATCH/$SAMPLE.root -chrom \$chr  -tree $BAMFILE -unique
done

# 5kbp
$CNVNATOR -root $SCRATCH/$SAMPLE.root -his 5000 -d $SCRATCH
$CNVNATOR -root $SCRATCH/$SAMPLE.root -stat 5000 
$CNVNATOR -root $SCRATCH/$SAMPLE.root -partition 5000 
$CNVNATOR -root $SCRATCH/$SAMPLE.root -call 5000 > $WD/$SAMPLE-5kbp.tsv

" > $WD/$SAMPLE.qsub
    # qsub $SAMPLE.qsub -A bws-221-ae
done

## Merge 5kbp results
echo -e "sample\tCNV_type\tcoordinates\tCNV_size\tnormalized_RD\teval1\teval2\teval3\teval4\tq0" > cnvs-CNVnator-twin-5kbp.tsv
sed 1d ../../bam-samples.tsv | while read line
do
SAMP=`echo $line | awk '{print $1}'`
echo $SAMP
awk -v samp=$SAMP '{print samp"\t"$0}' $SAMP-5kbp.tsv >> cnvs-CNVnator-twin-5kbp.tsv
done
