WD=`pwd`

## 5Kbp
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l walltime=12:00:00
#PBS -o $WD/cnmops.out
#PBS -e $WD/cnmops.err
#PBS -V
#PBS -N cnMOPs
Rscript $WD/cnmops-rscript.R $WD/../../5kbp/bc-gcCor-all.tsv.bgz 0 $WD/res-cnmops-twin-5kbp.RData
" > cnmops.qsub
qsub cnmops.qsub -A bws-221-ae


## Loose threshold
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l walltime=12:00:00
#PBS -o $WD/cnmops-loose.out
#PBS -e $WD/cnmops-loose.err
#PBS -V
#PBS -N cnMOPs
Rscript $WD/cnmops-rscript.R $WD/../../5kbp/bc-gcCor-all.tsv.bgz 0.5 $WD/res-cnmops-twin-5kbp-loose.RData
" > cnmops-loose.qsub
qsub cnmops-loose.qsub -A bws-221-ae
