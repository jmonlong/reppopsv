WD=`pwd`

## 5kbp
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -o $WD/cnmops.out
#PBS -e $WD/cnmops.err
#PBS -V
#PBS -N cnMOPs
Rscript $WD/cnmops-rscript-ck.R $WD/../../5kbp/files.RData 0 $WD/res-cnmops-cagekid-5kbp.RData
" > cnmops.qsub
qsub cnmops.qsub -A bws-221-ae


## Loose threshold 
echo "#"'!'"/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l walltime=12:00:00
#PBS -o $WD/cnmops-loose.out
#PBS -e $WD/cnmops-loose.err
#PBS -V
#PBS -N cnMOPs
Rscript $WD/cnmops-rscript-ck.R $WD/../../5kbp/files.RData 0.5 $WD/res-cnmops-cagekid-5kbp-loose.RData
" > cnmops-loose.qsub
qsub cnmops-loose.qsub -A bws-221-ae
