# Other scripts

## Running the methods

Find how PopSV and the other methods were run in their respective folder.
The raw calls of each method were imported in R and formatted by `formatCalls-CNVnator-FREEC-LUMPY-cnMOPS.R`.

### PopSV 

[PopSV](http://jmonlong.github.io/PopSV/) is a R package and the pipeline used [BatchJobs package](https://github.com/tudo-r/BatchJobs) to send jobs to our High Performance Computing resources (PBS/Torque).

In `PopSV` folder, find the scripts used to call CNVs in the Twins study (bin 5 Kbp), CageKid dataset (500 bp and 5 Kbp) and our epilepsy cohort (5 Kbp).

At the time the Twins and Cagekid were ran by manually writing each BatchJobs commands, while the epilepsy run used the functions written in the `automatedPipeline.R` file. More information on BatchJobs and PopSV are available at [http://jmonlong.github.io/PopSV/](http://jmonlong.github.io/PopSV//2-ClusterManagement.md/)

### Other methods

Find the scripts in their respective folder. The versions used for the analysis:

+ [cn.MOPS](http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html) v1.12.0
+ [Control-FREEC](http://boevalab.com/FREEC/) v7.2
+ [CNVnator](https://github.com/abyzovlab/CNVnator) v0.3.2
+ [LUMPY](https://github.com/arq5x/lumpy-sv) v0.2.13. (and [YAHA](https://github.com/GregoryFaust/yaha) v0.1.83)


