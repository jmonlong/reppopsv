## devtools::install_github("jmonlong/PopSV", ref="forPaper", args="-l /home/jmonlong/R/PopSVforPaper") ## Install
library(BatchJobs)
library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")

## Bam files and bin definition
files.df = read.table("../bam-samples.tsv",as.is=TRUE, header=TRUE)
bin.size = 5e3
##
files.df = init.filenames(files.df, code="5kbp")
save(files.df, file="files.RData")
bins.df = fragment.genome.hg19(bin.size)
save(bins.df, file="bins.RData")
rm(bins.df)


#### Get GC content in each bin
## system("rm -rf getGC-files")
getGC.reg <- makeRegistry(id="getGC")
getGC.f <- function(imF){
  load(imF)
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  bins.df = getGC.hg19(bins.df)
  save(bins.df, file=imF)
}
batchMap(getGC.reg, getGC.f,"bins.RData")
submitJobs(getGC.reg, 1, resources=list(walltime="1:0:0", nodes="1", cores="1",supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(getGC.reg)


#### Get bin counts in each sample and correct for GC bias
## system("rm -rf getBC-files")
getBC.reg <- makeRegistry(id="getBC")
getBC.f <- function(file.i, bins.f, files.df){
  load(bins.f)
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  bb.o = bin.bam(files.df$bam[file.i], bins.df, files.df$bc[file.i])
  correct.GC(files.df$bc.gz[file.i], bins.df, files.df$bc.gc[file.i])
  bb.o
}
batchMap(getBC.reg, getBC.f,1:nrow(files.df), more.args=list(bins.f="bins.RData", files.df=files.df))
submitJobs(getBC.reg, findNotDone(getBC.reg), resources=list(walltime="10:0:0", nodes="1", cores="1",supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(getBC.reg)

## OPTIONAL QC: check the total number of reads counted in all samples
library(ggplot2)
pdf("qc-nbReads.pdf")
qplot(reduceResultsVector(getBC.reg, fun=function(job, res)res$nb.reads)) + xlab("total number of reads") + ylab("number of samples") + theme_bw()
dev.off()
##



#### Sample QC and reference sample definition
### All samples are used as reference.
## system("rm -rf sampQC-files")
sampQC.reg <- makeRegistry(id="sampQC")
sampQC.f <- function(bins.f, files.df){
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  load(bins.f)
  pdf("sampQC.pdf")
  qc.o= qc.samples(files.df, bins.df, outfile.prefix="bc-gcCor-all.tsv", nb.cores=3)
  dev.off()
  qc.o
}
batchMap(sampQC.reg, sampQC.f, "bins.RData", more.args=list(files.df=files.df))
submitJobs(sampQC.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="3",queue="sw",supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(sampQC.reg)
samp.qc.o = loadResult(sampQC.reg, 1)


#### Normalize and compute Z scores
## system("rm -rf bcNormTN-files")
bcNormTN.reg <- makeRegistry(id="bcNormTN")
load("bins.RData")
bins.df = chunk.bin(bins.df, bg.chunk.size=1e5, sm.chunk.size=1e4, large.chr.chunks=TRUE)
save(bins.df, file="bins.RData")
bcNormTN.f <- function(chunk.id, file.bc, file.bin, cont.sample){
  load(file.bin)
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  bc.df = read.bedix(file.bc, subset(bins.df, bg.chunk==subset(bins.df, sm.chunk==chunk.id)$bg.chunk[1]))
  tn.norm(bc.df, cont.sample, bins=subset(bins.df, sm.chunk==chunk.id)$bin)
}
batchMap(bcNormTN.reg, bcNormTN.f,unique(bins.df$sm.chunk), more.args=list(file.bc=samp.qc.o$bc, file.bin="bins.RData",cont.sample=samp.qc.o$cont.sample))
submitJobs(bcNormTN.reg, findNotDone(bcNormTN.reg) , resources=list(walltime="12:0:0", nodes="1", cores="1",supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(bcNormTN.reg)

#### Write normalized bin counts and reference metrics
out.files = paste("ref-twin-kbp", c("bc-norm.tsv", "msd.tsv"), sep="-")
file.remove(out.files)
tmp = reduceResultsList(bcNormTN.reg, fun=function(res, job){
  write.table(res$bc.norm, file=out.files[1], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[1]), col.names=!file.exists(out.files[1]))
  write.table(res$norm.stats[,1:7], file=out.files[2], sep="\t", row.names=FALSE, quote=FALSE, append=file.exists(out.files[2]), col.names=!file.exists(out.files[2]))
})

#### Compute Z-scores in reference samples
## system("rm -rf zRef-files")
zRef.reg <- makeRegistry(id="zRef")
zRef.f <- function(bc.f, files.df){
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  z.comp(bc.f=bc.f, files.df=files.df, nb.cores=3, z.poisson=TRUE, chunk.size=1e4)
}
batchMap(zRef.reg, zRef.f,out.files[1], more.args=list(files.df=files.df))
submitJobs(zRef.reg, 1, resources=list(walltime="6:0:0", nodes="1", cores="3",queue="sw", supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(zRef.reg)

#### Abnormal bin calling
## system("rm -rf abCovCallCases-files")
abCovCallCases.reg <- makeRegistry(id="abCovCallCases")
abCovCallCases.f <- function(samp, files.df, norm.stats.f, bins.f,  bin.size){
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  load(bins.f)
  call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"/",samp,"-sdest-abCovCall.pdf"), FDR.th=.001, merge.cons.bins="stitch", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=bin.size*2+1, gc.df = bins.df)
}
batchMap(abCovCallCases.reg, abCovCallCases.f, ref.samples, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f="bins.RData", bin.size=bin.size))
submitJobs(abCovCallCases.reg, findNotDone(abCovCallCases.reg) , resources=list(walltime="1:0:0", nodes="1", cores="1",queue="sw", supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(abCovCallCases.reg)

res.df = do.call(rbind, reduceResultsList(abCovCallCases.reg))
save(res.df, file="cnvs-twin-5kbp-FDR001-mergeStitch-thSdest.RData")

#### Abnormal bin calling Loose
## system("rm -rf abCovCallCases05-files")
abCovCallCases05.reg <- makeRegistry(id="abCovCallCases05")
abCovCallCases05.f <- function(samp, files.df, norm.stats.f, bins.f,  bin.size){
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  load(bins.f)
  call.abnormal.cov(files.df=files.df, samp=samp, out.pdf=paste0(samp,"/",samp,"-sdest-abCovCall.pdf"), FDR.th=.05, merge.cons.bins="stitch", z.th="sdest", norm.stats=norm.stats.f, stitch.dist=bin.size*2+1, gc.df = bins.df)
}
batchMap(abCovCallCases05.reg, abCovCallCases05.f, ref.samples, more.args=list(files.df=files.df, norm.stats.f=out.files[2], bins.f="bins.RData", bin.size=bin.size))
submitJobs(abCovCallCases05.reg, findNotDone(abCovCallCases05.reg) , resources=list(walltime="1:0:0", nodes="1", cores="1",queue="sw", supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(abCovCallCases05.reg)

res.df = do.call(rbind,reduceResultsList(abCovCallCases05.reg))
save(res.df, file="cnvs-twin-5kbp-FDR05-mergeStitch-thSdest.RData")


######
## EXTRA : For examples, graphs
######

## Order, compress and index the Z-scores for each sample
## system("rm -rf indexExtra-files")
indexExtra.reg <- makeRegistry(id="indexExtra")
indexExtra.f <- function(dump, bc.f, msd.f){
  library(PopSV, lib.loc="/home/jmonlong/R/PopSVforPaper")
  comp.index.files(bc.f, rm.input=FALSE, reorder=TRUE)
  comp.index.files(msd.f, rm.input=FALSE, reorder=TRUE)
}
batchMap(indexExtra.reg, indexExtra.f, 1, more.args=list(bc.f=out.files[1], msd.f="ref-msd.tsv"))
submitJobs(indexExtra.reg, 1, resources=list(walltime="3:0:0", nodes="1", cores="1",queue="sw", supervisor.group="bws-221-ae"), wait=function(retries) 100, max.retries=10)
showStatus(indexExtra.reg)
