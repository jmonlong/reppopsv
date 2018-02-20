library(BatchJobs)
library(PopSV)

## The samples for EpiPopSV
controls = scan("controls.samples", "a")
affected = scan("affected.samples","a")

## Init file names
bam.files = data.frame(project="control", sample=controls, stringsAsFactors=FALSE)
bam.files = rbind(bam.files, data.frame(project="affected", sample=affected, stringsAsFactors=FALSE))
bam.files$bam = "noneed"
files.df = init.filenames(bam.files, code="EpiPopSV")

## Find already computed bin counts
path = "/RQexec/GROUP/pacoss/COMMUN/SV.analysis/wgs/PopSV/primary.analysis/5kb_bins/"
samples = list.files(path,"^S")
samples = grep(".pdf", samples, invert=TRUE, value=TRUE)
bins.f = paste0(path, samples, "/", samples, "-5kb_bins-bc-gcCor.tsv.bgz")
bins.f = bins.f[file.exists(bins.f)]
names(bins.f) = gsub(".*/(S.+)-5kb_bins-bc-gcCor.tsv.bgz", "\\1", bins.f)

## Update files.df
files.df$bc.gc.gz = NA
files.df$bc.gc.gz = bins.f[files.df$sample]
all(!is.na(files.df$bc.gc.gz))
save(files.df, file="files.RData")

## Get bins coordinates
bins.df = read.table(files.df$bc.gc.gz[1], as.is=TRUE, header=TRUE)
bins.df = bins.df[,-4]
save(bins.df, file="bins.RData")

source("automatedPipeline.R")

## Get GC content per bin (for later)
getGC = autoGCcounts("files.RData", "bins.RData", skip=2)

## Merge, normalize and call everything
res.df = autoNormTest("files.RData", "bins.RData", ref.samples=controls, FDR.th=.01, step.walltime=c(10,15,20,6,1,1), step.cores=c(6,6,12,6,1,1), skip=2)

## Write output file
res.df = merge(res.df, files.df[,c("sample","project")])
write.table(res.df, file="cnvs-5kb_bins-EpiPopSV-FDR01.tsv", sep="\t", quote=FALSE, row.names=FALSE)


## Quick count
bc.df = autoExtra("files.RData", "bins.RData", do=1)
save(bc.df, file="quickcount-1k.RData")
## Split BC norm for ref samples
splitRes= autoExtra("files.RData", "bins.RData", do=2)
## Quick count
bc.df = autoExtra("files.RData", "bins.RData", do=1, col.files="bc.gc.norm.gz")
save(bc.df, file="quickcount-1k-afterTN.RData")
## Order/compress/index BC norm for ref samples
ociRes= autoExtra("files.RData", "bins.RData", do=3)
