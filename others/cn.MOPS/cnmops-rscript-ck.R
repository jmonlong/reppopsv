args = commandArgs(TRUE)
files.f = args[1]  ## The tsv file with GC-corrected bin counts for all samples.
th.arg = as.numeric(args[2]) ## the calling threshold for copy number, 0 ~ the default/recommended.
outF = args[3] ## Output RData file

library(cn.mops)
library(parallel)

message("Read one sample...")

load(files.f)
bc1 = read.table(files.df$bc.gc.gz[1], as.is=TRUE, header=TRUE)

bc = lapply(files.df$bc.gc.gz, function(ff){
  message(ff)
  bc.f = read.table(ff, as.is=TRUE, header=TRUE)
  if(!all(bc.f$chr == bc1$chr & bc.f$start==bc1$start)){
    stop("Different order in bin count files, careful merging necessary.")
  }
  bc.f$bc
})

bc = do.call(cbind, bc)
colnames(bc) = files.df$sample

if(any(is.na(bc))){
    bc[is.na(bc)] = 0
}

bc.gr = GRanges(seq=bc1$chr,IRanges(start=bc1$start,end=bc1$end))
mcols(bc.gr) = bc

uTh = log2((3-th.arg)/2)
lTh = log2((1+th.arg)/2)

message("cn.MOPS...")

system.time({res.mops = cn.mops(bc.gr,upperThreshold=uTh,lowerThreshold=lTh)}) 
system.time({res.mops = calcIntegerCopyNumbers(res.mops)}) 

res.cnvs = cnvs(res.mops)

save(res.cnvs,file=outF)

sessionInfo()
