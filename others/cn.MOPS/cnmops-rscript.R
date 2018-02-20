args = commandArgs(TRUE)
bcF = args[1]  ## The tsv file with GC-corrected bin counts for all samples.
th.arg = as.numeric(args[2]) ## the calling threshold for copy number, 0 ~ the default/recommended.
outF = args[3] ## Output RData file

library(cn.mops)

message("Prepare data...")

bc = read.table(bcF,as.is=TRUE,header=TRUE)
bc.coord = bc[,1:3]
bc = as.matrix(bc[,-(1:3)])

if(any(is.na(bc))){
    bc[is.na(bc)] = 0
}

bc.gr = GRanges(seq=bc.coord$chr,IRanges(start=bc.coord$start,end=bc.coord$end))
mcols(bc.gr) = bc

head(bc.gr)

uTh = log2((3-th.arg)/2)
lTh = log2((1+th.arg)/2)

message("cn.MOPS...")

system.time({res.mops = cn.mops(bc.gr,upperThreshold=uTh,lowerThreshold=lTh)}) 
system.time({res.mops = calcIntegerCopyNumbers(res.mops)}) 

res.cnvs = cnvs(res.mops)

save(res.cnvs,file=outF)

sessionInfo()
