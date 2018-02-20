## This is the function that was used to import and merge the calls from other methods. The arguments are the path to the relevant files:
## others.f: the output file with the merged calls.
## files.f: the file with the normal/tumor information.

library(dplyr)
library(magrittr)

formatOtherCalls <- function(others.f, files.f, freec.s=NULL,freec.l=NULL,cnmops.s=NULL,cnmops.l=NULL, cnvnator.f=NULL, lumpy.f=NULL){
  others.df = NULL

  load(files.f)
  stringents = looses = files.df$sample
  if(any(colnames(files.df)=="status")){
    stringents = subset(files.df, status=="normal")$sample
    looses = subset(files.df, status!="normal")$sample
  }

  if(!is.null(freec.s)){
    res.m = read.table(freec.s, header=TRUE, as.is=TRUE, sep="\t")
    res.m = subset(res.m, sample %in% stringents)
    res.m$start = res.m$start + 1
    res.m$sig = NA
    others.df = rbind(others.df, data.frame(method="FREEC", set="stringent", res.m[, c("sample","chr","start","end","cn","sig")]))
  }

  if(!is.null(freec.l)){
    res.m = read.table(freec.l, header=TRUE, as.is=TRUE, sep="\t")
    res.m = subset(res.m, sample %in% looses)
    res.m$start = res.m$start + 1
    res.m$sig = NA
    others.df = rbind(others.df, data.frame(method="FREEC", set="loose", res.m[, c("sample","chr","start","end","cn","sig")]))
  }

  if(!is.null(cnmops.s)){
    res.m = load(cnmops.s)
    res.m = as.data.frame(eval(as.name(res.m[1])))[,c(1:3,6,9)]
    colnames(res.m) = c("chr","start","end","sample","cn")
    res.m = subset(res.m, sample %in% stringents)
    res.m$cn = as.numeric(gsub("CN","",res.m$cn))
    res.m$sig = NA
    others.df = rbind(others.df, data.frame(method="cn.MOPS", set="stringent", res.m[, c("sample","chr","start","end","cn","sig")]))
  }
  if(!is.null(cnmops.l)){
    res.m = load(cnmops.l)
    res.m = as.data.frame(eval(as.name(res.m[1])))[,c(1:3,6,9)]
    colnames(res.m) = c("chr","start","end","sample","cn")
    res.m = subset(res.m, sample %in% looses)
    res.m$cn = as.numeric(gsub("CN","",res.m$cn))
    res.m$sig = NA
    others.df = rbind(others.df, data.frame(method="cn.MOPS", set="loose", res.m[, c("sample","chr","start","end","cn","sig")]))
  }

  if(!is.null(cnvnator.f)){
    res.m = read.table(cnvnator.f, header=TRUE, as.is=TRUE, sep="\t")
    coord.l = strsplit(res.m$coordinates, ":")
    res.m$chr = unlist(lapply(coord.l, "[",1))
    pos.l = strsplit(unlist(lapply(coord.l, "[", 2)), "-")
    res.m$start = as.numeric(unlist(lapply(pos.l, "[",1)))
    res.m$end = as.numeric(unlist(lapply(pos.l, "[",2)))
    res.m$cn = 2 * res.m$normalized_RD
    res.m$sig = ifelse(res.m$eval1<res.m$eval2, res.m$eval1, res.m$eval2)
    others.df = rbind(others.df, data.frame(method="CNVnator", set="loose", subset(res.m, sample %in% looses)[, c("sample","chr","start","end","cn","sig")]))
    CNVNATOR.TH = 1e-5
    res.m = subset(res.m, sample %in% stringents & sig < CNVNATOR.TH)
    others.df = rbind(others.df, data.frame(method="CNVnator", set="stringent", res.m[, c("sample","chr","start","end","cn","sig")]))
  }

  if(!is.null(lumpy.f)){
    res.m = read.table(lumpy.f, header=TRUE, as.is=TRUE, sep="\t")
    colnames(res.m)[2:3] = c("chr","start")
    res.m$end = gsub(".*;END=([^;]+).*", "\\1", res.m$INFO)
    res.m$type = gsub(".*SVTYPE=([^;]+).*", "\\1", res.m$INFO)
    res.m$su = as.numeric(gsub(".*;SU=([^;]+).*", "\\1", res.m$INFO))
    ## DEL/DUP
    res.m.dd = subset(res.m, type!="BND")
    res.m.dd$end = as.numeric(res.m.dd$end)
    ## BND: if intra-chr, could be CNVs
    res.m.bnd = subset(res.m, type=="BND")
    res.m.bnd$event = gsub(".*EVENT=([^;]+).*", "\\1", res.m.bnd$INFO)
    res.m.bnd %<>% group_by(sample, type, event) %>% filter(length(unique(chr))==1) %>% summarize(chr=chr[1], end=max(start), start=min(start), su=min(su)) %>% ungroup %>% filter(end-start<4e6)
    ## Merge back
    res.m = rbind(select(res.m.bnd, sample, chr, start, end, type, su), select(res.m.dd, sample, chr, start, end, type, su))
    res.m = subset(res.m, chr %in% 1:22 & end-start>300)
    res.m$cn = ifelse(res.m$type=="DEL",0,3)
    res.m$sig = -res.m$su
    others.df = rbind(others.df, data.frame(method="LUMPY", set="loose", subset(res.m, sample %in% looses)[, c("sample","chr","start","end","cn","sig")]))
    LUMPY.TH = 4
    res.m = subset(res.m, su>LUMPY.TH & sample %in% stringents)
    others.df = rbind(others.df, data.frame(method="LUMPY", set="stringent", res.m[, c("sample","chr","start","end","cn","sig")]))
  }

  save(others.df, file=others.f)
}
