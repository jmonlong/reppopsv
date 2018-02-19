Download and format public annotations
======================================

Load packages, functions and data
---------------------------------

``` r
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)
library(knitr)
```

Download the data
-----------------

``` r
files.todownload = c("rmsk.txt.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz", 
    "gencode.v19.annotation.gtf.gz", "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz", 
    "cytoBandIdeo.txt.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz", 
    "gap.txt.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz", 
    "simpleRepeat.txt.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz", 
    "genomicSuperDups.txt.gz", "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz", 
    "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", 
    "chaisson_STR_expansions.bed", "http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/STR_expansions.bed", 
    "chaisson_insertions.bed", "http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/insertions.bed", 
    "chaisson_inversions.bed", "http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/inversions.bed", 
    "chaisson_deletions.bed", "http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/deletions.bed")
files.todownload %<>% matrix(ncol = 2, byrow = TRUE) %>% as.data.frame(stringsAsFactors = FALSE)
colnames(files.todownload) = c("out.filename", "url")
kable(files.todownload)
```

| out.filename                                      | url                                                                                                                    |
|:--------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------|
| rmsk.txt.gz                                       | <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz>                                                  |
| gencode.v19.annotation.gtf.gz                     | <ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz>                            |
| cytoBandIdeo.txt.gz                               | <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz>                                          |
| gap.txt.gz                                        | <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gap.txt.gz>                                                   |
| simpleRepeat.txt.gz                               | <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz>                                          |
| genomicSuperDups.txt.gz                           | <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz>                                      |
| ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz | <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz> |
| chaisson\_STR\_expansions.bed                     | <http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/STR_expansions.bed>            |
| chaisson\_insertions.bed                          | <http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/insertions.bed>                |
| chaisson\_inversions.bed                          | <http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/inversions.bed>                |
| chaisson\_deletions.bed                           | <http://eichlerlab.gs.washington.edu/publications/chm1-structural-variation/data/GRCh37/deletions.bed>                 |

``` r
## Download if not there
temp = lapply(1:nrow(files.todownload), function(ii) {
    if (!file.exists(paste0("../data/", files.todownload$out.filename[ii]))) {
        download.file(files.todownload$url[ii], paste0("../data/", files.todownload$out.filename[ii]))
    }
})
```

Repeat masker
-------------

``` r
rm = read.table("../data/rmsk.txt.gz", colClasses = c(rep("NULL", 5), "character", 
    "integer", "integer", "NULL", "NULL", rep("character", 3)))
colnames(rm) = c("chr", "start", "end", "repName", "repClass", "repFamily")
rm %<>% mutate(chr = gsub("chr", "", chr)) %>% filter(chr %in% 1:22) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
save(rm, file = "../data/rm.RData")
rm %>% as.data.frame %>% head %>% kable
```

| seqnames |  start|    end|  width| strand | repName   | repClass       | repFamily      |
|:---------|------:|------:|------:|:-------|:----------|:---------------|:---------------|
| 1        |  10000|  10468|    469| \*     | (CCCTAA)n | Simple\_repeat | Simple\_repeat |
| 1        |  10468|  11447|    980| \*     | TAR1      | Satellite      | telo           |
| 1        |  11503|  11675|    173| \*     | L1MC      | LINE           | L1             |
| 1        |  11677|  11780|    104| \*     | MER5B     | DNA            | hAT-Charlie    |
| 1        |  15264|  15355|     92| \*     | MIR3      | SINE           | MIR            |
| 1        |  16712|  16749|     38| \*     | (TGG)n    | Simple\_repeat | Simple\_repeat |

Gencode v19
-----------

``` r
gc.gtf = import("../data/gencode.v19.annotation.gtf.gz", "GTF")
seqlevels(gc.gtf) = gsub("chr", "", seqlevels(gc.gtf))
gc.gtf = subset(gc.gtf, gene_type == "protein_coding")
mcols(gc.gtf) = mcols(gc.gtf)[, c("type", "gene_name", "gene_id")]
genes = subset(gc.gtf, type == "gene")
exons = subset(gc.gtf, type == "exon")
save(genes, exons, file = "../data/gencodev19-proteincoding-genes-exons.RData")
genes %>% as.data.frame %>% head %>% kable
```

| seqnames |   start|     end|  width| strand | type | gene\_name | gene\_id          |
|:---------|-------:|-------:|------:|:-------|:-----|:-----------|:------------------|
| 1        |   69091|   70008|    918| +      | gene | OR4F5      | ENSG00000186092.4 |
| 1        |  134901|  139379|   4479| -      | gene | AL627309.1 | ENSG00000237683.5 |
| 1        |  367640|  368634|    995| +      | gene | OR4F29     | ENSG00000235249.1 |
| 1        |  621059|  622053|    995| -      | gene | OR4F16     | ENSG00000185097.2 |
| 1        |  738532|  739137|    606| -      | gene | AL669831.1 | ENSG00000269831.1 |
| 1        |  818043|  819983|   1941| +      | gene | AL645608.2 | ENSG00000269308.1 |

``` r
exons %>% as.data.frame %>% head %>% kable
```

| seqnames |   start|     end|  width| strand | type | gene\_name | gene\_id          |
|:---------|-------:|-------:|------:|:-------|:-----|:-----------|:------------------|
| 1        |   69091|   70008|    918| +      | exon | OR4F5      | ENSG00000186092.4 |
| 1        |  137621|  139379|   1759| -      | exon | AL627309.1 | ENSG00000237683.5 |
| 1        |  134901|  135802|    902| -      | exon | AL627309.1 | ENSG00000237683.5 |
| 1        |  367640|  368634|    995| +      | exon | OR4F29     | ENSG00000235249.1 |
| 1        |  621059|  622053|    995| -      | exon | OR4F16     | ENSG00000185097.2 |
| 1        |  739121|  739137|     17| -      | exon | AL669831.1 | ENSG00000269831.1 |

Centromere, telomere, gap annotation
------------------------------------

``` r
chr.band = read.table("../data/cytoBandIdeo.txt.gz", sep = "\t", as.is = TRUE)
colnames(chr.band) = c("chr", "start", "end", "band", "type")
centelgap = read.table("../data/gap.txt.gz", sep = "\t", as.is = TRUE)
centelgap = centelgap[, c(2:4, 8)]
colnames(centelgap) = c("chr", "start", "end", "type")
## Adding telomeres as 10kbp from each chromosomal end
centelgap = rbind(centelgap, chr.band %>% group_by(chr) %>% summarize(start = min(start), 
    end = 10000) %>% mutate(type = "telomere"))
centelgap = rbind(centelgap, chr.band %>% group_by(chr) %>% summarize(start = max(end) - 
    10000, end = max(end)) %>% mutate(type = "telomere"))
## 
centelgap %<>% mutate(chr = gsub("chr", "", chr)) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
save(centelgap, file = "../data/centelgap.RData")
centelgap %>% as.data.frame %>% head %>% kable
```

| seqnames |    start|      end|   width| strand | type     |
|:---------|--------:|--------:|-------:|:-------|:---------|
| 1        |        0|    10000|   10001| \*     | telomere |
| 1        |   177417|   227417|   50001| \*     | clone    |
| 1        |   267719|   317719|   50001| \*     | contig   |
| 1        |   471368|   521368|   50001| \*     | contig   |
| 1        |  2634220|  2684220|   50001| \*     | clone    |
| 1        |  3845268|  3995268|  150001| \*     | contig   |

Simple repeats
--------------

``` r
simprep = read.table("../data/simpleRepeat.txt.gz", as.is = TRUE, sep = "\t")
simprep = simprep[, c(2:4, 17)]
colnames(simprep) = c("chr", "start", "end", "sequence")
simprep %<>% mutate(chr = gsub("chr", "", chr)) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
save(simprep, file = "../data/simprep.RData")
simprep %>% as.data.frame %>% head %>% kable
```

| seqnames |  start|    end|  width| strand | sequence                                                                                                                  |
|:---------|------:|------:|------:|:-------|:--------------------------------------------------------------------------------------------------------------------------|
| 1        |  10000|  10468|    469| \*     | TAACCC                                                                                                                    |
| 1        |  10627|  10800|    174| \*     | AGGCGCGCCGCGCCGGCGCAGGCGCAGAG                                                                                             |
| 1        |  10757|  10997|    241| \*     | GGCGCAGGCGCAGAGAGGCGCGCCGCGCCGGCGCAGGCGCAGAGACACATGCTAGCGCGTCCAGGGGTGGAGGCGT                                              |
| 1        |  11225|  11447|    223| \*     | CGCCCCCTGCTGGCGACTAGGGCAACTGCAGGGTCCTCTTGCTCAAGGTGAGTGGCAGACGCCCACCTGCTGGCAGCCGGGGACACTGCAGGGCCCTCTTGCTTACTGTATAGTGGTGGCA |
| 1        |  11271|  11448|    178| \*     | AGTGGTGGCACGCCACCTGCTGGCAGCTAGGGACACTGCAGGGCCCTCTTGCTCAAGGTAT                                                             |
| 1        |  11283|  11448|    166| \*     | CGCCCCCTGCTGGCAGCTGGGGACACTGCAGGGCCCTCTTGCTCAAGGTATAGTGGCAGCA                                                             |

Segmental Duplications
----------------------

``` r
segdup = read.table("../data/genomicSuperDups.txt.gz", sep = "\t")
segdup = segdup[, c(2:4, 27)]
colnames(segdup) = c("chr", "start", "end", "fracMatch")
segdup %<>% mutate(chr = gsub("chr", "", chr)) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
save(segdup, file = "../data/segdup.RData")
segdup %>% as.data.frame %>% head %>% kable
```

| seqnames |  start|    end|  width| strand |  fracMatch|
|:---------|------:|------:|------:|:-------|----------:|
| 1        |  10000|  87112|  77113| \*     |   0.992904|
| 1        |  10000|  20818|  10819| \*     |   0.981582|
| 1        |  10000|  19844|   9845| \*     |   0.982898|
| 1        |  10000|  19844|   9845| \*     |   0.982898|
| 1        |  10464|  40733|  30270| \*     |   0.987856|
| 1        |  10485|  19844|   9360| \*     |   0.986907|

1000 Genomes Project - CNV catalog
----------------------------------

``` r
tgp = readVcf("../data/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz", "hg19")
tgp.df = cbind(tgp %>% info %>% as.data.frame %>% mutate(type = SVTYPE, end = END) %>% 
    dplyr::select(type, end), rowRanges(tgp) %>% as.data.frame %>% dplyr::select(seqnames, 
    start), geno(tgp)$GT %>% as.data.frame)
rownames(tgp.df) = NULL
nb.samps = ncol(geno(tgp)$GT)
tgp.df = tgp.df %>% mutate(type = ifelse(grepl("DEL", type), "DEL", type)) %>% 
    filter(seqnames != "X", type %in% c("CNV", "DEL", "DUP")) %>% gather(sample, 
    geno, 5:ncol(tgp.df)) %>% filter(geno != "0|0", geno != ".", geno != "0") %>% 
    group_by(seqnames, start, end) %>% mutate(prop = n()/nb.samps, svsize = end - 
    start)
tgp = makeGRangesFromDataFrame(tgp.df, keep.extra.columns = TRUE)
save(tgp, file = "../data/tgp-CNV-noX.RData")
tgp %>% as.data.frame %>% head %>% mutate(geno = gsub("\\|", "\\/", geno)) %>% 
    kable
```

| seqnames |     start|       end|  width| strand | type | sample  | geno |       prop|  svsize|
|:---------|---------:|---------:|------:|:-------|:-----|:--------|:-----|----------:|-------:|
| 1        |   4204667|   4204717|     51| \*     | DEL  | HG00096 | 1/0  |  0.9604633|      50|
| 1        |   6434482|   6445321|  10840| \*     | DEL  | HG00096 | 1/0  |  0.0087859|   10839|
| 1        |   6438160|   6445897|   7738| \*     | DEL  | HG00096 | 1/0  |  0.0119808|    7737|
| 1        |   7570074|   7571521|   1448| \*     | DEL  | HG00096 | 0/1  |  0.2919329|    1447|
| 1        |   8200707|   8211256|  10550| \*     | DUP  | HG00096 | 0/1  |  0.0163738|   10549|
| 1        |  10482499|  10483790|   1292| \*     | DEL  | HG00096 | 1/0  |  0.5371406|    1291|

The catalog is reformatted to have one row per sample/variant. We save the coordinates and SV type information.

To focus on CNV we kept the following *SVTYPE*s: *DEL*, *DUP* and *CNV*. To compare with PopSV we also removed variants smaller than 300bp, in chromosome X, or affecting more than 80% of the samples.

PacBio (Chaisson et al) SV catalog
----------------------------------

``` r
sv.chaisson = rbind(read.table("../data/chaisson_STR_expansions.bed", as.is = TRUE) %>% 
    mutate(type = "str"), read.table("../data/chaisson_inversions.bed", as.is = TRUE) %>% 
    mutate(type = "inversion"), read.table("../data/chaisson_deletions.bed", 
    as.is = TRUE, fill = TRUE, sep = "\t", flush = TRUE) %>% dplyr::select(V1, 
    V2, V3) %>% mutate(type = "deletion"), read.table("../data/chaisson_insertions.bed", 
    as.is = TRUE, fill = TRUE) %>% dplyr::select(V1, V2, V3) %>% mutate(type = "insertion"))
colnames(sv.chaisson)[1:3] = c("chr", "start", "end")
sv.chaisson$chr %<>% gsub("chr", "", .)
save(sv.chaisson, file = "../data/sv.chaisson.RData")
sv.chaisson %>% sample_n(6) %>% kable(row.names = FALSE)
```

| chr |      start|        end| type      |
|:----|----------:|----------:|:----------|
| 14  |   75520940|   75525521| deletion  |
| 16  |     804733|     804734| str       |
| 7   |  116366199|  116366303| insertion |
| 18  |   44528002|   44528003| str       |
| 15  |   30567170|   30567171| str       |
| 8   |     942889|     942995| deletion  |

OMIM disease genes
------------------

``` r
library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", 
    host = "Aug2017.archive.ensembl.org")
omim = getBM(attributes = c("hgnc_symbol", "mim_morbid_description"), mart = ensembl)
mergeDiseaseDesc <- function(desc) {
    desc %>% gsub(";;", ";", .) %>% strsplit(";") %>% unlist %>% unique %>% 
        paste(collapse = ";")
}
omim %<>% filter(mim_morbid_description != "") %>% group_by(hgnc_symbol) %>% 
    summarize(mim_morbid_description = mergeDiseaseDesc(mim_morbid_description))
save(omim, file = "../data/omim-genes.RData")
```

There are 3438 OMIM genes.
