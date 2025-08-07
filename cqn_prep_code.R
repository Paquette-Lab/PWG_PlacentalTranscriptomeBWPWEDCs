library(AnnotationHub)
hub <- AnnotationHub()
ensdb <- hub[["AH78783"]]
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
## get exon positions
ex <- exonsBy(ensdb, "gene")
## remove haplotypes and unplaced scaffolds
ex <- keepStandardChromosomes(ex, pruning.mode = "coarse")
## The exons can be overlapping, and we don't want to double count
ex <- reduce(ex)
## change to UCSC style (e.g., chr1 vs 1) to match BSgenome
seqlevelsStyle(ex) <- "UCSC"
## load("candle_20200226.Rdata") ## originally written by Jim for CANDLE WGCNA
load(here::here("Resubmission_July2024/NEW_RawData/gapps_20240325.Rdata"))
## order to match CANDLE data
## update to match GAPPS data
ex <- ex[gsub("\\.[0-9]+$", "", row.names(counts$counts))]
## get sequences
seq <- getSeq(Hsapiens, ex)

## these sequences are in a DNAStringSetList, which is a list
## of DNAStringSet objects. The alphabetFrequency function won't
## dispatch on DNAStringSetList objects, so we have to get the
## alphabet frequency for over 60k DNAStringSet objects, so we parallelize

## first a function
fun <- function(x) {
    tmp <- alphabetFrequency(x)
    ## some genes only have one exon, so alphabetFrequency returns a one-row
    ## matrix and colSums doesn't like that1
    if(nrow(tmp) == 1L) return(tmp[1,1:4]) else return(colSums(tmp[,1:4]))
}

## now parallelize using bplapply
## I have 20 cores available, so I use them all
## update 20240430: mariana changed MulticoreParam(20) to 1 
## (because only 1 core on local computer)
seq2 <- do.call(rbind, bplapply(seq, fun, BPPARAM = MulticoreParam(1)))

## now a dumb function to get the length and GC content
fun2 <- function(x) {
    rs <- rowSums(x)
    return(cbind(rs, rowSums(x[,c("C","G")])/rs))
}

forcqn <- fun2(seq2)
colnames(forcqn) <- c("length","gccontent")
saveRDS(forcqn, here::here("Resubmission_July2024/intermediateData/forcqn.Rds"))
