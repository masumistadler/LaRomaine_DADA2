# add this to every script
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.0") )
# otherwise dependencies are not found

############
# Packages #
############
library("dada2"); packageVersion("dada2")
library("tidyverse")
library("plyr")
library(data.table)
library(doMC); library(parallel)
# submit job from DADA2 folder

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# Set seed for session and reproducibility of permutations
set.seed(3)


pathF <- "./Renamed_2015-2018/cutadapt/forward" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "./Renamed_2015-2018/cutadapt/reverse" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered")

# create a vector with existing samples
sort.files <- function(x, pattern){
  sort(list.files(x, pattern = pattern, full.names = TRUE))
}

fnFs <- lapply(pathF, sort.files, pattern = "_R1.fastq")
baseFs <- unlist(sapply(fnFs, basename))
sample.names <- sapply(sapply(baseFs, strsplit, split = "_"), "[[", 1)
# Sanity check
any(duplicated(sample.names))

splitdf<-readRDS("~/projects/def-pauldel/mstadler/DADA2/Objects/splitdf_new.rds")
#splitdf<-readRDS("~/projects/def-pauldel/mstadler/DADA2/Objects/redo_filt_splitdf.rds")
#splitdf <- readRDS("./Objects/missing_split.rds")
#splitdf<-readRDS("./Objects/splitdf_new.rds")

deep <- splitdf %>% dplyr::filter(seq_depth == "Deep")
df18 <- splitdf %>% dplyr::filter(year == 2018)
splitdf <- rbind(deep, df18)

# create new folders to store data
dir.create(filtpathF)
dir.create(filtpathR)

# Separate data frame content into lists.
pathFs <- as.list(splitdf$pathFs)
pathRs <- as.list(splitdf$pathRs)
filtpathFs <- as.list(splitdf$filtpathFs)
filtpathRs <- as.list(splitdf$filtpathRs)
trimF <- as.list(splitdf$TrimF)
trimR <- as.list(splitdf$TrimR)


length(pathFs) == length(trimF)
# actual filter and trim function of DADA2
mapply(filterAndTrim, fwd = pathFs, filt = filtpathFs, rev = pathRs, filt.rev = filtpathRs,
                   truncLen = c(trimF,trimR), maxEE = c(5,5), truncQ = 2, maxN = 0, rm.phix=TRUE, compress=TRUE,
                   verbose=TRUE, multithread = T)
print("Done!")
quit(save = "no")
# Done
