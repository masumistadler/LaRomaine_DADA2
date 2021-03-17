############
# Packages #
############
library("dada2"); packageVersion("dada2")
library("tidyverse")
library("plyr")
library(data.table)
library(doMC); library(parallel)

#####################
# Quality filtering #
#####################

### Parallel set-up -----------------------------------------------------------------------
# register number of cores for parallel computing with 'apply' family
detectCores() # 24 (private PC: 4), we do not have 24 cores. It's 12 cores, 24 threads.
registerDoMC()
getDoParWorkers() # 12 (private PC: 2)

# prepare for parallel processing (with 'parallel') for 'vegan' functions
numCores <- detectCores()
cl <- makeCluster(numCores, type = "FORK") # using forking

# Set seed for session and reproducibility of permutations
set.seed(3)


# Some pre-processing

# in terminal
# @carrbas:~$ cd ./cutadapt
# mkdir forward
# mkdir reverse
# mv *_R1.fastq ./forward
# mv *_R1.fastq ./reverse


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

# read in master file for splitting by ID
#splitdf <- read.csv("./Meta/master_seqnames_splitID.csv", sep = ",", dec = ".", stringsAsFactors = F)
splitdf <- read.csv("./Meta/LaRomaine_molecular_naming_unified.csv", sep = ",", dec = ".", stringsAsFactors = F)

# any sample name not in the data frame?
any(splitdf$final_name %in% sample.names == F)
any(sample.names %in% splitdf$final_name == F)

# drop paired-end rows, keep one row per sample
splitdf <- splitdf %>% 
  select(final_name, dr_match_name, replicate, plate.id, year, dna_type, seq_depth) %>%
  distinct()

# Sanity check
nrow(splitdf) == length(sample.names)

# Add new paths
splitdf$pathFs <- file.path(pathF, paste0(splitdf$final_name, "_R1.fastq.gz"))
splitdf$pathRs <- file.path(pathR, paste0(splitdf$final_name, "_R2.fastq.gz"))
splitdf$filtpathFs <- file.path(pathF, "filtered", paste0(splitdf$final_name, "_R1.fastq.gz"))
splitdf$filtpathRs <- file.path(pathR, "filtered", paste0(splitdf$final_name, "_R2.fastq.gz"))

# add a split ID to run DADA2
# Read in meta data
meta <- read.csv("./Meta/main_bac_match_2021-02-22.csv", sep = ",", stringsAsFactors = F)
# select only necessary columns
meta <- meta %>% select(dna.match, sample.type.year, year, Season, sampling.date, lat, long)


# merge with splitdf
# sanity check
any(splitdf$dr_match_name %in% meta$dna.match == F)

setDT(meta); setDT(splitdf)
splitdf <- splitdf[meta, , on = c(dr_match_name = "dna.match")]

# any years don't match?
splitdf[year != i.year,] # all good
splitdf[, i.year := NULL] # remove
# export a portion
export <- splitdf %>% select(seq_name =  final_name,
                             dr_match_name, replicate, dna_type, seq_depth, plate.id,
                             year, season = Season, sampling.date, lat, long, sample.type.year)
write.table(export, "./Meta/sequence_metadata.csv", sep = ",", dec = ".", row.names = F)

# by plate, season and dna_type, except Deep sequences they run by plate
splitdf[seq_depth == "Deep", splitID := as.character(plate.id)]
splitdf[seq_depth == "Shallow", splitID := paste(plate.id, Season, dna_type, sep = "_")]

# pick a single sample for each splitID to determine where to cut the reads due to dropping quality
pickRandomRows <- function(df, numberOfRows = 1){
  df %>% slice(runif(numberOfRows, 0,  length(df[,1])))
}

# It seems that those samples with two genome quebec replicates had a low quality, and thus were run a second time
# In those cases, the second replicate is always of better quality
dupl <- splitdf[splitdf$dr_match_name %in% splitdf[grep("s2",splitdf$final_name),]$dr_match_name,]
qc <- plyr::alply(dupl[1:50,c("pathFs","pathRs")], .margins = 1, .fun = plotQualityProfile,
                  .parallel = TRUE)

# make two different groups for duplicates
splitdf[dr_match_name %in% unique(dupl$dr_match_name) & replicate == 1, splitID := paste(splitID,"LQ", sep = "_")]
# L330 is a true duplicate but not from GQ, correct
splitdf[final_name == "L330R",splitID := "7_Summer_cDNA"]
splitdf[final_name == "L330D",splitID := "5_Summer_DNA"]

# pick a random sample for each splitID
qualcheck <- splitdf %>% group_by(splitID) %>% 
  dplyr::slice(c(1)) %>%
  select(pathFs, pathRs, splitID)

# qualityPlot
qc <- plyr::alply(qualcheck[,c("pathFs","pathRs")], .margins = 1, .fun = plotQualityProfile,
                  .parallel = TRUE)
names(qc) <- qualcheck$splitID

# sanity check
names(qc) == sort(unique(splitdf$splitID))

# go through them
qc[[17]] ; names(qc)[17]
# plate 2, s2 seems better than s1

# Prepare for filtering
# add Trimming information, decided based on quality plots c(FWD, REV)
trim <- data.frame(splitID = sort(unique(splitdf$splitID)),
                   TrimF = c(225,150,225,150,225,
                             225,225,225,225,225,
                             225,225,150,160,160,
                             190,150,225,225,225,
                             225,225,225,225,225,
                             200,225),
                   TrimR = c(225,150,225,150,225,
                             225,225,225,225,225,
                             225,225,150,180,160,
                             190,180,225,225,225,
                             225,225,225,225,225,
                             225,225), stringsAsFactors = F)
# old
#trim <- data.frame(splitID = sort(unique(splitdf$splitID)),
#                   TrimF = c(225,225,200,210,210,
#                             220,210,190,225,225,
#                             220,225,225,225,225,
#                             225,225,225,225,225,
#                             225,225,225,225,
#                             210,225),
#                   TrimR = c(225,225,210,160,160,
#                             150,220,225,225,225,
#                             225,225,225,225,225,
#                             225,225,225,225,225,
#                             225,225,225,225,
#                             180, 180), stringsAsFactors = F)
splitdf <- left_join(splitdf, trim, by = "splitID")
saveRDS(splitdf, "./Objects/splitdf_new.rds")

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

## Filter
#`mapply` runs through the lists.

#system.time(for(i in 1:nrow(splitdf)){
#  filterAndTrim(fwd = splitdf$pathFs[i], filt = splitdf$filtpathFs[i], rev = splitdf$pathRs[i],
#                filt.rev = splitdf$filtpathRs[i], truncLen = c(splitdf$TrimF[i],splitdf$TrimR[i]),
#                maxEE = 2, truncQ = 2, maxN = 0, rm.phix=TRUE, compress=TRUE, verbose=TRUE,
#                multithread=TRUE)
#})

length(pathFs) == length(trimF)

system.time(mapply(filterAndTrim, fwd = pathFs, filt = filtpathFs, rev = pathRs, filt.rev = filtpathRs,
                   truncLen = c(trimF,trimR), maxEE = 2, truncQ = 2, maxN = 0, rm.phix=TRUE, compress=TRUE,
                   verbose=TRUE, multithread = T))
#     user    system   elapsed 
#59498.305  1256.957 60722.193
# around 16.8 hours

