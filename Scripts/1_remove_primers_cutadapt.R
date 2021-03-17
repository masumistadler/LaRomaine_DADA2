############
# Packages #
############
library("dada2"); packageVersion("dada2")
library("plyr"); library("dplyr")
library("ShortRead"); packageVersion("ShortRead")
library("Biostrings"); packageVersion("Biostrings")

##################
# Remove primers #
##################

path <- "./Renamed_2015-2018"  ## CHANGE ME to the directory containing the fastq files.

# Multi directory
#path <- list("/home/bioinf/data/Bioinf.LaRomaine/Raw/2015",
#             "/home/bioinf/data/Bioinf.LaRomaine/Raw/2016",
#             "/home/bioinf/data/Bioinf.LaRomaine/Raw/2017") ## CHANGE ME to the directory containing the fastq files.
# path.list <- lapply(list.files, path, pattern = "*fastq.gz")

path.list <- list.files(path, pattern = "*fastq.gz")

# sort files according to their names
sort.files <- function(x, pattern){
  sort(list.files(x, pattern = pattern, full.names = TRUE))
}

# extract file names separately for forward and reverse reads
# Multi-dir
#fnFs <- lapply(path, sort.files, pattern = "_R1.fastq.gz")
#fnRs <- lapply(path, sort.files, pattern = "_R2.fastq.gz")

fnFs <- sort.files(path, pattern = "_R1.fastq.gz")
fnRs <- sort.files(path, pattern = "_R2.fastq.gz")

# Example
fnFs[1] ; fnRs[1]
# Sanity check, same length?
length(fnFs) == length(fnRs)

# identify unique sequence of used primer
FWD <- "GTGCCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...

# identify all potential orientations of primer
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# identify file names without full path
# basename() takes the string after the last separator
# Multi dir
#baseFs <- lapply(fnFs, basename)
#baseRs <- lapply(fnRs, basename)

baseFs <- basename(fnFs)
baseRs <- basename(fnRs)

# Example
baseFs[2] ; baseRs[2]

# create a "filtN" folder where all filtered read files are going to be stored
# mapply() for multivariate objects as file.path and baseFs both are lists
# We will create a new folder in this directory
dir.create(file.path(path,"filtN"), showWarnings = F)

# Multi dir
#fnFs.filtN <- mapply(file.path, path, "filtN", baseFs) # Put N-filterd files in filtN/ subdirectory
#fnRs.filtN <- mapply(file.path, path, "filtN", baseRs)

fnFs.filtN <- file.path(path, "filtN",baseFs)
fnRs.filtN <- file.path(path, "filtN",baseRs)

# Example
fnFs.filtN[2] ; fnRs.filtN[2]

# actually filter and save reads
mapply(filterAndTrim, fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

# how often does the primer or reverse complement appear in our reads?
# we only inspect the first sample of each year
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Blank
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
# Sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[3]]))

#Multi-directory
#ff <- sapply(fnFs.filtN, "[[", 1) # extracts the first element in filtN folder
#fr <- sapply(fnRs.filtN, "[[", 1)
#for(i in 1:3){
#  x <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = ff[i]),
#             FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fr[i]), 
#             REV.ForwardReads = sapply(REV.orients, primerHits, fn = ff[i]), 
#             REV.ReverseReads = sapply(REV.orients, primerHits, fn = fr[i]))
#  print(x)
#}

## Remove primers ##

# load cutadapt
cutadapt <- "/home/bioinf/.local/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(command = cutadapt, args = "--version") # Run shell commands from R

# create folder cutadapt where the primer free sequences will be stored
#path.cut <- lapply(path, file.path, "cutadapt")
#lapply(path.cut, dir.create) # creates cutadapt folder
path.cut <- file.path(path, 'cutadapt')
dir.create(path.cut)

# preparation for cutting the primer
fnFs.cut <- file.path(path.cut, baseFs)
fnRs.cut <- file.path(path.cut, baseRs)
#Multi-dir
#fnFs.cut <- mapply(file.path, path.cut, baseFs) # combine the path and the sample names for FWD
#fnRs.cut <- mapply(file.path, path.cut, baseRs) # same for REV
FWD.RC <- dada2:::rc(FWD) # extract seqeuence of reverse complement of FWD
REV.RC <- dada2:::rc(REV) # extract seqeuence of reverse complement of REV

# create strings for the arguments fed to cutadapt
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g",FWD,"-a",REV.RC, sep = " ") 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G",REV,"-A",FWD.RC, sep = " ") 
# seems like the arg = argument of system2() needs the command arguments as elements of a vector
# meaning, if we paste it together as below it wont work
# I couldn't figure out how to make c() work with mapply so we do a work around
# first, paste together all elements into a string with an unique pattern '#' in between for later splitting

# -m argument is included as in 2015 samples, we have a lot of short reads
# it removes all reads that are shorter than 125 bases, I chose 125 as the minimum bases needed
# to have an overlap when merging paired-end reads is 127 for the V4 region of 16S rRNA

# -j 0 is added to allow multithread processing
# make sure to install pigz before if you work with compressed files
cut.arg <- paste(R1.flags, R2.flags, "-n", "2", "-m", "125", "-j", "0", "-o",
                 fnFs.cut, "-p", fnRs.cut, fnFs.filtN, fnRs.filtN, sep = "#")
# Multi-dir
#cut.arg <- mapply(paste, R1.flags, R2.flags, "-n", "2", "-m", "125", "-j", "0", "-o", fnFs.cut, "-p", fnRs.cut, fnFs.filtN, fnRs.filtN, sep = "#")  #
#names(cut.arg) <- c("y2015", "y2016", "y2017") # rename bin names for easier overview

# then we split by # into elements to combine later into a vector for cutadapt 
cut.arg <- strsplit(cut.arg, split = "[#]")
# Multi-dir
#cut.arg <- mapply(strsplit, cut.arg, split = "[#]") # create nested list

# could not get
# mapply(system2, command = cutadapt, args = cut.arg) to work...
# so we go the old fashioned way of looping

# to go through each file, we use two loops to access the bins inside bins
# NOTE: cutadapt prints a lot of text including an Unicode error that is not problematic at all
# Multi-dir
#for(i in 1:length(cut.arg)){
#  for(j in 1:length(cut.arg[[i]])){
#    system2(cutadapt, args = c(unlist(cut.arg[[i]][[j]])))
#  }
#}

for(i in 1:length(cut.arg)){
    system2(cutadapt, args = c(unlist(cut.arg[[i]])))
}

# as a sanity check, we repeat the same excercise as above
# Blank
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
# Sample
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[3]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[3]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[3]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[3]]))

#multi-dir
#ff <- sapply(fnFs.cut, "[[", 1)
#fr <- sapply(fnRs.cut, "[[", 1)
#for(i in 1:3){
#  x <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = ff[i]),
#             FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fr[i]), 
#             REV.ForwardReads = sapply(REV.orients, primerHits, fn = ff[i]), 
#             REV.ReverseReads = sapply(REV.orients, primerHits, fn = fr[i]))
#  print(x)
#}

# All zeroes? You've removed the primers successfully!