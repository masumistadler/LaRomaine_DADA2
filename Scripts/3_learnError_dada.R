
############
# Packages #
############
library("dada2"); packageVersion("dada2")
library("tidyverse")
library("plyr")
library("doParallel")

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

#Some pre-processing. Read-in `splitdf.rds` created in `4_quality_filtering.Rmd` for pooling.
#splitdf <- readRDS("./Objects/splitdf.rds")
splitdf <- readRDS( "./Objects/splitdf_new.rds")

# La Romaine 2018
#splitdf <- readRDS("./Objects/splitdf_2018.rds")
# LCare
#splitdf <- readRDS("./Objects/splitdf_LCare.rds")
#colnames(splitdf)[1:2] <- c("DadaNames","splitID")

# Extract only sample names.
baseFs <- unlist(sapply(splitdf$filtpathFs, basename))
sample.names <- sapply(sapply(baseFs, strsplit, split = "_"), "[[", 1)

# split samples into list bins to pool them by splitID
#Fs <- daply(splitdf[,c("num.splitID","filtpathFs")], .(num.splitID), .fun = list)
Fs <- daply(splitdf[,c("splitID","filtpathFs")], .(splitID), .fun = list)
Fs <- sapply(Fs, "[[", 2) # extract as vectors

#Rs <- daply(splitdf[,c("num.splitID","filtpathRs")], .(num.splitID), .fun = list)
Rs <- daply(splitdf[,c("splitID","filtpathRs")], .(splitID), .fun = list)
Rs <- sapply(Rs, "[[", 2)


# Run loop by bins
system.time(for(i in 1:length(Fs)){
  
  sample.names <- sapply(strsplit(basename(Fs[[i]]),"_"),`[`,1)
  sample.namesR <- sapply(strsplit(basename(Rs[[i]]),"_"),`[`,1)
  if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
  names(Fs[[i]]) <- sample.names
  names(Rs[[i]]) <- sample.namesR
  
  set.seed(3)
  
  # Learn forward error rates, save error plot and as R object
  errF <- learnErrors(Fs[[i]], nbases=1e8, multithread=TRUE, MAX_CONSIST = 20)
  saveRDS(errF, file = paste0("./Objects/",i,"_errF.rds"))
  ggsave(paste0("./Figures/Errors/",i,"_errF.png"), plot = plotErrors(errF, nominalQ = TRUE),
        units = "cm", width = 40, height = 40)
  
  #saveRDS(errF, file = paste0("../LCare/Objects/",i,"_errF.rds"))
  #ggsave(paste0("../LCare/Figures/",i,"_errF.png"), plot = plotErrors(errF, nominalQ = TRUE),
  #       units = "cm", width = 40, height = 40)
  
  # Learn reverse error rates, save error plot and as R object
  errR <- learnErrors(Rs[[i]], nbases=1e8, multithread=TRUE, MAX_CONSIST = 20)
  ggsave(paste0("./Figures/Errors/",i,"_errR.png"), plot = plotErrors(errR, nominalQ = TRUE),
         units = "cm", width = 40, height = 40)
  saveRDS(errR, file = paste0("./Objects/",i,"_errR.rds"))
  #ggsave(paste0("../LCare/Figures/",i,"_errR.png"), plot = plotErrors(errR, nominalQ = TRUE),
  #       units = "cm", width = 40, height = 40)
  #saveRDS(errR, file = paste0("../LCare/Objects/",i,"_errR.rds"))
  
  # Dereplicate
  # Name the derep-class objects by the sample names
  derepFs <- derepFastq(Fs[[i]], verbose=TRUE)
  names(derepFs) <- sample.names
  saveRDS(derepFs, file = paste0("./Objects/",i,"_derepFs.rds"))
  #saveRDS(derepFs, file = paste0("../LCare/Objects/",i,"_derepFs.rds"))
  
  derepRs <- derepFastq(Rs[[i]], verbose=TRUE)
  names(derepRs) <- sample.namesR
  saveRDS(derepRs, file = paste0("./Objects/",i,"_derepRs.rds"))
  #saveRDS(derepFs, file = paste0("../LCare/Objects/",i,"_derepRs.rds"))
  
  # run dada on pooled samples
  poolFs <- dada(derepFs, err=errF, pool=TRUE, multithread = TRUE)
  saveRDS(poolFs, file = paste0("./Objects/",i,"_poolFs.rds"))
  #saveRDS(poolFs, file = paste0("../LCare/Objects/",i,"_poolFs.rds"))
  
  poolRs <- dada(derepRs, err=errR, pool=TRUE, multithread = TRUE)
  saveRDS(poolRs, file = paste0("./Objects/",i,"_poolRs.rds"))
  #saveRDS(poolFs, file = paste0("../LCare/Objects/",i,"_poolRs.rds"))
  
  # merge paired-end reads
  mergers <- mergePairs(poolFs, derepFs, poolRs, derepRs)
  saveRDS(mergers, file = paste0("./Objects/",i,"_mergers.rds"))
  #saveRDS(mergers, file = paste0("../LCare/Objects/",i,"_mergers.rds"))
  
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab, file = paste0("./Objects/",i,"_seqtab.rds"))
  #saveRDS(seqtab, file = paste0("../LCare/Objects/",i,"_seqtab.rds"))
  print(paste("Finished... '", names(Fs)[i], "'. Number ", i, " in queue of ", length(Fs), ".", sep = ""))
})

#      user     system    elapsed 
# 5474158.69   31696.96 1242346.18  
# elapsed: 14.3 days
