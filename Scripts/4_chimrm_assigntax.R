############
# Packages #
############
library("dada2")
library("DECIPHER"); packageVersion("DECIPHER")
library("plyr")
library("doParallel"); library(doMC)

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

# Start analysis -------------------------------------------------------------------------------------------
# Merge multiple runs
merger <- list()
list.files("./Objects")
# we have 27 categories of Plate_Year_Season
for(i in 1:37){
  merger[[i]] <- readRDS(paste0("./Objects/",i,"_seqtab.rds"))
}

#seqtab <- readRDS("./Objects/25_seqtab.rds")
# any without rownames?
lapply(merger, row.names) # yes, 22

# check sample name
splitdf <- readRDS( "./Objects/splitdf_final.rds")
setDT(splitdf)

# merge all into one
st.all <- mergeSequenceTables(merger[[1]],merger[[2]],merger[[3]],merger[[4]],merger[[5]],
                              merger[[6]],merger[[7]],merger[[8]],merger[[9]],merger[[10]],
                              merger[[11]],merger[[12]],merger[[13]],merger[[14]],merger[[15]],
                              merger[[16]],merger[[17]],merger[[18]],merger[[19]],merger[[20]],
                              merger[[21]],merger[[22]],merger[[23]],merger[[24]], merger[[25]],
                              merger[[26]], merger[[27]], merger[[28]], merger[[29]], merger[[30]],
                              merger[[31]], merger[[32]], merger[[33]],  merger[[34]], merger[[35]],
                              merger[[36]], merger[[37]])
saveRDS(st.all, "./Objects/all_seqtab_2015-18.rds")

dim(st.all) # 1079 samples, 281263 ASVs
nrow(splitdf) # of 1080 samples
# one sample did not pass filtering

table(nchar(getSequences(st.all)))
# a few ASVs have weirdly long sequences, and a few very short
# Sequences that are much longer or shorter than expected may be the result of non-specific priming

# remove non-target-length sequences from sequence table
tar.seqtab <- st.all[,nchar(colnames(st.all)) %in% 250:258] # keep only reads that are 250 to 258 bp long
dim(tar.seqtab) # 271521 ASVs
table(nchar(getSequences(tar.seqtab)))

saveRDS(tar.seqtab, "./Objects/target_seqtab_2015-18.rds")
tar.seqtab <- readRDS("./Objects/target_seqtab_2015-18.rds")
# collapse ASVs with identical sequence
#col.st <- collapseNoMismatch(st.all, minOverlap = 20, verbose = T)
#col.st <- collapseNoMismatch(seqtab, minOverlap = 20, verbose = T)
#saveRDS(col.st, "./Objects/collapsed_seqtab_2018.rds")
# took 4 days for La Romaine 2015-2017 shallow + deep

# remove Chimeras
#col.st <- readRDS("./Objects/collapsed_seqtab_2018.rds")
seqtab.nochim <- removeBimeraDenovo(tar.seqtab, method="consensus", multithread=TRUE, verbose = T)
saveRDS(seqtab.nochim, "./Objects/nochim_seqtab_2015-18.rds")
seqtab.nochim <- readRDS('./Objects/nochim_seqtab_2015-18.rds')
# nochim2 <- readRDS('./Objects/nochim_seqtab_2015-18.rds')
# row.names(nochim2)[!(row.names(nochim2) %in% row.names(seqtab.nochim))]

dim(seqtab.nochim) # 170501 ASVs retained
sum(seqtab.nochim)/sum(tar.seqtab) # chimeras account for 7% of sequence reads - which is fine

# assign taxonomy
#tax <- assignTaxonomy(seqtab.nochim, "./DADA2/DB/silva_nr_v138_train_set.fa.gz", multithread = TRUE)
#saveRDS(tax, "./Objects/taxtab_silva_v138_2018.rds")

tax <- assignTaxonomy(seqtab.nochim, "./DB/rdp_train_set_18.fa.gz", multithread = TRUE)
saveRDS(tax, "./Objects/taxtab_rdp_v18_2015-2018.rds")

# assign taxa with DECIPHER
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("./DB/GTDB_r95-mod_August2020.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest, there is no species
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
#saveRDS(tax, "./Objects/taxtabspecies_silva_v138_2018.rds")
saveRDS(taxid, "./Objects/taxtab_gtdb_r95_2015-18.rds")

dim(seqtab.nochim)
dim(taxid); colnames(taxid)
head(taxid)
any(!is.na(taxid[,7])) #none assigned until species, there is no species column in trainingSet

# If we want until species, dada2 version of assign taxonomy:
taxa <- assignTaxonomy(seqtab.nochim, "./DB/GTDB_bac-arc_ssu_r86.fa.gz", multithread=TRUE)
taxsp <- addSpecies(taxa, "./DB/GTDB_dada2_assignment_species.fa.gz")

saveRDS(taxa, "./Objects/taxtab_dada_gtdb_r86_2015-18.rds")
saveRDS(taxsp, "./Objects/taxtabspecies_dada_gtdb_r86_2015-18.rds")

print("DONE!")
