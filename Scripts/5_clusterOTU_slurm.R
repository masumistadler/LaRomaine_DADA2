# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
# install packages in a separate run
#install.packages("data.table", lib = "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.0",
#                 repos = "https://cloud.r-project.org/")
#install.packages("BiocManager", lib = "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.0",
#                 repos = "https://cloud.r-project.org/")
#BiocManager::install("DECIPHER", lib = "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.0")

# add this to every script
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.0") )
# otherwise dependencies are not found

pckgs <- list("DECIPHER", "Biostrings", # OTU clustering
              "data.table")
### Load
invisible(lapply(pckgs, require,
                 character.only = T))

rm(pckgs)

# Set seed for session and reproducibility of permutations
# (just for consistency, no random iteration in this script)
set.seed(3)

# 2. Read and prepare data ---------------------------------------------------------------

# Read in dada2 output
seqtab <- readRDS("./Objects/nochim_seqtab_paper1.rds") # ASV table with raw sequences
tax <-
  readRDS("./Objects/taxtab_paper1_gtdb_r95_2018.rds") # assigned taxonomy

# extract taxonomy
tax <- as.data.frame(tax, stringsAsFactors = F)

# all at once
asv_sequences <- colnames(seqtab)
sample_names <- rownames(seqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

## Find clusters of ASVs to form the new OTUs
aln <- DECIPHER::AlignSeqs(dna, processors = NULL)
saveRDS(aln,"./Objects/aln_decipher.rds")

d <- DECIPHER::DistanceMatrix(aln, processors = NULL)
saveRDS(d,"./Objects/distmat_decipher.rds")

rm(aln)

clusters <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.01,
  # use `cutoff = 0.01` for a 99% OTU ; 0.003 for 99.7%
  # cutoff 0.02 for 98%
  # cutoff 0.03 for 97%
  processors = NULL)
saveRDS(clusters, "./Objects/OTU_clusters_99_decipher_paper1.rds")
