### PHYLOGENETIC TREE ###
# Slurm set-up -------------------------------------------------------------------------------------------------------
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.2") )
# otherwise dependencies are not found

# R set-up -----------------------------------------------------------------------------------------------------------
pckgs <- list("DECIPHER", "Biostrings", "dada2", # OTU clustering
              "data.table", "phangorn")
### Load
invisible(lapply(pckgs, require,
                 character.only = T))

rm(pckgs)

# Set seed for session and reproducibility of permutations
# (just for consistency, no random iteration in this script)
set.seed(3)

# Construct a phylogenetic tree ----------------------------------------------------------------------------------------
# Phylogenetic relatedness is commonly used to inform downstream analyses, 
# especially the calculation of phylogeny-aware distances between microbial communities. 
# The DADA2 sequence inference method is reference-free, so we must construct the phylogenetic tree
# relating the inferred sequence variants de novo. We begin by performing a multiple-alignment using 
# the DECIPHER R package. This section is based on an online workflow.


# Read in data -------------------------------------------------------------------------------------------------------------
# seqtab <- readRDS("./Objects/seqtab_otu99.rds")
# tax <- readRDS("./Objects/taxtabspecies_otu99.rds")

# Sequence table of ASVs clustered into OTUs
# Rows are sequences of merged OTUs
seqtab <- readRDS("./Objects/OTU_99_gtdb_2015-18_table_w.seq.rds")
tax <- readRDS("./Objects/OTU_99_gtdb_2015-18_taxonomy.rds")

seqs <- getSequences(seqtab)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
saveRDS(alignment,"./Objects/alignment_OTU_2015-18.rds")

## Find clusters of ASVs to form the new OTUs
# aln <- DECIPHER::AlignSeqs(dna, processors = NULL)
# saveRDS(aln,"/home/mstadler/scratch/aln_decipher.rds")

# Once the code above has run, there is no need to re-run everytime the script is launched
#alignment <- readRDS("./Objects/OTU_2015-18/alignment_OTU_2015-18.rds")

# The phangorn R package is then used to construct a phylogenetic tree. 
# Here we first construct a neighbour-joining tree, and then fit a GTR+G+I 
# (Generalized time-reversible with Gamma rate variation) maximum likelihood tree. 
# We are using the neighbour-joining tree as a starting point in order to obtain a 
# good-fitting model. A rough tree optimises later estimations. NJ trees are based on 
# distance and thus faster. Use of the as(alignment,"matrix") option means that an
# error will be thrown if sequences are of unequal length (which means they may not be aligned).

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align) # compute pairwise distances from sequences
treeNJ <- NJ(dm) # neighbor-joining tree estimation, un-rooted tree
# Note, tip order != sequence order

# Plot unrooted tree
#plot(treeNJ, "unrooted", cex = 0.5, show.tip.label = FALSE)

# Which model is the best? -------------------------------------------------------------------------------
# To identify the best model, out of all the options there are for constructing phylogenetic trees,
# we will use the modeltest function. This function can take a long time, thus we will use the 
# multicore option if several processors are available.

phy.modeltest <- modelTest(phang.align, tree = treeNJ, 
                           model = c("JC", "F81", "K80", "HKY", "SYM", "GTR"), G = TRUE, I = TRUE, 
                           k = 4, control = pml.control(epsilon = 1e-08, maxit = 3, trace = 1), multicore = TRUE, mc.cores = 24)
saveRDS(phy.modeltest, "./Objects/phy.models_OTU_2015-18.rds")

# To find the best model we will use a custom function that searches for the lowest AIC or BIC.

findBestModel <- function(mt, ic = "AIC") {
  if (ic == "AIC") {
    minic <- min(mt$AIC)
    bestmodel <- mt[mt$AIC == minic, ]
  } else if (ic == "BIC") {
    minic <- min(mt$BIC)
    bestmodel <- mt[mt$BIC == minic, ]
  }
  bestmodel
}

phy.bestmodel <- findBestModel(phy.modeltest)
#phy.bestmodel # best model is GT+G+I

saveRDS(phy.bestmodel, "./Objects/phy.bestmodel_OTU_2015-18.rds")

# The default pml model estimates a Jukes-Cantor model. 
# To change the model to a GTR+G+I model, we update certain arguments with the update function. 
# And then we optimise some model parameters with optim.pml. For larger trees the NNI (nearest-neighbor interchanges) 
# rearrangements often get stuck in local maxima. rearrangement="stochastic" performs stochastic rearrangements
# similar as, which makes random NNI permuation to the tree, which than gets optimised to escape local optima. 
# This algorithm may find better trees but they will also take more time.
# 
# fit <- pml(treeNJ, data=phang.align) # computes the likelihood of a phylogenetic tree
# # uses sequence alignment and a model
# ## negative edges length changed to 0!
# saveRDS(fit, "./Objects/phy.fit_OTU_2015-18.rds")
# 
# fitGTR <- update(fit, k=4, inv=0.2)
# # k = Number of intervals of the discrete gamma distribution
# # inv = proportion of invariable sites
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
# 
# rerooted <- midpoint(fitGTR$tree)
# saveRDS(fitGTR, "./Objects/phy.tree_OTU_2015-18.rds")
# saveRDS(rerooted, "./Objects/phy.tree.rerooted_OTU_2015-18.rds")
# detach("package:phangorn", unload=TRUE)
# 
# fitGTR <- readRDS("./Objects/phy.tree.rds")
# rerooted <- readRDS("./Objects/phy.tree.rerooted.rds")
# plot(fitGTR, cex = 0.5, show.tip.label = FALSE)
# 
# plot(rerooted, cex = 0.5, show.tip.label = FALSE)
# 
# 
# # concstruct phyloseq object from dada2 outputs
# ps <- phyloseq(otu_table(seqtab, taxa_are_rows = F),
#                sample_data(meta),
#                tax_table(tax),
#                phy_tree(phytree))