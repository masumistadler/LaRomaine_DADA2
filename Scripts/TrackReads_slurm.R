#add this to every script
.libPaths( c( .libPaths(), "/home/mstadler/projects/def-pauldel/R/x86_64-pc-linux-gnu-library/4.0") )
# otherwise dependencies are not found


library("dada2")
library("ShortRead"); packageVersion("ShortRead")
library("data.table")
library("plyr")

#setwd("/home/bioinf/data/Molecular/Masumi/DADA2")
ncores <- detectCores()

# # raw number of reads
# #path <- list("/home/bioinf/data/Bioinf.LaRomaine/Raw/2015/",
# #             "/home/bioinf/data/Bioinf.LaRomaine/Raw/2016/",
# #             "/home/bioinf/data/Bioinf.LaRomaine/Raw/2017/",
# #             "/home/bioinf/data/Bioinf.LaRomaine/Raw/2018/") ## CHANGE ME to the directory containing the fastq files.
# #path.list <- lapply(path, list.files, pattern = "*fastq.gz")
# #path.list <- mapply(paste0, path, path.list)
# #path.list <- unlist(path.list)
# 
# path <- "./Renamed_2015-2018"
# path.list <- list.files(path, pattern = "*fastq.gz")
# path.list <- paste(path, path.list, sep = "/")
# 
# 
# #path <- "/home/bioinf/data/LCare/Raw/"
# #path.list <- list.files(path, pattern = "*fastq.gz")
# #path.list <- paste0(path, path.list)
# 
# qa.sum <- qa(path.list, type = "fastq")
# raw <- data.frame(samples = gsub(".fastq.gz", "", rownames(qa.sum[["readCounts"]])),
#                   qa.sum[["readCounts"]][1],
#                   row.names = NULL, stringsAsFactors = F)
# colnames(raw)[2] <- "input"
# raw$samples <- sapply(strsplit(raw$samples, "_"), "[[",1)
# raw <- raw[duplicated(raw$samples),]
# saveRDS(raw, "./Objects/raw_trackreads.rds")
# 
# # number of reads after cutadapt
# path <- list("./Renamed_2015-2018/cutadapt/forward",
#              "./Renamed_2015-2018/cutadapt/reverse")
# path.list <- lapply(path, list.files, pattern = "*fastq.gz")
# path.list <- mapply(paste, path, path.list, sep = "/", SIMPLIFY = FALSE)
# path.list <- unlist(path.list)
# # path <- list("/home/bioinf/data/Bioinf.LaRomaine/Raw/2015/cutadapt/",
# #              "/home/bioinf/data/Bioinf.LaRomaine/Raw/2016/cutadapt/",
# #              "/home/bioinf/data/Bioinf.LaRomaine/Raw/2017/cutadapt/",
# #              "/home/bioinf/data/Bioinf.LaRomaine/Raw/2018/cutadapt/") ## CHANGE ME to the directory containing the fastq files.
# # path.list <- lapply(path, list.files, pattern = "*fastq.gz")
# # path.list <- mapply(paste0, path, path.list)
# # path.list <- unlist(path.list)
# 
# qa.sum <- qa(path.list, type = "fastq")
# cutadapt <- data.frame(samples = gsub(".fastq.gz", "", rownames(qa.sum[["readCounts"]])),
#                        qa.sum[["readCounts"]][1],
#                        row.names = NULL, stringsAsFactors = F)
# colnames(cutadapt)[2] <- "cutadapt"
# cutadapt$samples <- sapply(strsplit(cutadapt$samples, "_"), "[[",1)
# cutadapt <- cutadapt[duplicated(cutadapt$samples),]
# saveRDS(cutadapt,"./Objects/cutadapt_trackreads.rds")

# number of reads retained after filterAndTrim()
path <- list("./Renamed_2015-2018/cutadapt/forward/filtered",
             "./Renamed_2015-2018/cutadapt/reverse/filtered") ## CHANGE ME to the directory containing the fastq files.
path.list <- lapply(path, list.files, pattern = "*fastq.gz")
path.list <- mapply(paste, path, path.list, sep = "/",SIMPLIFY = FALSE)
path.list <- unlist(path.list)

#LCare
# path <- list("/home/bioinf/data/LCare/Raw/cutadapt/forward/filtered/",
#              "/home/bioinf/data/LCare/Raw/cutadapt/reverse/filtered/") ## CHANGE ME to the directory containing the fastq files.
# path.list <- lapply(path, list.files, pattern = "*fastq.gz")
# path.list <- mapply(paste0, path, path.list, SIMPLIFY = FALSE)
# path.list <- unlist(path.list)

qa.sum <- qa(path.list, type = "fastq")
filtered <- data.frame(samples = gsub(".fastq.gz", "", rownames(qa.sum[["readCounts"]])),
                       qa.sum[["readCounts"]][1],
                       row.names = NULL, stringsAsFactors = F)
colnames(filtered)[2] <- "filt"
filtered$samples <- sapply(strsplit(filtered$samples, "_"), "[[",1)
filtered <- filtered[duplicated(filtered$samples),]
saveRDS(filtered, "./Objects/filtered_trackreads.rds")

print("Done!")
quit(save = "no")

# filtered <- readRDS("./Objects/filtered_trackreads.rds")
# 
# # number of reads retained after dadaF()
# getN <- function(x) sum(getUniques(x))
# 
# #21 only has one sample and doesn't work
# dadaFs <- data.frame()
# for(i in c(1:21,23:27)){
#   pseudo <- readRDS(paste0("./Objects/",i,"_poolFs.rds"))
#   dadaFs <- rbind(dadaFs,data.frame(samples = names(sapply(pseudo, getN)),
#                                     dadaF = sapply(pseudo, getN),
#                                     row.names = NULL, stringsAsFactors = F))
# }
# 
# pos <- which(filtered$samples %in% dadaFs$samples == F)
# name <- filtered[!(filtered$samples %in% dadaFs$samples),]$samples
# pseudo <- readRDS(paste0("./Objects/",22,"_poolFs.rds"))
# 
# dadaFs<- rbind(dadaFs[1:pos-1,],data.frame(samples = name, dadaF = getN(pseudo)),
#                dadaFs[pos:nrow(dadaFs),])
# 
# saveRDS(dadaFs, "./Objects/dadaF_trackreads.rds")
# 
# # number of reads retained after dadaR()
# dadaRs <- data.frame()
# for(i in c(1:21,23:27)){
#   pseudo <- readRDS(paste0("./Objects/",i,"_poolRs.rds"))
#   dadaRs <- rbind(dadaRs,data.frame(samples = names(sapply(pseudo, getN)),
#                                     dadaR = sapply(pseudo, getN),
#                                     row.names = NULL, stringsAsFactors = F))
# }
# 
# pos <- which(filtered$samples %in% dadaRs$samples == F)
# name <- filtered[!(filtered$samples %in% dadaRs$samples),]$samples
# pseudo <- readRDS(paste0("./Objects/",22,"_poolRs.rds"))
# 
# dadaRs<- rbind(dadaRs[1:pos-1,],data.frame(samples = name, dadaR = getN(pseudo)),
#                dadaRs[pos:nrow(dadaRs),])
# 
# saveRDS(dadaRs,"./Objects/dadaR_trackreads.rds")
# 
# # number of reads retained after merging
# mergers <- data.frame()
# for(i in c(1:21,23:27)){
#   mer <- readRDS(paste0("./Objects/",i,"_mergers.rds"))
#   mergers <- rbind(mergers,data.frame(samples = names(sapply(mer, getN)),
#                                       merged = sapply(mer, getN),
#                                       row.names = NULL, stringsAsFactors = F))
# }
# 
# pos <- which(filtered$samples %in% mergers$samples == F)
# name <- filtered[!(filtered$samples %in% mergers$samples),]$samples
# pseudo <- readRDS(paste0("./Objects/",22,"_mergers.rds"))
# 
# mergers<- rbind(mergers[1:pos-1,],data.frame(samples = name, merged = getN(pseudo)),
#                 mergers[pos:nrow(mergers),])
# 
# saveRDS(mergers,"./Objects/merged_trackreads.rds")
# 
# # number of reads retained after chimera removal
# 
# nochim <- data.frame(samples = rownames(readRDS("./Objects/nochim_seqtab_2015-2018.rds")),
#                      noChim =rowSums(readRDS("./Objects/nochim_seqtab_2015-2018.rds")),
#                      stringsAsFactors = F)
# 
# raw <- readRDS("./Objects/raw_trackreads.rds")
# cutadapt <- readRDS("./Objects/cutadapt_trackreads.rds")
# filtered <- readRDS("./Objects/filtered_trackreads.rds")
# dadaFs <- readRDS("./Objects/dadaF_trackreads.rds")
# dadaRs <- readRDS("./Objects/dadaR_trackreads.rds")
# mergers <- readRDS("./Objects/merged_trackreads.rds")
# 
# track <- join_all(list(raw, cutadapt, filtered, dadaFs, dadaRs, mergers, nochim), by = 'samples', type = 'full')
# write.table(track,"./Objects/final_trackreads.csv", row.names = FALSE, sep = ",", dec = ".")
# 
# # calculate percentage retention
# setDT(track)
# 
# track[, c("filt.perc", "merg.perc","nochim.perc") := 
#         list(filt *100 / input,
#              merged *100/input,
#              noChim * 100 / input)]
# range(track$raw.nochim.perc)
# 
# View(track[raw.nochim.perc < 60,])
# 
# splitdf <- readRDS("./Objects/splitdf_new.rds")
# setDT(splitdf)
# 
# # problems at filt step, and merge
# splitdf[track, "perc.kept" := i.raw.nochim.perc, on = c("final_name" = "samples")]
# 
# View(splitdf[perc.kept < 60,])
