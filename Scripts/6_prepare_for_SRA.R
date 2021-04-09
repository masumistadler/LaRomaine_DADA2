# 1. R set-up ------------------------------------------------------------------------------
### Packages -------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(readxl)

### Functions -----------------------------------------------------------------------------
source("./Functions/custom_fun.R")

fin.df <- select_newest("./Output", "201520162017_fin_css_otu99_phantomcor_paper1_")
fin.df <- as.matrix(read.csv(
  paste0("./Output/", fin.df),
  sep = ";",
  dec = ".",
  row.names = 1,
  stringsAsFactors = F
))

splitdf <- readRDS("./Objects/splitdf_new.rds")

# Extract splitdf rows, that are represented in this data set
sub <- splitdf[splitdf$final_name %in% row.names(fin.df),]
sub <- sub[sub$final_name != "THW1D", ]
# Some downriver samples were wrongly assigned to "Downriver"
sub[grep("SWLR", sub$final_name),]$sample.type.year <- "Soilwater"
sub[grep("SLR", sub$final_name),]$sample.type.year <- "Soil"

# Save a sub splitdf to upload on Zenodo
write.table(sub, "./Output/dada2_filt_splitdf.csv", sep = ",", dec = ".", row.names = F)

# Read-in SRA file format
sra <- read.csv("./Output/Metagenome.environmental.1.0.tsv", sep = "\t", stringsAsFactors = F, skip = 10)
# must are sample_name, organism, isolation source, collection date, geographic location,
# latitude, longitude
# Unique identifiers such as nucleic_acid_type, replicate, sampling_depth(m)

setDT(sub)
sub[is.na(lat),]
sub[, lat_lon := paste(paste(lat, "N", sep = " "), paste(abs(long), "W", sep = " "), sep = " ")]
sub[, collection_date := sampling.date]
sub[, geo_loc_name := "Canada:Quebec"]
sub[, ecosystem := sample.type.year]
sub$ecosystem <- factor(sub$sample.type.year, levels = c("Soil","Sediment",
                                                         "Soilwater","Hyporheic", 
                                                         "Wellwater","Stream", "Tributary",
                                                         "HeadwaterLakes", "PRLake", "Lake", "IslandLake",
                                                         "Upriver", "RO3", "RO2", "RO1","Deep", "Downriver",
                                                         "Marine"),
                        labels = c("Soil","Sediment",
                                   "Soilwater","Hyporheic", 
                                   "Groundwater","Stream", "Tributary",
                                   "Lake", "Pond", "Lake", "Lake",
                                   "River",# "RO3", "RO2", "RO1","Deep",
                                   "Reservoir","Reservoir", "Reservoir","Reservoir",
                                   "River",
                                   "Estuary"))
sub[, organism := paste(ecosystem, "metagenome", sep = " ")]
sub[, isolation_source := ecosystem]
sub[ecosystem == "Stream" | ecosystem == "Tributary" |
      ecosystem == "Lake" | ecosystem == "Pond" | ecosystem == "River" |
      ecosystem == "Reservoir" | ecosystem == "Estuary", sampling_depth_m := 0.3]

# get sampling depths of hypolimnion
hypo <- read.csv("./MotherData/hypolim_depths.csv"); setDT(hypo)
sub[hypo, sampling_depth_m := i.sampling_depth, on = c("final_name" = "sample_name")]

# rename RNA to cDNA
sub[, nucleic_acid_type := dna_type][nucleic_acid_type == "RNA", nucleic_acid_type := "cDNA"]

# keep DNA-RNA match name as sample title
sub[, sample_title := dr_match_name]
sub[, sample_name := final_name]
# correct replicate number of sediment cores
sub[grep("L330.C1", sample_name), lat_lon := "51.18157 N 63.92587 W"]

# correct replicate number of sediment cores
sub[grep("L330CM", sample_name), replicate := str_sub(sample_title, start = -1)]
# duplicate soilwater
sub[grep("SWLR02", sample_name), replicate := str_sub(sample_title, start = -1)]
sub[grep("SWLR03", sample_name), replicate := str_sub(sample_title, start = -1)]

# Merge with SRA
sra <- read.csv("./Output/Metagenome.environmental.1.0.tsv", sep = "\t", skip = 10, header = T)
colnames(sra)
# must are:
# sample_name
# organism riverine metagenome
# isolation source = habitat type
# collection date "YYYY-MM-DD"
# geographic location
# latitude and longitude #d.dd N\S d.dd W\E e.g. 38.98 N 77.11 W
setDT(sra)
colnames(sra)<-c("sample_name","sample_title","bioproject_accession","organism",
                 "host","isolation_source","collection_date","geo_loc_name",
                 "lat_lon","ref_biomaterial","rel_to_oxygen","samp_collect_device",
                 "samp_mat_process","samp_size","source_material_id","description")
sra.out <- data.frame(sample_name = sub$sample_name,
                      sample_title = sub$sample_title,
                      bioproject_accession = "PRJNA693020",
                      organism = sub$organism,
                      host = NA,
                      isolation_source = sub$isolation_source,
                      collection_date = sub$collection_date,
                      geo_loc_name = sub$geo_loc_name,
                      lat_lon = sub$lat_lon,
                      ref_biomaterial = NA,
                      rel_to_oxygen = NA,
                      samp_collect_device = NA,
                      samp_mat_process = NA,
                      samp_size = NA,
                      source_material_id = NA,
                      description = NA,
                      nucleic_acid_type = sub$nucleic_acid_type,
                      replicate = sub$replicate,
                      sampling_depth_m = sub$sampling_depth_m)

View(sra.out[sra.out$sample_name %in% non.uni,])

write.table(sra.out, "./Output/sra_lrpaper1_metagenome.tsv",  quote=FALSE, sep='\t',
            row.names = F)

# Meta data ------------------------------------------------------------------------------

meta.ov <- read_excel("./Output/SRA_metadata.xlsx", sheet = 2)
# 2 rows per sample = forward and reverse
colnames(meta.ov)

meta <- data.frame(sample_name = sra.out$sample_name,
                   library_ID = sra.out$sample_name,
                   type = sra.out$isolation_source,
                   title = "unique",
                   library_strategy = "AMPLICON",
                   library_source = "METAGENOMIC",
                   library_selection = "PCR",
                   library_layout = "paired",
                   platform = "ILLUMINA",
                   instrument_model = "Illumina MiSeq",
                   design_description = "Conducted by Genome Quebec",
                   filetype = "fastq",
                   filename = "for",
                   filename2 = "rev")

setDT(meta)

meta[, title := paste("16S rRNA (V4) amplicon sequencing of", type, "bacterioplankton community", sep = " ")]
meta[, filename := paste0(sample_name, "_R1.fastq.gz")]
meta[, filename2 := paste0(sample_name, "_R2.fastq.gz")]

meta[, type := NULL]

write.table(meta, "./Output/sra_lrpaper1_metadata.tsv",  quote=FALSE, sep='\t',
            row.names = F)

# Copy files to be uploaded into a new folder
meta <- read.csv("./Output/sra_lrpaper1_metadata.tsv", sep = "\t", header = T, stringsAsFactors = F)

# path to original files
path <- "./Renamed_2015-2018/"

all.fastq <- list.files(path,pattern = ".fastq.gz")
# get those in meta
sub <- all.fastq[all.fastq %in% c(meta$filename, meta$filename2)]
# sanity check
length(sub) == (length(meta$filename) + length(meta$filename2)) # TRUE

# copy files
copy.path <- "./SRA_submission/SUB8904525/"

file.copy(paste0(path, sub), copy.path)

# Upload to SRA
# $ ftp -i
# $ open ftp-private.ncbi.nlm.nih.gov
# enter user
# enter password
# $ cd user_folder
# $ mkdir submissionfolder
# $ cd submissionfolder
# $ mput *.fastq.gz
