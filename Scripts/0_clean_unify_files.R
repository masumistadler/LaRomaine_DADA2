library(tidyverse)
library(data.table)
library(filesstrings)
library(stringr)

# Unify data files ---------------------------------------------------------------------------------------------------------
# It seems like there are duplicates (created by Genome Quebec), especially in 2015.
# There are more files than we have submitted to GQ.

# Let's unify the naming strategy across all years
# make a data frame with each file name as a row across all years

df <- 
  data.frame(stringsAsFactors = F)

folders <- list.dirs("./Raw")
folders <- folders[grep("LaRomaine", folders)]

for(i in 1:length(folders)){
  # extract all file names
  files <- list.files(path = folders[i], pattern = "*fastq.gz")
  if(length(files) != 0L){
    first.split <- strsplit(files, split = "[.]") # split to get sample name
    ele.no <- sapply(first.split, length) # some misc. files from GQ are inside
    unit <- file.names <- paste(sapply(first.split[which(ele.no >= 5)],"[[",2),
                                sapply(first.split[which(ele.no >= 5)],"[[",4), sep = "_")
    file.names <- sapply(strsplit(sapply(first.split[which(ele.no >= 5)],"[[",5), split = "_"), "[[", 1)
    if(i == 7 & any(grep("Hy", sapply(strsplit(sapply(first.split[which(ele.no >= 5)],"[[",5), split = "_"), "[[", 2)))){
      w.hy <- grep("Hy", sapply(strsplit(sapply(first.split[which(ele.no >= 5)],"[[",5), split = "_"), "[[", 2))
      file.names[w.hy] <- paste(file.names[w.hy],
                                sapply(strsplit(sapply(first.split[which(ele.no >= 5)],"[[",5), split = "_"), "[[", 2)[w.hy],
                                sep = "_")
    }
    
    # extract some metadata
    submission <- sapply(strsplit(folders[i], "_"), "[[", 3)
    year <- sapply(strsplit(folders[i], "_"), "[[", 2)
    #combine to data frame
    tmp <- data.frame(folder.path = folders[i],
                      year = year, submission = submission,
                      files = files[which(ele.no >= 5)], unit = unit, sample_name = file.names)
    df <- rbind(df, tmp)
  }
}
setDT(df)
df[, merge.id := paste(sample_name, year, submission, sep = "_")]

# read in plate IDs
plate <- read.csv("./Meta/laromaine_nanuq.csv", sep = ",")
setDT(plate)
plate[,sample_name := str_replace(sample_name, "[.]", "-")]
plate[year == 2018,sample_name := str_replace(sample_name, "_replacement", "")]
#df[!(df$sample_name %in% plate$sample_name),]
plate[, merge.id := paste(sample_name, year, submission, sep = "_")]
df[plate, plate.id := i.plate, on = .(merge.id)]

# check if we are dealing with unique rows
df %>%
  distinct() %>%
  nrow() == nrow(df)
# yes

# Check for duplicates
df[duplicated(df$sample_name),] # there are a bunch

# Check if there are duplicates across years?
df[duplicated(df$sample_name),] %>%
  group_by(sample_name) %>%
  mutate(no.years = length(levels(factor(year)))) %>%
  filter(no.years > 1)
# there are a few that overlap between 2015 and 2017

# make a repcliate column (later)
# before we add some useful info
# first extract those that are true replicates, not because we have different years
setDT(df)
df[, pe_pairs := as.numeric(factor(paste(df$unit,df$sample_name)))]
df <- df[grep("_R2.fastq.gz",files), pairend := "RVS"]
df <- df[grep("_R1.fastq.gz",files), pairend := "FWD"]

# add sequencing depth info
df[, seq_depth := "Shallow"]
df[str_detect(df$sample_name, "d$"), seq_depth := "Deep"]
# correct typo
df[, new_name := str_replace(df$sample_name, "DSEq", "DSeq")]
df[str_detect(df$new_name, "DSeq$"), seq_depth := "Deep"]

# remove d and DSeq's to match with file names
df[, id_name := str_remove(new_name, "d$")]
df[, id_name := str_remove(id_name, "DSeq$")]

# add DNA type info
df[year == "2015", dna_type := "DNA"]
df[year == "2016", dna_type := "DNA"]
# For 2016, naming strategy depends on plate
df[(plate.id == 6 | plate.id == 7) & str_detect(df$new_name, "R$"), dna_type := "cDNA"]
# For 2017 DNA and RNA are different plates
df[plate.id == 11, dna_type := "DNA"]
df[plate.id == 12, dna_type := "cDNA"]
# For 2018, everything is indicated by the sample name
df[year == "2018" & str_detect(df$new_name, "D$"), dna_type := "DNA"]
df[year == "2018" & str_detect(df$new_name, "R$"), dna_type := "cDNA"]

#df[is.na(dna_type),]

# import actual 2015 sample_names
y2015.id <- read.csv("./Meta/2015_sampleID.csv", sep = ";")
# replace spaces
y2015.id$sample.names <- str_replace_all(y2015.id$sample.names, " ", "")

# merge with 2015 df
setDT(y2015.id)
y2015.id[!(y2015.id$DNA_ID %in% df$id_name),] #LR108 missing
df[year == "2015",]$id_name %in% y2015.id$DNA_ID # all there
  
only.2015 <- df[year == "2015",][y2015.id, ID.2015 := i.sample.names, on = c(id_name = "DNA_ID")]

# merge back with other years
df <- rbind(only.2015, df[year != "2015", ][,ID.2015 := NA])

# merge names into one column
df[, final_name := new_name]
df[!is.na(ID.2015), final_name := ID.2015]

# correct hypolimnion samples
# first, add splitter to samples with only "Hy"
other <- df[plate.id <= 13,]
other[,final_name := str_replace(final_name, "Hy", ".Hypo.")]
# extract samples with "_" separator
hypos <- sapply(strsplit(other$files, "_"), length)
ext <- sapply(strsplit(other$files, "_")[hypos > 3], "[[", 3)
other[hypos > 3,final_name := paste0(final_name,".", ext, ".")]

sub18 <- df[plate.id >=14,]
sub18[,final_name := str_replace(final_name, "_Hy", ".Hypo.")]

# in 2018 we have "true replicates" that were sent to sequencing twice.
# remove the "-2" identifier for duplicate as it will be added later
sub18[,final_name := str_replace(final_name, "-2", "")]

df <- rbind(other, sub18)



# add replicate column
df[, replicate := 1:nrow(.SD), by = .(final_name, seq_depth, pairend, dna_type)]

# remove hyphens to unfiy naming
df[, final_name := str_remove(final_name, "-")]

# unify hypolimnion indication
#df[, final_name := str_replace(final_name, "Hypo", "Hy")]
# 2016 plate 5, "D", does not indicate DNA but "deep" ??
# Not sure, leave it as replicates
#df[plate.id == 5, final_name := str_replace(final_name, "D$", "Hy")]

# one depth is wrong
df[, final_name := str_replace(final_name, "60m", "90m")]

# remove all D and R's and unify naming
df[, final_name := str_remove(final_name, "D$")]
df[, final_name := str_remove(final_name, "R$")]

# add replicate ID
df[replicate == 2 & year == 2018, final_name :=paste0(final_name, "-t2-")]
# "s" for sequencing replicate, vs "t" for true replicate
# some replicates are not "true, meaning that we did not send replicates to be sequenced
# genome quebec gave us more files back for the same sample = "s"
df[replicate == 2 & dna_type != "cDNA" & year != 2018, final_name :=paste0(final_name, "-s2-")]

# add deep Sequence identifier
df[seq_depth == "Deep", final_name := str_remove(final_name, "d$")]
df[seq_depth == "Deep", final_name := str_remove(final_name, "DSeq$")]
df[seq_depth == "Deep", final_name := paste0(final_name, "DSeq")]

# add Dna Type identifier
df[seq_depth == "Shallow" & dna_type == "DNA", final_name := paste0(final_name, "D")]
df[seq_depth == "Shallow" & dna_type == "cDNA", final_name := paste0(final_name, "R")]

# duplicates ?
dups <- df[duplicated(final_name) & pairend == "FWD",]$final_name
head(df[final_name %in% dups,])
# these are actual duplicates, that we sent for sequencing twice

# update replicate column for those
df[final_name %in% dups, replicate := 1:nrow(.SD), by = .(final_name, pairend)]
# add to final sample name, 't' for true replicate
df[final_name %in% dups & replicate == 2, final_name :=paste0(paste0(str_remove(final_name, "D$"), "-t2-"), "D")]

# recheck, if there are duplicates
dups <- df[duplicated(final_name) & pairend == "FWD",]$final_name
head(df[final_name %in% dups,])

# no duplicates!

# correct a few wrong samples
df[final_name == "RO2R52R", final_name := "RO252R"]
df[final_name == "SWR34R", final_name := "SW34R"]
df[final_name == "RO236pR", final_name := "RO236R"]
df[final_name == "RO236pD", final_name := "RO236D"]
df[final_name == "RO230DPR", final_name := "RO230R"]
df[final_name == "RO31.Hypo.D", final_name := "RO301.Hypo.D"]
df[final_name == "RO31D", final_name := "RO301D"]
df[final_name == "RO34D", final_name := "RO304D"]
df[final_name == "RO37D", final_name := "RO307D"]
df[final_name == "L230R", final_name := "L330.C1R"]
df[final_name == "L230R", final_name := "L330.C1R"]
df[final_name == "RO2328D", final_name := "RO328D"]

# remove Ds and Rs to match DNA and RNA
df[dna_type == "DNA", dr_match_name := str_replace(df[dna_type == "DNA",]$final_name, "D$", "")]
df[dna_type == "cDNA", dr_match_name := str_replace(df[dna_type == "cDNA",]$final_name, "R$", "")]

# and remove replicate ID
df[, dr_match_name := sapply(strsplit(df$dr_match_name, "-"), "[[", 1)]
# remove DSeq
df[, dr_match_name := str_remove(dr_match_name, "DSeq")]

# Now that we have a unified, coherent naming strategy, we will copy the raw files into a new folder and
# rename the samples so that each file is unique
# add R1 and R2 identifier into data frame
df[pairend == "FWD", ext := "_R1.fastq.gz"]
df[pairend == "RVS", ext := "_R2.fastq.gz"]

# only 2018 submission 2
#df <- df[year == 2018 & submission == 2,]

# make folder to store
dir.create("./Renamed_2015-2018/")
raw.files <- paste(df$folder.path, df$files, sep = "/")
copy.as <- paste0("./Renamed_2015-2018/", df$final_name, df$ext)

# last sanity check, any duplicates?
copy.as[duplicated(copy.as)] # no

head(data.frame(raw.files, copy.as))
tail(data.frame(raw.files, copy.as))

# lets copy and rename
#file.copy(raw.files, copy.as)

# Check if all files were copied.
nrow(df)
new.files <- list.files("./Renamed_2015-2018/")
length(new.files) == nrow(df) # yes

# save naming meta data
# select necessary columns
df[,dr_match_name := str_replace(dr_match_name, ".Hypo.", ".Hypo")]
df[,dr_match_name := str_replace(dr_match_name, "m.", "m")]

# PCR blanks have NA in dr_match_name
df[is.na(dr_match_name), dr_match_name := "pcrblank18"]

out <- df %>%
  select(folder.path, file.path = files,
         original.name = sample_name,
         ID.2015, pairend,
         final_name, dr_match_name,
         replicate, plate.id,
         year, dna_type, seq_depth,
         ext)



write.table(out, "./Meta/LaRomaine_molecular_naming_unified.csv", sep = ",",
            row.names = F)
# D = DNA
# R = RNA/cDNA
# DSeq = Deeply sequenced
# -s2- : second (file) replicate provided by Genome Quebec (did not send replicates to GQ, they gave us two files)
# -t2- : second TRUE replicate, we did send repelicates to GQ
# Hypo : Hypolimnion, for specific depths, xx m is given (e.g. 90 m)

samples <- read.csv("./Meta/LaRomaine_molecular_naming_unified.csv", sep = ",",
                       stringsAsFactors = F)
sample.names <- unique(samples$dr_match_name)

# Match names to Database names (aka "MasterFile")
library(data.table)
main <- read.csv("./Meta/LaRomaine_MasterFile_withDist_2021-02-22.csv", sep = ",", dec= '.', stringsAsFactors = F)
setDT(main)

# remove all rows with .1 in reservoirs
ignore <- str_detect(main[sample.type != "Lake" & sample.type != "Soilwater",]$sample.name,"[.]")
no.lake <- main[sample.type != "Lake" & sample.type != "Soilwater",]
ignore <- no.lake[ignore,]$sample.name
sub.points <- c(main[!(sample.name %in% ignore),]$sample.name,"LR09.1")
sub.main <- main[sample.name %in% sub.points,]

# replace "hypo" with "Hypo" and manually correct a few samples
sub.main[, dna.match := str_replace(sample.name, "hypo", ".Hypo")]
sub.main[, dna.match := str_replace(dna.match, "[_]", ".")]
sub.main[, dna.match := str_replace(dna.match, "-4m", ".4m")]
sub.main[sample.type == "Lake", dna.match := str_replace(dna.match, "[.]","")]

# remove some white space in a few sample names
sub.main[, dna.match := str_trim(dna.match)]

# remove all "-"
sub.main[, dna.match := str_replace(dna.match, "[-]","")]

# check if there are samples with no match
sample.names[!(sample.names %in% sub.main$dna.match)]

# manually correct a few samples
sub.main[dna.match == "L330C1", dna.match := "L330.C1"]
sub.main[dna.match == "L330C2", dna.match := "L330"]
#sub.main[sample.type == "Lake", dna.match := str_replace(sub.main[sample.type == "Lake",]$dna.match, "[.]","")]
sub.main[dna.match == "T49.C2", dna.match := "TR49"]
# could be spring or summer..., no other tributary was samples in spring, take summer
sub.main[dna.match == "T56", dna.match := "TR56"]
sub.main[dna.match == "T57", dna.match := "TR57"]
sub.main[dna.match == "RO2111.Hypo", dna.match := "RO2111.90m"]

sub.main[dna.match == "SWPR02.1", dna.match := "SWPR02"]
sub.main[dna.match == "SPR4", dna.match := "SPR04"]
sub.main[dna.match == "LR09.1", dna.match := "LR09"]
#sub.main[dna.match == "RO328", dna.match := "RO2328"]

sub.main[dna.match == "SN47", dna.match := "Snow47"]
sub.main[dna.match == "Snow47", c("sample.type.year","sample.type") := list("Snow","Snow")]

# there are no entries for S33, S34, S35, S36, S37, S39, S40, S42
# take campaign and year from SW with same number
# 44,45,46,47,48,49,50,51,52,53,54,5556,57,58,59
extra.soil <- sub.main[dna.match %in% c("SW33", "SW34", "SW35", "SW36", "SW37", "SW39", "SW40", "SW42",
                                        "SW44","SW45","SW46","ST47","ST48","SW49","ST50",
                                        "ST51","ST52","SW53","SW54","SW55","ST57","HW58","HW59"),] %>%
  select(sample.name, year, campaign,julian.day, sampling.date, lat, long, distance.from.mouth, catchment.area)
setDT(extra.soil)
extra.soil[, c("sample.type","sample.type.year") := list("Soil","Soil")]
extra.soil[, dna.match := c("S33", "S34", "S35", "S36", "S37", "S39", "S40", "S42",
                            "S47","S48","S44","S45","S46","S50","S51",
                            "S52","S49","S53","S58","S59","S57","S54","S5556")]
# merge with main data

all <- bind_rows(sub.main, extra.soil)

# there are no entries for sediments extract from surface samples
sed <- data.frame(original = c(sample.names[str_detect(sample.names, "SED")],
                               sample.names[str_detect(sample.names, "CM")],
                               "L331S","L332S"))
setDT(sed)
sed[, match := str_replace(original, "SED","")]
sed[str_detect(original, "L330"), match := "L330_C2"]
sed[, match := str_replace(match, "S","")]

# some samples are missing an "L"
sed[, correct := match]
sed[str_detect(match,"LR") != T, correct := paste0("L",str_replace(match, "L",""))]
# extract corresponding meta data
sed.meta <- sub.main[sample.name %in% sed$correct,] %>%
  select(sample.name, year, campaign,julian.day, sampling.date, lat, long, distance.from.mouth, catchment.area)
setDT(sed.meta)
sed.meta[, c("sample.type","sample.type.year") := list("Sediment","Sediment")]
sed <- sed[sed.meta, , on = c("correct" = "sample.name")]
sed[, c("match","correct") := list(NULL,NULL)]
colnames(sed)[1] <- "dna.match"

all <- bind_rows(all, sed)

# there are no entries for bioassays extract from surface samples
bio <- data.frame(original = sample.names[str_detect(sample.names, "BIO")])
setDT(bio)
bio[, match := str_replace(original, "A1BIO", "")]
bio[, match := str_replace(match, "B1BIO", "")]
# extract corresponding meta data
bio.meta <- sub.main[dna.match %in% bio$match,] %>%
  select(sample.name, year, campaign,julian.day, sampling.date, lat, long, distance.from.mouth, catchment.area, dna.match)
setDT(bio.meta)
bio.meta[, c("sample.type","sample.type.year") := list("Bioassay","Bioassay")]
bio <- bio[bio.meta, , on = c("match" = "dna.match")]
bio[, c("match") := list(NULL)]
colnames(bio)[1] <- "dna.match"

all <- bind_rows(all, bio)

# Some miscallaneaous samples
misc <- data.frame(dna.match = c("THW1","LR17blankc2","Blank","pcrblank18","Blank19"),
                   campaign = c(1,2,1,1,1),
                   year = c(2016, 2017, 2016,2018,2018),
                   sample.type.year = c("Unknown","Blank","Blank","Blank","Blank"))

all <- bind_rows(all, misc)

all <- all[all$dna.match %in% sample.names,]
# sanity check
sample.names[!(sample.names %in% all$dna.match)]
nrow(all)
length(sample.names)

# add Season column
all[campaign == 1, Season := "Spring"]
all[campaign == 2, Season := "Summer"]
all[campaign == 3, Season := "Autumn"]

# export
write.table(all, "./Meta/main_bac_match_2021-05-11.csv", sep = ",", row.names = F)
