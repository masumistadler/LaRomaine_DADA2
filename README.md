# DADA2 pipeline used for the La Romaine project

This repository stores the scripts and files used to process the 16S rRNA sequencing reads of the La Romaine project. The project is part of the Industrial Research Chair in Carbon Biogeochemistry in Boreal Aquatic systems (CarBBAS Chair) led by Paul A. del Giorgio.

## Workflow overview

Samples were sequenced with a Illumina (MiSeq and/or HiSeq) paired-end sequencing platform at Genome Quebec.
In the given workflow, we process reads sequenced for the V4 region of the 16S rRNA.

In brief, the workflow was largely inspired by the [official DADA2 tutorials](https://benjjneb.github.io/dada2/tutorial.html). Our pipeline is mainly designed for big-data, however, it is similarly applicable to smaller datasets. 

We start by renaming sample file names, to simplify the names returned by our sequencing service. Followed by removing the primers using `cutadapt`.
After quality checking, the samples are run through the `dada` step on a plate basis with 'pseudo-pooling' enabled to ensure the retention of singleton reads within samples. Pooling only removes singletons within sampling campaigns. We cluster the dervied ASVs by a 99% similarity threshold into Operational taxonomic units (OTUs), to improve DNA-RNA matching. The level of resolution of ASVs are believed to be too high to be able to match potential DNA and RNA reads of the same taxon.

Samples in 2015 showed lower quality and hence, many reads are being lost through the pipeline. Hence, `maxEE` was increased to `maxEE = c(5,5)` in the `filterAndTrim` step. This improved the read retention for most samples.

---

The original workflow was designed for samples that were split in several folders by years (e.g. 2015, 2016, 2017) and thus makes use of lists frequently. This is now deprecated and all samples are within a single folder. The list `apply` version is still inside the scripts, but commented out with `#`.

Steps that required higher computing power and did not successfully finish on our lab computer with 60 GB RAM and 12 cores were run on Compute Candada's Cedar supercomputer. The equivalent scripts are provided and indicated as `*_slurm.R` and the corresponding `.sh` files are available as well.

---

Raw sequences can be found on SRA under the Bioproject number: PRJNA693020.
Intermediate processing files are stored on Zenodo: 10.5281/zenodo.4611421.

Files are being uploaded as manuscripts are published.

Currently available files:

- 2015-2017: 16S rRNA gene and transcripts (DNA and cDNA) in spring, summer, autumn (shallow sequencing)
   - Part of the manuscript: Stadler M, del Giorgio PA. Terrestrial connectivity , upstream aquatic history and seasonality shape bacterial community assembly within a large boreal aquatic network. ISME J 2022; 16: 937–947. 
   - All files part of this paper are named with: '_paper1_'

---

Final scripts used:

- Scripts/0_clean_unify_files.R
- Scripts/1_remove_primers_cutadapt.R
- Scripts/2_quality_filtering.R OR 2_quality_filt_slurm.R
- Scripts/3_learnError_dada.R OR 3_learnError_dada_slurm.R
- Scripts/4_chimrm_assigntax.R
- Scripts/5_clusterOTU_slurm.R
- Scripts/6_prepare_for_SRA.R

Final data files of all La Romaine data run together:

- Cleaned sequences, assigned ASVs (run through DADA2): Objects/nochim_seqtab_2015-2018.rds
- Assigned taxonomy using GTDB database v.95: Objects/taxtav_gtdb_r95_2015-18.rds
- Distance matrix calculated as part of 99% OTU clustering: Objects/distmat_decipher_2015-18.rds
- Sequences aligned for 99% OTU clustering: Objects/aln_decipher_2015-18.rds

---

## Attribution:

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.


Please cite the Zenodo DOI if you use either the data or the scripts on this Github repository.
