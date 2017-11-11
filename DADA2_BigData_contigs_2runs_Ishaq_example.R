##----------------------------------------------
## This is the workflow I compiled off of the DADA2 tutorial, DADA2 BigData Tutorial, and numerous comments online.
## This particular project used Illumina MiSeq 16S rRNA paired-end reads, from two sequencing runs.
## This is my first time using DADA2, and I thought it would be helpful to others trying to learn to have another full walk-through.
##----------------------------------------------

## install all the tools we'll need, including some we'll use later for diversity analysis with phyloseq.
	library(devtools)
	install_github('benjjneb/dada2')
	source('https://bioconductor.org/biocLite.R')
	install.packages("GenomeInfoDB")
	biocLite("dada2")
	biocLite('ShortRead')
	biocLite('phyloseq')
	biocLite('GenomeInfoDbData') 
	library(ape); packageVersion("ape")
	library(phyloseq); packageVersion("phyloseq")
	library(ggplot2); packageVersion("ggplot2")

## load libraries and set wd
	setwd('/directory/')
	library(dada2); packageVersion('dada2')
	library(ShortRead)
	library(ggplot2)
	library(dplyr)
	library(phyloseq)
	library(GenomeInfoDBData) 
	set.seed(777)

##---------------------------
## Start with file parsing, filtering and trim. Repeat for each cohort (sequencing run), forward and reverse.  
##---------------------------
# Specify the path to the unfiltered sequences, I happen to have them separate by run, and by forward and reverse.
	path1F <- '/directory/Cohort1_forward/raw'
	path1R <- '/directory/Cohort1_reverse/raw'
	path2F <- '/directory/Cohort2_forward/raw'
	path2R <- '/directory/Cohort2_reverse/raw'

# Specify the file path to place filtered files, which I put in a subdirectory.
	filtpath1F <- '/directory/Cohort1_forward/filteredF_cohort1_contigs'
	filtpath1R <- '/directory/Cohort1_reverse/filteredR_cohort1_contigs'
	filtpath2F <- '/directory/Cohort2_forward/filteredF_cohort2_contigs'
	filtpath2R <- '/directory/Cohort2_reverse/filteredR_cohort2_contigs'

# Specify the file designation to look for, this looks for zipped gz files to grab the file names.
	fns1F <- list.files(path1F)
	fastqs1F <- fns1F[grepl('.gz$', fns1F)] # CHANGE if different file extensions

	fns1R <- list.files(path1R)
	fastqs1R <- fns1R[grepl('.gz$', fns1R)] # CHANGE if different file extensions

	fns2F <- list.files(path2F)
	fastqs2F <- fns2F[grepl('.gz$', fns2F)] # CHANGE if different file extensions

	fns2R <- list.files(path2R)
	fastqs2R <- fns2R[grepl('.gz$', fns2R)] # CHANGE if different file extensions

# Fully specify the path, and file names for the reads.  
	fns1F <- file.path(path1F, fastqs1F)
	fns1R <- file.path(path1R, fastqs1R)
	fns2F <- file.path(path2F, fastqs2F)
	fns2R <- file.path(path2R, fastqs2R)

# Pull the sequence names using the forward sequenes and take the sequencing name chunk.
	sample.names1 <- sapply(strsplit(basename(fastqs1F), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
	names(fastqs1F) <- sample.names1 ## 339 raw samples

	sample.names2 <- sapply(strsplit(basename(fastqs2F), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
	names(fastqs2F) <- sample.names2 ## 344 raw samples

# Plot quality profiles, on a randomly selected sample number.
	plotQualityProfile(fns1F[[136]])

# Filter so that only matching sequences pass. fastqPairedFilter wants one sample pair at a time, so you need to create a loop to do all of them.  
# A 250-cylce kit was used, so seq lengths should all be 250 b
# The maxEEparameter sets the maximum number of “expected errors” allowed in a read.  I relaxed this from maxEE=2 to increase the number of reads that are retained downstream.  

# Start with Cohort 1, set up a null variable to catch the reads.in and filtered.reads.out stats.
	filtoutput1 <- NULL; 

# Filter using a loop to reduce memory usage.
for(sam in sample.names1){ # for each sample in the sample.name list
	fastqF <- fns1F[grepl(paste0(sam,'*'),fns1F)] #search for the F file with that name
	fastqR <- fns1R[grepl(paste0(sam,'*'),fns1R)] #search for the R file with that name
	  
# Filter that sample pair and output as new variable names.	  
tmp <- fastqPairedFilter(
	  	c(fastqF, fastqR), 
	  	c(c(file.path(filtpath1F, paste0(sam, "_F_filt.fastq.gz"))), c(file.path(filtpath1R, paste0(sam, "_R_filt.fastq.gz")))), 
	  	trimLeft=c(10, 10), 
	  	truncLen=c(240, 230), 
	  	maxEE=c(2,3), 
	  	maxN= 0, 
	  	rm.phix=TRUE, verbose=TRUE)
	  	filtoutput1 <-rbind(filtoutput1, tmp) #send the output to that variable for each pair
} 

# Save the output to reload as needed.
	saveRDS(filtoutput1, "filtered_cohort1_output_contigs.RDS")


# Do the same with cohort 2, set up a null variable to catch the reads.in and filtered.reads.out stats.
	filtoutput2 <- NULL; 
# Filter using a loop to reduce memory usage.
for(sam in sample.names2){ #for each item in the sample name list
	fastqF <- fns2F[grepl(paste0(sam,'*'),fns2F)] #search for the F file with that name
	fastqR <- fns2R[grepl(paste0(sam,'*'),fns2R)]	#search for the R file with that name
	  
# Filter that sample pair and output as new same names. 
tmp <- fastqPairedFilter(
	  c(fastqF, fastqR), 
	  c(c(file.path(filtpath2F, paste0(sam, "_F_filt.fastq.gz"))), c(file.path(filtpath2R, paste0(sam, "_R_filt.fastq.gz")))), 
	  trimLeft=c(10, 10), 
	  truncLen=c(240, 230), 
	  maxEE=c(2,3), 
	  maxN= 0, 
	  rm.phix=TRUE, verbose=TRUE)
	  filtoutput2 <-rbind(filtoutput2, tmp) 
}

# Save the output to reload as needed.
	saveRDS(filtoutput2, "filtered_cohort2_output_contigs.RDS")


##---------------------------
## Learn error rates, needs to be done for each cohort/run separately, and forward and reverse separately.
## DADA2 pipeline recommends using a variety of quality samples, so I included all of my samples- including the negative controls.
## I also increased the number of reads that were used to calculate errors.
## Based on a number of run-throughs, I found this gave the most realistic number of SVs in my samples, negative controls, and mock community.
##---------------------------
set.seed(1234)

# Set the file names for cohort 1 F.
	fns1F <- list.files(filtpath1F, full.names = TRUE)
	filts1F <- fns1F[grepl('.gz$', fns1F)] # CHANGE if different file extensions
	sample.names1F <- sapply(strsplit(basename(filts1F), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq
	names(filts1F) <- sample.names1F 

# Set the file names for cohort 1 R.
	fns1R <- list.files(filtpath1R, full.names = TRUE)
	filts1R <- fns1R[grepl('.gz$', fns1R)] # CHANGE if different file extensions
	sample.names1R <- sapply(strsplit(basename(filts1R), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq
	names(filts1R) <- sample.names1R 

# Set the file names for cohort 2 F.
	fns2F <- list.files(filtpath2F, full.names = TRUE)
	filts2F <- fns2F[grepl('.gz$', fns2F)] # CHANGE if different file extensions
	sample.names2F <- sapply(strsplit(basename(filts2F), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq
	names(filts2F) <- sample.names2F 

# Set the file names for cohort 2 R.
	fns2R <- list.files(filtpath2R, full.names = TRUE)
	filts2R <- fns2R[grepl('.gz$', fns2R)] # CHANGE if different file extensions
	sample.names2R <- sapply(strsplit(basename(filts2R), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq
	names(filts2R) <- sample.names2R 

# Learn the error rates from cohort 1 and 2, F/R.  When using smaller datasets, DADA2 recommended 1e6 reads, but for the big data workflow, it recommended 2e6.
	err1F <- learnErrors(filts1F, nreads = 2e6, multithread=TRUE, randomize=TRUE)
	err1R <- learnErrors(filts1R, nreads = 2e6, multithread=TRUE, randomize=TRUE)
	err2F <- learnErrors(filts2F, nreads = 2e6, multithread=TRUE, randomize=TRUE)
	err2R <- learnErrors(filts2R, nreads = 2e6, multithread=TRUE, randomize=TRUE)

# Save error profiles to output files as needed.
	saveRDS(err1F, 'error-profile-cohort1F.RDS')
	saveRDS(err1R, 'error-profile-cohort1R.RDS')
	saveRDS(err2F, 'error-profile-cohort2F.RDS')
	saveRDS(err2R, 'error-profile-cohort2R.RDS')

##---------------------------
## Sample inference; picking the sequence variants.  
## Do not pool different sequencing runs as they have their own inherent error rates. 
## Can set the error rate as a matrix that has already been learned, or learn now as selfConsist=TRUE.
##---------------------------
# Dereplication, sample inference, and merger of paired-end reads for Cohort 1, in a loop to save memory. 
	mergers1 <- vector("list", length(sample.names1F))
	names(mergers1) <- sample.names1F
	for(sam in sample.names1F) {
	  cat("Processing:", sam, "\n")
	  	derep1F <- derepFastq(filts1F[[sam]])
	    dada1Fs <- dada(derep1F, err=err1F, multithread=TRUE) # sample inference on the forwards
	    derep1R <- derepFastq(filts1R[[sam]])
	    dada1Rs <- dada(derep1R, err=err1R, multithread=TRUE) # sample inference on the reverse
	    merger1 <- mergePairs(dada1Fs, derep1F, dada1Rs, derep1R) # merge the forwards and reverses
	    mergers1[[sam]] <- merger1
	}

# Dereplication, sample inference, and merger of paired-end reads for Cohort 1, in a loop to save memory. 
	mergers2 <- vector("list", length(sample.names2F))
	names(mergers2) <- sample.names2F
	for(sam in sample.names2F) {
	  cat("Processing:", sam, "\n")
	    derep2F <- derepFastq(filts2F[[sam]])
	    dada2Fs <- dada(derep2F, err=err2F, multithread=TRUE) # sample inference on the forwards
	    derep2R <- derepFastq(filts2R[[sam]])
	    dada2Rs <- dada(derep2R, err=err2R, multithread=TRUE) # sample inference on the reverse
	    merger2 <- mergePairs(dada2Fs, derep2F, dada2Rs, derep2R) # merge the forwards and reverses
	    mergers2[[sam]] <- merger2
	}

##---------------------------
## Make a sequence table and remove chimeras, again per cohort/run.
##---------------------------
# Construct sequence table.
	seqtab1 <- makeSequenceTable(mergers1)
# Remove bimeras.	
	seqtab1.nochim <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=TRUE)  
# Save the output for later.
	saveRDS(seqtab1.nochim, '/directory/seqtab-cohort1-contigs-nochim.rds')
	
# Construct sequence table.
	seqtab2 <- makeSequenceTable(mergers2)
# Remove bimeras.	
	seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE) 
# Save the output for later.
	saveRDS(seqtab2.nochim, '/directory/seqtab-cohort2-contigs-nochim.rds')
	

##---------------------------
## Workflow verification steps
##---------------------------
# Track Cohort 1 seqs through the analysis, by Treatment group (Controls are a treatment).
	getN <- function(x) sum(getUniques(x))
# Load a metadata file with a column "Treatment".
	map_cohort1 <- read.csv('/directory/map_cohort1.csv', header = TRUE, row.names = 1)
# Paste together the filtering output, and the output from removing chimeras.
	track1 <- cbind(filtoutput1, rowSums(seqtab1.nochim), mapall1$Treat)
	colnames(track1) <- c("reads.in","filtered", "nonchimeras", "Treatment")
	rownames(track1) <- sample.names1
	
## Plot by Treatment level and relabel x-axis by Treatment group (alphabetical order).
	plot(track1[,"Treatment"], track1[,"reads.in"], xaxt="n", main="Raw reads for cohort 1") + axis(1, at=1:8, labels=c("Group1", "Group2", ...))
	plot(track1[,"Treatment"], track1[,"filtered"], xaxt="n", main="Filtered reads for cohort 1") + axis(1, at=1:8, labels=c("Group1", "Group2", ...))
	plot(track1[,"Treatment"], track1[,"nonchimeras"], xaxt="n", main="Chimera-free contigs for cohort 1") + axis(1, at=1:8, labels=c("Group1", "Group2", ...))

# Track Cohort 2 seqs through the analysis, by Treatment group (Controls are a treatment).
	getN <- function(x) sum(getUniques(x))
# Load a metadata file with a column "Treatment".
	map_cohort2 <- read.csv('/directory/map_cohort2.csv', header = TRUE, row.names = 1)
# Paste together the filtering output, and the output from removing chimeras.
	track2 <- cbind(filtoutput2, rowSums(seqtab2.nochim), mapall2$Treat)
	colnames(track2) <- c("reads.in","filtered", "nonchimeras", "Treatment")
	rownames(track2) <- sample.names2

## Plot by Treatment level and relabel x-axis by Treatment group (alphabetical order).
	plot(track2[,"Treatment"], track2[,"reads.in"], xaxt="n", main="Raw reads for cohort 2") + axis(1, at=1:11, cex.axis=.7, labels=c("Group1", "Group2", ...))
	plot(track2[,"Treatment"], track2[,"filtered"], xaxt="n", main="Filtered reads for cohort 2") + axis(1, at=1:11, cex.axis=.7, labels=c("Group1", "Group2", ...))
	plot(track2[,"Treatment"], track2[,"nonchimeras"], xaxt="n", main="Chimera-free contigs for cohort 2") + axis(1, at=1:11, cex.axis=.7, labels=c("Group1", "Group2", ...))


# Check the SVs for the mock sample community of 10 species: ZymoBiomics Microbial Community Standard.
	unqs.mock <- seqtab2.nochim["lane1-s344-index-CAACACCT-CATTGACG-PCRpos4-mock",]
	unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
	cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n") 
	
# How populated are the SVs in the control?
	plot(unqs.mock, ylab="Number of seqs in SV, Mock", xlab="Number of SVs, Mock")

# Check that the mock samples classify appropriately
	mock.ref <- getSequences(file.path(path, "ZymoBiomics_16S_mock.fasta"))
	match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
	cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")


##---------------------------
## Merge seq tables here between cohorts/seq runs
##---------------------------
# Merge both cohorts/runs.
	seqtab.all.nochim <- mergeSequenceTables(seqtab1.nochim, seqtab2.nochim, orderBy= "abundance")
# Save the output.
	saveRDS(seqtab.all.nochim, '/directory/seqtab-all-nochim.rds') 

# Assign taxonomy.
	all.taxa <- assignTaxonomy(seqtab.all.nochim, '/directory/silva_nr_v128_train_set.fa.gz', minBoot = 75)
# Save the output.
	saveRDS(all.taxa, '/directory/taxa-silva-all-contigs-nochim.rds') 