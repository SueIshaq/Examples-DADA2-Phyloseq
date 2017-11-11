#######------------------------------
## This script is for checking through negative controls using phyloseq, but not including diversity analysis.
## I compiled this from phyloseq tutorials, lots of online comments, and creating some parts myself.
## This is my first time using phyloseq, and I thought other people learning it would appreciate another example.

## These papers are great references for the problem with batch effects in sequencing data: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-014-0087-z, http://onlinelibrary.wiley.com/doi/10.1111/nph.12923/abstract
###### -----------------------------

## My dataset contains ~600 samples and ~ 80 negative controls, from various georaphic locations, times, and treatments, as well as different sequencing kits and DNA extraction batches.  They are also spread across 2 sequence runs.  Thus, there are a lot of opportunities to introduce technical error or biological contamination that could be conflated with treatment effects.  Importantly, each of my actual samples is somewhat unique, therefore I don't want to remove things on the basis of abundance alone.  It's an active project, so I am not sharing the data yet.

## I have three types of negative controls: kit controls where a blank extraction was run for each batch of extractions, blank DNA sewabs which were extracted and run, and PCR negative controls.  But, I wasn't sure how I wanted to deal with each type of control.  I came up with three approaches, and then looked at unweighted, non-rarified PCoA to watch how my axes changed. 


# I used DADA2 to quality trim my raw data and pick Sequence Variants (SVs), and then passed it to phyloseqas in the DADA2 tutorial.
 
# Load the sequence table from DADA2 as needed.  Can also be a sequence tables where columns are SVs/OTUs/taxa, and rows are samples.
	seqtab.all.nochim	<- readRDS('/directory/seqtab-all-nochim.rds') 

# Load the taxonomy from DADA2 as needed.
	all.taxa <- readRDS('/directory/taxa-silva-all-contigs-nochim.rds')

# Load metadata table as needed. Columns are SVs, and rows are samples.
	metadata <- read.csv("metadata_all.csv", header=T, row.names=1)

## create a phyloseq object with all samples
	ps_ex <- phyloseq(otu_table(seqtab.all.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(all.taxa))

# Check the object to verify the number of samples and SVs
	ps_ex

## Plot the taxa sums to see how abundant they all are.  I limited the y-axis to better see the long tail of rare taxa.
	plot(sort(taxa_sums(ps_ex), TRUE), type="h", ylim=c(0, 20))

# See if obscure taxon (choose one) is prevalent across real samples or just controls. Plot by group, using the Treatment column in my metadata file.
	ps_ex_Chro = subset_taxa(ps_ex, Genus=="Chroococcidiopsis")
	plot_bar(ps_ex_Chro, "Treatment", fill = "Genus", title="Chroococcidiopsis across samples")

# See if common taxon is prevalent across real samples or just controls.
	ps_ex_Staph = subset_taxa(ps_ex, Genus=="Staphylococcus")
	plot_bar(ps_ex_Staph, "Treatment", fill = "Genus", title="Staphylococcus across samples")


# What does my ordination look like?
	C1.ord <- ordinate(ps_ex, method ="PCoA", "jaccard")
	plot_ordination(ps_ex, C1.ord, type="samples", color="Treatment")
	
	# For mine, I see definite clustering, but depending on which variable I use to color the graph, I can get cluseters by Treatment, dna_extraction_batch, dna_extraction_kit, etc. Because samples were not completely randomized by DNA extraction batch or sequencing run, it's not clear if there is a batch effect or just a treatment effect.  Either way, that potential batch effect is driving my primary axis.



###### -----------------------------
## Look negative controls to give a better idea of total seqs and total taxa present.
###### -----------------------------
# Subset out just the negative controls. My negative controls are all labelled as such in the Treatment column of my metadata.
	ps_ex_controls = subset_samples(ps_ex, Treatment == "NegCon") 

# Clean out taxa/SV columns that are no longer present.
	ps_ex_controls <- prune_taxa(taxa_sums(ps_ex_PCR_controls) > 0, ps_ex_PCR_controls)
	
# How abundant are the SVs in the negative controls?	
	plot(sort(taxa_sums(ps_ex_PCR_controls), TRUE), type="h")
	
# What's in them?	
	plot_bar(ps_ex_PCR_controls, "Sample_ID", fill = "Genus")
	
# Find taxa that are shared among a certain number of negative controls, here I use 5 samples.  This might indicate if there is a common contaminant that might be coming from a common reagent, often water.
	controls_shared_taxa = filter_taxa(ps_ex_PCR_controls, function(x) sum(x >= 1) == (5), TRUE)
	plot_bar(controls_shared_taxa, "Sample_ID", fill = "Genus")
	

###### -------------------------------
## Removing effects of negative controls.  
###### -------------------------------


########### Approach 1- subtraction + outright removal

## This approach is pertinant for working with OTUs, or situations where you wouldn't want to remove the whole OTU, just substract out a certain number sequences from specific columns. There is currently no way to do that just from phyloseq, so I made a work-around.
## If you used DADA2, and have SVs, you can skip this section because SVs are more exact and I found it's better to remove them.  I ended up going with Approach 3, as I tried this method and wasn't satifyied with my "cleaned" PCoA. 
#####--------------------------------

# First, I need to subset out my runs because my phyloseq object contains two cohorts/sequencing runs.
	cohort1 = subset_samples(ps_ex, Cohort == "1") ## 339 samples


# Then I need to subset by dna-extraction batch. I have one negative kit control for each DNA extraction, so I need to remove those sequences per each extraction batch. 
	cohort1_batch1 = subset_samples(cohort1, dna_extraction_batch == "1") 

# Subset out the negative control for that extraction batch.
	cohort1_batch1_kit = subset_samples(cohort1_batch1, Treatment == "NegCon_kit")
	
# Then caluclate the column sums, which are the SV/OTU sums.  My sample's ID is KC01	
	KC01_sums <- colSums(otu_table(cohort1_batch1_kit))
	
## move the kit SV sums to a vector so it can be subtracted.
	KC01_sums_vec <- as.vector(KC01_sums)
	
## move cohort1_batch1 to a dataframe to subtract things.
	C1B1 = as(otu_table(cohort1_batch1), "matrix")
	C1B1df = as.data.frame(C1B1)
	
## subtract SVs totals from KC01 from other samples using sweep. So, if you have 10 seqs in your first Negative Control SV, it'll subtract 10 from each sample row in the first column.
	C1B1df[,1:56132] <- sweep(C1B1df[,1:56132],2,KC01_sums_vec)
	C1B1df[,1] # check that it worked
	
## replace all negative numbers with 0, since we can't have negative seqs.
	C1B1df <- replace(C1B1df, C1B1df < 0, 0)
	C1B1df[,1] # check that it worked

# Do those 7 steps again for each extraction batch to remove the contamination from the extraction kit negative control.  I had 15 batches, and I over-cautiously ran that chunk of code for each of them, changing the names each time so that I had 15 different dataframes at the end.  But, you could probably just pass a list to a vector and make a loop of it.

## I also have PCR neg controls which didn't have a DNA extraction, but I wanted to be able to knit them back in so I just made them into a dataframe for consistency.
	C1BPCRneg = as(otu_table(cohort1_PCRneg), "matrix")
	C1BPCRnegdf = as.data.frame(C1BPCRneg)


# Now, we can merge all our new OTU tables back together for Cohort 1, and reload as a phyloseq object along with the previous taxa table and metedata, since we haven't removed any samples rows or SVs/taxa columns.  Once in phyloseq, we can remove columns/Svs that are now 0.
	cohort1_kitcleaned <- rbind(C1B1df, C1B2df, C1B3df, ..., C1BPCRnegdf)
	dim(cohort1_kitcleaned) # should rows and columns of match what went in

# Because my sequencing runs also had different extraction batches, I repeated this entire process up to this point for my 2nd cohort/run, and after combining cohort 2 back together, I combined both cohorts.
	cohorts_kitcleaned <- rbind(cohort1_kitcleaned, cohort2_kitcleaned)


# Create a new phyloseq object with all samples again.
	ps_ex_kitcleaned <- phyloseq(otu_table(cohorts_kitcleaned, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(all.taxa))
               
# Check the object to make sure all the samples and the SVs came back.  The columns numbers will be the same, even if they are empty now.
	ps_ex_kitcleaned

# Clean out taxa/SV columns that are no longer present.
	ps_ex_kitcleaned <- prune_taxa(taxa_sums(ps_ex_kitcleaned) > 0, ps_ex_kitcleaned)

# Clean out samples that are now empty (NegCon_kit)
	ps_ex_kitcleaned <- prune_samples(sample_sums(ps_ex_kitcleaned) > 0, ps_ex_kitcleaned)
	
# Plot again to see if that long tail of rare taxa has shortened.	
	plot(sort(sample_sums(ps_ex_kitcleaned), TRUE), type="h")
	

# Next I needed to take care of the PCR negative controls and the swab negative controls.  Since they are specific to cohort/run, I repeated what I did above twice: once for cohort 1 and once for cohort 2.

# Subset by run.  I am starting with cohort/run 1.
	ps_ex_kitcleaned_cohort1 = subset_samples(ps_ex_kitcleaned, Cohort == "1") 
	
# Subset the controls and prune to be only those taxa.
	ps_ex_controls_cohort1 = subset_samples(ps_ex_kitcleaned_cohort1, Treatment == "NegCon_swab" | Treatment == "NegCon_PCR")

# Remove SVs/columns/taxa that are empty.
	ps_ex_controls_cohort1 <- prune_taxa(taxa_sums(ps_ex_controls_cohort1) > 0, ps_ex_controls_cohort1)
	
# First, make the taxa names/OTU names/SVs into a vector so you can remove them. Make one for the controls, and one from the whole of cohort 1.
	control_vec <- as.vector(taxa_names(ps_ex_controls_cohort1))
	cohort1_vec <- as.vector(taxa_names(ps_ex_cohort1))

# Subtract the controls vector from the cohort 1 vector to make a list of SV names to keep.	
	keep_cohort1 <- setdiff(cohort1_vec, control_vec)

# Then, use the keep vector for the prune taxa argument in phyloseq, because it wants the argument to be true (matching).  It will not accept non-matching (!=)
	ps_ex_cohort1_clean <- prune_taxa(keep_cohort1, ps_ex_kitcleaned_cohort1)
	

# Repeat these 6 steps for cohort 2.  


# Did it work?	
	C1.ord <- ordinate(ps_ex_cohort1_clean, method ="PCoA", "jaccard")
	plot_ordination(ps_ex_cohort1_clean, C1.ord, type="samples", color="dna_extraction_batch")

	C2.ord <- ordinate(ps_ex_cohort2_clean, method ="PCoA", "jaccard")
	plot_ordination(ps_ex_cohort2_clean, C2.ord, type="samples", color="dna_extraction_batch")


## Yes, Approach 1 works, and now that primary axis of potential batch effect is now my secondary axis. When I recolor by different variables, there is much more clustering by Treatment than by any batch effects. However, that second axis is also one of my time variables, so don't want to get rid of all of the variation on that axis.  

## If it looks good, go ahead and put the cohorts back together.  Since we stayed in phyloseq for that, we can use phyloseq to recomine those objects.
	ps_ex_totally_clean <- merge_phyloseq(ps_ex_cohort1_clean, ps_ex_cohort1_clean)



###### ----------------------------------
### Approach 2- subtract out all SVs that are found in any of the negative controls. I did this for each cohort, but I did not differentiate by DNA ecxtraction batch this time.

# Subset out by cohort/run.
	ps_ex_cohort1 = subset_samples(ps_ex, Cohort == "1") 
	ps_ex_cohort2 = subset_samples(ps_ex, Cohort == "2") 
	
# Subset the controls and prune to be only those taxa.
	ps_ex_controls_cohort1 = subset_samples(ps_ex_cohort1, Treatment == "NegCon_swab" | Treatment == "NegCon_PCR" | Treatment == "NegCon_kit")
	ps_ex_controls_cohort1 <- prune_taxa(taxa_sums(ps_ex_controls_cohort1) > 0, ps_ex_controls_cohort1)
	
	ps_ex_controls_cohort2 = subset_samples(ps_ex_cohort2, Treatment == "NegCon_swab" | Treatment == "NegCon_PCR" | Treatment == "NegCon_kit")
	ps_ex_controls_cohort2 <- prune_taxa(taxa_sums(ps_ex_controls_cohort2) > 0, ps_ex_controls_cohort2)	
	

 # Make the taxa names into a vector, for the controls and your samples, so you can remove them.
	control1_vec <- as.vector(taxa_names(ps_ex_controls_cohort1))
	cohort1_vec <- as.vector(taxa_names(ps_ex_cohort1))
	keep_cohort1 <- setdiff(cohort1_vec, control1_vec)


	control2_vec <- as.vector(taxa_names(ps_ex_controls_cohort2))
	cohort2_vec <- as.vector(taxa_names(ps_ex_cohort2))
	keep_cohort2 <- setdiff(cohort2_vec, control2_vec)
	
	
# Then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching)
	ps_ex_cohort1_clean <- prune_taxa(keep_cohort1, ps_ex_cohort1)

	ps_ex_cohort2_clean <- prune_taxa(keep_cohort2, ps_ex_cohort2)

# Clean out taxa/SV columns that are no longer present.
	ps_ex_cohort1_clean <- prune_taxa(taxa_sums(ps_ex_cohort1_clean) > 0, ps_ex_cohort1_clean)
	ps_ex_cohort1_clean <- prune_samples(sample_sums(ps_ex_cohort1_clean) > 0, ps_ex_cohort1_clean)

	ps_ex_cohort2_clean <- prune_taxa(taxa_sums(ps_ex_cohort2_clean) > 0, ps_ex_cohort2_clean)
	ps_ex_cohort2_clean <- prune_samples(sample_sums(ps_ex_cohort2_clean) > 0, ps_ex_cohort2_clean)

# Did it work?	
	C1.ord <- ordinate(ps_ex_cohort1_clean, method ="PCoA", "jaccard")
	plot_ordination(ps_ex_cohort1_clean, C1.ord, type="samples", color="dna_extraction_batch")

	C2.ord <- ordinate(ps_ex_cohort2_clean, method ="PCoA", "jaccard")
	plot_ordination(ps_ex_cohort2_clean, C2.ord, type="samples", color="dna_extraction_batch")


## Yes, Approach 2 works, and now that primary axis of potential batch effect is now my secondary axis. When I recolor by different variables, there is much more clustering by Treatment than by any batch effects. However, that second axis is also one of my time variables, so don't want to get rid of all of the variation on that axis. But, since my negative kit controls showed a lot of variation in number and types of taxa, I don't want to remove everything found there from all samples indiscriminatly.  

## If it looks good, go ahead and put the cohorts back together.  Since we stayed in phyloseq for that, we can use phyloseq to recomine those objects.
	ps_ex_totally_clean <- merge_phyloseq(ps_ex_cohort1_clean, ps_ex_cohort1_clean)



### ------------------------------------
### Approach 3- remove PCR and swab contamint SVs fully from each cohort, respectively, and remove kit contaminants fully from each dna_extraction_batch, respectively.


###### Cohort 1-----------
# Subset out by cohort/run.
	ps_ex_cohort1 = subset_samples(ps_ex, Cohort == "1") 
	
# Subset the controls and prune to be only those taxa	
	ps_ex_controls_cohort1 = subset_samples(ps_ex_cohort1, Treatment == "NegCon_swab" | Treatment == "NegCon_PCR")
	ps_ex_controls_cohort1 <- prune_taxa(taxa_sums(ps_ex_controls_cohort1) > 0, ps_ex_controls_cohort1)

# Make the taxa names into a vector so you can remove them. 
	control_vec <- as.vector(taxa_names(ps_ex_controls_cohort1))
	cohort1_vec <- as.vector(taxa_names(ps_ex_cohort1))
	keep_cohort1 <- setdiff(cohort1_vec, control_vec)

# Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
	ps_ex_cohort1_clean <- prune_taxa(keep_cohort1, ps_ex_cohort1)
	
# Then remove the sampels which are now empty, namely the NegCon_PCR and NegCon_swab.
	ps_ex_cohort1_clean <- prune_samples(sample_sums(ps_ex_cohort1_clean) > 0, ps_ex_cohort1_clean)
	

# Subset out each dna_extraction_batch
	cohort1_batch1 = subset_samples(ps_ex_cohort1_clean, dna_extraction_batch == "1") 
	cohort1_batch2 = subset_samples(ps_ex_cohort1_clean, dna_extraction_batch == "2") 
	cohort1_batch3 = subset_samples(ps_ex_cohort1_clean, dna_extraction_batch == "3")
	...


# Subset the controls and prune to be only those taxa.  Some of my kit controls had no taxa after removing PCR contamination, and I got an error message that it was empty. I noted those and skipped the next steps for those batches.
	cohort1_batch1_kit = subset_samples(cohort1_batch1, Treatment == "NegCon_kit")
	cohort1_batch2_kit = subset_samples(cohort1_batch2, Treatment == "NegCon_kit")
	cohort1_batch3_kit = subset_samples(cohort1_batch3, Treatment == "NegCon_kit")
	...

# Prune taxa to remove empty columns.
	cohort1_batch1_kit <- prune_taxa(taxa_sums(cohort1_batch1_kit) > 0, cohort1_batch1_kit)
	cohort1_batch2_kit <- prune_taxa(taxa_sums(cohort1_batch2_kit) > 0, cohort1_batch2_kit)
	cohort1_batch3_kit <- prune_taxa(taxa_sums(cohort1_batch3_kit) > 0, cohort1_batch3_kit)
	...

# Make the taxa names into a vector so you can remove them. 
	b1_control_vec <- as.vector(taxa_names(cohort1_batch1_kit))
	b1_cohort1_vec <- as.vector(taxa_names(cohort1_batch1))
	b1_keep_cohort1 <- setdiff(b1_cohort1_vec, b1_control_vec)
# Then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
	b1_cohort1_clean <- prune_taxa(b1_keep_cohort1, cohort1_batch1)

# Make the taxa names into a vector so you can remove them. 
	b2_control_vec <- as.vector(taxa_names(cohort1_batch2_kit))
	b2_cohort1_vec <- as.vector(taxa_names(cohort1_batch2))
	b2_keep_cohort1 <- setdiff(b2_cohort1_vec, b2_control_vec)
# Then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
	b2_cohort1_clean <- prune_taxa(b2_keep_cohort1, cohort1_batch2)

# Make the taxa names into a vector so you can remove them. 
	b3_control_vec <- as.vector(taxa_names(cohort1_batch3_kit))
	b3_cohort1_vec <- as.vector(taxa_names(cohort1_batch3))
	b3_keep_cohort1 <- setdiff(b3_cohort1_vec, b3_control_vec)
# Then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
	b3_cohort1_clean <- prune_taxa(b3_keep_cohort1, cohort1_batch3)

	...


# Merge the phyloseq objects back together, then remove any blank taxa or samples.  For any batches that DIDN'T have a kit control removed because it was empty, I simply put back in the unmodified subset.
	ps_ex_cohort1_clean <- merge_phyloseq(b1_cohort1_clean, b2_cohort1_clean, b3_cohort1_clean, b4_cohort1_clean, cohort1_batch5, cohort1_batch6, b7_cohort1_clean, cohort1_batch8, b9_cohort1_clean, b10_cohort1_clean, b11_cohort1_clean, b12_cohort1_clean, b13_cohort1_clean, b14_cohort1_clean, b15_cohort1_clean)
	
	
# Clean out taxa/SV columns that are no longer present.
	ps_ex_cohort1_clean <- prune_taxa(taxa_sums(ps_ex_cohort1_clean) > 0, ps_ex_cohort1_clean)
	ps_ex_cohort1_clean <- prune_samples(sample_sums(ps_ex_cohort1_clean) > 0, ps_ex_cohort1_clean)


# Repeat all of the Approach 3 steps for cohort 2.  


# Did it work?	
	C1.ord <- ordinate(ps_ex_cohort1_clean, method ="PCoA", "jaccard")
	plot_ordination(ps_ex_cohort1_clean, C1.ord, type="samples", color="dna_extraction_batch")

	C2.ord <- ordinate(ps_ex_cohort2_clean, method ="PCoA", "jaccard")
	plot_ordination(ps_ex_cohort2_clean, C2.ord, type="samples", color="dna_extraction_batch")


## Yes, Approach 3 works, and now that primary axis of potential batch effect is my secondary axis. When I recolor by different variables, there is much more clustering by Treatment and almost none by batch effects. I say almost none, because some of my DNA extraction batches also happen to be Treatment batches, as they represent a subset of samples from a different location.  Thus, I can't tell if those samples cluster separately solely because of location or also because of batch effect.  

## If it looks good, go ahead and put the cohorts back together.  We can use phyloseq to recomine those objects.
	ps_ex_totally_clean <- merge_phyloseq(ps_ex_cohort1_clean, ps_ex_cohort1_clean)

## Now I can FINALLY start analyzing my actual samples!