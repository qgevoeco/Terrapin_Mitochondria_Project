# Terrapin Mitochondria Project
Sequence assembly, variant calling, and analyses for diamondback terrapin mitochondrial genome

All work conducted by Samuel Weaver, Tonia Schwartz, and Matthew Wolak at Auburn University

### 1) First, we aligned raw reads to the reference genome using bwa. 

Code to do this is in the `MitochondrialAlignment.sh` script.
These `bams` were used to extract consensus sequences and as input to `genious` to generate `phylip`, `fasta`, and `nexus`-formatted alignments.

We used the GUI interface of `PopArt` to generate Median-Joining Haplotype Networks.

### 2) Following alignment, we then called variants using bcftools. 
The code to do this can be found in the `Mpileup_Variant_calling.sh` script.

This vcf was used to generate `phylip` and `nexus` alignments (`invariant_*.phylip`).

### 3) Once we had a VCF, we estimated summary statistics for all Alabama samples
Code for this is in the the `Mito_Genome_Summaries.R` script.

### 4) Make maps of sample locations and their mitochondrial haplotype
Information about the unique sample locations is found in `wgsSampleLocations.csv` for our own sequence data.

Information about sample locations for the sequence data from Parham et al. (2008 _Biology Letters_ doi:10.1098/rsbl.2007.0599) was extracted from their supporting information and is contained in `Parham2008Locations.csv`.

The `R` code to generate maps is contained in `mapCode.R`
