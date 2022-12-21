# Terrapin Mitochondria Project
Sequence assembly, variant calling, and analyses for diamondback terrapin mitochondrial genome
All work conducted with Tonia Schwartz and Matthew Wolak at Auburn University

### 1) First, we aligned raw reads to the reference genome using bwa. 

Code to do this is in the Mitochondrial alignment.sh script.
These bams were used to extract consensus sequences and as input to genious to generate phylip, fasta, and nexus-formatted alignments.
We used the GUI interface of PopArt to generate Median-Joining Haplotype Networks.

### 2) Following alignment, we then called variants using bcftools. 
The code to do this can be found in the Mpileup_Variant_calling.sh script.
This vcf was used to generate phylip and nexus alignments (invariant_*.phylip)

### 3) Once we had a VCF, we estimated summary statistics for all Alabama samples
Code for this is in the the Mito_Genome_Summaries.R script
