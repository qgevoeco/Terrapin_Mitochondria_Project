#Let's try variant calling with mpileup
#Don't need to submit anything, just run it as an interactive job
srun -N1 -n4 --pty /bin/bash

FinalAlignmentDirectory="/scratch/sbw0033/Mitochondria_Project/Data/Alignments/FinalAlignments"
ReferenceGenome="/scratch/sbw0033/Mitochondria_Project/Data/Terrapin_Mitochondrial_Reference.fasta"
BamList="/scratch/sbw0033/Mitochondria_Project/Data/bamlist.txt"

#Load modules using specific version to ensure future compatibility
module load bwa/0.7.17
module load samtools/1.11
module load picard/2.23.9
module load java/15.0.1
module load bcftools
module load vcftools

cd "${FinalAlignmentDirectory}"

#Make invariant vcf
bcftools mpileup -f $ReferenceGenome -b $BamList --max-depth 5000 -P 1.1e-5 | bcftools call -c --ploidy 1 -O v -o Invariant_Raw_Mito_Variants.vcf
#Make variant vcf
bcftools mpileup -f $ReferenceGenome -b $BamList --max-depth 5000 -P 1.1e-5 | bcftools call -c --variants-only --ploidy 1 -O v -o SNPs_Only_Raw_Mito_Variants.vcf

#Transfer to my computer
#scp -r sbw0033@easley.auburn.edu:/scratch/sbw0033/Mitochondria_Project/Data/Alignments/FinalAlignments/Invariant_Raw_Mito_Variants.vcf /Users/samweaver/Docs/TerrapinRProject/Data
#scp -r sbw0033@easley.auburn.edu:/scratch/sbw0033/Mitochondria_Project/Data/Alignments/FinalAlignments/SNPs_Only_Raw_Mito_Variants.vcf /Users/samweaver/Docs/TerrapinRProject/Data

#Summarize depth of coverage, site quality scores, and missingness across all SNPs
#vcftools --vcf /scratch/sbw0033/Mitochondria_Project/Data/Alignments/FinalAlignments/SNPs_Only_Raw_Mito_Variants.vcf --freq2 --out AlleleFrequency.txt --max-alleles 2
#vcftools --vcf /scratch/sbw0033/TerrapinGenomics/Data/SNPs_Only_Raw_Mito_Variants.vcf --depth --out CoveragePerIndividual.txt --max-alleles 2
#vcftools --vcf /scratch/sbw0033/TerrapinGenomics/Data/SNPs_Only_Raw_Mito_Variants.vcf --site-mean-depth --out CoveragePerSite.txt --max-alleles 2
#vcftools --vcf /scratch/sbw0033/TerrapinGenomics/Data/SNPs_Only_Raw_Mito_Variants.vcf --site-quality --out SiteQuality.txt --max-alleles 2
#vcftools --vcf /scratch/sbw0033/TerrapinGenomics/Data/SNPs_Only_Raw_Mito_Variants.vcf --missing-indiv --out Indiv_Missing.txt --max-alleles 2
#vcftools --vcf /scratch/sbw0033/TerrapinGenomics/Data/SNPs_Only_Raw_Mito_Variants.vcf --missing-site --out Site_Missing.txt --max-alleles 2

#Use VCFtools to summarize a few things across the mitochondrial genome
cd /Users/samweaver/Docs/TerrapinRProject/Data
#First, need to get only diploid sites and sites with no missing data
vcftools --vcf Invariant_Raw_Mito_Variants.vcf --remove-indels --max-alleles 2 --out Filtered_Mito_VCF --recode
