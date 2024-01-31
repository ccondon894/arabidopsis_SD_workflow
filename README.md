# arabidopsis_SD_workflow
This Snakemake workflow phases F1 Arabidopsis hybrids using parental short read info, and uses a maximum likelihood framework to estimate segregation distortion effect size and loci along each chromosome.

This workflow includes a conda environment that loads all tools necessary to complete the workflow.

# How to Use:
Update the config.yml file for the parameters 'fastq_folder', 'bam_folder', 'vcf_path', 'vcf_file', species_folder', and 'ref_genome' according to the arabidopsis species and data location. 

Run the Snakemake workflow specifying the number of cores needed:

`snakemake --cores 8 --use-conda`
