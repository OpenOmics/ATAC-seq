#ATAC-seq
Running SF's ATAC-seq pipeline on Biowulf:

1. Copy the contents of this folder to a new working directory.

2. Create a folder within the working directory called "fastq" containing the raw fastq files or their symlinks. All files must use the following nomenclature: {samplename}_R1.fastq.gz and {samplename}_R2.fastq.gz where samplename should not contain any "." or "_" symbols. "-" symbols are allowed.

3. Create a text file called "config.py". An example can be found within this folder. It must have at least two rows in the format listed here:
   analysis="{workdirectory}"
   ref="{referencegenome}"
where workdirectory is the full path to the working directory and referencegenome is one of "hg38", "hg19", "mm10, "mm9", and "RheMac10".

4. To run a dry run:
   module load snakemake/5.5.2
   snakemake -npr

5. To run the pipeline, edit the header of submit.sbatch to your email address and then:
    sbatch submit.sbatch



Some notes:
- This is currently set up to treat each ATAC-seq sample independently. There is code from SF to run replicates as a set through Genrich, but it has not be incorporated yet.

Group names may not contain any of the following symbols: "-" "." "_"
