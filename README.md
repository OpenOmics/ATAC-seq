#ATAC-seq
Running SF's ATAC-seq pipeline on Biowulf:

1. Copy the contents of this folder to a new working directory.

2. Create a folder within the working directory called "fastq" containing the raw fastq files or their symlinks. All files must use the following nomenclature: {samplename}\_R1.fastq.gz and {samplename}\_R2.fastq.gz where "samplename" should not contain any "." symbols. "\_" and "-" symbols are allowed.

3. Create a text file called "config.json". An example can be found within this folder. 
"ref" must be defined, but the other variables are only for differential analysis at this time.
Options are: "hg38", "hg19", "mm10, "mm9", and "rheMac10".

i.e.
   
{ "ref": "hg38",
  "projectID": "NCBR-test",
}

4. To run a dry run:  
   module load snakemake/5.5.2  
   snakemake -npr  

5. To run the pipeline, edit the header of submit.sbatch to your email address and then:  
    sbatch submit.sbatch  

For Differential Analyses (once the first half is complete):

1. Update config.json to include three more variables: projectID, groups, and contrasts following the example here.
Group names may not contain any of the following symbols: "-", ".", "_". 
For groups, each sample should be listed using "samplename" as defined in step 2 above. Samples can be placed in more than one group. Contrasts must each be given a unique key for organization in the json (no matter how many contrasts there are), but they can be anything you like and will not be incorporated into the analysis in any way.

i.e.
   
{ "ref": "hg38",
  "projectID": "NCBR-test",
  "groups": { "EV": ["EVp1", "EVp3", "EVp5"],
              "OE": ["OEp2", "OEp4", "OEp6"]
            },
  "contrasts": { "1": ["OE", "EV"]
               }
}

2. To run a dry run:  
   module load snakemake/5.5.2  
   snakemake -npr -s Snakefile.Differential  

3. To run the pipeline, edit the header of submit.sbatch to your email address and then:  
    sbatch submit_Differential.sbatch  
