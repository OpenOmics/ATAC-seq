import os,re
from os import listdir
import glob
import json

import program
import reference

wdir = os.getcwd()
sourcedir = wdir + "/fastq/"

with open(wdir+ '/cluster.json') as file:
    clusterConfig = json.load(file)

with open(wdir+ '/config.json') as file:
    config = json.load(file)


# this script currently assumes that naming is _R1.fastq.gz
sample = [os.path.basename(file).split('.')[0] for file in glob.glob(sourcedir+'/*')]

samps = []
for item in sample:
        newvar = item.split("_R1")
        othervar = item.split("_R2")
        samps.append(newvar[0])
new = []
for item in samps:
        if '_R2' not in item:
                new.append(item)
samples = new


#if not os.path.isfile(config.analysis + '/metadata.txt'):
rule all:
        input:  expand(wdir + "/QC/{sample}_trimmed_R1_fastqc.html", sample=samples),
                expand(wdir + "/QC/{sample}_R1_fastqc.html", sample=samples),
                expand(wdir + "/bam/{sample}.sorted.markdup.bam.bai", sample=samples),
                expand(wdir + "/Genrich/{sample}.genrich.RPM.bw", sample=samples),
                wdir + "/multiqc.html",
                wdir + "/Deeptools/PCA.png",
                #peak_result,

rule cutadapt:
        input: R1 = sourcedir + "{sample}_R1.fastq.gz", 
               R2 = sourcedir + "{sample}_R2.fastq.gz"
        output: R1 = temp(wdir + "/fastq2/{sample}_trimmed_R1.fastq.gz"), 
                R2 = temp(wdir + "/fastq2/{sample}_trimmed_R2.fastq.gz"), 
        shell: """
module load {program.cutadapt}
cutadapt -j {clusterConfig[cutadapt][threads]} -b file:{program.adapters} -B file:{program.adapters} --trim-n -m 20 -o {output.R1} -p {output.R2} {input.R1} {input.R2}
"""

rule fastqc:
        input: R1 = wdir + "/fastq2/{sample}_trimmed_R1.fastq.gz", 
               R2 = wdir + "/fastq2/{sample}_trimmed_R2.fastq.gz"
        output: forward = wdir + "/QC/{sample}_trimmed_R1_fastqc.html", 
                reverse = wdir + "/QC/{sample}_trimmed_R2_fastqc.html", 
        shell:  """
module load {program.fastQC} 
fastqc -o "QC" --noextract -k 5 -t {clusterConfig[fastqc][threads]} -f fastq {input.R1} {input.R2}
"""

rule fastqc1:
        input: R1 = sourcedir + "{sample}_R1.fastq.gz",
               R2 = sourcedir + "{sample}_R2.fastq.gz"
        output: forward = wdir + "/QC/{sample}_R1_fastqc.html", 
                reverse = wdir + "/QC/{sample}_R2_fastqc.html", 
        shell: """
module load {program.fastQC} 
fastqc -o "QC" --noextract -k 5 -t {clusterConfig[fastqc][threads]} -f fastq {input.R1} {input.R2}
"""

rule seqtk:
        input: R1 = wdir + "/fastq2/{sample}_trimmed_R1.fastq.gz", 
               R2 = wdir + "/fastq2/{sample}_trimmed_R2.fastq.gz"
        output: R1 = temp(wdir + "/fastq2/{sample}_trimmed_R1.sub.fastq"), 
                R2 = temp(wdir + "/fastq2/{sample}_trimmed_R2.sub.fastq")
        shell: """
module load {program.seqtk} 
seqtk sample -s100 {input.R1} 5000000 >{output.R1} && seqtk sample -s100 {input.R2} 5000000 >{output.R2}
"""

rule fastqscreen:
        input: R1 = wdir + "/fastq2/{sample}_trimmed_R1.sub.fastq", 
               R2 = wdir + "/fastq2/{sample}_trimmed_R2.sub.fastq"
        output: one = wdir + "/QC/{sample}_trimmed_R1.sub_screen.png", 
                two = wdir + "/QC/{sample}_trimmed_R2.sub_screen.png", 
        shell: """
module load {program.fastq_screen}
module load {program.bowtie2}
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
cd /lscratch/$SLURM_JOBID
fastq_screen --outdir . --threads {clusterConfig[fastqscreen][threads]} --nohits --conf {program.conf_chip} --subset 0 --aligner bowtie2 {input.R1} {input.R2}
mv /lscratch/$SLURM_JOBID/{wildcards.sample}*png {wdir}/QC
mv /lscratch/$SLURM_JOBID/{wildcards.sample}*txt {wdir}/QC
mv /lscratch/$SLURM_JOBID/{wildcards.sample}*html {wdir}/QC
"""

rule bowtie2:
        input: R1 = wdir + "/fastq2/{sample}_trimmed_R1.fastq.gz", 
               R2 = wdir + "/fastq2/{sample}_trimmed_R2.fastq.gz"
        output: final = temp(wdir + "/bam/{sample}.sorted.bam")
        shell: """
module load {program.samtools}
module load {program.bowtie2} 
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
cd /lscratch/$SLURM_JOBID
bowtie2 -p {clusterConfig[bowtie2][threads]} -k 10 -x {reference.refgen} -1 {input.R1} -2 {input.R2} | samtools view -bS - | samtools sort -m 2G -@ {clusterConfig[bowtie2][threads]} - > {output.final}
"""

rule kraken2:
        input: R1 = wdir + "/fastq2/{sample}_trimmed_R1.sub.fastq",
               R2 = wdir + "/fastq2/{sample}_trimmed_R2.sub.fastq",
        output: result = temp(wdir + "/QC/{sample}.kraken"), 
                report = wdir + "/QC/{sample}.kraken2.report.txt", 
                krona = wdir + "/QC/{sample}.kraken.krona"
        params: prefix = wdir + "/QC/{sample}.kraken"
        shell: """
module load {program.kraken2} 
kraken2 --threads {clusterConfig[kraken2][threads]} --db {program.kraken2db} --output {params.prefix} --report {output.report} --paired {input.R1} {input.R2} && cut -f2,3 {output.result} > {output.krona}
"""

rule krona:
    input: wdir + "/QC/{sample}.kraken.krona"
    output: out = wdir + "/QC/{sample}.kraken.krona.html"
    shell: "module load {program.krona}; krona {input} -c -o {output.out}"

rule sortByRead:
    input: wdir + "/bam/{sample}.sorted.bam"
    output: final = temp(wdir + "/bam/{sample}.sortedByRead.bam"), 
    shell: """
module load {program.samtools}
samtools sort {input} -n -m 2G -@ {clusterConfig[bowtie2][threads]} -o {output.final}
"""

rule index:
    input: wdir + "/bam/{sample}.sorted.bam"
    output: temp(wdir + "/bam/{sample}.sorted.bam.bai"), 
    shell: "module load {program.samtools}; samtools index -@ {clusterConfig[index][threads]} {input}"

rule markdup:
    input: wdir + "/bam/{sample}.sorted.bam"
    output: out = wdir +"/bam/{sample}.sorted.markdup.bam", 
            metric = wdir + "/bam/{sample}_MARKEDUPmetrics.txt", 
    shell: """
module load {program.markedup}
#module load {program.samtools}
if [ ! -e /lscratch/$SLURM_JOBID ]; then mkdir /lscratch/$SLURM_JOBID ;fi
java -Xmx{clusterConfig[markdup][mem]} -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT={input} OUTPUT={output.out} METRICS_FILE={output.metric} ASSUME_SORTED=true CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=/lscratch/$SLURM_JOBID
"""

rule markdupIndex:
    input:wdir +"/bam/{sample}.sorted.markdup.bam"
    output:wdir+ "/bam/{sample}.sorted.markdup.bam.bai"
    shell: "module load {program.samtools}; samtools index {input}"

if config["ref"] == "rheMac10":
    rule genrich:
        input: wdir+"/bam/{sample}.sortedByRead.bam"
        output: genrich = wdir+"/Genrich/{sample}.narrowPeak.tmp", 
                bdg = temp(wdir+"/Genrich/{sample}.bdg"), 
                bedgraph = temp(wdir+"/Genrich/{sample}.bedgraph"), 
                bed = temp(wdir+"/Genrich/{sample}.bed"), 
        shell: """
module load {program.genrich} 
{program.genrichPath} -t {input} -o {output.genrich} -j -y -r -v -d 150 -m 5 -e chrM,chrY -f {output.bdg} -b {output.bed}
cut -f 1,2,3,4 {output.bdg} | grep -v 'chrUn' | grep -v 'NW' | tail -n +2 > {output.bedgraph}
"""
else:
    rule genrich:
        input: wdir+"/bam/{sample}.sortedByRead.bam"
        output: genrich = temp(wdir+"/Genrich/{sample}.narrowPeak.tmp"), 
                bdg = temp(wdir+"/Genrich/{sample}.bdg"), 
                bedgraph = temp(wdir+"/Genrich/{sample}.bedgraph"), 
                bed = temp(wdir+"/Genrich/{sample}.bed"), 
        shell: """
module load {program.genrich} 
{program.genrichPath} -t {input} -o {output.genrich} -j -y -r -v -d 150 -m 5 -e chrM,chrY -E {reference.blacklist} -f {output.bdg} -b {output.bed}
cut -f 1,2,3,4 {output.bdg} | grep -v 'chrUn' | grep -v 'NW' | tail -n +2 > {output.bedgraph}
            """

rule remove_chrUn:
    input: wdir+"/Genrich/{sample}.narrowPeak.tmp"
    output: wdir+"/Genrich/{sample}.narrowPeak"
    shell:"""
        set +e
        grep -v 'chrUn' {input} | grep -v 'NW' > {output}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            touch {output}
        else
            exit 0
        fi
        """

def file_len(fname):
    with open(fname) as f:
        for i,l in enumerate(f,1):
            pass
    return i

rule NormBdg:
    input: bedgraph = wdir+"/Genrich/{sample}.bedgraph", 
           bed = wdir+"/Genrich/{sample}.bed"
    output: normBdg = temp(wdir + "/Genrich/{sample}.genrich.RPM.bedgraph")
    params: rname = "pl:NormBdg"
    run:
        normFactor = file_len(input.bed)/1e6
        f2 = open(output.normBdg,"w")
        with open(input.bedgraph) as f1:
            for cnt, line in enumerate(f1):
                if 'NA' not in line:
                    line2 = line.strip().split('\t')
                    line2[3] = str( round( float(line2[3]) / normFactor, 6 ) )
                    f2.write( "\t".join(line2) + "\n" )
        f2.close()

rule Bdg2Bw:
    input: bedgraph = wdir + "/Genrich/{sample}.genrich.RPM.bedgraph"
    output: bigwig = wdir +"/Genrich/{sample}.genrich.RPM.bw", 
            sorted = temp(wdir + "/Genrich/{sample}.tmp.bedgraph")
    shell: """
module load ucsc
sort -k1,1 -k2,2n -T ./ --parallel=8 {input.bedgraph} > {output.sorted}
bedGraphToBigWig {output.sorted} {reference.sizes} {output.bigwig}
 """

def plotFingerprintInput():
    inputfiles  =  ' '.join(["bam/%s.sorted.markdup.bam" % i for i in samples])
    bigwig =  ' '.join(["Genrich/%s.genrich.RPM.bw" % i for i in samples])
    labels = ' '.join(samples)
    return inputfiles, labels, bigwig

rule plotFingerprint:
    input: lambda wildcards: expand(wdir + "/bam/{sample}.sorted.markdup.bam", sample=samples),
           lambda wildcards: expand(wdir + "/Genrich/{sample}.genrich.RPM.bw", sample=samples)
    params: bam = plotFingerprintInput()[0], 
            labels = plotFingerprintInput()[1], 
            bigwig = plotFingerprintInput()[2]
    output: fp = wdir + "/Deeptools/plotFingerQualityMetrics.png", 
            tab = wdir + "/Deeptools/plotFingerQualityMetrics.tab"
    shell: """
module load {program.deeptools}
plotFingerprint --numberOfProcessors 8 -b {params.bam} --labels {params.labels} --minMappingQuality 5 --skipZeros --ignoreDuplicates --outQualityMetrics plotFingerprint_QC_metrics.txt -T "Fingerprints of ATAC-seq samples" --plotFile {output.fp} --outQualityMetrics {output.tab} --outRawCounts Deeptools/RawCounts.tab
"""

rule bigwig_summary:
    input: lambda wildcards: expand(wdir + "/bam/{sample}.sorted.markdup.bam", sample=samples),
           lambda wildcards: expand(wdir + "/Genrich/{sample}.genrich.RPM.bw", sample=samples)
    output: bw = wdir + "/Deeptools/bigwig_summary.npz"
    params: bigwig = plotFingerprintInput()[2], 
            labels = plotFingerprintInput()[1]
    shell: """
module load {program.deeptools}
multiBigwigSummary BED-file -p 8 -b {params.bigwig} --labels {params.labels} -o {output.bw} --BED {reference.deeptools_bed}
"""

rule plotCorrelation:
    input: lambda wildcards: expand(wdir + "/bam/{sample}.sorted.markdup.bam",sample=samples), 
           bw = wdir + "/Deeptools/bigwig_summary.npz"
    output: spearman = wdir +"/Deeptools/Spearman_Correlation_of_Samples_Heatmap.png", 
            pearson = wdir + "/Deeptools/Pearson_Correlation_of_Samples_Heatmap.png",  
    params: bigwig = plotFingerprintInput()[2], labels = plotFingerprintInput()[1]
    shell: """
module load {program.deeptools}
plotCorrelation -in {input.bw} --corMethod spearman --skipZeros --whatToPlot heatmap -T "Spearman Correlation of ATAC-seq samples" -o {output.spearman} --outFileCorMatrix Deeptools/SpearmanCor_bigwigScores.tab
plotCorrelation -in {input.bw} --corMethod pearson --skipZeros --whatToPlot heatmap -T "Pearson Correlation of ATAC-seq samples" -o {output.pearson} --outFileCorMatrix Deeptools/PearsonCor_bigwigScores.tab
"""

rule plotGeneHeat:
    input: lambda wildcards: expand(wdir + "/bam/{sample}.sorted.markdup.bam",sample=samples), 
           bw = wdir + "/Deeptools/bigwig_summary.npz"
    output: computeMat = temp(wdir + "/Deeptools/compute_matrix.gz"), 
            heatmap = wdir+"/Deeptools/Heatmap_of_Gene_Regions_mqc.png"
    params: bigwig = plotFingerprintInput()[2], labels = plotFingerprintInput()[1]
    shell: """
module load {program.deeptools}
computeMatrix scale-regions -p 16 -b 5000 -a 5000 -m 8000 -S {params.bigwig} -R {reference.deeptools_bed} --skipZeros -o {output.computeMat} --samplesLabel {params.labels}
plotHeatmap -m {output.computeMat} -out {output.heatmap} --dpi 100
"""

if len(samples)> 2:
    rule plotPCA:
        input: bw = wdir + "/Deeptools/bigwig_summary.npz"
        output: wdir + "/Deeptools/PCA.png"
        shell:"""
module load {program.deeptools}
plotPCA -in {input} -o {output} --outFileNameData Deeptools/PCA.tab
            """

rule alignmentSieve:
    input: sortBam = wdir + "/bam/{sample}.sorted.markdup.bam"
    output: bam = temp(wdir + "/bam/{sample}.sorted.markdup.filtered.bam"), 
           sort = temp(wdir + "/bam/{sample}.sorted.markdup.filtered.sorted.bam")
    params: smp = "{sample}"
    shell: """
module load {program.samtools}
module load {program.deeptools}
alignmentSieve -p 8 -b {input.sortBam} -o {output.bam} -l {params.smp} --ATACshift --ignoreDuplicates --minMappingQuality 5 --minFragmentLength 20
samtools sort -@24 {output.bam} -o {output.sort}
samtools index {output.sort}
        """

def bamPEFragmentSizeInput():
    inputfiles  =  ' '.join(["bam/%s.sorted.markdup.filtered.sorted.bam" % i for i in samples])
    labels = ' '.join(samples)
    return inputfiles, labels


rule bamPEFragmentSize:
    input :lambda wildcards: expand(wdir + "/bam/{sample}.sorted.markdup.filtered.sorted.bam", sample=samples)
    output: wdir + "/Deeptools/outRawFragmentLengths.txt"
    params: bamFiles = bamPEFragmentSizeInput()[0], labels = bamPEFragmentSizeInput()[1]
    shell: """
module load {program.deeptools}
bamPEFragmentSize -T "Fragment size of ATAC-seq data" --maxFragmentLength 1000 -b {params.bamFiles} --samplesLabel {params.labels} --outRawFragmentLengths Deeptools/outRawFragmentLengths.txt --table Deeptools/table.txt
        """

rule modPeak:
    input: wdir + "/Genrich/{sample}.narrowPeak"
    output: temp(wdir + "/Genrich/{sample}.narrowPeak.saf")
    shell: """
awk 'BEGIN {{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}} {{print $4, $1, $2+1, $3, "."}}' {input} > {output}
        """

rule featureCounts:
    input: saf = wdir + "/Genrich/{sample}.narrowPeak.saf", 
           bam = wdir + "/bam/{sample}.sortedByRead.bam"
    output: fc = temp(wdir + "/Genrich/{sample}.genrich.fc.txt"), 
            err = wdir + "/Genrich/{sample}.genrich.fc.err", 
            log = temp(wdir + "/Genrich/{sample}.genrich.fc.log")
    shell:"""
module load {program.featureCounts}
featureCounts -p -T 8 -a {input.saf} -F SAF -o {output.fc} {input.bam}  2>{output.err} 1>{output.log}
        """

rule fripScore:
    input: wdir + "/Genrich/{sample}.genrich.fc.err"
    output: wdir + "/Genrich/{sample}.frip.score.txt"
    run:
        myfile = open(str(output), 'w')
        with open(str(input)) as f:
            for line in f:
                if 'Process BAM file' in line:
                    sample = line.split(' ')[4][:-3]
                if 'Successfully assigned alignments' in line:
                    score = line.split("(")[1].split(')')[0]
                    myfile.write(str(sample) + '\t' + str(score) + '\n')
            myfile.close()

rule ChIPseeker:
    input: expand(wdir + "/Genrich/{sample}.narrowPeak", sample=samples)
    output: wdir + "/Genomic_Annotation_among_different_ATACseq_data_mqc.png"
    shell: "module load {program.rver}; Rscript {program.chipseeker} {config[ref]} {input}"

rule multiqc:
    input: expand(wdir + "/bam/{sample}.sorted.markdup.bam", sample=samples), 
           expand(wdir + "/QC/{sample}.kraken2.report.txt", sample=samples), 
           expand(wdir+ "/QC/{sample}_trimmed_R1.sub_screen.png", sample=samples), 
           expand(wdir +"/QC/{sample}_trimmed_R2.sub_screen.png", sample=samples),
           expand(wdir + "/Genrich/{sample}.frip.score.txt", sample=samples), 
           wdir + "/Deeptools/bigwig_summary.npz",
           wdir + "/Deeptools/Heatmap_of_Gene_Regions_mqc.png", 
           wdir + "/Deeptools/plotFingerQualityMetrics.tab", 
           wdir + "/Deeptools/outRawFragmentLengths.txt", 
           wdir + "/Deeptools/Pearson_Correlation_of_Samples_Heatmap.png", 
           wdir + "/Deeptools/PCA.png",
           wdir + "/Genomic_Annotation_among_different_ATACseq_data_mqc.png"
    output: wdir + "/multiqc.html"
    shell: """
module load {program.multiqc} 
multiqc -f ./ -n {output} -c {program.multiqcYaml}
"""

