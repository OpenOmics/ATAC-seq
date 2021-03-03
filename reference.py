import config

if config.ref == "hg38":
    refgen = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/indexes/hg38"
    blacklist = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/indexes/hg38.blacklist.bed.gz"
    sizes = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/indexes/hg38.fa.sizes"
    deeptools_bed = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/genes.bed"

if config.ref == "hg19":
    refgen = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/indexes/hg19"
    blacklist = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/indexes/hg19.blacklist.bed.gz"
    sizes = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/indexes/hg19.fa.sizes"
    deeptools_bed = "/data/CCBR_Pipeliner/db/PipeDB/Indices/hg19_basic/genes.bed"

if config.ref == "mm10":
    refgen = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/mm10"
    blacklist = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/mm10.blacklist.bed.gz"
    sizes = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/indexes/mm10.fa.sizes"
    deeptools_bed = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/genes.bed"

if config.ref == "mm9":
    refgen = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm9_basic/indexes/mm9"
    blacklist = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm9_basic/indexes/mm9.blacklist.bed.gz"
    sizes = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm9_basic/indexes/mm9.fa.sizes"
    deeptools_bed = "/data/CCBR_Pipeliner/db/PipeDB/Indices/mm9_basic/genes.bed"

if config.ref == "rheMac10":
    refgen = "/data/NCBR/references/RheMac10/bowtie2/rheMac10"
    blacklist = ""
    sizes = "/data/NCBR/references/RheMac10/rheMac10.chrom.sizes"
    deeptools_bed = "/data/NCBR/references/RheMac10/rheMac10.ncbiRefSeq.deeptools.sorted.merged.bed"


