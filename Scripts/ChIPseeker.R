library(ChIPseeker)
library(ggplot2)
#library(clusterProfiler)
#library(ReactomePA)
args <- commandArgs(trailingOnly = TRUE)

#args <- list("hg38", "Genrich/GH3564.narrowPeak", "Genrich/GH3565.narrowPeak", "Genrich/GH3566.narrowPeak", "Genrich/GH3568.narrowPeak", "Genrich/GH3562.narrowPeak", "Genrich/GH3560.narrowPeak", "Genrich/GH3563.narrowPeak", "Genrich/GH3570.narrowPeak", "Genrich/GH3561.narrowPeak", "Genrich/GH3571.narrowPeak", "Genrich/GH3567.narrowPeak", "Genrich/GH3569.narrowPeak")
genome <- args[1]

if (genome == "mm10") {
  library(org.Mm.eg.db)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  orgdb = "org.Mm.eg.db"
  organism = "mmu"
}else if(genome == "hg38" ){
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  orgdb = "org.Hs.eg.db"
  organism = "hsa"
}else if(genome == "hg19" ){
    library(org.Hs.eg.db)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    orgdb = "org.Hs.eg.db"
    organism = "hsa"
}else if(genome == "rheMac10" ){
  if (!requireNamespace("org.Mmu.eg.db", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
    BiocManager::install("org.Mmu.eg.db")
  }
  if (!requireNamespace("TxDb.Mmulatta.UCSC.rheMac10.refGene", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
    BiocManager::install("TxDb.Mmulatta.UCSC.rheMac10.refGene")
  }
  library(org.Mmu.eg.db)
  library(TxDb.Mmulatta.UCSC.rheMac10.refGene)
  txdb <- TxDb.Mmulatta.UCSC.rheMac10.refGene
  orgdb = "org.Mmu.eg.db"
  organism = "Mmu"
}

input <- args[2:length(args)]
input1 <- input[sapply(input, file.size) > 0]
nm <- lapply(gsub('.narrowPeak', '' , lapply(gsub('.*/Genrich/', '', input1), function(x) (x))), function(y) (y))

print (input1)
names(input1) <- nm
input1 <- input1[order(names(input1))]

#########################################################
## Profile of ChIP peaks binding to TSS regions
#########################################################
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList_all <- lapply(input1, getTagMatrix, windows=promoter)
tagMatrixList <- tagMatrixList_all[sapply(tagMatrixList_all, nrow) > 0]

if (is.list(tagMatrixList) & length(tagMatrixList) != 0) {
  if (length(input1) <= 10) { rel_size=1; ncolLegend=1 
  } else if (length(input1) <= 40) {rel_size=0.8; ncolLegend=1 
  } else {rel_size=0.4; ncolLegend=2}
  png('Average_Profile_of_ATAC_peaks_binding_to_TSS_region_mqc.png',res = 200, width = 6, height = 4, units = "in")
  p <- plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Frequency")
  p <- p + theme(legend.text=element_text(size=rel(rel_size)), legend.key.size = unit(0.2, "cm")) +
           guides(color=guide_legend(ncol=ncolLegend))
  print(p)
  dev.off()
}

#########################################################
## Peak Annotation
#########################################################
png('Genomic_Annotation_among_different_ATACseq_data_mqc.png',res = 200, width = 6, height = 4, units = "in")
peakAnnoList <- lapply(input1, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

if (length(input1) <= 10) { rel_size=1 
} else if (length(input1) <= 40) {rel_size=0.8
} else {rel_size=0.4}
p <- plotAnnoBar(peakAnnoList)
p <- p + theme(axis.text.y=element_text(size=rel(rel_size)))
print(p)
dev.off()

#########################################################
## Visualize distribution of TF-binding loci relative to TSS
#########################################################

png('Distribution_of_Binding_Sites_among_different_ATACseq_data_mqc.png', res = 200, width = 8, height = 4, units = "in")
if (length(input1) <= 10) { rel_size=1 
} else if (length(input1) <= 40) {rel_size=0.8 
} else {rel_size=0.4}
p <- plotDistToTSS(peakAnnoList)
p <- p + theme(axis.text.y=element_text(size=rel(rel_size)))
print(p)
dev.off()

