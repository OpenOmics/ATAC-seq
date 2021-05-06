args = commandArgs(trailingOnly=TRUE)
DiffBind <- read.table(args[1],header=T,stringsAsFactors=F)
Homer <- read.delim(args[2],stringsAsFactors=F)
DiffBind <- data.frame(DiffBind,Peak_ID=paste0("Peak_",1:nrow(DiffBind)),stringsAsFactors=F)
names(DiffBind)[1:3] <- c("Chr","Start","End")
Homer$Start <- Homer$Start -1
names(Homer)[1] <- "Peak_ID"
alldata <- merge(DiffBind,Homer)
alldata <- alldata[,c(4,1,2,3,8:12,16:18,20,23:25,27)]
alldata <- alldata[order(alldata$FDR, alldata$Peak_ID),]
root <- unlist(strsplit(args[2],split=".txt"))
outFile <- paste0(root,".merged.txt")
write.table(alldata, outFile, quote=F, sep="\t", row.names=F)
