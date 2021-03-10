args = commandArgs(trailingOnly=TRUE)
DiffBind <- read.table(args[1],header=T,stringsAsFactors=F)
Uropa <- read.table(args[2],header=T,stringsAsFactors=F)
DiffBind <- data.frame(DiffBind,peak_id=paste0("Peak_",1:nrow(DiffBind)),stringsAsFactors=F)
names(DiffBind)[1:3] <- c("peak_chr","peak_start","peak_end")
alldata <- merge(DiffBind,Uropa)
alldata <- alldata[,c(4,1,2,13,3,8:12,14:26)]
alldata <- alldata[order(alldata$FDR, alldata$peak_id),]
root <- unlist(strsplit(args[2],split=".txt"))
outFile <- paste0(root,".merged.txt")
write.table(alldata, outFile, quote=F, sep="\t", row.names=F)


