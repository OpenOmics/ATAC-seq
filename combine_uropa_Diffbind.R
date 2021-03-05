# I grabbed this code from a script I was doing for a slightly different purpose so the naming isn't quite right, but I have uused it for what we want

# this works for "starts" and finalhits files

merge_starts_DiffBind <- function(startsFile,DiffBindFile) {
        starts <- read.table(startsFile,header=T,stringsAsFactors=F)
        DiffBind <- read.table(DiffBindFile,header=T,stringsAsFactors=F)
        # pick one of the next two lines depending on whether the peak numbering started at 0 or 1
#       DiffBind <- data.frame(DiffBind,peak_id=paste0("Peak",1:nrow(DiffBind)),stringsAsFactors=F)
        DiffBind <- data.frame(DiffBind,peak_id=paste0("Peak",0:(nrow(DiffBind)-1)),stringsAsFactors=F)
        names(DiffBind)[1:3] <- c("peak_chr","peak_start","peak_end")
        alldata <- merge(DiffBind,starts)
        alldata <- alldata[,c(4,1,2,13,3,8:12,14:26)]
        alldata <- alldata[order(alldata$FDR, alldata$peak_id),]
        root <- unlist(strsplit(startsFile,split=".txt"))
        outFile <- paste0(root,".score.txt")
        write.table(alldata, outFile, quote=F, sep="\t", row.names=F)
}