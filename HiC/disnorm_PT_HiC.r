#pos_raw <- read.table("/scratch/sh8tv/Project/TCF1_Tcell/Data/process/HiC/CD8_ctrl/rawData_5000_abs.bed")
#mat2_raw <- read.table("/scratch/sh8tv/Project/TCF1_Tcell/Data/process/HiC/CD8_ctrl/rawData_5000.matrix")
#mat1_raw <- read.table("/scratch/sh8tv/Project/TCF1_Tcell/Data/process/HiC/CD8_Tcf1Lef1dko/rawData_5000.matrix")
a <- commandArgs(T)
pos_raw <- read.table(a[1])
mat1_raw <- read.table(a[2])
mat2_raw <- read.table(a[3])
useChrm <- a[4]
outname <- a[5]

total1 <- sum(as.numeric(mat1_raw[,3]))
total2 <- sum(as.numeric(mat2_raw[,3]))
seq_depth_factor <- total2/total1

binsize <- pos_raw[1,3] - pos_raw[1,2]
dis_range <- 1e6
COLN <- (2*dis_range / binsize) + 1
midline <- dis_range/binsize

use_range <- 2e5
usebinnum <- use_range/binsize
if(useChrm == "all"){
    usechrbin <- 1:nrow(pow_raw)
}else{
    usechrbin <- which(pos_raw[,1]==useChrm)
}

load_mat <- function(posMat, countMat){
    tmpmat <- matrix(rep(0,nrow(posMat)*COLN),ncol=COLN)
    colnames(tmpmat) <- seq(-midline,midline)*binsize
    usemat <- countMat[which(countMat[,2]-countMat[,1] <= midline & countMat[,1] %in% usechrbin),]
    for(i in 1:nrow(usemat)){
        tmpmat[usemat[i,1], midline+1+usemat[i,2]-usemat[i,1]] <- usemat[i,3]
        tmpmat[usemat[i,2], midline+1-(usemat[i,2]-usemat[i,1])] <- usemat[i,3]
    }
    return(tmpmat)# / as.numeric(apply(tmpmat,2,mean))) )
}
PT <- function(inline,binnum){
    L <- length(inline)
    binindex <- c((200-binnum+1):200, 202:(202+binnum-1))
    #binindex <- c(-201)
    a <- inline[1:(L/2)][binindex]
    b <- inline[((L/2)+1):L][binindex]

    #LFC <- mean(a)-mean(b)
    LFC <- mean(log2(a+0.5) - log2(b+0.5))
    if(LFC > 0){
        P <- t.test(a,b,alternative="greater",paired=T)$p.val
    }else{
        P <- t.test(a,b,alternative="less",paired=T)$p.val
    }
    return(c(LFC,-log10(P)))
}
PTdown <- function(inline,binnum){
    L <- length(inline)
    binindex <- c(202:(202+binnum-1))
    a <- inline[1:(L/2)][binindex]
    b <- inline[((L/2)+1):L][binindex]
    #LFC <- mean(a)-mean(b)
    LFC <- mean(log2(a+0.5) - log2(b+0.5))
    if(LFC > 0){
        P <- wilcox.test(a,b,alternative="greater",paired=T)$p.val
    }else{
        P <- wilcox.test(a,b,alternative="less",paired=T)$p.val
    }
    return(c(LFC,-log10(P)))
}
raw1 <- load_mat(pos_raw,mat1_raw)*seq_depth_factor
raw2 <- load_mat(pos_raw,mat2_raw)
norm1 <- scale(raw1)
norm2 <- scale(raw2)

useidx <- which(apply(raw1,1,max)>0 & apply(raw2,1,max)>0 )
#usedata <- cbind(norm1,norm2)[useidx,]
usedata <- cbind(raw1,raw2)[useidx,]

PVFC <- t(apply(usedata,1,PT,usebinnum))
colnames(PVFC) <- c("log2FC","logP")

PVFCdown <- t(apply(usedata,1,PTdown,usebinnum))
colnames(PVFCdown) <- c("log2FC","logP")

PVFC[3579,]
PVFCdown[3579,]

PT(usedata[3579,],usebinnum)
PTdown(usedata[3579,],usebinnum)


usePos <- pos_raw[useidx,]
usePos[,4] <- paste0("b",usePos[,4])
outdata <- cbind(usePos, PVFC)
write.table(outdata,file=outname,row.names=F,col.names=F,sep="\t",quote=F)

#write.table(outdata,file="test_Tcf1ko_diffHiC_scale_chr10.bed",row.names=F,col.names=F,sep="\t",quote=F)


#tmpmat1 <- matrix(rep(0,nrow(pos_raw)*COLN),ncol=COLN)
#colnames(tmpmat1) <- seq(-midline,midline)*binsize
#
#tmpmat2 <- matrix(rep(0,nrow(pos_raw)*COLN),ncol=COLN)
#colnames(tmpmat2) <- seq(-midline,midline)*binsize
#
#
#
##usemat1 <- mat1_raw[which(mat1_raw[,2]-mat1_raw[,1] <= midline),]
##usemat2 <- mat2_raw[which(mat2_raw[,2]-mat2_raw[,1] <= midline),]
#usemat1 <- mat1_raw[which(mat1_raw[,2]-mat1_raw[,1] <= midline & mat1_raw[,1] %in% usechrbin),]
#usemat2 <- mat2_raw[which(mat2_raw[,2]-mat2_raw[,1] <= midline & mat2_raw[,1] %in% usechrbin),]
#
##for(i in 1:nrow(usemat1)){
#ptm <- proc.time()
#for(i in 1:nrow(usemat1)){
#    tmpmat1[usemat1[i,1], midline+1+usemat1[i,2]-usemat1[i,1]] <- usemat1[i,3]
#    tmpmat1[usemat1[i,2], midline+1-(usemat1[i,2]-usemat1[i,1])] <- usemat1[i,3]
#    if(i %% 100000 == 0){
#        print(c(i, (proc.time() - ptm)[3]))
#        ptm <- proc.time()
#    }
#}
#
#ptm <- proc.time()
#for(i in 1:nrow(usemat2)){
#    tmpmat2[usemat2[i,1], midline+1+usemat2[i,2]-usemat2[i,1]] <- usemat2[i,3]
#    tmpmat2[usemat2[i,2], midline+1-(usemat2[i,2]-usemat2[i,1])] <- usemat2[i,3]
#    if(i %% 100000 == 0){
#        print(c(i, (proc.time() - ptm)[3]))
#        ptm <- proc.time()
#    }
#}
##mean1 <- apply(tmpmat1,2,mean)
#norm1 <- t(t(tmpmat1) / as.numeric(apply(tmpmat1,2,mean)))#  %*% as.matrix(1/as.numeric(apply(tmpmat1,2,mean)))
#norm2 <- t(t(tmpmat2) / as.numeric(apply(tmpmat2,2,mean)))#as.matrix(tmpmat2) %*% as.matrix(1/as.numeric(apply(tmpmat2,2,mean)))
#useidx <- which(apply(norm1,1,max)>0 & apply(norm2,1,max)>0 )
#
#PT <- function(inline){
#    L <- length(inline)
#    a <- inline[1:(L/2)]
#    b <- inline[((L/2)+1):L]
#    LFC <- log2(mean(a)/mean(b))
#    if(LFC > 0){
#        P <- t.test(a,b,alternative="greater",paired=T)$p.val
#    }else{
#        P <- t.test(a,b,alternative="less",paired=T)$p.val
#    }
#    return(c(LFC,P))
#}
#
#usedata <- cbind(norm2,norm1)[useidx,]
#
#PVFC <- apply(usedata,1,PT)
#
#
#
#
#
#
#
#
#
#
#
#