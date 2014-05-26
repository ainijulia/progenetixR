##Author: prath
##Part 1

samples.read <- read.table("/Volumes/pgweb/Desktop/Prisni/CervCa/gsm-all", sep="\t", header=T, row.names=NULL)

samples.read <- as.matrix(samples.read)
samples.read.filter <- samples.read[-1,]
samples.read.filter <- samples.read.filter[,-1]
#ncol(samples.read.filter)-17 = 3092
samples.read.filter <- samples.read.filter[,1:3092]
samples.filter <- t(samples.read.filter[,-1])
#nrow(samples.read.filter)

#find sd

mean.profile <- c()
sd.profile <- c()
for(i in 1:nrow(samples.filter)){
  mean.profile[i] <- mean(as.numeric(samples.filter[i,]), na.rm=TRUE)
  sd.profile[i] <- sd(as.numeric(samples.filter[i,]), na.rm=TRUE)
}

## Part2

data.CervCa<- read.table("/Volumes/pgweb/Desktop/Prisni/CervCa/gsm-all", sep="\t", header=T, row.names=NULL)
data.CervCa <- as.matrix(data.CervCa)

scale.param <- c()
for(i in 1:length(sd.profile)){
  if(sd.profile[i]==0){
    scale.param[i] <- FALSE
  }
  else {scale.param[i] <- TRUE}
}


data.CervCa.scale <- matrix(NA, nrow=nrow(data.CervCa), ncol=ncol(data.CervCa))
for( i in 3:3093){
  data.CervCa.scale[2:nrow(data.CervCa),i] <- scale(as.numeric(unlist(data.CervCa[2:nrow(data.CervCa),i])), center=TRUE, scale=scale.param[i-2])
}
#sd(only.scale.CervCa[,1:3091]
only.scale.CervCa <- data.CervCa.scale[2:nrow(data.CervCa), 3:3093]

scaled.CervCa <- only.scale.CervCa[,-3056:-3091]
scaled.CervCa <- scaled.CervCa[,-3034:-3036]
scaled.CervCa <- scaled.CervCa[,-2829:-2842]
scaled.CervCa <- scaled.CervCa[,-2782:-2790]
scaled.CervCa <- scaled.CervCa[,-2309:-2326]
scaled.CervCa <- scaled.CervCa[,-2202:-2219]
scaled.CervCa <- scaled.CervCa[,-2087:-2103]

corr.cervCa.mat <- cor(scaled.CervCa)

library(ggplot2)
library(reshape)


corr.cervCa.m <- melt(corr.cervCa.mat)

ggplot(corr.cervCa.m, aes(X1, X2, fill=value)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "yellow")

## Part 3
## find indices for columns with corr coeff >= 0.5

col.names1 <- as.matrix(colnames(data.CervCa[,3:3093]))
col.names2 <- as.matrix(col.names1[-3056:-3091,])
col.names2  <- as.matrix(col.names2 [-3034:-3036,])
col.names2  <- as.matrix(col.names2 [-2829:-2842,])
col.names2  <- as.matrix(col.names2 [-2782:-2790,])
col.names2  <- as.matrix(col.names2 [-2309:-2326,])
col.names2  <- as.matrix(col.names2 [-2202:-2219,])
col.names2  <- as.matrix(col.names2 [-2087:-2103,])

colnames(corr.cervCa.mat) <- col.names2[,1]
colnames(scaled.CervCa) <- col.names2[,1]

#mat.cutoff3 <- corr.cervCa.mat[rowSums(corr.cervCa.mat) >= 0.5, colSums(corr.cervCa.mat) >= 0.5, drop=FALSE]
mat.cutoff4 <- ifelse(corr.cervCa.mat >= 0.5, NA, corr.cervCa.mat)
colnames(mat.cutoff4) <- col.names2[,1]

index <- which(is.na(mat.cutoff4[1,]))
highlycorr.matrix <- corr.cervCa.mat[1,index]
colnames.sample <- col.names2[index,1]
for(i in 2:nrow(mat.cutoff4)){
  index <- which(is.na(mat.cutoff4[i,]))
  highlycorr.matrix <- rbind(highlycorr.matrix, corr.cervCa.mat[i,index])
  colnames.sample <- rbind(colnames.sample, col.names2[index,1])
}
index.matrix <- as.matrix(index)

## index.matrix gives the indices for higly correlated 

sumary <- summary(mat.cutoff4)


# in the scaled matrix remove off all the columns corresponding to the correlation cut off



highlyCor <- findCorrelation(corMatMy, 0.70)
#Apply correlation filter at 0.70,
#then we remove all the variable correlated with more 0.7.
datMyFiltered.scale <- datMy.scale[,-highlyCor]
corMatMy <- cor(datMyFiltered.scale)
corrplot(corMatMy, order = "hclust")
