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
sd(only.scale.CervCa[,1:3091]
 only.scale.CervCa <- data.CervCa.scale[2:nrow(data.CervCa), 3:3093]
 
corr.cervCa.mat <- cor(scaled.CervCa)
corrplot(corr.cervCa.mat, order="hclust")

filter.data <- function(mat){
	index <- c()
	c <- 0
	for(i in 1:ncol(mat)){
		if(sd(mat[,i])==0){
			c <- c+1
			index[c] <- i
		}
	}
	return(index)
}

index.sd.0 <- filter.data(only.scale.CervCa)

scaled.CervCa <- only.scale.CervCa[,-3056:-3091]
scaled.CervCa <- scaled.CervCa[,-3034:-3036]
scaled.CervCa <- scaled.CervCa[,-2829:-2842]
scaled.CervCa <- scaled.CervCa[,-2782:-2790]
scaled.CervCa <- scaled.CervCa[,-2309:-2326]
scaled.CervCa <- scaled.CervCa[,-2202:-2219]
scaled.CervCa <- scaled.CervCa[,-2087:-2103]



data.CervCa.scale.1 <- scale(as.numeric(data.CervCa[2:nrow(data.CervCa),3]), center=TRUE, scale=TRUE)
sd(data.CervCa.scale.1)
data.CervCa.scale.2 <- scale(as.numeric(data.CervCa[2:nrow(data.CervCa),4]), center=TRUE, scale=TRUE)
sd(data.CervCa.scale.2)
data.CervCa.scale.3 <- scale(as.numeric(data.CervCa[2:nrow(data.CervCa),5]), center=TRUE, scale=TRUE)
sd(data.CervCa.scale.3)
corr.cervCa.mat <- cor(data.CervCa.scale)
corrplot(corr.cervCa.mat, order="hclust")
datMy.scale<- scale(datMy[2:ncol(datMy)],center=TRUE,scale=TRUE);
#scale all the features (from feature 2 bacause feature 1 is the predictor outp
corMatMy <- cor(datMy.scale)
#compute the correlation matrix
corrplot(corMatMy, order = "hclust")
