samples.read <- read.table("/Volumes/pgweb/Desktop/Prisni/OvaCa/gsm-all", sep="\t", header=T, row.names=NULL)

samples.read <- as.matrix(samples.read)
samples.read.filter <- samples.read[-1,]
samples.read.filter <- samples.read.filter[,-1]
#ncol(samples.read.filter)-17 = 3092
samples.read.filter <- samples.read.filter[,1:3092]
samples.filter <- t(samples.read.filter[,-1])
#nrow(samples.read.filter)

##mean
mean.profile <- c()
sd.profile <- c()
for(i in 1:nrow(samples.filter)){
	mean.profile[i] <- mean(as.numeric(samples.filter[i,]), na.rm=TRUE)
	sd.profile[i] <- sd(as.numeric(samples.filter[i,]), na.rm=TRUE)
}
# y <- as.matrix(samples.read[1,3:3093])

#cbind(y[,1], sd.profile)

plot(sd.profile, xaxt = 'n', main="standard deviation plot accross genomic regions in 220 ova cancer samples", xlab="segments", ylab="standard deviation", cex=.5, col="red", ylim=c(min(sd.profile), max(sd.profile)) )
#par(new=TRUE)
#points(mean.profile, ylab="mean", cex=.5, col="green", ylim=c(min(mean.profile), max(mean.profile)))
axis(1, at=1:length(sd.profile), labels=sd.label[,1])


plot(mean.profile, xaxt = 'n', main="mean plot accross genomic regions in 220 ova cancer samples", xlab="segments", ylab="mean", cex=.5 ,  col="green", ylim=c(min(mean.profile), max(mean.profile)))
axis(1, at=1:length(mean.profile), labels=sd.label[,1])

#df <- data.frame(y[,1], sd.profile, mean.profile)
#require(ggplot2)

#ggplot(df, aes(y[,1])) +                    # basic graphical object
#geom_point(aes(y=sd.profile), colour="red") +  # first layer
#geom_point(aes(y=mean.profile), colour="green" +
#ylab("segment and mean profiles") + xlab("segments")
#)  # second layer

#####


mean.count <- table(mean.profile)
sd.count <- table(sd.profile)
sd.count[names(sd.count)==0] #115
sd.count[names(sd.count)<=0.2]
sum(sd.count[names(sd.count)<=0.2])
#  115         7            6                10                 4                35
# 177 segment regions

new.samples.append <- cbind(samples.filter, sd.profile)
new.samples.append <- cbind(new.samples.append, mean.profile)
region.names <- rownames(new.samples.append)
new.samples.append <- cbind(new.samples.append, region.names)
sample.names <- samples.read.filter[,1]
sample.names <- c(sample.names, "sd.profile", "mean.profile", "regions.names")
colnames(new.samples.append) <- sample.names

#samples.new.append.list <- data.frame(new.samples.append)

#samples.shortlist = samples.new.append.list[as.numeric(samples.new.append.list$sd.profile)<=0.2, ]

#dataOnBoth = data[data$value_2 > 0,]


filter.by.sd  <- function(new.samples.append){
	sample.profiles.shortlist <- matrix(NA, nrow=1, ncol=ncol(new.samples.append))
	#index <- c()
	c <- 0
	colnames(sample.profiles.shortlist) <- colnames(new.samples.append)
	for(i in 1:nrow(new.samples.append)){
		if(as.numeric(new.samples.append[i,length(samples.read.filter[,1])+1])<=0.2){
			#c <- c+1
			sample.profiles.shortlist <- rbind(sample.profiles.shortlist, new.samples.append[i,1:ncol(new.samples.append)]) 
			#index[c] <- i
		}
	}
	return(sample.profiles.shortlist)
	#return(index)
}

#seg.def <- rownames(new.samples.append)



shortlist <- filter.by.sd(new.samples.append)
shortlist.ovaCa <- shortlist[-1,]

write.table(shortlist.ovaCa, file="/Users/rath/PrisniProjects/Results/shortlist.profiles/OvaCa.shortlist.tab", sep="\t", row.names=T, col.names=T, quote=F)




 