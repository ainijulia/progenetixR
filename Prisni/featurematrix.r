readFile <- function(f) {
	newsegfile <- read.table(f, sep="\t")
	newsegfile_mat <- as.matrix(newsegfile)
	return(newsegfile_mat)
}

featurematrix <- function(f, data.samples){
	data.sample1 <- readFile(f)
	data.samples <- rbind(data.sample, data.sample1)
	return(data.samples)
}

setwd("~/Desktop/Prisni/GSE11960")
dir_list <- list.dirs()
dir_list <- as.matrix(dir_list)
file.list <- dir_list[-1,]
data.samples <- readFile(file.list[1,1])
for(i in 2:nrow(file.list)){
	feature.matrix <- feature.matrix(file.list[i,1], data.samples) 
}
data.sample1 <- readFile("GSM302687/normalized.segments.tab", sep="\t")
data.sample2 <- readFile("GSM302688/normalized.segments.tab", sep="\t")
data.samples <- rbind(data.sample1, data.sample2)
