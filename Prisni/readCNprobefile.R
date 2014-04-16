readCNprobefile <- function(file) {
	print("hello0")
	lstall <- read.table(file, sep="\t")
	print("hello1")
	newmat <- as.matrix(lstall)
	print("hello2")
 combine <- matrix(NA, ncol=6, nrow=nrow(lstall)-1)
 for (i in 1:nrow(lstall)){
 	combine[i,1] <- newmat[i+1,2]
 	combine[i,2] <- newmat[i+1,3]
 	combine[i,3] <- newmat[i+1,4]
 	combine[i,4] <- newmat[i+1,5]
 	combine[i,5] <- newmat[i+1,6]
 	combine[i,6] <- newmat[i+1,1]
 }
 
 combine[,1] <- paste("chr", combine[,1], sep="")
 return(combine)
}
