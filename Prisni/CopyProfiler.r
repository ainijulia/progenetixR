copy.profile.matrix <- function(xmat){
		x <- matrix(NA,nrow=nrow(xmat), ncol=7)
		x[,1:6] <- xmat[,1:6]
		for(i in 2:nrow(x)){
			if(as.numeric(x[i,5]>as.numeric(0.15))){
				x[i,7] <- 1 
			}
			else if(as.numeric(x[i,5]<as.numeric(-0.15))){
				x[i,7] <- -1  
			}
			else {
				x[i,7] <- 0
			}
		}
		return(x)
	}
	
probe.count <- function(mat){
	count <- 0
	for(i in 2:nrow(mat)){
		count <- count+as.numeric(mat[i,6])
	}
	return(count)
}
	
args <- commandArgs(TRUE)
setwd("~/Desktop/GSE11960/")
dir_list <- list.dirs()
segfilename <- c()
#read all new_segments file. The new_segments files are the one formatted to BED file format structure.
for(i in 2:length(dir_list)){
	segfilename[i-1] <- paste(dir_list[i], "segments.tab", sep="/")
	##1
	segfile_mat <- readFile(segfilename[i-1])
	new.seg.matrix <- copy.profile.matrix(segfile_mat)
	count[i-1] <- probe.count(new.seg.matrix)
	}