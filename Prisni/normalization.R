args <- commandArgs(TRUE)
setwd("~/Desktop/GSE11960/")
dir_list <- list.dirs()
segfilename <- c()
chr <- c(1:23)
#read all new_segments file. The new_segments files are the one formatted to BED file format structure.
for(i in 2:length(dir_list)){
	segfilename[i-1] <- paste(dir_list[i], "new_segments.tab", sep="/")
	##1
	segfile_mat <- readFile(segfilename[i-1])
	
	}

######### 


segfile_sample1_copy_matrix

	counter <- countCopyProfile(segfile_sample1, gain, chr[1])
	segment_profile <- sampleChrCopyProfiler(segfile_sample1, chr[1], counter)

##1
readFile <- function(f) {
	newsegfile <- read.table(f, sep="\t")
	newsegfile_mat <- as.matrix(newsegfile)
	return(newsegfile_mat)
}

##2
# for each segment file the following function reads column number 2 for chromosome.  and returns a list of chrmosome number. chr {1 .. 23}

getChromosomeSet <- function(mat){
	return(intersect(mat[,2], mat[,2]))
}

##3
# for each chromosome from the "chr" list, find the segment size

getAllSegSize <- function(mat, chr){
	for(i in 2:nrow(mat)){
		if(as.numeric(chr)==as.numeric(mat[i,2])){
			print("segment size calculating")
			segment_size <- as.numeric(mat[i,4]) - as.numeric(mat[i,3])
		}
	}
	return(sum(segment_size))
}

##3.1

countCopyProfile <- function(mat, copy_status, chr){
	count <- 0
	for(i in 2:nrow(mat)){
		#print(i)
		if(copy_status == as.numeric(mat[i,7])){
				if(as.numeric(chr)==as.numeric(mat[i,2])){
				count <- count+1
				}
		}
	}
	return(count) 
} 

##3.2 Per chromosome: table of following contents:
## sample	start stop	probes	log2	copy_status
sample_chr_profile <- c()
sampleChrCopyProfiler <- function(mat, chr, count, sample_chr_profile){
	
for(i in 2:nrow(mat)){
	if(as.numeric(chr)==as.numeric(mat[i,2])){
		sample_chr_profile[[1]][i] <- mat[i,1] 
		sample_chr_profile[[2]][i] <- mat[i,2] 
		sample_chr_profile[[3]][i] <- mat[i,3] 
		sample_chr_profile[[4]][i] <- mat[i,4] 
		sample_chr_profile[[5]][i] <- mat[i,5]
		sample_chr_profile[[6]][i] <- mat[i,6]
		sample_chr_profile[[7]][i] <- mat[i,7]   
		}
	}
	
	return(sample_chr_profile)
}

##4
# Create a list of Nmax_local for all samples. Nmax_global = max(Namx_local_list)


##5

getStartPosition <- function(mat, chr){
	startpos <- c()
	for(i in 2:nrow(mat)){
		if(as.numeric(chr[[2]])==as.numeric(mat[i,2])){
			print("startpos")
			startpos[[i-1]] <- mat[i,3]
		}
	}
	return(startpos)
}

getEndPosition <- function(mat, chr){
	endpos <- c()
	for(i in 2:nrow(mat)){
		if(as.numeric(chr[[2]])==as.numeric(mat[i,2])){
			print("endpos")
			endpos[[i-1]] <- mat[i,4]
		}
	}
	return(endpos)
	}
	
##5
# recreate segment files with new start and end positions in BED format

getBEDFormatMatrix <- function(chro, start, end, Nscale){
	bed_matrix <- matrix(NA, nrow=Nscale, ncol=3)
	for(i in 1:Nscale){
		bed_matrix[i,1] <- paste("chr",chro, sep="")
		bed_matrix[i,2] <- start[[i]]
		bed_matrix[i,3] <- end[[i]]
	}
	return(bed_matrix)
}










