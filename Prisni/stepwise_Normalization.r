

# read file as matrix
readFile <- function(f) {
	newsegfile <- read.table(f, sep="\t")
	newsegfile_mat <- as.matrix(newsegfile)
	return(newsegfile_mat)
}


# find chromosome matrix

findChromosomeWiseData <- function(chr, mat){
	profile_chr <- c()
	c <- 0
	for(i in 2:nrow(mat)){
		if(as.numeric(chr)==mat[i,2]){
			print(i)
			c <- c+1
			profile_chr[[c]] <- c(mat[i,1],mat[i,2],mat[i,3],mat[i,4],mat[i,5],mat[i,6],mat[i,7]) 
		}
	}
	return(profile_chr)
}

create.sample.chr.matrix <- function(profile.list){
	sample.chr.mat <- rbind(profile.list[[1]])
	for(i in 2:length(profile.list)){
	sample.chr.mat <- rbind(sample.chr.mat, profile.list[[i]])
	}
	
	return(sample.chr.mat)
}
# from this matrix separate normal, gain and loss segment  profiles

# count gain/loss/normal profiles
countCopyProfile <- function(data, copy_status){
	count <- 0
	for(i in 1:length(data)){
		#print(i)
		if(as.numeric(copy_status)==as.numeric(data[[i]][7])){
		count <- count+1
		}
	}
	return(count) 
}

##usage
## gives number of loss, normal and amplified segments per chromosome
countTotalCopyProfiles <- c(countCopyProfile(chr5_sample1, -1),countCopyProfile(chr5_sample1, 0), countCopyProfile(chr5_sample1, 1))

# separate copy status based profiles

#find gain/loss/normal matrix
segmentProfilesPerChromosome <- function(data, copy_status, copy_count){
	segment_profile_matrix <- matrix(NA, nrow=copy_count, ncol=7)
	c <- 0
for(i in 1:length(data)){
	if(as.numeric(copy_status)==as.numeric(data[[i]][7])){
		c <- c+1
		segment_profile_matrix[c,1] <- data[[i]][1]
		segment_profile_matrix[c,2] <- data[[i]][2]
		segment_profile_matrix[c,3] <- data[[i]][3]
		segment_profile_matrix[c,4] <- data[[i]][4]
		segment_profile_matrix[c,5] <- data[[i]][5]
		segment_profile_matrix[c,6] <- data[[i]][6]
		segment_profile_matrix[c,7] <- data[[i]][7]
	}
}
return(segment_profile_matrix)
}

##usage
chr5_sample1_loss <- segmentProfilesPerChromosome(chr5_sample1, -1, countTotalCopyProfiles[1])
chr5_sample1_normal <- segmentProfilesPerChromosome(chr5_sample1, 0, countTotalCopyProfiles[2])
chr5_sample1_gain <- segmentProfilesPerChromosome(chr5_sample1, 1, countTotalCopyProfiles[3])

# convert to BED format

getBEDFormatMatrix <- function(sample_matrix){
	bed_matrix <- matrix(NA, nrow=nrow(sample_matrix), ncol=7)
	for(i in 1:nrow(sample_matrix)){
		bed_matrix[i,1] <- paste("chr",sample_matrix[i,2], sep="")
		bed_matrix[i,2] <- sample_matrix[i,3]
		bed_matrix[i,3] <- sample_matrix[i,4]
		bed_matrix[i,4] <- sample_matrix[i,5]
		bed_matrix[i,5] <- sample_matrix[i,6]
		bed_matrix[i,6] <- sample_matrix[i,7]
		bed_matrix[i,1] <- sample_matrix[i,1]
	}
	return(bed_matrix)
}


# check total number of probes. For samples belonging to same platform IDs have same count of probe marker points
probeCount <- function(mat){
	count <- 0
	for(i in 2:nrow(mat)){
		count <- count+as.numeric(mat[i,6])
	}
	return(count)
}

# 
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# Normalized Matrix
normalize.matrix <- function(nmat){
	for(i in 1:nrow(nmat)){
		non.norm.matrix <- matrix(exmat[i,1:ncol(nmat)], nrow=1)
		norm.matrix <- rep.row(non.norm.matrix[1,1:ncol(nmat)],as.numeric(non.norm.matrix[1,6]))
	}
	norm.matrix <- rbind(norm.matrix, norm.matrix)
	return(norm.matrix)
}

# Final Normalized Matrix (merge norm.matrix for copy gain, loss and neutral)
norm.matrix.gain <- normalize.matrix(nmat)
norm.matrix.loss <- normalize.matrix(nmat)
norm.matrix.neutral <- normalize.matrix(nmat)

