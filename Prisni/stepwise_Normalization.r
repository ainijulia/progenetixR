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


