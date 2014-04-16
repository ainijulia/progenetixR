args <- commandArgs(TRUE)
setwd(args[1])
dir_list <- list.dirs()
filename <- c()
for(i in 2:length(dir_list)){
	filename[i-1] <- paste(dir_list[i], "segments.tab", sep="/")
	processfile(filename[i-1])
}

processfile <- function(f){
	seg_read <- read.table(f, sep="\t")
	seg_read_mat <- as.matrix(seg_read)
		new_seg_read_mat <- matrix(NA, nrow=nrow(seg_read_mat), ncol=ncol(seg_read_mat)+1)

			for( i in 2:nrow(seg_read_mat)){
				if(as.numeric(seg_read_mat[i,5]) > 0.15){
					new_seg_read_mat[i,7] <- 1
					}
				else if(as.integer(seg_read_mat[i,5]) < -0.15){
						new_seg_read_mat[i,7] <- -1

						} else {
								new_seg_read_mat[i,7] <- 0
								}
				}
			for( i in 1:nrow(seg_read_mat)){
				new_seg_read_mat[i,1] <- seg_read_mat[i,1]
				new_seg_read_mat[i,2] <- seg_read_mat[i,2]
				new_seg_read_mat[i,3] <- seg_read_mat[i,3]
				new_seg_read_mat[i,4] <- seg_read_mat[i,4]
				new_seg_read_mat[i,5] <- seg_read_mat[i,5]
				new_seg_read_mat[i,6] <- seg_read_mat[i,6]
				}
	new_seg_read_mat[1,7] <- "copy_status"
write.table(new_seg_read_mat, file=gsub("segments","new_segments",f), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
}
