################################
##This script converts hg18 to hg19 build formats##
################################
#if you want see commands in output file
#options(echo=TRUE)
#args <- commandArgs(trailingOnly = TRUE)
#print(args)
#setwd(args)
################################
file <- c()
setwd("~/Desktop/GSE11960") #set ~/Desktop/GSE11960 as working directory
dirs <- list.dirs()

################################
##This function reads the CNprobe format file and converts it to a BED standard format.
##The new BED format file is saved in the same directory as CNprobe file
cnprobeToBed <- function(file){
	lstall <- read.table(file, skip=1, sep="\t", fill=TRUE, nrow=length(count.fields(file))-1)
	newmat <- as.matrix(lstall)

	combine <- matrix(NA, ncol=5, nrow=nrow(lstall))
		for (i in 1:nrow(lstall)){
			print(i)
 			combine[i,1] <- newmat[i,2]
 			combine[i,2] <- newmat[i,3]
 			combine[i,3] <- as.numeric(newmat[i,3])+1
 			combine[i,4] <- newmat[i,4]
 			combine[i,5] <- newmat[i,1]
 			print("done")
 	 	}
	combine[,1] <- paste("chr", combine[,1], sep="")
write.table(combine, file=paste(gsub("CNprobes.tab", "CNprobes18.tab", file), sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
print("written to CNprobes18 file")
}

###################################
##This function reads in the BED format file and runs liftOver Command utility tool
##Two new files will be created in the same folder called CNprobes19.mapped.tab and CNprobes19_u.tab

stepLiftOver <- function(file){
	inputfile <- gsub("CNprobes", "CNprobes18", file)
	outputmappedfile <- gsub("CNprobes.tab", "CNprobes19.mapped.tab", file)
	outputunmappedfile <- gsub("CNprobes.tab", "CNprobes19_u.tab", file)
	#Syntax: liftOver <inputfile> <location of chain file> <outputfile> <unmappedfile>
	#NB: the chain file location has to be hardcoded and correct chain file has to be pasted to ensure functioning of the this tool
	command <- paste("~/Desktop/bin/liftOver", inputfile, "~/Desktop/hg18ToHg19.over.chain", outputmappedfile, outputunmappedfile, sep=" ") 
	system(command)
	print("successful conversion")
}

#################################
##This fucntion reads in the mapped format converted files in BED format and generates a CNprobe format file.
##The new file is saved as CNprobes19.tab file in the same folder as the input file
bedToCNprobe <- function(file){
	cnprobefile <- gsub("CNprobes.tab", "CNprobes19.mapped.tab", file)
	newlist <- read.table(cnprobefile, skip=1, sep="\t", fill=TRUE)
	print("read CNprobes 19 mapped file in BED format")
	listmat <- as.matrix(newlist)
	newcombine <- matrix(NA, ncol=4, nrow=nrow(newlist))
	newcombine[1,1] <- "ID"
	newcombine[1,2] <- "chro"
	newcombine[1,3] <- "pos"
	newcombine[1,4] <- "log2"

	for(i in 2:nrow(newlist)){ #for loop no. 2 start
	
	chro <- listmat[i,1]
	startpos <- listmat[i,2]
	endpos <- listmat[i,3]
	log2 <- listmat[i,4]
	ID <- listmat[i,5]
	
	pos <- round(mean(startpos:endpos), digit=0)
	
	newcombine[i,1] <- ID
	newcombine[i,2] <- chro
	newcombine[i,3] <- pos
	newcombine[i,4] <- log2
	
	} #for loop no. 2 end
		#write final file matrix as a tab file
	newcombine[,2] <- gsub("chr", "", newcombine[,2])
	newcombine[1,2] <- gsub("o", "chro", newcombine[1,2])
	write.table(newcombine, file=paste(gsub("CNprobes.tab", "CNprobes19.tab", file), sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
	print("CNprobes19 file written in CNprobes format")
	}

##################################
##execute the functions for all files in the working directory.

for(i in 2:length(dirs)){
file[i-1] <- paste(dirs[i], "CNprobes.tab", sep="/")
cnprobeToBed(file[i-1])
stepLiftOver(file[i-1])
bedToCNprobe(file[i-1])
}


