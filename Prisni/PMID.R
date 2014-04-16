series_ids <- read.table("~/Downloads/seriesids.tab")
series_PMIDs <- as.matrix(series_ids)
url <- "http://www.ncbi.nlm.nih.gov/geo/query/"

new_url <- c()
PMID <- c()

for(i in 1:nrow(series_ids)){
new_url[i] <- paste(url, series_ids[[1]][i], sep="acc.cgi?acc=") ## new url contruction
PMID[[i]] <- extractPMID(new_url[i])
}

series_PMIDs[,2] <- PMID

## Function to extract PMIDs

extractPMID <- function(new_url){
pg <- readLines(new_url)
for(i in 1:length(pg)){
	if(grepl("citation", pg[i], ignore.case=TRUE)){
		cit <- pg[i+1]
		new_cit <- strsplit(cit, " ")[[1]]
		id <- substr(new_cit[3],5,(nchar(new_cit[3])-4))
		
	}
}
return(id)
}


###

for(i in 1:length(gsm)){
	if(grepl("Characteristics", gsm[i], ignore.case=TRUE)){
		descr <- gsm[i+1]
		split_descr <- strsplit(characteristics, "<br>")[[1]]
		grade <- substr(split_descr[2], 12,nchar(split_descr[2]))
		age <- 
		
	}
}