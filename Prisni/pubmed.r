series_ids <- read.table("~/Downloads/seriesids.tab")
series_PMIDs <- as.matrix(series_ids)
url <- "http://www.ncbi.nlm.nih.gov/geo/query/"

new_url <- c()
PMID <- c()
# 3699
for(i in 1:length(series_ids)){
print(series_ids[i])
new_url[i] <- paste(url, series_ids[i], sep="acc.cgi?acc=") ## new url contruction
PMID[[i]] <- extractPMID(new_url[i])
print("PMID extracted")
}

series_PMIDs <- PMID

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


author_line <- c()
for(i in 1:length(lines)){
if(grepl("author information", lines[i], ignore.case=TRUE)){
author_line <- strsplit(lines[i], "><")

}
}
for(i in 1:length(author_line)){
if(grepl("author information", author_line[[1]][i])){
	print(author_line[[1]][i])
}
}