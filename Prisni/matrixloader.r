c.sample <- read.table("~/Desktop/Prisni/RenalCa/renal-carcinoma-UID.txt")
uid <- c()
url <- c()
for( i in 1:nrow(c.sample)){
url[i] <- paste("http://www.progenetix.org/cgi-bin/api.cgi?db=arraymap&uid_m=", c.sample[[1]][i],"&output=matrix", sep="")
}

pg.matrix.loader <- function(url){
	pg.frame <- read.table(url(url), header=T, sep="\t", na="NA")
	pg.matrix <- as.matrix(pg.frame)
	return(pg.matrix)
} 

write.to.file <- function(mat, gsm){
write.table(mat, file=paste("~/Desktop/Prisni/RenalCa", gsm, sep="/"), quote=FALSE, row.names=T, col.names=T, sep="\t")
print("written to file...")
print(gsm)
}


pg.frame.final <- colnames(pg.matrix.loader(url[1]))
for( i in 1:nrow(c.sample)){
url <- paste("http://www.progenetix.org/cgi-bin/api.cgi?db=arraymap&uid_m=", c.sample[[1]][i],"&output=matrix", sep="")
pg.frame.matrix <- pg.matrix.loader(url)
pg.frame.final <- rbind(pg.frame.final, pg.frame.matrix)
write.to.file(pg.frame.matrix, c.sample[[1]][i])
}

write.to.file(pg.frame.final, "gsm-all")

pgframe1 <- read.table(url("http://www.progenetix.org/cgi-bin/api.cgi?db=arraymap&uid_m=GSM302781&output=matrix"), header=T, sep="\t", na="NA")