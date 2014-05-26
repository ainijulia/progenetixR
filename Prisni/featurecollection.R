samples.read1 <- as.matrix(read.table("/Volumes/pgweb/Desktop/Prisni/CervCa/gsm-all", sep="\t", header=T, row.names=NULL))
samples.read2 <- as.matrix(read.table("/Volumes/pgweb/Desktop/Prisni/OvaCa/gsm-all", sep="\t", header=T, row.names=NULL))
samples.read3 <- as.matrix(read.table("/Volumes/pgweb/Desktop/Prisni/Melanoma/gsm-all", sep="\t", header=T, row.names=NULL))
samples.read4 <- as.matrix(read.table("/Volumes/pgweb/Desktop/Prisni/RenalCa/gsm-all", sep="\t", header=T, row.names=NULL))
samples.matrix <- rbind(samples.read1[-1,-1], samples.read2[-1,-1])
samples.matrix <- rbind(samples.matrix, samples.read3[-1,-1])
samples.matrix <- rbind(samples.matrix, samples.read4[-1,-1])

feature.names <- colnames(samples.matrix)
feature.names # column names of samples.matrix .. 3109

feature.names[3101] # diagnosis text

samples.names <- samples.matrix[,1] # all sample ids .. 1163
class.names <- samples.matrix[,3101]
#assigning numeric values to class.names
## class.id == 00001 -> CervCa
## class.id == 00010 -> OvaCa
## class.id == 00100 -> Melanoma
## class.id == 01000 -> RenCa
## class.id == 10000 -> not assigned any of the above classes

class.id <- c()
for(i in 1:length(class.names)){
if(grepl(class.names[i], "Cervical Cancer", ignore.case=TRUE)) {
  class.id[[i]] <- unlist(c(0,0,0,0,1))
}
else if(grepl(class.names[i], "Ovarian Carcinoma", ignore.case=TRUE)){
  class.id[[i]] <- unlist(c(0,0,0,1,0))
}
else if(grepl(class.names[i], "Melanoma", ignore.case=TRUE)){
  class.id[[i]] <- unlist(c(0,0,1,0,0))
}
else if(grepl(class.names[i], "Renal Carcinoma", ignore.case=TRUE)){
  class.id[[i]] <- unlist(c(0,1,0,0,0))
}
else {
  class.id[[i]] <- unlist(c(1,0,0,0,0))
}
}

class.id.matrix <- as.matrix(class.id[[1]])
for(j in 2:length(class.id)){
  class.id.matrix <- cbind(class.id.matrix, as.matrix(class.id[[j]]))
}

# ncol(class.id.matrix) = number of features 
# nrow(class.id.matrix) = state values to identify a cancer class type)


or.feature.matrix <- samples.matrix[,-1] # original samples vs segments matrix

or.feature.matrix <- or.feature.matrix[,1:3091] # 1163X3091 matrix; 3091 features for 1163 samples

#keep diagnosis text in the original sample feature matrix

new.class.feature.matrix <- cbind(class.id, or.feature.matrix) # class vs feature signature matrix

#compute feature versus class correlations .. pairwise correlations between each feature and class type

## cor(x,y) --> x <- feature(i) && y <- class.id.matrix[i,j] ( i -> one of the 5 vectors; j -> feature columns)
# number of features (total) = ncol(or.feature.matrix)
## X <- feature(i)
## Y <- class.id.matrix[i,j] ; i <- row number (class id locator), j <- column number (state swtich)
## such that corr.list[[i]] <- {c.1, c.2, c.3, c.4, c.5} 

corr.list1 <- c()
corr.list2 <- c()
corr.list3 <- c()
corr.list4 <- c()
corr.list5 <- c()
for(i in 1:ncol(or.feature.matrix)){
  
  X <- as.numeric(or.feature.matrix[,i])
  Y1 <- class.id.matrix[1,]
  Y2 <- class.id.matrix[2,]
  Y3 <- class.id.matrix[3,]
  Y4 <- class.id.matrix[4,]
  Y5 <- class.id.matrix[5,]
  
  corr.list1[[i]] <- cor(X,Y1) 
  corr.list2[[i]] <- cor(X,Y2) 
  corr.list3[[i]] <- cor(X,Y3) 
  corr.list4[[i]] <- cor(X,Y4) 
  corr.list5[[i]] <- cor(X,Y5) 
}

corr.list <- cbind(corr.list1, corr.list2)
corr.list <- cbind(corr.list, corr.list3)
corr.list <- cbind(corr.list, corr.list4)
corr.list <- cbind(corr.list, corr.list5)

# corr.list is a 5 column matrix repsenting correlations between features and the id switches
# assigning a threshold 'd'
## say, d == 0.15 (preserving 278 features)

abs.corr.list1 <- abs(corr.list1)
corr.extr1 <- abs.corr.list[abs.corr.list1>=0.15]

abs.corr.list2 <- abs(corr.list2)
corr.extr2 <- abs.corr.list[abs.corr.list2>=0.15]

abs.corr.list3 <- abs(corr.list3)
corr.extr3 <- abs.corr.list[abs.corr.list3>=0.15]

abs.corr.list4 <- abs(corr.list4)
corr.extr4 <- abs.corr.list[abs.corr.list4>=0.15]

abs.corr.list5 <- abs(corr.list5)
corr.extr5 <- abs.corr.list[abs.corr.list5>=0.15]


corr.extr <- cbind(corr.extr1, corr.extr2)
corr.extr <- cbind(corr.extr, corr.extr3)
corr.extr <- cbind(corr.extr, corr.extr4)
corr.extr <- cbind(corr.extr, corr.extr5)

# 'corr.feature.class' is the filtered list of correlation coefficients for pairwise correlations between features and classes

order.corr.list <- sort(corr.feature.class, decreasing=TRUE) # order the filtered correlation coefficients in descending order



 