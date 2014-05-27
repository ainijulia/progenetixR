
# @author prath  ----------------------------------------------------------

# @date 24.05.2014 --------------------------------------------------------

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

feature.names <- feature.names[2:3092] # contains feature names to be selected

feature.id <- c(1:3091)

feature.map <- cbind(feature.id, feature.names)

library(hash)

feature.hash <- hash(feature.id, feature.names)
# usage feature.hash$"1" returns "seg_chr1.0.999999"


samples.names <- samples.matrix[,1] # all sample ids .. 1163
samples.hash <- hash(1:1163, samples.names) #usage: samples.hash$"1" returns "GSM253202"

class.names <- samples.matrix[,3101]
class.hash <- hash(1:1163, class.names)

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
X <- as.numeric(or.feature.matrix[,1])
Y1 <- class.id.matrix[1,]
Y2 <- class.id.matrix[2,]
Y3 <- class.id.matrix[3,]
Y4 <- class.id.matrix[4,]
Y5 <- class.id.matrix[5,]

corr.list1 <- cor(X,Y1) 
corr.list2 <- cor(X,Y2) 
corr.list3 <- cor(X,Y3) 
corr.list4 <- cor(X,Y4) 
corr.list5 <- cor(X,Y5)

for(i in 2:ncol(or.feature.matrix)){
  
  X <- as.numeric(or.feature.matrix[,i])
  
  corr.list1 <- append(corr.list1, cor(X,Y1))
  corr.list2 <- append(corr.list2, cor(X,Y2)) 
  corr.list3 <- append(corr.list3, cor(X,Y3)) 
  corr.list4 <- append(corr.list4, cor(X,Y4)) 
  corr.list5 <- append(corr.list5, cor(X,Y5)) 
}

corr1.feature <- hash(feature.id, corr.list1)
corr2.feature <- hash(feature.id, corr.list2)
corr3.feature <- hash(feature.id, corr.list3)
corr4.feature <- hash(feature.id, corr.list4)
corr5.feature <- hash(feature.id, corr.list5)

corr.list <- cbind(corr.list1, corr.list2)
corr.list <- cbind(corr.list, corr.list3)
corr.list <- cbind(corr.list, corr.list4)
corr.list <- cbind(corr.list, corr.list5)

##########################
####> dim(corr.list)
####[1] 3091    5
####> c(max(corr.list1), min(corr.list1))
####[1]  0.3700273 -0.2576420
####> c(max(corr.list2), min(corr.list2))
####[1]  0.1008075 -0.1178245
####> c(max(corr.list3), min(corr.list3))
####[1]  0.4240914 -0.4393086
####> c(max(corr.list4), min(corr.list4))
####[1]  0.3634608 -0.3824388
####> c(max(corr.list5), min(corr.list5))
####[1]  0.2518153 -0.2524150
##########################

par(mfrow=c(2,3))
corrplot(corr.list)
plot(corr.list1)
plot(corr.list2)
plot(corr.list3)
plot(corr.list4)
plot(corr.list5)

# corr.list is a 5 column matrix repsenting correlations between features and the id switches
# assigning a threshold 'd'
## say, d == 0.10

abs.corr.list1 <- abs(corr.list1)
abs.corr.list2 <- abs(corr.list2)
abs.corr.list3 <- abs(corr.list3)
abs.corr.list4 <- abs(corr.list4)
abs.corr.list5 <- abs(corr.list5) ## SU(i,c)

abs.corr.list1 <- as.list(abs.corr.list1)
abs.corr.list2 <- as.list(abs.corr.list2)
abs.corr.list3 <- as.list(abs.corr.list3)
abs.corr.list4 <- as.list(abs.corr.list4)
abs.corr.list5 <- as.list(abs.corr.list5)


names(abs.corr.list1) <- feature.id
names(abs.corr.list2) <- feature.id
names(abs.corr.list3) <- feature.id
names(abs.corr.list4) <- feature.id
names(abs.corr.list5) <- feature.id

abs.corr1.feature <- hash(feature.id, abs.corr.list1)
abs.corr2.feature <- hash(feature.id, abs.corr.list2)
abs.corr3.feature <- hash(feature.id, abs.corr.list3)
abs.corr4.feature <- hash(feature.id, abs.corr.list4)
abs.corr5.feature <- hash(feature.id, abs.corr.list5)

corr.extr1 <- abs.corr.list[abs.corr.list1>=0.05] ## for cancer type 5 .. Cervical Cancer .. S_list 
corr.extr2 <- abs.corr.list[abs.corr.list2>=0.05] ## for cancer type 4 .. Ovarian Carcinoma 
corr.extr3 <- abs.corr.list[abs.corr.list3>=0.05] ## for cancer type 3 ..Melanoma 
corr.extr4 <- abs.corr.list[abs.corr.list4>=0.05] ## for cancer type 2 .. Renal Carcinoma 
corr.extr5 <- abs.corr.list[abs.corr.list5>=0.05] ## for cancer type 1 .. Unknown


corr.extr1 <- abs.corr.list1[abs.corr.list1>=0.10] ## for cancer type 1 .. Cervical Cancer .. S_list 
corr.extr2 <- abs.corr.list2[abs.corr.list2>=0.10] ## for cancer type 2 .. Ovarian Carcinoma  
corr.extr3 <- abs.corr.list3[abs.corr.list3>=0.10] ## for cancer type 3 .. Melanoma 
corr.extr4 <- abs.corr.list4[abs.corr.list4>=0.10] ## for cancer type 4 .. Renal Carcinoma  
corr.extr5 <- abs.corr.list5[abs.corr.list5>=0.10] ## for cancer type 5 .. Unknown 

selected.features1 <- names(corr.extr1)
selected.features2 <- names(corr.extr2)
selected.features3 <- names(corr.extr3)
selected.features4 <- names(corr.extr4)
selected.features5 <- names(corr.extr5)


c(max(abs.corr.list1), min(abs.corr.list1))
c(max(abs.corr.list2), min(abs.corr.list2))
c(max(abs.corr.list3), min(abs.corr.list3))
c(max(abs.corr.list4), min(abs.corr.list4))
c(max(abs.corr.list5), min(abs.corr.list5))

##> c(max(abs.corr.list1), min(abs.corr.list1))
##[1] 0.3700273 0.0000000 :::: 0.05
##> c(max(abs.corr.list2), min(abs.corr.list2))
##[1] 0.1178245 0.0000000 ::: 0.05
##> c(max(abs.corr.list3), min(abs.corr.list3))
##[1] 0.4393086375 0.0000327509 ::: 0.05
##> c(max(abs.corr.list4), min(abs.corr.list4))
##[1] 0.3824388 0.0000000 ::: 0.05
##> c(max(abs.corr.list5), min(abs.corr.list5))
##[1] 0.252415 0.000000 ::: 0.05

par(mfrow=c(2,3))
plot(corr.extr1)
plot(corr.extr2)
plot(corr.extr3)
plot(corr.extr4)
plot(corr.extr5)


corr.extr <- cbind(corr.extr1, corr.extr2)
corr.extr <- cbind(corr.extr, corr.extr3)
corr.extr <- cbind(corr.extr, corr.extr4)
corr.extr <- cbind(corr.extr, corr.extr5) 

# 'corr.feature.class' is the filtered list of correlation coefficients for pairwise correlations between features and classes

order.corr.list1 <- sort.int(as.numeric(corr.extr1), decreasing=TRUE, index.return=TRUE)
order.corr.list2 <- sort.int(as.numeric(corr.extr2), decreasing=TRUE, index.return=TRUE)
order.corr.list3 <- sort.int(as.numeric(corr.extr3), decreasing=TRUE, index.return=TRUE)
order.corr.list4 <- sort.int(as.numeric(corr.extr4), decreasing=TRUE, index.return=TRUE)
order.corr.list5 <- sort.int(as.numeric(corr.extr5), decreasing=TRUE, index.return=TRUE) # order the filtered correlation coefficients in descending order

index.list1 <- hash(1:length(selected.features1), selected.features1)
index.list2 <- hash(1:length(selected.features2), selected.features2)
index.list3 <- hash(1:length(selected.features3), selected.features3)
index.list4 <- hash(1:length(selected.features4), selected.features4)
index.list5 <- hash(1:length(selected.features5), selected.features5)

feature.id.names.list <- as.list(feature.hash)
feature.id.names.list <- as.matrix(feature.id.names.list)
feature.index <- rownames(feature.id.names.list)
feature.index <- cbind(feature.index, feature.id.names.list[,1])

index.list1 <- as.list(index.list1)
index.list1 <- as.matrix(index.list1)
row.l1 <- rownames(index.list1)
index.list1 <- cbind(row.l1, index.list1[,1]) 

feature.i1 <- c()
for(i in 1:nrow(index.list1)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(index.list1[i,2]) == as.numeric(feature.index[j,1])){
      feature.i1[i] <- feature.index[j,2] 
    }
  }
}
feature.i1 <- as.matrix(feature.i1) # features for list 1

index.feature.list1 <- cbind(index.list1, feature.i1[,1])


index.list2 <- as.list(index.list2)
index.list2 <- as.matrix(index.list2)
row.l2 <- rownames(index.list2)
index.list2 <- cbind(row.l2, index.list2[,1]) 
feature.i2 <- c()
for(i in 1:nrow(index.list2)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(index.list2[i,2]) == as.numeric(feature.index[j,1])){
      feature.i2[i] <- feature.index[j,2] 
    }
  }
}
feature.i2 <- as.matrix(feature.i2) # features for list 2

index.feature.list2 <- cbind(index.list2, feature.i2[,1])


index.list3 <- as.list(index.list3)
index.list3 <- as.matrix(index.list3)
row.l3 <- rownames(index.list3)
index.list3 <- cbind(row.l3, index.list3[,1]) 
feature.i3 <- c()
for(i in 1:nrow(index.list3)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(index.list3[i,2]) == as.numeric(feature.index[j,1])){
      feature.i3[i] <- feature.index[j,2] 
    }
  }
}
feature.i3 <- as.matrix(feature.i3) # features for list 3

index.feature.list3 <- cbind(index.list3, feature.i3[,1])


index.list4 <- as.list(index.list4)
index.list4 <- as.matrix(index.list4)
row.l4 <- rownames(index.list4)
index.list4 <- cbind(row.l4, index.list4[,1]) 
feature.i4 <- c()
for(i in 1:nrow(index.list4)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(index.list4[i,2]) == as.numeric(feature.index[j,1])){
      feature.i4[i] <- feature.index[j,2] 
    }
  }
}
feature.i4 <- as.matrix(feature.i4) # features for list 4

index.feature.list4 <- cbind(index.list4, feature.i4[,1])

index.list5 <- as.list(index.list5)
index.list5 <- as.matrix(index.list5)
row.l5 <- rownames(index.list5)
index.list5 <- cbind(row.l5, index.list5[,1]) 
feature.i5 <- c()
for(i in 1:nrow(index.list5)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(index.list5[i,2]) == as.numeric(feature.index[j,1])){
      feature.i5[i] <- feature.index[j,2] 
    }
  }
}
feature.i5 <- as.matrix(feature.i5) # features for list 5

index.feature.list5 <- cbind(index.list5, feature.i5[,1])


keys.list1 <- as.numeric(unlist(order.corr.list1[2]))
keys.list2 <- as.numeric(unlist(order.corr.list2[2]))
keys.list3 <- as.numeric(unlist(order.corr.list3[2]))
keys.list4 <- as.numeric(unlist(order.corr.list4[2]))
keys.list5 <- as.numeric(unlist(order.corr.list5[2]))

keys.list1 <- as.matrix(keys.list1)

ordered.feature1 <- c()
for(i in 1:nrow(keys.list1)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(keys.list1[i,1]) == as.numeric(feature.index[j,1])){
      ordered.feature1[i] <- feature.index[j,2] 
    }
  }
}
ordered.feature1 <- as.matrix(ordered.feature1) # ordered features for list 1 in descending order
ordered.key.feature1 <- cbind(keys.list1[,1], ordered.feature1[,1])

keys.list2 <- as.matrix(keys.list2)

ordered.feature2 <- c()
for(i in 1:nrow(keys.list2)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(keys.list2[i,1]) == as.numeric(feature.index[j,1])){
      ordered.feature2[i] <- feature.index[j,2] 
    }
  }
}
ordered.feature2 <- as.matrix(ordered.feature2) # ordered features for list 2 in descending order
ordered.key.feature2 <- cbind(keys.list2[,1], ordered.feature2[,1]) 

keys.list3 <- as.matrix(keys.list3)

ordered.feature3 <- c()
for(i in 1:nrow(keys.list3)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(keys.list3[i,1]) == as.numeric(feature.index[j,1])){
      ordered.feature3[i] <- feature.index[j,2] 
    }
  }
}
ordered.feature3 <- as.matrix(ordered.feature3) # ordered features for list 3 in descending order
ordered.key.feature3 <- cbind(keys.list3[,1], ordered.feature3[,1])

keys.list4 <- as.matrix(keys.list4)

ordered.feature4 <- c()
for(i in 1:nrow(keys.list4)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(keys.list4[i,1]) == as.numeric(feature.index[j,1])){
      ordered.feature4[i] <- feature.index[j,2] 
    }
  }
}
ordered.feature4 <- as.matrix(ordered.feature4) # ordered features for list 4 in descending order
ordered.key.feature4 <- cbind(keys.list4[,1], ordered.feature4[,1])

keys.list5 <- as.matrix(keys.list5)

ordered.feature5 <- c()
for(i in 1:nrow(keys.list5)){
  for(j in 1:nrow(feature.index)){
    if(as.numeric(keys.list5[i,1]) == as.numeric(feature.index[j,1])){
      ordered.feature5[i] <- feature.index[j,2] 
    }
  }
}
ordered.feature5 <- as.matrix(ordered.feature5) # ordered features for list 5 in descending order
ordered.key.feature5 <- cbind(keys.list5[,1], ordered.feature5[,1])
