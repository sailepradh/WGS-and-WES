########################################################################
################## Both BWAvsBWT Analysis ##############################
## (BWA vs BWT) summary report #########################################

setwd("/Users/salendrapradh/Documents/BWAvsBWT/Both/")
temp =list.files (pattern ="*both_BWAvsBWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96)  
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)
number = numeric(length=96)

## test with 2 for BWT and 3 for BWA 
for (i in 1:96){
  int_BWT =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",1))
  CD_ALT_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",2))
  CD_ALT_2_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",3))
  CD_ALT_2_BWT[is.na(CD_ALT_2_BWT)] <- 0
  
  int_BWA =  sapply(strsplit((myfiles[[i]][[3]]), ":"),"[",2)
  CD_Ref_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",1))
  CD_ALT_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",2))
  CD_ALT_2_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",3))
  CD_ALT_2_BWA[is.na(CD_ALT_2_BWA)] <- 0
  
  MRR_1_BWT <- round((CD_ALT_BWT+CD_ALT_2_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_2_BWT <- round((CD_Ref_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_BWT <- pmin(MRR_1_BWT,MRR_2_BWT)
  
  MRR_1_BWA <- round((CD_ALT_BWA+CD_ALT_2_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_2_BWA <- round((CD_Ref_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_BWA <- pmin(MRR_1_BWA,MRR_2_BWA)
  
  mean_MRR[i] = paste(round(mean(as.numeric(MRR_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(MRR_BWA),na.rm=TRUE), digit =2), sep =":")
  
  CD_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  CD_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",3))
  mean_CD [i] <- paste(round(mean(as.numeric(CD_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(CD_BWA),na.rm=TRUE), digit =2), sep =":")
  
  GQ_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  GQ_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",4))
  mean_GQ [i] <- paste(round(mean(as.numeric(GQ_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(GQ_BWA),na.rm=TRUE), digit =2), sep =":")

  number[i] <- length(CD_ALT_BWT)
  test[[i]][[1]] = cbind (CD_Ref_BWT,CD_ALT_BWT,MRR_BWT,CD_BWT,GQ_BWT)
  test[[i]][[2]] = cbind (CD_Ref_BWA,CD_ALT_BWA,MRR_BWA,CD_BWA,GQ_BWA)
}


Summary_Both_BWABWT = cbind (sample,number, mean_CD, mean_GQ, mean_MRR)
write.table(Summary_Both_BWABWT,"/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Summary_Both_BWABWT.tsv",
            sep="\t", quote =F , row.names = F)


################## Reference (0/0) Homo Analysis #############################
setwd("/Users/salendrapradh/Documents/BWAvsBWT/Both/Con_Dis/")
temp =list.files (pattern ="*whomo_both_BWAvsBWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_Ref = vector ("list", length(sample))
names(test_Ref) <- sample
names(myfiles) <- sample

mean_CD= numeric(length=96)
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)
number = numeric(length=96)


## test with 2 for BWT and 3 for BWA 
for (i in 1:96){
  int_BWT =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",1))
  CD_ALT_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",2))
  CD_ALT_2_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",3))
  CD_ALT_2_BWT[is.na(CD_ALT_2_BWT)] <- 0
  
  int_BWA =  sapply(strsplit((myfiles[[i]][[3]]), ":"),"[",2)
  CD_Ref_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",1))
  CD_ALT_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",2))
  CD_ALT_2_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",3))
  CD_ALT_2_BWA[is.na(CD_ALT_2_BWA)] <- 0
  
  MRR_1_BWT <- round((CD_ALT_BWT+CD_ALT_2_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_2_BWT <- round((CD_Ref_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_BWT <- pmin(MRR_1_BWT,MRR_2_BWT)
  
  MRR_1_BWA <- round((CD_ALT_BWA+CD_ALT_2_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_2_BWA <- round((CD_Ref_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_BWA <- pmin(MRR_1_BWA,MRR_2_BWA)
  
  mean_MRR[i] = paste(round(mean(as.numeric(MRR_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(MRR_BWA),na.rm=TRUE), digit =2), sep =":")
  
  CD_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  CD_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",3))
  mean_CD [i] <- paste(round(mean(as.numeric(CD_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(CD_BWA),na.rm=TRUE), digit =2), sep =":")
  
  GQ_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  GQ_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",4))
  mean_GQ [i] <- paste(round(mean(as.numeric(GQ_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(GQ_BWA),na.rm=TRUE), digit =2), sep =":")
  
  number[i] <- length(CD_ALT_BWT)
  test_Ref[[i]][[1]] = cbind (CD_Ref_BWT,CD_ALT_BWT,MRR_BWT,CD_BWT,GQ_BWT)
  test_Ref[[i]][[2]] = cbind (CD_Ref_BWA,CD_ALT_BWA,MRR_BWA,CD_BWA,GQ_BWA)
}

# par (mfrow =c(2,2))
# hist (mean_CD, breaks = 20)
# #hist(as.numeric(Reference_Homo_BWAvsBWT_BWA_Summary[,2]), breaks=20)
# hist (mean_GQ, breaks = 20)
# hist (mean_MRR, breaks =20)
# 
# Int_1 <- myfiles$S0156[,3][which(test$S0156[,4] > 150)]
# Int_1
# #Int_2<-S0156_whomo_g[which(GQ ==80)]
# #Int_3 <- S0156_whomo_g[which(is.nan(MRR))]
# #Int_4 <- S0156_whomo_g[which(MRR > 0.00)]
# 

Reference_BothBWAvsBWT_Summary = cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Reference_BothBWAvsBWT_Summary,
            "/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Reference_BothBWAvsBWT_Summary.tsv",
            sep="\t", quote =F , row.names = F)


##################################################################################
############ Alternate (0/1) Heterozugous Genome Analysis ########################

setwd("/Users/salendrapradh/Documents/BWAvsBWT/Both/Con_Dis/")
temp =list.files (pattern ="*het_both_BWAvsBWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_het = vector ("list", length(sample))
names(test_het) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96)
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)
number = numeric(length=96)

## test with 2 for BWT and 3 for BWA 

for (i in 1:96){
  int_BWT =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",1))
  CD_ALT_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",2))
  CD_ALT_2_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",3))
  CD_ALT_2_BWT[is.na(CD_ALT_2_BWT)] <- 0
  
  int_BWA =  sapply(strsplit((myfiles[[i]][[3]]), ":"),"[",2)
  CD_Ref_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",1))
  CD_ALT_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",2))
  CD_ALT_2_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",3))
  CD_ALT_2_BWA[is.na(CD_ALT_2_BWA)] <- 0
  
  MRR_1_BWT <- round((CD_ALT_BWT+CD_ALT_2_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_2_BWT <- round((CD_Ref_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_BWT <- pmin(MRR_1_BWT,MRR_2_BWT)
  
  MRR_1_BWA <- round((CD_ALT_BWA+CD_ALT_2_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_2_BWA <- round((CD_Ref_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_BWA <- pmin(MRR_1_BWA,MRR_2_BWA)
  
  mean_MRR[i] = paste(round(mean(as.numeric(MRR_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(MRR_BWA),na.rm=TRUE), digit =2), sep =":")
  
  CD_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  CD_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",3))
  mean_CD [i] <- paste(round(mean(as.numeric(CD_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(CD_BWA),na.rm=TRUE), digit =2), sep =":")
  
  GQ_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  GQ_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",4))
  mean_GQ [i] <- paste(round(mean(as.numeric(GQ_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(GQ_BWA),na.rm=TRUE), digit =2), sep =":")
  
  number[i] <- length(CD_ALT_BWT)
  test_het[[i]][[1]] = cbind (CD_Ref_BWT,CD_BWT,MRR_BWT,CD_BWT,GQ_BWT)
  test_het[[i]][[2]] = cbind (CD_Ref_BWA,CD_BWA,MRR_BWA,CD_BWA,GQ_BWA)
}


Het_BothExoVsExo_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Het_BothExoVsExo_Summary,
            "/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Het_BothExoVsExo_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant (1/1) Homozygous Genome Analysis ########################

setwd("/Users/salendrapradh/Documents/BWAvsBWT/Both/Con_Dis/")
temp =list.files (pattern ="*mhomo_both_BWAvsBWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_homo = vector ("list", length(sample))
names(test_homo) <- sample
names(myfiles) <- sample

## test with 2 for BWT and 3 for BWA 

for (i in 1:96){
  int_BWT =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",1))
  CD_ALT_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",2))
  CD_ALT_2_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",3))
  CD_ALT_2_BWT[is.na(CD_ALT_2_BWT)] <- 0
  
  int_BWA =  sapply(strsplit((myfiles[[i]][[3]]), ":"),"[",2)
  CD_Ref_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",1))
  CD_ALT_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",2))
  CD_ALT_2_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",3))
  CD_ALT_2_BWA[is.na(CD_ALT_2_BWA)] <- 0
  
  MRR_1_BWT <- round((CD_ALT_BWT+CD_ALT_2_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_2_BWT <- round((CD_Ref_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
  MRR_BWT <- pmin(MRR_1_BWT,MRR_2_BWT)
  
  MRR_1_BWA <- round((CD_ALT_BWA+CD_ALT_2_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_2_BWA <- round((CD_Ref_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
  MRR_BWA <- pmin(MRR_1_BWA,MRR_2_BWA)
  
  mean_MRR[i] = paste(round(mean(as.numeric(MRR_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(MRR_BWA),na.rm=TRUE), digit =2), sep =":")
  
  CD_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  CD_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",3))
  mean_CD [i] <- paste(round(mean(as.numeric(CD_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(CD_BWA),na.rm=TRUE), digit =2), sep =":")
  
  GQ_BWT <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  GQ_BWA <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[3]]), ":"),"[",4))
  mean_GQ [i] <- paste(round(mean(as.numeric(GQ_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(GQ_BWA),na.rm=TRUE), digit =2), sep =":")
  
  number[i] <- length(CD_ALT_BWT)
  test_homo[[i]][[1]] = cbind (CD_Ref_BWT,CD_BWT,MRR_BWT,CD_BWT,GQ_BWT)
  test_homo[[i]][[2]] = cbind (CD_Ref_BWA,CD_BWA,MRR_BWA,CD_BWA,GQ_BWA)
}

Mut_BothExoVsExo_Summary=cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Mut_BothExoVsExo_Summary,
            "/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Mut_BothExoVsExo_Summary.tsv",
            sep="\t", quote =F , row.names = F)


##############Discordant Variant Analysis##############################
#######################################################################
setwd("~/Desktop/New_analysis/WES_BWA-WES_BWT/Both/Con_Dis/")
temp =list.files (pattern ="*discordant_both_BWAvsBWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_dis = vector ("list", length(sample))
names(test_dis) <- sample
names(myfiles) <- sample

for (i in 1:96){
  test_dis[[i]][[1]] <- myfiles[[i]][(grepl("0/0", myfiles[[i]][2]$V2) & (grepl("0/1", myfiles[[i]][3]$V3))),]
  test_dis[[i]][[2]] <- myfiles[[i]][(grepl("0/0", myfiles[[i]][2]$V2) & (grepl("1/1", myfiles[[i]][3]$V3))),]
  test_dis[[i]][[3]] <- myfiles[[i]][(grepl("0/1", myfiles[[i]][2]$V2) & (grepl("0/0", myfiles[[i]][3]$V3))),]
  test_dis[[i]][[4]] <- myfiles[[i]][(grepl("0/1", myfiles[[i]][2]$V2) & (grepl("1/1", myfiles[[i]][3]$V3))),]
  test_dis[[i]][[5]] <- myfiles[[i]][(grepl("1/1", myfiles[[i]][2]$V2) & (grepl("0/0", myfiles[[i]][3]$V3))),]
  test_dis[[i]][[6]] <- myfiles[[i]][(grepl("1/1", myfiles[[i]][2]$V2) & (grepl("0/1", myfiles[[i]][3]$V3))),]
}


## test with 2 for BWT and 3 for BWA 

test_BWTBWA = vector ("list", length(sample))
names(test_BWTBWA ) <- sample

mean_CD = data.frame(matrix(NA, nrow=96, ncol=6))
rownames(mean_CD)<-sample
colnames(mean_CD)<- c("Ref_BWT:Het_BWA", "Ref_BWT:Hom_BWA","Het_BWT:Ref_BWA","Het_BWT:Homo_BWA","Homo_BWT:Ref_BWA","Homo_BWT:Het_BWA") 


mean_GQ = data.frame(matrix(NA, nrow=96, ncol=6))
rownames(mean_GQ)<-sample
colnames(mean_GQ)<- c("Ref_BWT:Het_BWA", "Ref_BWT:Hom_BWA","Het_BWT:Ref_BWA","Het_BWT:Homo_BWA","Homo_BWT:Ref_BWA","Homo_BWT:Het_BWA") 

mean_MRR = data.frame(matrix(NA, nrow=96, ncol=6))
rownames(mean_MRR)<-sample
colnames(mean_MRR)<- c("Ref_BWT:Het_BWA", "Ref_BWT:Hom_BWA","Het_BWT:Ref_BWA","Het_BWT:Homo_BWA","Homo_BWT:Ref_BWA","Homo_BWT:Het_BWA") 

number = data.frame(matrix(NA, nrow=96, ncol=6))
rownames(number)<-sample
colnames(number)<- c("Ref_BWT:Het_BWA", "Ref_BWT:Hom_BWA","Het_BWT:Ref_BWA","Het_BWT:Homo_BWA","Homo_BWT:Ref_BWA","Homo_BWT:Het_BWA") 


for (i in 1:96){
    for (j in 1:6){
      
      int_BWT =  as.character(sapply(strsplit((test_dis[[i]][[j]]$V2), ":"),"[",2))
      CD_Ref_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",1))
      CD_ALT_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",2))
      CD_ALT_2_BWT <- as.numeric(sapply(strsplit(int_BWT, ","), "[",3))
      CD_ALT_2_BWT[is.na(CD_ALT_2_BWT)] <- 0
      
      
      int_BWA =  as.character(sapply(strsplit((test_dis[[i]][[j]]$V3), ":"),"[",2))
      CD_Ref_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",1))
      CD_ALT_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",2))
      CD_ALT_2_BWA <- as.numeric(sapply(strsplit(int_BWA, ","), "[",3))
      CD_ALT_2_BWA[is.na(CD_ALT_2_BWA)] <- 0
      
      MRR_1_BWT <- round((CD_ALT_BWT+CD_ALT_2_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
      MRR_2_BWT <- round((CD_Ref_BWT)/(CD_Ref_BWT+CD_ALT_BWT+CD_ALT_2_BWT), digits =2)
      MRR_BWT <- pmin(MRR_1_BWT,MRR_2_BWT)
      
      MRR_1_BWA <- round((CD_ALT_BWA+CD_ALT_2_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
      MRR_2_BWA <- round((CD_Ref_BWA)/(CD_Ref_BWA+CD_ALT_BWA+CD_ALT_2_BWA), digits =2)
      MRR_BWA <- pmin(MRR_1_BWA,MRR_2_BWA)
      
      mean_MRR[i] = paste(round(mean(as.numeric(MRR_BWT),na.rm=TRUE), digit =2), round(mean(as.numeric(MRR_BWA),na.rm=TRUE), digit =2), sep =":")
      
      CD_BWT <- as.numeric(sapply(strsplit(as.character (test_dis[[i]][[j]]$V2), ":"),"[",3))
      CD_BWA <- as.numeric(sapply(strsplit(as.character (test_dis[[i]][[j]]$V3), ":"),"[",3))
      mean_CD [i,j] <- paste(round( mean (as.numeric(CD_BWT),na.rm=TRUE), digits =2),round( mean (as.numeric(CD_BWA),na.rm=TRUE), digits =2)  ,  sep =":")
      
      GQ_BWT <- as.numeric(sapply(strsplit(as.character (test_dis[[i]][[j]]$V2), ":"),"[",4))
      GQ_BWA <- as.numeric(sapply(strsplit(as.character (test_dis[[i]][[j]]$V3), ":"),"[",4))
      mean_GQ [i,j] <- paste(round(mean(as.numeric(GQ_BWT),na.rm=TRUE), digits =2), round(mean(as.numeric(GQ_BWA),na.rm=TRUE), digits =2), sep =":")
      
      number[i,j] <- length(CD_ALT_BWA)
      test_BWTBWA[[i]][[1]] = cbind (CD_Ref_BWT,CD_Ref_BWA,CD_ALT_BWT,CD_ALT_BWA,MRR_BWT,MRR_BWA,CD_BWT,CD_BWA,GQ_BWT,GQ_BWA)
    }
}

Summary_discordant = data.frame(matrix(NA, nrow=96, ncol=6))
rownames(Summary_discordant)<-sample
colnames(Summary_discordant)<- c("Ref_BWT:Het_BWA", "Ref_BWT:Hom_BWA","Het_BWT:Ref_BWA","Het_BWT:Homo_BWA","Homo_BWT:Ref_BWA","Homo_BWT:Het_BWA") 
Summary_discordant$`Ref_BWT:Het_BWA`<- paste(number$`Ref_BWT:Het_BWA`,mean_CD$`Ref_BWT:Het_BWA`, mean_GQ$`Ref_BWT:Het_BWA`,mean_MRR$`Ref_BWT:Het_BWA`, sep=",")
Summary_discordant$`Ref_BWT:Hom_BWA`<- paste(number$`Ref_BWT:Hom_BWA`,mean_CD$`Ref_BWT:Hom_BWA`, mean_GQ$`Ref_BWT:Hom_BWA`,mean_MRR$`Ref_BWT:Hom_BWA`, sep=",")
Summary_discordant$`Het_BWT:Ref_BWA`<- paste(number$`Het_BWT:Ref_BWA`,mean_CD$`Het_BWT:Ref_BWA`, mean_GQ$`Het_BWT:Ref_BWA`,mean_MRR$`Het_BWT:Ref_BWA`, sep=",")
Summary_discordant$`Het_BWT:Homo_BWA`<- paste(number$`Het_BWT:Homo_BWA`,mean_CD$`Het_BWT:Homo_BWA`, mean_GQ$`Het_BWT:Homo_BWA`,mean_MRR$`Het_BWT:Homo_BWA`, sep=",")
Summary_discordant$`Homo_BWT:Ref_BWA`<- paste(number$`Homo_BWT:Ref_BWA`,mean_CD$`Homo_BWT:Ref_BWA`, mean_GQ$`Homo_BWT:Ref_BWA`,mean_MRR$`Homo_BWT:Ref_BWA`, sep=",")
Summary_discordant$`Homo_BWT:Het_BWA`<- paste(number$`Homo_BWT:Het_BWA`,mean_CD$`Homo_BWT:Het_BWA`, mean_GQ$`Homo_BWT:Het_BWA`,mean_MRR$`Homo_BWT:Het_BWA`, sep=",")

write.table(Summary_discordant,
            "~/Desktop/New_analysis/WES_BWA-WES_BWT/Homo_Exome_BothBWAvsBWT_BWTBWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)



# par (mfrow =c(2,2))
# hist (mean_CD, breaks = 20)
# #hist(as.numeric(Reference_Homo_BWAvsBWT_BWA_Summary[,2]), breaks=20)
# hist (mean_GQ, breaks = 20)
# hist (mean_MRR, breaks =20)
# 
# 
# 
# hist(as.numeric(Summary_BWA_Exome[,2]), breaks=20, 
#      xlab = "mean_CD", ylab = "Counts", main = NULL,
#      xlim =c(0,150), ylim = c(0,30), col=rgb(1,1,0,0.25)) #rgb(0.25,0.2,0,0.5))
# 
# hist (as.numeric(Reference_Homo_BWAvsBWT_BWA_Summary[,2]), breaks=10, col=rgb(0,1,1,0.25),add = T) #rgb(0,0.25,0.05)
# 
# hist (as.numeric(Het_Exome_BWAvsBWT_BWT_Summary[,2]), breaks=20, col=rgb(1,0,1,0.25),add = T) #rgb(0,1,0,0.5)
# hist (as.numeric(Homo_Exome_BWAvsBWT_BWA_Summary[,2]), breaks=20, col=rgb(1,0.5,0,0.25),add = T)
# 
# 
# library("ggplot2")
# qplot(as.data.frame(Het_Exome_BWAvsBWT_BWT_Summary[,2]), data =as.data.frame(Het_Exome_BWAvsBWT_BWT_Summary), geom = 'histogram')





###########################Not done############################
##################################################################################
############ Mutant homo second allele(2/2) case Genome Analysis ########################
setwd("~/Desktop/New_analysis/individual_sample_exome_BWA/Variant_list")
temp =list.files (pattern ="*mhomo_second_allele_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=95)  #Since one file did't have this kind of combination so 95 samples were taken
mean_GQ = numeric(length=95)
mean_MRR = numeric(length=95)

for (i in 1:95){
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  MRR <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Mut_Second_Allele_Exome_BWA_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Mut_Second_Allele_Exome_BWA_Summary,"~/Desktop/Compared/New_analysis/Mut_Second_Allele_Exome_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant homo second allele(0/2) case Genome Analysis ########################

temp =list.files (pattern ="*mhomo_second_allele_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96)  
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)

for (i in 1:96){
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  MRR <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Mut_Second_Allele_Exome_BWA_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Mut_Second_Allele_Exome_BWA_Summary,"~/Desktop/Compared/New_analysis/Mut_Second_Allele_Exome_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant het second allele(0/2) case Genome Analysis ##################
setwd("~/Desktop/New_analysis/individual_sample_exome_BWA/Variant_list")
temp =list.files (pattern ="*het_second_allele_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96) 
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)

for (i in 1:96){
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  MRR <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Het_Second_Allele_Exome_BWA_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Het_Second_Allele_Exome_BWA_Summary ,"~/Desktop/Compared/New_analysis/Het_Second_Allele_Exome_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant het second allele(1/2) case Genome Analysis ##################
setwd("~/Desktop/New_analysis/individual_sample_exome_BWA/Variant_list")
temp =list.files (pattern ="*het_first_second_allele_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96) 
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)

for (i in 1:96){
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  MRR <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Het_First_Second_Allele_Exome_BWA_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Het_First_Second_Allele_Exome_BWA_Summary,"~/Desktop/Compared/New_analysis/Het_First_Second_Allele_Exome_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)

###########################################################







par (mfrow =c(2,2))
hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),1]), breaks=50,prob=TRUE,border="white", 
     ylim = c(0, 0.25),main ="Read Depth Exome Sequencing", xlab ="Approximate Read Depth", ylab="Density",
     cex.lab = 0.95, cex.main = 0.95)
for(i in 1:96){
  lines(density(as.numeric(test[[i]][,4]),na.rm=TRUE))}

hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),1]), breaks=50,prob=TRUE,border="white", ylim = c(0, 0.075),
     main ="Read Depth Genome Sequencing", xlab ="Approximate Read Depth", ylab="Density",
     cex.lab = 0.95, cex.main = 0.95)
#hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),2]), breaks=50,prob=TRUE,border="white")
for(i in 1:96){
  lines(density(as.numeric(test[[i]][,5]),na.rm=TRUE))}

hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),3]), breaks=50,prob=TRUE,border="white", ylim = c(0,20)
     ,main ="Read Depth Genome Sequencing", xlab ="Approximate Read Depth", ylab="Density",
     cex.lab = 0.95, cex.main = 0.95)
#hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),2]), breaks=50,prob=TRUE,border="white")
for(i in 1:96){
  lines(density(as.numeric(test[[i]][,3]),na.rm=TRUE))}
