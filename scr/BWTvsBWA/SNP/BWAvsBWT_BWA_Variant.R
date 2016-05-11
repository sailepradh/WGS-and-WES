########################################################################
################## Exome - BWA Analysis ################################ 
################## Reference (0/0) Homo Analysis #############################

setwd("/Users/salendrapradh/Documents/BWAvsBWT/WES_BWA/Variant_list_WES_BWA/")
temp =list.files (pattern ="*whomo_WES_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_Ref = vector ("list", length(sample))
names(test_Ref) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96)
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)
number = numeric(length=96)

for (i in 1:96){
  int =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  CD_ALT_2[is.na(CD_ALT_2)] <- 0
  
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_ALT)
  test_Ref[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

par (mfrow =c(2,2))
hist (mean_CD, breaks = 20)
#hist(as.numeric(Reference_Homo_BWAvsBWT_BWA_Summary[,2]), breaks=20)
hist (mean_GQ, breaks = 20)
hist (mean_MRR, breaks =20)

Int_1 <- myfiles$S0156[,3][which(test$S0156[,4] > 150)]
Int_1
#Int_2<-S0156_whomo_g[which(GQ ==80)]
#Int_3 <- S0156_whomo_g[which(is.nan(MRR))]
#Int_4 <- S0156_whomo_g[which(MRR > 0.00)]

Reference_ExoVsExo_BWA_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Reference_ExoVsExo_BWA_Summary,
            "/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Reference_ExoVsExo_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Alternate (0/1) Heterozugous Genome Analysis ########################

setwd("/Users/salendrapradh/Documents/BWAvsBWT/WES_BWA/Variant_list_WES_BWA//")
temp =list.files (pattern ="*het_WES_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_het = vector ("list", length(sample))
names(test_het) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96)
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)
number = numeric(length=96)

for (i in 1:96){
  int =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  CD_ALT_2[is.na(CD_ALT_2)] <- 0
  
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_ALT)
  test_het[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Het_ExoVsExo_BWA_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Het_ExoVsExo_BWA_Summary,
            "/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Het_ExoVsExo_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant (1/1) Homozygous Genome Analysis ############################

setwd("/Users/salendrapradh/Documents/BWAvsBWT/WES_BWA/Variant_list_WES_BWA/")
temp =list.files (pattern ="*mhomo_WES_BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test_homo = vector ("list", length(sample))
names(test_homo) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=96)
mean_GQ = numeric(length=96)
mean_MRR = numeric(length=96)
number = numeric(length=96)

for (i in 1:96){
  int =  sapply(strsplit((myfiles[[i]][[2]]), ":"),"[",2)
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  CD_ALT_2[is.na(CD_ALT_2)] <- 0
  
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[2]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_ALT)
  test_homo[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Mut_ExoVsExo_BWA_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Mut_ExoVsExo_BWA_Summary,
            "/Users/salendrapradh/WGS-and-WES/result/BWAvsBWT/Mut_ExoVsExo_BWA_Summary.tsv",
            sep="\t", quote =F , row.names = F)


par (mfrow =c(2,2))
hist (mean_CD, breaks = 20)
#hist(as.numeric(Reference_Homo_BWAvsBWT_BWA_Summary[,2]), breaks=20)
hist (mean_GQ, breaks = 20)
hist (mean_MRR, breaks =20)



hist(as.numeric(Summary_BWA_Exome[,2]), breaks=20, 
     xlab = "mean_CD", ylab = "Counts", main = NULL,
     xlim =c(0,150), ylim = c(0,30), col=rgb(1,1,0,0.25)) #rgb(0.25,0.2,0,0.5))

hist (as.numeric(Reference_Homo_BWAvsBWT_BWA_Summary[,2]), breaks=10, col=rgb(0,1,1,0.25),add = T) #rgb(0,0.25,0.05)

hist (as.numeric(Het_Exome_BWAvsBWT_BWT_Summary[,2]), breaks=20, col=rgb(1,0,1,0.25),add = T) #rgb(0,1,0,0.5)
hist (as.numeric(Homo_Exome_BWAvsBWT_BWA_Summary[,2]), breaks=20, col=rgb(1,0.5,0,0.25),add = T)





library("ggplot2")
qplot(as.data.frame(Het_Exome_BWAvsBWT_BWT_Summary[,2]), data =as.data.frame(Het_Exome_BWAvsBWT_BWT_Summary), geom = 'histogram')




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
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
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