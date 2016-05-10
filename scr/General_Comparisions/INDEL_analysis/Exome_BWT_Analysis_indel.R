########################################################################
################## Exome - BWT Analysis ################################ 
################## Reference (0/0) Homo Analysis #############################

setwd("~/Desktop/New_analysis/INDEL_analysis/individual_exomes_BWT/Variant_list")
temp =list.files (pattern ="*whomo_BWT.txt")
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
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  CD_ALT_2[is.na(CD_ALT_2)] <- 0
  CD_ALT_3 <- as.numeric(sapply(strsplit(int, ","), "[",4))
  CD_ALT_3[is.na(CD_ALT_3)] <- 0
  CD_ALT_4 <- as.numeric(sapply(strsplit(int, ","), "[",5))
  CD_ALT_4[is.na(CD_ALT_4)] <- 0
  CD_ALT_5 <- as.numeric(sapply(strsplit(int, ","), "[",6))
  CD_ALT_5[is.na(CD_ALT_5)] <- 0
  CD_ALT_6 <- as.numeric(sapply(strsplit(int, ","), "[",7))
  CD_ALT_6[is.na(CD_ALT_6)] <- 0
  
  
  MRR_1 <- round((CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6)/(CD_Ref+CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_Ref)
  # test_homo[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ))
  test_Ref[[i]] = cbind (CD_Ref,CD_ALT,CD_ALT_2,CD_ALT_3,CD_ALT_4,CD_ALT_5,CD_ALT_6,MRR,CD,GQ)
}

#hist (CD, breaks = 100)
#hist (GQ, breaks = 100)
#hist (MRR, breaks =100)

#Int_1 <- S0156_whomo_g[which(CD > 60)]
#Int_2<-S0156_whomo_g[which(GQ ==80)]
#Int_3 <- S0156_whomo_g[which(is.nan(MRR))]
#Int_4 <- S0156_whomo_g[which(MRR > 0.00)]

Ref_Exo_BWT_Indel_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Ref_Exo_BWT_Indel_Summary,"/Users/salendrapradh/WGS-and-WES/result/General Comparison/INDEL/Ref_Exo_BWT_Indel_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Alternate (0/1) Heterozugous Genome Analysis ########################

setwd("~/Desktop/New_analysis/INDEL_analysis/individual_exomes_BWT/Variant_list")
temp =list.files (pattern ="*het_BWT.txt")
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
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  CD_ALT_2[is.na(CD_ALT_2)] <- 0
  CD_ALT_3 <- as.numeric(sapply(strsplit(int, ","), "[",4))
  CD_ALT_3[is.na(CD_ALT_3)] <- 0
  CD_ALT_4 <- as.numeric(sapply(strsplit(int, ","), "[",5))
  CD_ALT_4[is.na(CD_ALT_4)] <- 0
  CD_ALT_5 <- as.numeric(sapply(strsplit(int, ","), "[",6))
  CD_ALT_5[is.na(CD_ALT_5)] <- 0
  CD_ALT_6 <- as.numeric(sapply(strsplit(int, ","), "[",7))
  CD_ALT_6[is.na(CD_ALT_6)] <- 0
  
  
  MRR_1 <- round((CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6)/(CD_Ref+CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_Ref)
  test_het[[i]] =  cbind (CD_Ref,CD_ALT,CD_ALT_2,CD_ALT_3,CD_ALT_4,CD_ALT_5,CD_ALT_6,MRR,CD,GQ)
  #test_het[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Het_Exo_BWT_Indel_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Het_Exo_BWT_Indel_Summary,"/Users/salendrapradh/WGS-and-WES/result/General Comparison/INDEL/Het_Exo_BWT_Indel_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant (1/1) Homozygous Genome Analysis ########################

setwd("~/Desktop/New_analysis/INDEL_analysis/individual_exomes_BWT/Variant_list")
temp =list.files (pattern ="*mhomo_BWT.txt")
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
  int =  sapply(strsplit((myfiles[[i]][[6]]), ":"),"[",2)
  
  CD_Ref <- as.numeric(sapply(strsplit(int, ","), "[",1))
  CD_ALT <- as.numeric(sapply(strsplit(int, ","), "[",2))
  CD_ALT_2 <- as.numeric(sapply(strsplit(int, ","), "[",3))
  CD_ALT_2[is.na(CD_ALT_2)] <- 0
  CD_ALT_3 <- as.numeric(sapply(strsplit(int, ","), "[",4))
  CD_ALT_3[is.na(CD_ALT_3)] <- 0
  CD_ALT_4 <- as.numeric(sapply(strsplit(int, ","), "[",5))
  CD_ALT_4[is.na(CD_ALT_4)] <- 0
  CD_ALT_5 <- as.numeric(sapply(strsplit(int, ","), "[",6))
  CD_ALT_5[is.na(CD_ALT_5)] <- 0
  CD_ALT_6 <- as.numeric(sapply(strsplit(int, ","), "[",7))
  CD_ALT_6[is.na(CD_ALT_6)] <- 0
  
  
  MRR_1 <- round((CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6)/(CD_Ref+CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2+CD_ALT_3+CD_ALT_4+CD_ALT_5+CD_ALT_6), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_Ref)
  test_het[[i]] =  cbind (CD_Ref,CD_ALT,CD_ALT_2,CD_ALT_3,CD_ALT_4,CD_ALT_5,CD_ALT_6,MRR,CD,GQ)
  #test_het[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Homo_Exo_BWT_Indel_Summary =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Homo_Exo_BWT_Indel_Summary ,"/Users/salendrapradh/WGS-and-WES/result/General Comparison/INDEL/Homo_Exo_BWT_Indel_Summary.tsv",
            sep="\t", quote =F , row.names = F)





##################################################################################
############ Mutant homo second allele(2/2) case Genome Analysis #################

setwd("~/Desktop/New_analysis/individual_sample_exome_BWT/Variant_list")
temp =list.files (pattern ="*mhomo_second_allele_BWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=92)  #Since four file did't have this kind of combination so 95 samples were taken
mean_GQ = numeric(length=92)
mean_MRR = numeric(length=92)

for (i in 1:92){
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

Mut_Second_Allele_Exome_BWT_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Mut_Second_Allele_Exome_BWT_Summary,"~/Desktop/Compared/New_analysis/Mut_Second_Allele_Exome_BWT_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant het second allele(0/2) case Genome Analysis ##################
setwd("~/Desktop/New_analysis/individual_sample_exome_BWT/Variant_list")
temp =list.files (pattern ="*het_second_allele_BWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
names(myfiles) <- sample
mean_CD= numeric(length=95)  # one file missing
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

Het_Second_Allele_Exome_BWT_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Het_Second_Allele_Exome_BWT_Summary ,"~/Desktop/Compared/New_analysis/Het_Second_Allele_Exome_BWT_Summary.tsv",
            sep="\t", quote =F , row.names = F)

##################################################################################
############ Mutant het second allele(1/2) case Genome Analysis ##################
setwd("~/Desktop/New_analysis/individual_sample_exome_BWT/Variant_list")
temp =list.files (pattern ="*het_first_second_allele_BWT.txt")
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

Het_First_Second_Allele_Exome_BWT_Summary =cbind(sample,mean_CD,mean_GQ,mean_MRR)
write.table(Het_First_Second_Allele_Exome_BWT_Summary,"~/Desktop/Compared/New_analysis/Het_First_Second_Allele_Exome_BWT_Summary.tsv",
            sep="\t", quote =F , row.names = F)

###########################################################