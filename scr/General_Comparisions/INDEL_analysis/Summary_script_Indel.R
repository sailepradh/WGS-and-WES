##############################################################################
## Genome BWA summary report #################################################
##############################################################################

setwd("~/Desktop/New_analysis/INDEL_analysis/individual_sample_genomes/")
temp =list.files (pattern ="*g.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
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
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Summary_Gen_Indel =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Summary_Gen_Indel,"/Users/salendrapradh/WGS-and-WES/result/General Comparison/INDEL/Summary_Gen_indel.tsv",
    sep="\t", quote =F , row.names = F)


#############################################################################
## Exome-BWA summary report #################################################
#############################################################################

setwd("~/Desktop/New_analysis/INDEL_analysis/individual_exomes_BWA/")
temp =list.files (pattern ="*BWA.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
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
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Summary_Exo_BWA_Indel =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Summary_Exo_BWA_Indel,"/Users/salendrapradh/WGS-and-WES/result/General Comparison/INDEL/Summary_Exo_BWA_Indel.tsv",
            sep="\t", quote =F , row.names = F)

#########################################################################################
## Exome-BWT summary report ########################################################
###################################################################################
setwd("~/Desktop/New_analysis/INDEL_analysis/individual_exomes_BWT/")
temp =list.files (pattern ="*BWT.txt")
myfiles = lapply (temp,  read.table, sep = " ", stringsAsFactor=FALSE) 
sample <- sapply(strsplit(temp, "_"), "[",1)
test = vector ("list", length(sample))
names(test) <- sample
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
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Summary_Exo_BWT_Indel =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Summary_Exo_BWT_Indel,"/Users/salendrapradh/WGS-and-WES/result/General Comparison/INDEL/Summary_Exo_BWT_Indel.tsv",
            sep="\t", quote =F , row.names = F)




