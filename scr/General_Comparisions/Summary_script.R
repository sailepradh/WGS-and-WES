#########################################################################################
## Genome BWA summary report ########################################################

setwd("~/Documents/GenvsExo/individual_sample_genomes/")
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
  
  
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_ALT)
  #test[[i]] = cbind (CD_Ref,CD_ALT,CD_ALT_2,MRR,CD,GQ)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Summary_Genome =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Summary_Genome,"/Users/salendrapradh/WGS-and-WES/result/Summary_Genome.tsv",
            sep="\t", quote =F , row.names = F)

################# Visulaization stuffs ##############################
hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),4]),
     breaks=1000, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     ylim = c(0,18000), xlim = c(0,150),
     main ="Read Depth Genome Sequencing", xlab ="Approximate Read Depth", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 


for(i in 1:96){ 
  hist(as.numeric(test[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(0,0,1,1))}

for(i in 1:96){ 
  hist(as.numeric(test_Ref[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(0,1,0,1))}

for(i in 1:96){ 
  hist(as.numeric(test_het[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(1,0,0,1))}

for(i in 1:96){ 
  hist(as.numeric(test_homo[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='', ylab = 'Number of Variants',
       main ="", border=rgb(0,1,1,1))}

leg.txt <- c("All variant", "Reference Homozygous(A/A)",
             "Heterozygous(a/A)", "Mutant Homozygous (a/a)")
col =c(rgb(0,0,1,1),
       rgb(0,1,0,1),
       rgb(1,0,0,1),
       rgb(0,1,1,1))

legend("topright", col=col, leg.txt,cex=1, pch=20, pt.cex = 1,bty = "n") 

################################################################
##### MRR Ratio ########################################

par(mfrow=c(2,2))

hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 120000),
     main ="Minor Allele Ratio in all variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


hist(as.numeric(test_Ref$S0156[1:(dim(test_Ref$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 110000),
     main ="Minor Allele Ratio in Ref(A/A) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_Ref[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}

hist(as.numeric(test_het$S0156[1:(dim(test_het$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 2000),
     main ="Minor Allele Ratio in Het(a/A) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_het[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


hist(as.numeric(test_homo$S0156[1:(dim(test_homo$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 15000),
     main ="Minor Allele Ratio in Homo(a/a) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_homo[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}

########################################################################
par(mfrow=c(2,2))
hist(as.numeric(test_Ref$S0156[1:(dim(test_Ref$S0156)[1]),5]),
     breaks=100,prob=TRUE,border="white", ylim = c(0, 1.25), 
     main ="Genotyping Quality BWA Exome Sequencing Homo (a/a) variants", xlab ="Genotyping Quality", ylab="Density", cex.lab =
       0.95, cex.main = 0.95) 


for(i in 1:96){ 
  lines(density(as.numeric(test[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_Ref[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_het[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_homo[[i]][,5]),na.rm=TRUE),col="darkorange")}


#############################################################################
## Exome-BWA summary report #################################################
#############################################################################
setwd("~/Documents/GenvsExo/individual_exomes_BWA/")
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
  
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_ALT)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Summary_Exome_BWA =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Summary_Exome_BWA,"/Users/salendrapradh/WGS-and-WES/result/Summary_Exome_BWA.tsv",
            sep="\t", quote =F , row.names = F)

##############################################################################

hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),4]),
     breaks=1000, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     ylim = c(0,8500), xlim = c(0,700),
     main ="Read Depth BWA Exome Sequencing", xlab ="Read Depth", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 


for(i in 1:96){ 
  hist(as.numeric(test[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(0,1,1,1))}

for(i in 1:96){ 
  hist(as.numeric(test_Ref[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(1,1,0,1))}

for(i in 1:96){ 
  hist(as.numeric(test_het[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(1,0,0,1))}

for(i in 1:96){ 
  hist(as.numeric(test_homo[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='', ylab = 'Number of Variants',
       main ="", border=rgb(0,0.5,0.5,0.5))}

leg.txt <- c("All variant", "Reference Homozygous(A/A)",
             "Heterozygous(a/A)", "Mutant Homozygous (a/a)")
col =c(rgb(0,1,1,1),
       rgb(1,1,0,1),
       rgb(1,0,0,1),
       rgb(0,0.5,0.5,0.5))

legend("topright", col=col, leg.txt,cex=1, pch=20, pt.cex = 1,bty = "n") 

#######################################################################
par(mfrow=c(2,2))

hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.75),yaxt="n",ylim =c(0, 120000),
     main ="Minor Allele Ratio in all variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


hist(as.numeric(test_Ref$S0156[1:(dim(test_Ref$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.75),yaxt="n",ylim =c(0, 110000),
     main ="Minor Allele Ratio in Ref(A/A) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_Ref[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}

hist(as.numeric(test_het$S0156[1:(dim(test_het$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.75),yaxt="n",ylim =c(0, 3000),
     main ="Minor Allele Ratio in Het(a/A) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_het[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


hist(as.numeric(test_homo$S0156[1:(dim(test_homo$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.75),yaxt="n",ylim =c(0, 15000),
     main ="Minor Allele Ratio in Homo(a/a) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_homo[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


##############################################################
par(mfrow=c(2,2))
hist(as.numeric(test_Ref$S0156[1:(dim(test_Ref$S0156)[1]),5]),
     breaks=100,prob=TRUE,border="white", ylim = c(0, 1.25), 
     main ="Genotyping Quality BWA Exome Sequencing Homo (a/a) variants", xlab ="Genotyping Quality", ylab="Density", cex.lab =
       0.95, cex.main = 0.95) 


for(i in 1:96){ 
  lines(density(as.numeric(test[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_Ref[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_het[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_homo[[i]][,5]),na.rm=TRUE),col="darkorange")}

#########################################################################################
## Exome-BWT summary report ########################################################
###################################################################################
setwd("~/Documents/BWAvsBWT/individual_exomes_BWT/")
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
  
  MRR_1 <- round((CD_ALT+CD_ALT_2)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR_2 <- round((CD_Ref)/(CD_Ref+CD_ALT+CD_ALT_2), digits =2)
  MRR <- pmin(MRR_1,MRR_2)
  mean_MRR[i] = round(mean(as.numeric(MRR),na.rm=TRUE), digit =2)
  
  CD <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",3))
  mean_CD [i] <- round( mean (as.numeric(CD),na.rm=TRUE), digits =2)
  GQ <- as.numeric(sapply(strsplit(as.character (myfiles[[i]][[6]]), ":"),"[",4))
  mean_GQ [i] <- round(mean(as.numeric(GQ),na.rm=TRUE), digits =2)
  number[i] <- length(CD_ALT)
  test[[i]] = cbind (CD_Ref,CD_ALT,MRR,CD,GQ)
}

Summary_Exome_BWT =cbind(sample,number,mean_CD,mean_GQ,mean_MRR)
write.table(Summary_Exome_BWT,"/Users/salendrapradh/WGS-and-WES/result/Summary_Exome_BWT.tsv",
            sep="\t", quote =F , row.names = F)


hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),4]),
     breaks=1000, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     ylim = c(0,10000), xlim = c(0,600),
     main ="Read Depth BWT Exome Sequencing", xlab ="Read Depth", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 


for(i in 1:96){ 
  hist(as.numeric(test[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(0,1,1,1))}

for(i in 1:96){ 
  hist(as.numeric(test_Ref[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(1,1,0,1))}

for(i in 1:96){ 
  hist(as.numeric(test_het[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", border=rgb(1,0,0,1))}

for(i in 1:96){ 
  hist(as.numeric(test_homo[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
       xlab='', ylab = 'Number of Variants',
       main ="", border=rgb(0,0.5,0.5,0.5))}

leg.txt <- c("All variant", "Reference Homozygous(A/A)",
             "Heterozygous(a/A)", "Mutant Homozygous (a/a)")
col =c(rgb(0,1,1,1),
       rgb(1,1,0,1),
       rgb(1,0,0,1),
       rgb(0,0.5,0.5,0.5))

legend("topright", col=col, leg.txt,cex=1, pch=20, pt.cex = 1,bty = "n") 


################################################################
par(mfrow=c(2,2))

hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 100000),
     main ="Minor Allele Ratio in all variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


hist(as.numeric(test_Ref$S0156[1:(dim(test_Ref$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 110000),
     main ="Minor Allele Ratio in Ref(A/A) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_Ref[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}

hist(as.numeric(test_het$S0156[1:(dim(test_het$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 2000),
     main ="Minor Allele Ratio in Het(a/A) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_het[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}


hist(as.numeric(test_homo$S0156[1:(dim(test_homo$S0156)[1]),3]),
     breaks=100, col=rgb(0,0,0,0),border =rgb(0,0,0,0),
     xlim = c(0,0.6),yaxt="n",ylim =c(0, 15000),
     main ="Minor Allele Ratio in Homo(a/a) variants", xlab ="MRR ", ylab="Number of Variants")
#cex.lab = 0.95,cex.main = 0.95) 

axis(2)
for(i in 1:96){ 
  hist(as.numeric(test_homo[[i]][,3]),na.rm=TRUE, breaks=100,add=TRUE, 
       xlab='Approximate Read Depth', ylab = 'Number of Variants',
       main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,1,1))}

##################################################################
par(mfrow=c(2,2))
hist(as.numeric(test_Ref$S0156[1:(dim(test_Ref$S0156)[1]),5]),
     breaks=100,prob=TRUE,border="white", ylim = c(0, 1.25), 
     main ="Genotyping Quality BWT Exome Sequencing Homo(a/a) variants", xlab ="Genotyping Quality", ylab="Density", cex.lab =
       0.95, cex.main = 0.95) 


for(i in 1:96){ 
  lines(density(as.numeric(test[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_Ref[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_het[[i]][,5]),na.rm=TRUE),col="darkorange")}

for(i in 1:96){ 
  lines(density(as.numeric(test_homo[[i]][,5]),na.rm=TRUE),col="darkorange")}


##################################################################

# boxplot(as.numeric(test$S0156[1:(dim(test$S0156)[1]),1]))
# hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),1]), breaks =50)
# 
# 
# # par (mfrow =c(5,5)) 
# hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),4]),
#   breaks=1000, col=rgb(0,0,0,0),border =rgb(0,0,1,1),
#   ylim = c(0,15000),
#   main ="Read Depth Genome Sequencing", xlab ="Approximate Read Depth", ylab="Number of Variants")
#   #cex.lab = 0.95,cex.main = 0.95) 
# 
# for(i in 1:96){ 
#   hist(as.numeric(test[[i]][,4]),na.rm=TRUE, breaks=1000,add=TRUE, 
#        xlab='Approximate Read Depth', ylab = 'Number of Variants',
#         main ="Read Depth Genome Platform sample[1]", col=rgb(0,0,0,0), border=rgb(0,0,1,1))}
# 
# 
# par (mfrow =c(2,2))
# hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),1]), 
#      # breaks=100,
#      prob=TRUE,
#      # border="white", 
#      # ylim = c(0, 0.15),
#      main ="Read Depth Genome Sequencing", xlab ="Approximate Read Depth", ylab="Density")
#      # cex.lab = 0.95, cex.main = 0.95)
# for(i in 1:96){
#   lines(hist(as.numeric(test[[i]][,4]),na.rm=TRUE))}




hist.data = hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),5]), plot=F)
hist.data$counts = log(hist.data$counts)
par (mfrow =c(2,2))
hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),5]))
plot(hist.data,ylab='log10(Frequency)')





hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),3]),
breaks=50,prob=TRUE,border="white", ylim = c(0,20) ,main ="Read Depth Genome
Sequencing", xlab ="Approximate Read Depth", ylab="Density", cex.lab = 0.95,
cex.main = 0.95) #hist(as.numeric(test$S0156[1:(dim(test$S0156)[1]),2]),
breaks=50,prob=TRUE,border="white") for(i in 1:96){ 
lines(density(as.numeric(test[[i]][,3]),na.rm=TRUE))}


