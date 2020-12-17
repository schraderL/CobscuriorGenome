library(dplyr)
library(data.table)
library(tableHTML)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(reshape2)
library(zoo)
library(scales)
library("ape")
library("ggpubr")
library(ggrepel)
#Counting the total coverage of TE sequences in Genome.
read.delim("/global/homes/jg/r_stef04/Master/ResultFiles/Old/TEbed/Total/ObirTEsorted.csv.bed.bed", header=F, sep = "\t",  stringsAsFactors = FALSE) -> TEsorted
TEsorted$V5<-TEsorted$V3-TEsorted$V2
(sum(TEsorted$V5)/223876465)*100

#Counting the total coverage of classified TE sequences in Genome.
read.delim("/global/homes/jg/r_stef04/Master/ResultFiles/TEbed/Waur.TEclass.bed.MergedClassTE.bed", header=F, sep = "\t",  stringsAsFactors = FALSE) -> TEsorted
TEsorted$V5<-TEsorted$V3-TEsorted$V2
(sum(TEsorted$V5)/223876465)*100

#Counting the individual coverage of superfamilies (nested TEs included, do not add up to 100%). Do not forget to change the genome size at the specified lines in this skript!
read.delim("/global/homes/jg/r_stef04/Master/ResultFiles/TESorted_RMOutfiles/CfloTEsorted.csv", header=F, sep = "\t",  stringsAsFactors = FALSE) -> TEsorted

read.delim("/global/homes/jg/r_stef04/Master/ResultFiles/TERMout_csv/ObirTEsorted.csv", header=F, sep = "\t",  stringsAsFactors = FALSE) -> TEsorted

#read.delim("/home/r/r_stef04/Master/GenomesTE/Sinv2/Sinv2TEsorted.csv", header=F, sep = "\t",  stringsAsFactors = FALSE) -> TEsorted3
TEsorted$V9[TEsorted$V9=="DNA/TcMar-m44"|TEsorted$V9=="DNA/TcMar?"|TEsorted$V9=="DNA/TcMar-ISRm11"|TEsorted$V9=="DNA/TcMar-Pogo"|TEsorted$V9=="DNA/TcMar-Tc4?"|TEsorted$V9=="DNA/Tc1"|TEsorted$V9=="DNA/TcMar-Tc2"|TEsorted$V9=="DNA/TcMar-Fot1"|TEsorted$V9=="DNA/TcMar"|TEsorted$V9=="DNA/TcMar-Mariner"|TEsorted$V9=="DNA/TcMar-Tc1"|TEsorted$V9=="DNA/TcMar-Tc1?"|TEsorted$V9=="DNA/TcMar-Tigger"|TEsorted$V9=="DNA/TcMar-Tigger?"|TEsorted$V9=="DNA/TcMar-Tc4"] <- "DNA/TcMar"
TEsorted$V9[TEsorted$V9=="DNA/CMC-Chapaev"|TEsorted$V9=="DNA/CMC-Chapaev-3"|TEsorted$V9=="DNA/CMC-EnSpm"|TEsorted$V9=="DNA/CMC-EnSpm?"|TEsorted$V9=="DNA/CMC-Transib"] <- "DNA/CMC"
TEsorted$V9[TEsorted$V9=="DNA/hAT"|TEsorted$V9=="DNA/hAT-hobo"|TEsorted$V9=="DNA/hAT-Tag1"|TEsorted$V9=="DNA/hAT-hATx"|TEsorted$V9=="DNA/hAT-Ac"|TEsorted$V9=="DNA/hAT-hAT19"|TEsorted$V9=="DNA/hAT-hAT19?"|TEsorted$V9=="DNA/hAT-hAT5"|TEsorted$V9=="DNA?/hAT?"|TEsorted$V9=="DNA/hAT-Blackjack"|TEsorted$V9=="DNA/hAT-Charlie"|TEsorted$V9=="DNA/hAT?"|TEsorted$V9=="hAT-Charlie"|TEsorted$V9=="hAT-hAT19"|TEsorted$V9=="DNA/hAT-hATm"|TEsorted$V9=="DNA/hAT-Pegasus"|TEsorted$V9=="DNA/hAT-Tip100"] <- "DNA/hAT"
TEsorted$V9[TEsorted$V9=="DNA/Academ"|TEsorted$V9=="DNA/Academ-1"] <- "DNA/Academ"
TEsorted$V9[TEsorted$V9=="DNA/Crypton"|TEsorted$V9=="DNA/Crypton-V"|TEsorted$V9=="DNA/Crypton-I"] <- "DNA/Crypton"
TEsorted$V9[TEsorted$V9=="DNA/Helitron"|TEsorted$V9=="DNA/RC/Helitron"] <- "DNA/Helitron"
TEsorted$V9[TEsorted$V9=="DNA/Ginger"|TEsorted$V9=="DNA/Ginger-1"] <- "DNA/Helitron"
TEsorted$V9[TEsorted$V9=="DNA/MULE-NOF"|TEsorted$V9=="DNA/MuLE-MuDR"|TEsorted$V9=="DNA/MULE-MuDR"|TEsorted$V9=="DNA/MULE-MuDR?"] <- "DNA/MULE"
TEsorted$V9[TEsorted$V9=="DNA/Sola"|TEsorted$V9=="DNA/Sola-1"|TEsorted$V9=="DNA/Sola-2"|TEsorted$V9=="DNA/Sola-3"] <- "DNA/Sola"
TEsorted$V9[TEsorted$V9=="LINE/I-Jockey"|TEsorted$V9=="LINE/Jockey"|TEsorted$V9=="LINE/I"|TEsorted$V9=="LINE/I-Nimb"] <- "LINE/I_Jockey"
TEsorted$V9[TEsorted$V9=="LINE/L1"|TEsorted$V9=="LINE/L1-Tx1"] <- "LINE/L1"
TEsorted$V9[TEsorted$V9=="LINE/R1"|TEsorted$V9=="LINE/R1-LOA"] <- "LINE/R1"
TEsorted$V9[TEsorted$V9=="LINE/R2"|TEsorted$V9=="LINE/R2-NeSL"] <- "LINE/R2"
TEsorted$V9[TEsorted$V9=="LINE/CR1"|TEsorted$V9=="LINE/CR1-Zenon"] <- "LINE/CR1"
TEsorted$V9[TEsorted$V9=="LTR/Gypsy"|TEsorted$V9=="LTR/Gypsy-Cigr"] <- "LTR/Gypsy"
TEsorted$V9[TEsorted$V9=="LINE/RTE"|TEsorted$V9=="LINE/RTE-BovB"|TEsorted$V9=="LINE/RTE-X"|TEsorted$V9=="LINE/RTE-RTE"] <- "LINE/RTE"
TEsorted$V9[TEsorted$V9=="DNA/PIF-Harbinger"|TEsorted$V9=="DNA/PIF-ISL2EU"|TEsorted$V9=="DNA/PIF-ISL2EU?"|TEsorted$V9=="DNA/PIF-Spy"] <- "DNA/PIF"
TEsorted$V9[TEsorted$V9=="DNA/Kolobok"|TEsorted$V9=="DNA/Kolobok-Hydra"|TEsorted$V9=="DNA/Kolobok-T2"|TEsorted$V9=="DNA/Kolobok-T2?"|TEsorted$V9=="DNA/Kolobok-E"] <- "DNA/Kolobok"
TEsorted$V9[TEsorted$V9=="DNA"|TEsorted$V9=="DNA/DNA"] <- "DNA/Unclassified"
TEsorted$V9[TEsorted$V9=="DNA/PiggyBac"|TEsorted$V9=="DNA/PiggyBac?"] <- "DNA/PiggyBac"
TEsorted$V9[TEsorted$V9=="LINE"] <- "LINE/Unclassified"
TEsorted$V9[TEsorted$V9=="LTR"] <- "LTR/Unclassified"
TEsorted$V9[TEsorted$V9=="SINE"] <- "SINE/Unclassified"
TEsorted$V17<-TEsorted$V6-TEsorted$V5
#Genome size has to be changed for current genome!
TEsorted$V18<-100*(TEsorted$V17/223876465)
TEsorted2<-TEsorted[c("V9","V17","V18")]
split(TEsorted2, TEsorted$V9)-> splitTE
#Genome size has to be changed for current genome!
sum(TEsorted$V17)/223876465
sum(TEsorted$V18)
x=0
for (i in splitTE) {
  x=x+1}
df <- data.frame(matrix(NA, nrow = x, ncol = 0)) 
y=0
for (i in splitTE) {
  y=y+1
  print(i$V9[1])
  df$TE_Class[y] <- (i$V9[1])
  df$Perc_of_genome[y] <- sum(i$V18)
  df$Number_of_copies[y] <- nrow(i)
  df$Number_of_masked_bases[y] <- sum(i$V17)}
#df$Perc_of_genom<-100*(df$Number_of_masked_bases/317671980)
#sum(df$Perc_of_genome)
#sum(TEsorted2$V18)
df<-separate(data = df, col = TE_Class, into = c("TE_Class", "TE_Superfamily"), sep = "/")
df2<-as.data.frame(tapply(df$Perc_of_genome, df$TE_Class, FUN=sum))
sum(df2$`tapply(df$Perc_of_genome, df$TE_Class, FUN = sum)`)
write.table(df, file="/global/homes/jg/r_stef04/Master/ResultFiles/TERMout_csv/Obir_perc.csv", sep="\t", col.names=TRUE,row.names=FALSE)
row.names.remove <- c("LINE", "LTR", "SINE", "Retro")
df3 <- subset(df2, rownames(df2) %in% row.names.remove)
df2$sum_RNA<-sum(df3$`tapply(df$Perc_of_genome, df$TE_Class, FUN = sum)`)
write.table(df2, file="/global/homes/jg/r_stef04/Master/ResultFiles/TERMout_csv/Obir_perc2.csv", sep="\t", col.names=FALSE,row.names=TRUE)


df$TE_Class[df$TE_Class=="LINE"] <- "RNA/LINE"
df$TE_Class[df$TE_Class=="SINE"] <- "RNA/SINE"
df$TE_Class[df$TE_Class=="LTR"] <- "RNA/LTR"
df$log_Perc <- log(df$Perc_of_genome, 10)
df$log_Number <- log(df$Number_of_copies, 10)
ggplot(data=df, aes(x=log_Number, y=log_Perc, shape=TE_Class, size=Perc_of_genome))+ 
  labs(shape = "TE classification", size="Amount of genome in %", x = "log10-transformed number of TE copies", y = "log10-transformed proportion of genome (in %)")+
  #geom_label_repel(aes(label = TE_Superfamily),
  #                 size =2,
  #                 box.padding   = 1, 
  #                 point.padding = 0.01,
  #                 segment.color = 'grey50') +
  theme_classic()+
  geom_point(aes(color=TE_Class)) + labs(color = "TE classification")+ggtitle("Temnothorax curvispinosus genome TE content")


read.delim("/home/r/r_stef04/Master/Results/CompleteData.csv", header=T, comment.char="#", sep = "\t", stringsAsFactors=FALSE) -> N
N <- as.data.frame(N)
N <- N[order(N$TE_base_perc),]
N$log_N50 <- log(N$N50_cont)

  ggplot(data = N)+
  geom_point(aes(x=Genome_size,y=TE_base_perc,size=N50_cont,color=Most_common_TE))+
  geom_smooth(aes(x=Genome_size,y=TE_base_perc),method = 'lm')+
  labs(color = "Most common TE family", size="Quality of genome in N50", x="Genome size in bp", y="Amount of TE's in genome (in %)")+ggtitle("Correlation between total amount of TE bases and genome size")

  ggplot(data = N)+
  geom_point(aes(x=Genome_size,y=All_OR,size=N50_cont,color=Most_common_TE))+
  geom_smooth(aes(x=Genome_size,y=All_OR),method = 'lm')+
  labs(color = "Most common TE family", size="Quality of genome in N50", x="Genome size in bp", y="Number of OR gene models found")+ggtitle("Correlation between found OR gene models and genome size")


lm1<-lm(N$TE_base_perc~N$Genome_size+N$N50_cont)
summary(lm1)
lm2<-lm(N$All_OR~N$Genome_size+N$N50_cont)
summary(lm2)


