library(gtools)
library(tidyr)
read.delim("/global/homes/jg/r_stef04/Master/TEpipeline/RS/CobsTE/CobsTE.classif", header=F, sep = "\t") -> TE1
TElist <- list()

#For-Schleife über alle Reihen des data frames
for (i in 1:nrow(TE1)) {
  row <- TE1[i, ]
  if (grepl("coding",row$V8) == T){
    codin1<-gsub(".*coding=\\((.*?)\\)\\;.*","\\1", row$V8,perl=T)
    codin1<-gsub("TE_BLRtx:(.*?)$|TE_BLRtx:(.*?)\\;","\\1\\2",codin1,perl=T)
    codin1<-gsub("\\:"," ",codin1,perl=T)
    codin1<-gsub("\\,","\t",codin1,perl=T)
    codin1<-gsub("%","",codin1,perl=T)
    codin1<-gsub(";","",codin1,perl=T)
    codin1<-gsub("TE_BLRx ","",codin1,perl=T)
    c<-strsplit(x = codin1,"\t")
    c<-as.data.frame(c)
    colnames(c) <- "X"
    separate(c, X, c("a", "b", "c", "d", "e", "f", "g", NA), sep = " ", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn") -> c2
    Class1<- c2[c("b", "c", "d", "e", "g")]
    #Class1[complete.cases(Class1), ]
    #Class1 %>% rename(b = TE_name, c = TE_class, d = TE_family, e = TE_subfamily, g = TE_percentage)
    colnames(Class1)[colnames(Class1) == "b"] <- "TE_name"
    colnames(Class1)[colnames(Class1) == "c"] <- "TE_class"
    colnames(Class1)[colnames(Class1) == "d"] <- "TE_family"
    colnames(Class1)[colnames(Class1) == "e"] <- "TE_subfamily"
    colnames(Class1)[colnames(Class1) == "g"] <- "perc"
    Class1$perc <- as.numeric(Class1$perc)
    Class2<-Class1[order((Class1$perc), decreasing = TRUE),]
    Class2<-head(Class2[1, ])
  }else{
    Class2<-data.frame(OldClass=NA,TE_name=NA,Class=NA,TE_class=NA,Family=NA,TE_family=NA,TE_subfamily=NA)
  }
  
  Class2$OldClass <- row$V1
  Class2$Class <-row$V5
  Class2$Family <-row$V6
  Class3 <- Class2[c("OldClass","TE_name", "Class","TE_class","Family","TE_family","TE_subfamily")]
  TElist[[i]] <-Class3
}
FinalTE = do.call(rbind, TElist)

#New columns containing only actual classes, families and subfamilies (merging other data as well)
FinalTE$Subfamily <- NA
FinalTE$Family2 <- NA

FinalTE$Subfamily[FinalTE$TE_subfamily=="Bel-Pao"] <- "BEL"
FinalTE$Subfamily[FinalTE$TE_subfamily=="CACTA"] <- "CMC-EnSpm"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Chapaev"] <- "CMC-Chapaev"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Copia"] <- "Copia"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Crypton"] <- "Crypton"
FinalTE$Subfamily[FinalTE$TE_subfamily=="DIRS"] <- "DIRS"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Ginger1"] <- "Ginger-1"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Gypsy"] <- "Gypsy"
FinalTE$Subfamily[FinalTE$TE_subfamily=="hAT"] <- "hAT"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Helitron"] <- "Helitron"
FinalTE$Subfamily[FinalTE$TE_subfamily=="ISL2EU"] <- "PIF-ISL2EU"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Jockey"] <- "I-Jockey"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Kolobok"] <- "Kolobok"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Maverick"] <- "Maverick"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Merlin"] <- "Merlin"
FinalTE$Subfamily[FinalTE$TE_subfamily=="MuDR"] <- "MULE-MuDR"
FinalTE$Subfamily[FinalTE$TE_subfamily=="P"] <- "P"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Penelope"] <- "Penelope"
FinalTE$Subfamily[FinalTE$TE_subfamily=="PIF-Harbinger"] <- "PIF-Harbinger"
FinalTE$Subfamily[FinalTE$TE_subfamily=="PiggyBac"] <- "PiggyBac"
FinalTE$Subfamily[FinalTE$TE_subfamily=="R2"] <- "R2"
FinalTE$Subfamily[FinalTE$TE_subfamily=="RTE"] <- "RTE"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Sola"] <- "Sola"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Tc1-Mariner"] <- "TcMar-Tc1"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Transib"] <- "CMC-Transib"
FinalTE$Subfamily[FinalTE$TE_subfamily=="I"] <- "I"
FinalTE$Subfamily[FinalTE$TE_subfamily=="Academ"] <- "Academ"
FinalTE$Family2[FinalTE$Class=="II" & is.na(FinalTE$Family2)] <- "DNA"
FinalTE$Subfamily[FinalTE$Family=="DIRS"] <- "DIRS"
FinalTE$Family2[FinalTE$TE_family=="LTR"] <- "LTR"
FinalTE$Subfamily[FinalTE$TE_family=="LTR" & FinalTE$Subfamily==""] <- "LTR"
FinalTE$Family2[FinalTE$TE_family=="LINE"] <- "LINE"
FinalTE$Family2[FinalTE$Family=="LINE"] <- "LINE"
#FinalTE$Subfamily[FinalTE$TE_family=="LINE" & FinalTE$Subfamily==""] <- "LINE"
FinalTE$Subfamily[FinalTE$Family=="LINE" & FinalTE$Subfamily==""] <- "LINE"
FinalTE$Family2[FinalTE$TE_family=="SINE"] <- "SINE"
FinalTE$Family2[FinalTE$Family=="SINE"] <- "SINE"
#FinalTE$Subfamily[FinalTE$TE_family=="SINE" & FinalTE$Subfamily==""] <- "SINE"
#FinalTE$Subfamily[FinalTE$Family=="SINE" & FinalTE$Subfamily==""] <- "SINE"
FinalTE$Family2[FinalTE$Subfamily=="TcMar-Tc1"] <- "DNA"
#FinalTE$Family2[FinalTE$Subfamily=="Penelope"] <- "LINE"
FinalTE$Family2[FinalTE$Class=="I" & is.na(FinalTE$TE_family) & is.na(FinalTE$Family2)] <- ""
FinalTE$Family2[FinalTE$Class=="noCat"] <- "Unknown"
FinalTE$Family2[FinalTE$Family=="PotentialHostGene"] <- "PotentialHostGene"
FinalTE$Subfamily[is.na(FinalTE$Subfamily)] <- ""
FinalTE$Family2[is.na(FinalTE$Family2)] <- "Unknown"
FinalTE$Family2[FinalTE$Family=="SINE|TRIM"] <- "SINE"
FinalTE$Subfamily[FinalTE$Family=="SINE|TRIM"] <- ""
FinalTE$Family2[FinalTE$Family=="PLE|MITE" & FinalTE$Subfamily=="Penelope"] <- "LINE"
FinalTE$Family2[FinalTE$Family=="LTR|TIR" & FinalTE$Subfamily=="Gypsy"] <- "LTR"
FinalTE$Family2[FinalTE$Family=="PLE|TIR" & FinalTE$Subfamily=="Gypsy"] <- "LTR"
FinalTE$Family2[FinalTE$Family=="DIRS|TIR" & FinalTE$Subfamily=="BEL"] <- "LTR"
FinalTE$Family2[FinalTE$Subfamily=="Sola"] <- "DNA"
#FinalTE$Family2[FinalTE$Class=="I" & FinalTE$Family2 == "Unknown"] <- "RNA"
FinalTE$Subfamily[FinalTE$Family2=="Retroposon?"] <- "Retroposon?"
FinalTE$Family2[FinalTE$Family2 == "SINE" & FinalTE$Subfamily== "Retroposon?"] <- "Unknown"
FinalTE$Subfamily[FinalTE$Family=="LARD" & FinalTE$Subfamily==""] <- "LARD"
FinalTE$Subfamily[FinalTE$Family=="Maverick" & FinalTE$Subfamily==""] <- "Maverick"
FinalTE$Subfamily[FinalTE$Family=="Helitron" & FinalTE$Subfamily==""] <- "Helitron"
FinalTE$Subfamily[FinalTE$Family=="DIRS"] <- "DIRS"
FinalTE$Family2[FinalTE$Subfamily=="DIRS"] <- "LTR"
FinalTE$Subfamily[FinalTE$TE_subfamily=="5S"] <- "5S"
FinalTE$Family2[FinalTE$Family=="LINE" & FinalTE$Subfamily==""] <- "LINE"
FinalTE$Subfamily[FinalTE$Family=="TRIM" & FinalTE$Subfamily==""] <- "TRIM"
FinalTE$Subfamily[FinalTE$Family=="LTR" & FinalTE$Subfamily==""] <- "TRIM"
FinalTE$Family2[FinalTE$Subfamily=="TRIM"] <- "LTR"
FinalTE$Family2[FinalTE$Subfamily=="LARD"] <- "LTR"
FinalTE$Subfamily[FinalTE$Family2=="DNA" & FinalTE$Subfamily=="" & FinalTE$Family=="TIR"] <- "TIR"
FinalTE$Subfamily[FinalTE$Family2=="DNA" & FinalTE$Subfamily=="" & FinalTE$Family=="MITE"] <- "MITE"
FinalTE$Family2[FinalTE$Family2=="Unknown" & FinalTE$Subfamily=="Penelope"] <- "LINE"
FinalTE$Subfamily[FinalTE$Subfamily=="LTR" & FinalTE$Family2=="LTR"] <- ""

# To transform the fasta file into a data frame, use awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < (input file) > out.fa
#in the command line first

read.delim("/home/r/r_stef04/Master/TEpipeline/RS/CobsTE/out.fa", header=F, comment.char="#", sep = "\t") -> Lib1 

Lib2 <- separate(data = Lib1, col = V1, into = c("V0", "V1"), sep = ">")
Lib3 <- Lib2[c("V1","V2")]

#Creating new headers
TEline <- list()
TEline2 <- list()
for (i in 1:nrow(FinalTE)) {
  if (FinalTE$Subfamily[i]==""){
    newline <-paste(FinalTE[i,1],"#",FinalTE[i,9]," @Hymenoptera","  [S:] ",FinalTE[i,1],sep="")
  }else{
    newline <-paste(FinalTE[i,1],"#",FinalTE[i,9],"/",FinalTE[i,8]," @Hymenoptera","  [S:] ",FinalTE[i,1],sep="")
  }
  oldline <- paste(FinalTE[i,1])
  TEline[[i]] <- newline
  TEline2[[i]] <- oldline
}

#Combining new headers with sequences
FinalTE2 = do.call(rbind, TEline)
FinalTE4 = do.call(rbind, TEline2)

Lib2<-Lib3
Lib2$id<-gsub(".*_(R.[0-9]+$)","\\1",Lib2$V1,perl=T)
FTE3<-as.data.frame(FinalTE2)
FTE3$id<-gsub(".* (R\\.[0-9]+$)","\\1",FTE3$V1,perl=T)
Final<-merge(Lib2,FTE3,by="id",all.x=T,all.y=T)
Final <- Final[c("V1.y","V2")]

write.table(Final, "/home/r/r_stef04/Master/TEpipeline/RS/CobsTE/Cobs.txt", row.names = FALSE, sep="\t", quote = FALSE)

