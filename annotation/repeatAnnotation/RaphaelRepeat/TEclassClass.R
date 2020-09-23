library(gtools)
library(tidyr)
#1
read.delim("/global/homes/jg/r_stef04/Master/Library/FinalLib/UnknownClass/UC.fa", header=F, comment.char="", sep = "\t") -> Lib1 
Lib1<-Lib1[c("V1")]
write.table(Lib1, "/global/homes/jg/r_stef04/Master/Library/FinalLib/UnknownClass/ListHeader.csv", row.names = FALSE, col.names=FALSE, sep=" ", quote = FALSE)
#2
read.delim("/global/homes/jg/r_stef04/Master/Library/FinalLib/UnknownClass/ListHeader.csv", header=F, comment.char="", sep = " ") -> Lib2
Lib3<-Lib2[c("V1","V7","V2","V3","V4","V5")]
Lib4 <- separate(data = Lib3, col = V1, into = c("Class", "Seq"), sep = "#")
Lib5<-Lib4[c("Class","V7","V2","V3","V4","V5")]
Lib5$Fin<-paste(Lib5$Class,"#",Lib5$V7," ",Lib5$V2,"  ",Lib4$V4," ",Lib5$V5, sep = "")
Lib6<-Lib5[c("Fin")]
write.table(Lib6, "/global/homes/jg/r_stef04/Master/Library/FinalLib/UnknownClass/ListHeader1.csv", row.names = FALSE, col.names=FALSE, sep=" ", quote = FALSE)
#3
read.delim("/global/homes/jg/r_stef04/Master/Library/FinalLib/UnknownClass/AddUnCl2.csv", header=F, comment.char="", sep = " ") -> Liba2
Liba2$V2[Liba2$V2=="result:"|Liba2$V2=="RepeatModeler"] <- ""
Liba3 <- separate(data = Liba2, col = V1, into = c("Class", "Seq"), sep = "#")
Liba3$Fin<-paste(Liba3$Class,"#",Liba3$V3," ",Liba3$V2, sep = "")
Liba4<-Liba3[c("Fin")]
write.table(Liba4, "/global/homes/jg/r_stef04/Master/Library/FinalLib/UnknownClass/ListHeader2.csv", row.names = FALSE, col.names=FALSE, sep=" ", quote = FALSE)
