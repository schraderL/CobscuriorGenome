library("systemPipeR")
library("GenomicFeatures")
library("plyr")
library("tidyverse")
library("dplyr")
txdb <- makeTxDbFromGFF(file="/global/homes/jg/r_stef04/Master/ResultFiles/ORgff/Veme.OR.gff3", format="gff3", organism="Vollenhovia emeryi")
read.delim("/global/homes/jg/r_stef04/Master/ResultFiles/ORgff/ShortenedGff/Veme.OR.new.gff3", header=F, comment.char="#", sep = "\t") -> GeneOR
feat <- genFeatures(txdb, featuretype="all", reduce_ranges=TRUE)
G <- as.data.frame(feat$intergenic, stringsAsFactors=FALSE)
Gene <- as.data.frame(feat$gene, stringsAsFactors=FALSE)
#Take out all genes with intergenic content longer than 20000 bp
tmp <- subset(G, width <= 4000)

#Subsetting according to chromosomes/scaffolds 
x <- split(tmp,tmp$seqnames)


#Seperating gene description
#tmp2 <- x[[1]] %>% remove_rownames %>% column_to_rownames(var="feature_by") %>% as.data.frame()
DFEstablishingFunction<-function(listElement){listElement %>% remove_rownames %>% column_to_rownames(var="feature_by") %>% as.data.frame()}
DFEstablished<-lapply(x,DFEstablishingFunction)
DFEstablished<-DFEstablished[sapply(DFEstablished, nrow)>1]
SeperateGeneNames<- function(listElement){separate(data = listElement, col = featuretype_id, into = c("Gene1", "Gene2"), sep = "__")}
SeperatedGenes<- lapply(DFEstablished, SeperateGeneNames)
PutDown <- function(listElement){listElement %>% mutate(GeneNext=lag(Gene2))}
SeperatedGenes<- lapply(SeperatedGenes, PutDown)

FullOR <- do.call(rbind.data.frame, SeperatedGenes)

FullOR$Array<-NA

#for (i in 1:nrow(FullOR)) {
#  if (is.na(FullOR$GeneNext[i])){
#    FullOR$Array[i]="n"
#  } else{
#    FullOR$Array[i]="y"}
#}

FullOR$GeneNext[is.na(FullOR$GeneNext)] <- "None"


for (i in 1:nrow(FullOR)) {
  if (FullOR$Gene1[i]==FullOR$GeneNext[i]){
    FullOR$Array[i]="y"
  } else{
    FullOR$Array[i]="n"}
}

#tmp3$Array[is.na(tmp3$Array)] <- "n"

FullOR<-FullOR[!is.na(FullOR$Gene2),]

FullOR$Array2<-NA

for (i in 1:nrow(FullOR)) {
  if (FullOR$Array[i]=="n"){
    FullOR$Array2[i]="Start"
  } else {FullOR$Array2[i]="Continuation"}
}

n=0

for (i in 1:nrow(FullOR)) {
  if (FullOR$Array2[i]=="Continuation"){
    FullOR$Array3[i]=n
  } else {
    FullOR$Array3[i]=n+1
    n=n+1}
}

FullORArrays <- FullOR %>%
  group_by(Array3) %>%
  filter(n()>1)

n=0

for (i in 1:nrow(FullORArrays)) {
  if (FullORArrays$Array2[i]=="Continuation"){
    FullORArrays$Array3[i]=n
  } else {
    FullORArrays$Array3[i]=n+1
    n=n+1}
}

ArraysList10 <- split(FullORArrays,FullORArrays$Array3)


library("ape")
library("ggpubr")

#ggscatter(GeneOR, x = "Bp.distance.threshold", y = "log_Genes", ass = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", xlab = "BP", ylab = "Gene number")
#ggscatter(GeneOR, x = "Bp.distance.threshold", y = "log_Arrays", ass = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", xlab = "BP", ylab = "Number of arrays")
#ggscatter(GeneOR, x = "Genes.in.arrays", y = "Number.of.arrays", ass = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", xlab = "Gene number", ylab = "Number of arrays")

#hist(G$width, 5000, xlim=c(0,20000))

#tmp2<-vector()
#for(i in 1:length(ArraysList10)){tmp2[i]<-nrow(ArraysList10[[i]])+1}
#boxplot(tmp2)

#mean(FullORArrays$width)

#for (i in ArraysList10){print(mean(i$width))}

#Step 2: Fusing into the arrays

n=1
d<-setNames(data.frame(matrix(ncol = 5, nrow = length(ArraysList10))), c("scaffold", "start", "end", "gene1", "gene2"))
for (listElement in ArraysList10) {
#  d$scaffold[n]<-as.character(listElement$seqnames[1])
  print(levels(listElement$seqnames[1]))
  d$start[n]<-min(listElement$start)
  d$end[n]<-max(listElement$end)
  d$gene1[n]<-print(listElement$Gene1[1])
  d$gene2[n]<-print(listElement$Gene2[nrow(listElement)])
  n=n+1
}

#Step3: Giving the "real" start and endpoints of array gene couples (earlier start and end were of the intergenic content, NOT the array genes themselves)
GeneOR <- GeneOR[GeneOR$V3=="gene",]
GeneOR %>% separate(V9, c("Name", "Else"), sep = ";") -> GeneOR
GeneOR %>% separate(Name, c("Else", "Name"), sep = "=") -> GeneOR
GeneArray <- merge(GeneOR, FullORArrays, by.x="Name", by.y="Gene1")
GeneArray2 <- merge(GeneOR, GeneArray, by.x="Name", by.y="Gene2")
GeneArray3 <- GeneArray2[c("V1.x", "Name", "Name.y", "V4.y", "V5.x", "Array3")]
ArraysList <- split(GeneArray3,GeneArray3$Array3)
GeneArray3$sum<-GeneArray3$V5.x-GeneArray3$V4.y
GeneArray3 %>%
  group_by(GeneArray3$Array3) %>%
  dplyr::slice(c(which.max(V5.x))) -> GeneArray4
GeneArray3 %>%
  group_by(GeneArray3$Array3) %>%
  dplyr::slice(c(which.min(V4.y))) -> GeneArray5
GeneArray5$End<-GeneArray4$V5.x
GeneArray5$Begin<-GeneArray5$V4.y
GeneArray5$Begin<-GeneArray5$Begin-1
GeneArray5$End<-GeneArray5$End-1
GeneArray5$Array3<-paste("Array",GeneArray5$Array3)
BedArray<-GeneArray5[c("V1.x","Begin","End","Array3")]
BedArray2<-BedArray %>% 
  dplyr::rename(
    chrom = V1.x,
    chromStart = Begin,
    chromEnd = End,
    name = Array3
  )
as.data.frame(lengths(BedArray))->ArrayNumber
as.data.frame(nrow(GeneArray))->GeneNumber
GeneNumber$Arrays<-ArrayNumber$`lengths(BedArray)`[1]
GeneNumber$mean<-GeneNumber$`nrow(GeneArray)`/GeneNumber$Arrays
write.table(GeneNumber, "/global/homes/jg/r_stef04/Master/ResultFiles/Arrays/VemeStatsArrays.csv", row.names = FALSE, sep="\t", quote = FALSE)
write.table(BedArray2, "/global/homes/jg/r_stef04/Master/ResultFiles/Arrays/VemeArrays.bed", row.names = FALSE, sep="\t", quote = FALSE)
#read.delim("/global/homes/jg/r_stef04/Master/ResultFiles/GenomeFiles/Bed/Veme.fa.fai.bed", header=F, comment.char="#", sep = "\t") -> BedGenome
#BedGenome<-BedGenome %>% 
#  dplyr::rename(
#    chrom = V1
#  )
#BedAdd <- merge(BedGenome, BedArray2, by.x="chrom")
#BedAdd$DivEnd <- BedAdd$V3-BedAdd$chromEnd
#BedAdd$DivBegin <- BedAdd$chromStart
#CA <- data.frame(Scaffold=character(),
#                 Array=integer(),
#                 Begin=integer(),
#                 End=integer(),
#                 stringsAsFactors=FALSE)
#for (i in ArraysList) {
#  i[order((i$V4.y), decreasing = TRUE),]
#  #i[order((i$V5.x), decreasing = TRUE),]
#}
#ArraysList3 = do.call(rbind, ArraysList)
#Benchmark program
read.delim("/home/r/r_stef04/Master/Results/ORarrayBenchmark.csv", header=T, comment.char="#", sep = "\t") -> GeneOR2
GeneOR2$log_Genes <- log(GeneOR2$Genes.in.arrays)
#my_data$log_N50 <- log(my_data$N50)
p = ggplot() + 
  geom_line(data = GeneOR2, aes(x = Bp.distance.threshold, y = Genes.in.arrays), color = "blue") +
  geom_line(data = GeneOR2, aes(x = Bp.distance.threshold, y = Number.of.arrays*10), color = "red") + 
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Number of arrays (red)"))
print(p)