options(stringsAsFactors=F)

library(goseq) #load goseq package
source("http://bioconductor.org/biocLite.R")
biocLite("geneLenDataBase")
library("geneLenDataBase")

data(hg19.refGene.LENGTH) #load the matrix with length of gene, here you can modify with file containing number of HPAII site or probe 
list=hg19.refGene.LENGTH

a=hg19.refGene.LENGTH #or name of your file with number of observations by gene
#LOI= list of gene of interest, refseqID
a=cbind(a[,2],rep(0,nrow(a)))
b=which(a[,1]%in%LOI) #select gene overlapping between your LOI and refGene
a[b,2]<-1 #assign 1 to overlapping gene
x=as.vector(a[,2])
x=as.numeric(x)
names(x)<-list[,2]
x<-x[-which(duplicated(names(x)))]
pwf=nullp(x,'h19','refGene',bias.data=probe[-which(duplicated(probe[,2])),3])


#biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

#loading dataset for KEGG database

en2eg = as.list(org.Hs.egREFSEQ2EG) #loading gene name matching to refseq ID
eg2kegg = as.list(org.Hs.egPATH)
grepKEGG = function(id, mapkeys) {
unique(unlist(mapkeys[id], use.names = FALSE))
}
kegg = lapply(en2eg, grepKEGG, eg2kegg)
pathKEGG=goseq(pwf,"hg19","refGene",gene2cat=kegg,test.cats=c("KEGG")) #list of pathways based on KEGG database

###matching KEGG pathways names and gene names to the output table

xx <- as.list(org.Hs.egPATH2EG)

tableR=mat.or.vec(20,max(pathKEGG[1:20,4])) #crating a matrix based on the top20 pathways

for (i in c(1:20)){
print(i)
tableR[i,]=c(unique(CI[CI%in%list[which(list[,1]%in%xx[[pathKEGG[i,1]]]),2]]),(rep(0,(max(pathKEGG[1:20,4])-(length(unique(CI[CI%in%list[which(list[,1]%in%xx[[pathKEGG[i,1]]]),2]])))))))}

source("http://bioconductor.org/biocLite.R")
biocLite("KEGG.db")
library(KEGG.db) #loading pathways name matching with KEGG ID

pathR<- as.list(KEGGPATHNAME2ID)


testpathR=apply(as.matrix(pathKEGG[1:20,1]),1,function(each){return(names(pathR[pathR%in%each]))})
a<- data.frame(matrix(unlist(testpathR), nrow=20, byrow=T))


name=refseq[,5:6]
test=mat.or.vec(20,max(pathKEGG[1:20,4]))
for (i in c(1:20)){
print(i)
test[i,]=c(unique(name[(name[,1]%in%tableR[i,]),2]),rep(0,max(pathKEGG[1:20,4])-(length(unique(name[(name[,1]%in%tableR[i,]),2])))))}

final=data.frame(a,pathKEGG[1:20,],tableR,test)




##you can do the same using different database 

#GO database

en2eg = as.list(org.Hs.egREFSEQ2EG)
eg2GO = as.list(org.Hs.egGO)
grepGO = function(id, mapkeys) {
unique(unlist(mapkeys[id], use.names = FALSE))
}
GO = lapply(en2eg, grepGO, eg2GO)

pathGO=goseq(pwf,"hg19","refGene",gene2cat=GO,test.cats=c("GO:BP"))

xx <- as.list(org.Hs.egGO2EG)

pathR<- as.list(GOTERM)

#REACTOME database

library(reactome.db)
en2eg = as.list(org.Hs.egREFSEQ2EG)
eg2REACT = as.list(reactomeEXTID2PATHID)
grepREACT = function(id, mapkeys) {
unique(unlist(mapkeys[id], use.names = FALSE))
}
REACT = lapply(en2eg, grepREACT, eg2REACT)

pathREACT=goseq(pwf,"hg19","refGene",gene2cat=REACT)

xx <- as.list(reactomePATHID2EXTID)

pathR<- as.list(reactomePATHNAME2ID)