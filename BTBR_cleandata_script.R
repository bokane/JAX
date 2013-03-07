###########################################################
#  This is a script to demo the use of Attie BTBR eQTL data
#  January 8, 2013 - GAC
###########################################################
# 
#set working directory
setwd("C:/Users/Bo/Desktop/Jax_Sue")

#load qtl library
library(qtl)
library(ggplot2)
# for windows
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

###########################################################
### load the data 
load("BTBR_Data.RData")
load("annot.RData")
#
ls()

rm(list=c("batch.adipose","batch.gastroc","batch.hypo",     
  "batch.islet","batch.kidney","batch.liver"))
ls()

###########################################################
### synch the annotation with the gex data
names(annot)
dim(annot)
dim(liver.mlratio)
table(annot$probe_use)

annot.use <- subset(annot, probe_use==1)
dim(annot.use)

###transform data
rz.transform<-function(y) {
   rankY=rank(y, ties.method="average", na.last="keep")
   rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
   rzT
}

liver.rz <- apply(liver.mlratio[,annot.use$a_gene_id],2,rz.transform)
dim(liver.rz)

adipose.rz <- apply(adipose.mlratio[,annot.use$a_gene_id],2,rz.transform)
dim(adipose.rz)

gastroc.rz <- apply(gastroc.mlratio[,annot.use$a_gene_id],2,rz.transform)
dim(gastroc.rz)

hypo.rz <- apply(hypo.mlratio[,annot.use$a_gene_id],2,rz.transform)
dim(hypo.rz)

islet.rz <- apply(islet.mlratio[,annot.use$a_gene_id],2,rz.transform)
dim(islet.rz)

kidney.rz <- apply(kidney.mlratio[,annot.use$a_gene_id],2,rz.transform)
dim(kidney.rz)

rm(list=ls(pattern="*.mlratio"))
ls()

###########################################################
### clean up and rename phenotypes

### lipomics 

names(lipomics)
(long.names.lipomics <-names(lipomics)[-c(2:12,19,26,33,40,44)])
lipomics<- lipomics[,-c(2:12,19,26,33,40,44)]
#
names(lipomics) <- c("MouseNum","WEIGHT.4wk","LENGTH.4wk",            
	"GLU.4wk","INS.4wk","TRIG.4wk","HOMA.4wk","WT.6wk","LEN.6wk",            
	"GLU.6wk","INS.6wk","TRIG.6wk","HOMA.6wk","WT.8wk","LEN.8wk",            
	"GLU.8wk","INS.8wk","TRIG.8wk","HOMA.8wk","WT.10wk","LEN.10wk",
	"GLU.10wk","INS.10wk","TRIG.10wk","HOMA.10wk","insulin.10wk",
	"cpeptide.10wk","pepins.10wk","glucose.10wk")
#
cbind(names(lipomics),long.names.lipomics)
#
dim(lipomics)

###rbm
long.names.rbm <- names(rbm)
names(rbm) <- c("MouseNum", "ApoA1","Beta.2.Microglobulin","Calbindin","CD40","CD40.Ligand",                                                              
	"Clusterin","CRP","Cystatin.C","EGF","Endothelin.1","Eotaxin","Factor.VII",                                                               
	"FGF.9","FGF.basic","Fibrinogen","clotted","GCP.2","GM.CSF","Growth.Hormone",                                                           
	"GST.alpha","GST.Mu","Haptoglobin","IFN.gamma","IgA","IL.10","IL.11","IL.12p70",                                             
	"IL.17","IL.18","IL.1alpha","IL.1beta","IL.2","IL.3","IL.4","IL.5","IL.6","IL.7",                                                     
	"Insulin.rbm","IP.10","KC.GROalpha","Leptin","LIF","Lymphotactin","MCP.1","MCP.3",                               
	"MCP.5","M.CSF","MDC","MIP.1alpha","MIP.1beta","MIP.1gamma","MIP.2","MIP.3beta",                        
	"MMP.9","MPO","Myoglobin","NGAL","OSM","Osteopontin","RANTES","SAP","SCF","SGOT",                           
	"TIMP.1","Tissue.Factor","TNF.alpha","TPO","VCAM.1","VEGF","vWF","Fibrinogen") 
#	
cbind(names(rbm),long.names.rbm)

### cpl
long.names.cpl <- names(cpl) 
names(cpl) <- c("MouseNum","NEFA","LDL","HDL","CHOL")
#	
cbind(names(cpl),long.names.cpl)

###d2o
d2o <- d2o[,1:7]
names(d2o)

### liverTG
long.names.liverTG <- names(liverTG)[c(1,6,7,8,9)]
liverTG <- liverTG[,c(1,6,7,8,9)]
names(liverTG) <- c("MouseNum","TG","TG.homogenate","Protein","liver.TG")
#	
cbind(names(liverTG),long.names.liverTG)

### necropsy
long.names.necropsy <- names(necropsy)[-c(2,3,21)]
necropsy <- necropsy[,-c(2,3,21)]
names(necropsy) <- c("MouseNum","SVL.Length","Hypothalamus.wt","Brain.wt",                               
	"Liver.wt","R.Kidney.wt","R.Adipose.wt","Fat.wt","L.Adipose.wt","Liver.wt",                          
	"Rem.Liver.wt","Soleus.wt","Gastroc.wt","Spleen.wt","Heart.wt","n.islets",                        
	"Weight","Gastroc.wt")
#	
cbind(names(necropsy),long.names.necropsy)

### plasmaurine
long.names.plasmaurine <- names(plasmaurine)[-c(2,16)]
plasmaurine <- plasmaurine[,-c(2,16)]
names(plasmaurine) <- c("MouseNum","Plasma.Urea","Plasma.Creatinine","Plasma.Sodium",        	"Plasma.Potassium","Plasma.Chloride","Plasma.Date","Urine.Volume",           	"Urinary.Creatinine","Urinary.Sodium","Urinary.Potassium","Urinary.Calcium",    	"Urinary.Magnesium","Urinary.Protein")
#	
cbind(names(plasmaurine),long.names.plasmaurine)

# and in the darkness cbind them
phenotypes <- cbind(cpl, d2o[,-1],lipomics[,-1],liverTG[,-1],necropsy[,-1],
	plasmaurine[,-1],rbm[,-1])
class(phenotypes)
dim(phenotypes)	

tmp <- c(long.names.cpl, names(d2o)[-1],long.names.lipomics[-1],
	long.names.liverTG[-1],long.names.necropsy[-1],
	long.names.plasmaurine[-1],long.names.rbm[-1])
#
cbind(tmp,names(phenotypes))	
#
long.names<-tmp	
rm(list=ls(pattern="long.names.*"))
rm(cpl,d2o,lipomics,liverTG,necropsy,plasmaurine,rbm)
rm(tmp)
#
annot <- annot.use
rm(annot.use)
ls()
	
#
phenotypes.rz <- phenotypes
n.pheno <- dim(phenotypes)[2]-1
phenotypes.rz[,2:n.pheno+1]<-apply(phenotypes[,2:n.pheno+1],2,rz.transform)
#

#check row correspondence
cbind(phenotypes$MouseNum, rownames(liver.rz), rownames(adipose.rz))

ls()

save(file="BTBR.clean.data.Rdata")



