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
library(lattice)
# For Windows
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}
###########################################################
### load the data 
load("BTBR.clean.data.RData")
#load("annot.RData")
#
ls()

rm(list=c("batch.adipose","batch.gastroc","batch.hypo",     
  "batch.islet","batch.kidney","batch.liver","cpl","d2o",            
  "gastroc.mlratio","hypo.mlratio","islet.mlratio","kidney.mlratio", 
  "liverTG","rbm","necropsy","plasmaurine"))
ls()
# #########################################################
# # levelplot  # will use with fewer genotypes
# jpeg(file='levelplot.jpg',width = 15000, height = 15000,
#      pointsize = 12, quality = 75, bg = "white", res = 800)
# levelplot <- cor(phenotypes.rz[,-c(1,2,18,19,76)])
# require(lattice)
# levelplot(levelplot)
# dev.off()

#########################################################
# heat map
phenotypes1.rz <-phenotypes.rz[,-c(1,218,19,76)]
names(phenotypes1.rz)

# #quartz()
# jpeg(file='heatmap.jpg',width = 15000, height = 15000,
#      pointsize = 12, quality = 75, bg = "white", res = 800)
# heatmap(cor(phenotypes1.rz[-c(1,2)],use="complete.obs"))
# dev.off()

# heat map on heat spot_1
# phenotypes2.rz <-phenotypes.rz[,c(77,150,2,144,108,80,17,95,89,64,83,46,55,27,34,32,25,81,12,39,41,48,113,90,
#                                   6,53,91,152,58,21,79,105,104,159,50,29,22,28,26,96,36,13,42,43,35,24,23,121,122,123,117,126,
#                                   115,20,70,116,112,107,78,99,94,109,102,16,61,100,11,7,9,111)]
# names(phenotypes2.rz)
# #quartz()
# jpeg(file='heatmap_1.jpg',width = 8000, height = 8000,
#      pointsize = 12, quality = 75, bg = "white", res = 800)
# heatmap(cor(phenotypes2.rz,use="complete.obs"))
# dev.off()


###########################################################
###find my genes 
(my.gene.names <- grep('^Npy',annot$gene_symbol,value=TRUE))
#
#find the index
(my.gene.indx <- grep('^Npy',annot$gene_symbol)[4])
#
# grab the data columns for these genes
Npy.expr <- as.data.frame(cbind(
  hypo.rz[,annot[my.gene.indx,"a_gene_id"]],
  adipose.rz[,annot[my.gene.indx,"a_gene_id"]]
  ))
# 
#look
head(Npy.expr)
#
#rename vars
names(Npy.expr) <- c("Npy.hypo","Npy.adip")
#
# look at scatterplots (not that important)
quartz()  ## tissue - tissue gene expression correlation 
qplot(Npy.hypo, Npy.adip, data=Npy.expr) + geom_smooth(method="lm")
#
cor(Npy.expr)
#
#clean up
rm(my.gene.indx, my.gene.names) 
ls()

###########################################################
###look for correlates of Npy (gene)
#
#create a correlation "scan" function
# x is a matrix or dataframe with lots varaibles 
# y is a vector that we will assign outisde of the function
# "apply" this function to compute correlation of y with each element of x
mycorr <- function(x){
  cor(x,y,use="complete.obs")
  }

# compute correlations of Npy.hypo(y) with gene expression in hypothalamus tissue(x)
y <- Npy.expr$Npy.hypo
hypo.cor <- apply(hypo.rz,2,"mycorr")
#
quartz()

hist(hypo.cor, breaks=100) # Correlation coefficient not normally distributed. 
                           # Note: correlation cofficient is not transformed
#
quartz()
qqnorm(hypo.cor)
#
# what are the most highly correlated genes?
sum(abs(hypo.cor)>0.5)
sum(abs(hypo.cor)>0.6)
sum(abs(hypo.cor)>0.7)
sum(abs(hypo.cor)>0.85)
indx.hypo <- which(abs(hypo.cor)>0.85) 
annot[indx.hypo, c("a_gene_id","gene_symbol")] # returns a list of genes names (with gene IDs) of genes expressed in hypo that are
                                               # highly correlated with Npy gene expression in 
                                               # in hypothalamus tissue


#look at phenotype - gene correlations
y <- Npy.expr$Npy.hypo
hypo.phenotypes.cor <- apply(phenotypes.rz[,-1],2,"mycorr")
sort(hypo.phenotypes.cor)

# compute correlations of Npy expression in adipose tissue with other genes expressed in the tissue
y <- Npy.expr$Npy.adip
adip.cor <- apply(adipose.rz,2,"mycorr")
#
quartz()
hist(adip.cor, breaks=100)
#
quartz()
qqnorm(adip.cor)
#
sum(abs(adip.cor)>0.5)
sum(abs(adip.cor)>0.6)
sum(abs(adip.cor)>0.7)
indx.adip <- which(abs(adip.cor)>0.7)
annot[indx.adip, c("a_gene_id","gene_symbol")] # 10024398389.1

#look at phenotype correlations
y <- Npy.expr$Npy.adip
adip.phenotypes.cor <- apply(phenotypes.rz[,-1],2,"mycorr")
sort(adip.phenotypes.cor)

###########################################################
###create cross structure for gex mapping analysis 

#adds expression data to temporary variable
tmp <- as.data.frame(cbind(hypo.rz[,indx.hypo],adipose.rz[,indx.adip]))
names(tmp)  <- c(paste("hypo.",annot[indx.hypo, "gene_symbol"],sep=''),
  paste("adip.",annot[indx.adip, "gene_symbol"],sep=''))
#adds all phenotype data to temporary variable
#(del)select columns as you wish to limit amount of 
#imported data

tmp2 <-as.data.frame(phenotypes.rz[,-1])

f2g$pheno <- cbind(f2g$pheno[,c(2,6,7)])

f2g$pheno <- cbind(f2g$pheno, tmp, tmp2)

#double check columns
names(f2g$pheno)

#preview data before beginning more stuff,
#make sure everything is entered properly
summary(f2g)

#clear workspace; remove temp variables from above
rm(tmp2)
rm(tmp)
###########################################################
### pairwise scatterplots

# some useful plotting functions
# see documentation for "pairs" function 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor*0.5, col=c("gray60", "black")[(abs(r)>0.65)+1])
}
#
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2],0,1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
# look at phenotype correlations

# pheno.cor <- cor(phenotypes.rz)
# sort(pheno.cor)

# # 0.65 as threshold
# 
# pheno.number <- which(abs(pheno.cor) < 0.999 & abs(pheno.cor) > 0.65 , arr.ind=TRUE)
# 
# jpeg(file='all.pheno.cor.jpg',width = 15000, height = 15000,
#      pointsize = 12, quality = 75, bg = "white", res = 800)
# pairs(phenotypes1.rz[,c(2,3,4,21,22,23,24,26,43,45,47,54,55,57,58,64,65,101,117,124,125,126,127,128,130,131)], # Column number
#       upper.panel=panel.cor,diag.panel=panel.hist)
# dev.off()

# 0.90 as threshold

# pheno.number <- which(abs(pheno.cor) < 0.999 & abs(pheno.cor) > 0.91 , arr.ind=TRUE)
# 
# jpeg(file='0.90.pheno.cor.jpg',width = 5000, height = 5000,
#      pointsize = 12, quality = 75, bg = "white", res = 500)
# pairs(phenotypes1.rz[,c(3,4,43,45,47,54,55,57,58,64,65)], # Column number
#       upper.panel=panel.cor,diag.panel=panel.hist)
# dev.off()

###################################################
#
quartz()
pairs(f2g$pheno[,4:9], # Column number
  upper.panel=panel.cor,diag.panel=panel.hist)
quartz()
pairs(f2g$pheno[,10:14], # Column number
      upper.panel=panel.cor,diag.panel=panel.hist)

# #
# quartz()
# pairs(f2g$pheno[,10:15],
#    upper.panel=panel.cor,diag.panel=panel.hist)
# 
# 
# #
# quartz()
# pairs(f2g$pheno[,16:21],
#       upper.panel=panel.cor,diag.panel=panel.hist)
###########################################################
### scan for QTL
#
# setup for scanning
f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

#convenient to keep "sex" handy as a numeric variable
sex <- as.numeric(f2g$pheno$Sex)

###########################
# 2 dimensional QTL scans #
###########################

#keep commented out until run time can be improved

#my.scan2a <- scantwo(f2g, pheno.col=c(4:8), addcovar=sex, method="hk")

#plot(my.scan2a, verbose=FALSE)
#######################
# Regular  QTL  scans #
#######################

#scan with sex as an additive covariate
my.scan1a <- scanone(f2g, pheno.col=c(4:148), addcovar=sex, method="hk")
#
#run permutations
my.perm1a <-scanone(f2g,pheno.col=8,addcovar=sex,method="hk",n.perm=100,perm.Xsp=TRUE)

# > names(f2g$pheno)
# [1] "MouseNum"      "Sex"           "pgm"           "hypo.Tmem215"  "hypo.Dlx1"     "hypo.Lhx8"    
# [7] "hypo.Npy"      "hypo.Nts"      "hypo.Nmu"      "adip.Gpr88"    "adip.Npy"      "adip.Tubal3"  
# [13] "adip.AY512931" "adip.Klk11"

# plot scans
#quartz()
for(i in 4){  # But what is i? Seems to be pheno.col - 1 (hypo.Npy)
quartz()
plot(my.scan1a,lodcolumn=i)
add.threshold(my.scan1a, lodcolumn=1,
  perms=my.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(my.scan1a, lodcolumn=1,
  perms=my.perm1a, alpha=0.20,lty="dashed",lwd=2,col="green")
add.threshold(my.scan1a, lodcolumn=1,
              perms=my.perm1a, alpha=0.63,lty="dashed",lwd=2,col="blue")
}
for(i in 8){  # But what is i? Seems to be pheno.col - 1 (adip.Npy)
  quartz()
  plot(my.scan1a,lodcolumn=i)
  add.threshold(my.scan1a, lodcolumn=1,
                perms=my.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
  add.threshold(my.scan1a, lodcolumn=1,
                perms=my.perm1a, alpha=0.20,lty="dashed",lwd=2,col="green")
  add.threshold(my.scan1a, lodcolumn=1,
                perms=my.perm1a, alpha=0.63,lty="dashed",lwd=2,col="blue")
}
#report results
summary(my.perm1a)
# 
summary(my.scan1a, perms=my.perm1a, alpha=0.10, format="allpheno")# Shows all phenotypes; but eQTL only ran against the last phenotype

#
#effect plots
quartz()
effectplot(f2g, "hypo.Def8", mname2=find.marker(f2g,10,31.7), 
    main="Hypo Def8 @ Chr 10 31.7 cM",mname1="Sex")

#confidence interval & physical map pos 
(intQ2 <- bayesint(my.scan1a[c("chr","pos","hypo.Ogn")], chr=10, expandtomarkers=TRUE)) 
c(pmap[[10]][rownames(intQ2)[1]], pmap[[10]][rownames(intQ2)[3]])
    

#phenotypes.rz <-phenotypes.rz[,-c(40:55)]  
#phenotypes.rz1 <- rowMeans(is.na(phenotypes.rz)) 
#phenotypes.rz2 <- colMeans(is.na(phenotypes.rz))

#quartz()
#heatmap(cor(phenotypes.rz[,2:4],use="complete.obs"))

#names(phenotypes.rz)
