rm(list=ls())

library(chron)
library(colorspace)
library(mime)
library(munsell)
library(labeling)
library(rlang)
library(stringi)
library(evaluate)
library(highr)
library(markdown)
library(yaml)
library(backports)
library(jsonlite)
library(digest)
library(plyr)
library(reshape2)
library(scales)
library(tibble)
library(lazyeval)
library(RColorBrewer)
library(stringr)
library(knitr)
library(magrittr)
library(checkmate)
library(htmlwidgets)
library(viridisLite)
library(Rcpp)
library(Formula)
library(ggplot2)
library(latticeExtra)
library(acepack)
library(gtable)
library(data.table)
library(htmlTable)
library(viridis)
library(htmltools)
library(base64enc)
library(minqa)
library(RcppEigen)
library(lme4)
library(SparseM)
library(MatrixModels)
library(pbkrtest)
library(quantreg)
library(car)
library(htmlTable)
library(Hmisc)
library(survival)
library(foreign)
library(bitops)
library(caTools)
library(gplots)
library(ROCR)
library(mice)
library(writexl)
library(HardyWeinberg)
library(officer)
library(uuid)
library(compareGroups)
library(nlme)
library(vcd)
library(boot)
library(tibble)
library(haven)
library(icenReg)
library(MASS)
library(sandwich)   
library(lmtest)
library(gam)
library(smoothHR)
library(metafor)
library(DBI)
library(mitools)
library(RcppArmadillo)
library(miceadds)
library(dplyr)
library(estimatr)
library(lubridate)
library(snakecase)
library(janitor)
library(miceadds)
library(fmsb)

RutinesLocals<- "N:/data/durable/Syntax/ah/routines"
source(file.path(RutinesLocals,"carrega.llibreria.r"))
source(file.path(RutinesLocals,"merge2.r"))
source(file.path(RutinesLocals,"fix2.r"))
source(file.path(RutinesLocals,"table2.r"))
source(file.path(RutinesLocals,"subset2.r"))
source(file.path(RutinesLocals,"format2.r"))
source(file.path(RutinesLocals,"order2.r"))
source(file.path(RutinesLocals,"intervals.r"))


### GUAPAS ###
##############

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.0001,signif(x,1),
                   ifelse(abs(x)<0.001,signif(x,1),
                          ifelse(abs(x)<0.1,round(x,3),
                                 ifelse(abs(x)<1,round(x,2),signif(x,3)))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",round(x,3)))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

pval_ic_guapa <- function(ratio, lo, up) 
{
  z<-qnorm(1-0.05/2)
  x<-log(ratio)/((log(up)-log(lo))/(2*z))
  pval<-pval_guapa(exp(-0.717*x-(0.416*x^2)))
  return(pval)
}


header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

setwd("N:/data/durable/Projects/Magnus_MR_BMI/gwas")


###########################################
### GENERATION OF DATASETS FOR THE GWAS ###
###########################################

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")
dat<-subset2(dat,"dat$anthrop_na_mom==0 & !is.na(dat$bmi_grs_mom)")
dat<-dat[,c("sentrixid_mom","subfertility_12plus")]
dat<-rename.vars(dat,from=c("sentrixid_mom","subfertility_12plus"),to=c("IID","infertility"))
dat<-dat[order(dat$IID,-abs(dat$infertility)),]
dat<-dat[!duplicated(dat$IID),]
dat$X.FID<-dat$IID
dat<-dat[,c("X.FID","IID","infertility")]
write.table(dat,"./Maternal_Infertility.txt",sep="\t",row.names=FALSE, quote=FALSE)

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")
dat<-subset2(dat,"dat$anthrop_na_dad==0 & !is.na(dat$bmi_grs_dad)")
dat<-dat[,c("sentrixid_dad","subfertility_12plus")]
dat<-rename.vars(dat,from=c("sentrixid_dad","subfertility_12plus"),to=c("IID","infertility"))
dat<-dat[order(dat$IID,-abs(dat$infertility)),]
dat<-dat[!duplicated(dat$IID),]
dat$X.FID<-dat$IID
dat<-dat[,c("X.FID","IID","infertility")]
write.table(dat,"./Paternal_Infertility.txt",sep="\t",row.names=FALSE, quote=FALSE)


#####################
### GWAS IN PLINK ###
#####################


############################################
### MERGING GWAS RESULTS IN STUDY BADGES ###
############################################

### HARVEST - MOTHERS ###

h_m_chr1<-as.data.frame(read.delim("./MoBa_results/h_m_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr2<-as.data.frame(read.delim("./MoBa_results/h_m_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr3<-as.data.frame(read.delim("./MoBa_results/h_m_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr4<-as.data.frame(read.delim("./MoBa_results/h_m_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr5<-as.data.frame(read.delim("./MoBa_results/h_m_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr6<-as.data.frame(read.delim("./MoBa_results/h_m_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr7<-as.data.frame(read.delim("./MoBa_results/h_m_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr8<-as.data.frame(read.delim("./MoBa_results/h_m_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr9<-as.data.frame(read.delim("./MoBa_results/h_m_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr10<-as.data.frame(read.delim("./MoBa_results/h_m_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr11<-as.data.frame(read.delim("./MoBa_results/h_m_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr12<-as.data.frame(read.delim("./MoBa_results/h_m_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr13<-as.data.frame(read.delim("./MoBa_results/h_m_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr14<-as.data.frame(read.delim("./MoBa_results/h_m_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr15<-as.data.frame(read.delim("./MoBa_results/h_m_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr16<-as.data.frame(read.delim("./MoBa_results/h_m_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr17<-as.data.frame(read.delim("./MoBa_results/h_m_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr18<-as.data.frame(read.delim("./MoBa_results/h_m_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr19<-as.data.frame(read.delim("./MoBa_results/h_m_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr20<-as.data.frame(read.delim("./MoBa_results/h_m_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr21<-as.data.frame(read.delim("./MoBa_results/h_m_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
h_m_chr22<-as.data.frame(read.delim("./MoBa_results/h_m_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

harvest_mom<-rbind(h_m_chr1,h_m_chr2,h_m_chr3,h_m_chr4,h_m_chr5,h_m_chr6,h_m_chr7,h_m_chr8,h_m_chr9,h_m_chr10,h_m_chr11,
                   h_m_chr12,h_m_chr13,h_m_chr1,h_m_chr15,h_m_chr16,h_m_chr17,h_m_chr18,h_m_chr19,h_m_chr20,h_m_chr21,h_m_chr22)
harvest_mom$LOG_OR<-log(harvest_mom$OR)
write.table(harvest_mom,"./MoBa_results/harvest_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(harvest_mom,file="./MoBa_results/harvest_mom.RData")


### ROTTERDAM1 - MOTHERS ###

r1_m_chr1<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr2<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr3<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr4<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr5<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr6<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr7<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr8<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr9<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr10<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr11<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr12<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr13<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr14<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr15<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr16<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr17<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr18<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr19<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr20<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr21<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_m_chr22<-as.data.frame(read.delim("./MoBa_results/r1_m_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

rotterdam1_mom<-rbind(r1_m_chr1,r1_m_chr2,r1_m_chr3,r1_m_chr4,r1_m_chr5,r1_m_chr6,r1_m_chr7,r1_m_chr8,r1_m_chr9,r1_m_chr10,r1_m_chr11,
                      r1_m_chr12,r1_m_chr13,r1_m_chr1,r1_m_chr15,r1_m_chr16,r1_m_chr17,r1_m_chr18,r1_m_chr19,r1_m_chr20,r1_m_chr21,r1_m_chr22)
rotterdam1_mom$LOG_OR<-log(rotterdam1_mom$OR)
write.table(rotterdam1_mom,"./MoBa_results/rotterdam1_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(rotterdam1_mom,file="./MoBa_results/rotterdam1_mom.RData")


### ROTTERDAM2 - MOTHERS ###

r2_m_chr1<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr2<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr3<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr4<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr5<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr6<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr7<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr8<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr9<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr10<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr11<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr12<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr13<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr14<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr15<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr16<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr17<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr18<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr19<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr20<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr21<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_m_chr22<-as.data.frame(read.delim("./MoBa_results/r2_m_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

rotterdam2_mom<-rbind(r2_m_chr1,r2_m_chr2,r2_m_chr3,r2_m_chr4,r2_m_chr5,r2_m_chr6,r2_m_chr7,r2_m_chr8,r2_m_chr9,r2_m_chr10,r2_m_chr11,
                      r2_m_chr12,r2_m_chr13,r2_m_chr1,r2_m_chr15,r2_m_chr16,r2_m_chr17,r2_m_chr18,r2_m_chr19,r2_m_chr20,r2_m_chr21,r2_m_chr22)
rotterdam2_mom$LOG_OR<-log(rotterdam2_mom$OR)
write.table(rotterdam2_mom,"./MoBa_results/rotterdam2_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(rotterdam2_mom,file="./MoBa_results/rotterdam2_mom.RData")


### NORMENT_FEB18 - MOTHERS ###

nf18_m_chr1<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr2<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr3<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr4<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr5<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr6<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr7<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr8<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr9<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr10<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr11<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr12<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr13<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr14<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr15<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr16<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr17<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr18<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr19<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr20<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr21<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_m_chr22<-as.data.frame(read.delim("./MoBa_results/nf18_m_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

norment_feb18_mom<-rbind(nf18_m_chr1,nf18_m_chr2,nf18_m_chr3,nf18_m_chr4,nf18_m_chr5,nf18_m_chr6,nf18_m_chr7,nf18_m_chr8,nf18_m_chr9,nf18_m_chr10,nf18_m_chr11,
                         nf18_m_chr12,nf18_m_chr13,nf18_m_chr1,nf18_m_chr15,nf18_m_chr16,nf18_m_chr17,nf18_m_chr18,nf18_m_chr19,nf18_m_chr20,nf18_m_chr21,nf18_m_chr22)
norment_feb18_mom$LOG_OR<-log(norment_feb18_mom$OR)
write.table(norment_feb18_mom,"./MoBa_results/norment_feb18_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(norment_feb18_mom,file="./MoBa_results/norment_feb18_mom.RData")


### NORMENT_MAY16 - MOTHERS ###

nm16_m_chr1<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr2<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr3<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr4<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr5<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr6<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr7<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr8<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr9<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr10<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr11<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr12<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr13<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr14<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr15<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr16<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr17<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr18<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr19<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr20<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr21<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_m_chr22<-as.data.frame(read.delim("./MoBa_results/nm16_m_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

norment_may16_mom<-rbind(nm16_m_chr1,nm16_m_chr2,nm16_m_chr3,nm16_m_chr4,nm16_m_chr5,nm16_m_chr6,nm16_m_chr7,nm16_m_chr8,nm16_m_chr9,nm16_m_chr10,nm16_m_chr11,
                         nm16_m_chr12,nm16_m_chr13,nm16_m_chr1,nm16_m_chr15,nm16_m_chr16,nm16_m_chr17,nm16_m_chr18,nm16_m_chr19,nm16_m_chr20,nm16_m_chr21,nm16_m_chr22)
norment_may16_mom$LOG_OR<-log(norment_may16_mom$OR)
write.table(norment_may16_mom,"./MoBa_results/norment_may16_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(norment_may16_mom,file="./MoBa_results/norment_may16_mom.RData")


### HARVEST - FATHERS ###

h_p_chr1<-as.data.frame(read.delim("./MoBa_results/h_p_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr2<-as.data.frame(read.delim("./MoBa_results/h_p_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr3<-as.data.frame(read.delim("./MoBa_results/h_p_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr4<-as.data.frame(read.delim("./MoBa_results/h_p_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr5<-as.data.frame(read.delim("./MoBa_results/h_p_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr6<-as.data.frame(read.delim("./MoBa_results/h_p_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr7<-as.data.frame(read.delim("./MoBa_results/h_p_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr8<-as.data.frame(read.delim("./MoBa_results/h_p_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr9<-as.data.frame(read.delim("./MoBa_results/h_p_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr10<-as.data.frame(read.delim("./MoBa_results/h_p_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr11<-as.data.frame(read.delim("./MoBa_results/h_p_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr12<-as.data.frame(read.delim("./MoBa_results/h_p_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr13<-as.data.frame(read.delim("./MoBa_results/h_p_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr14<-as.data.frame(read.delim("./MoBa_results/h_p_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr15<-as.data.frame(read.delim("./MoBa_results/h_p_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr16<-as.data.frame(read.delim("./MoBa_results/h_p_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr17<-as.data.frame(read.delim("./MoBa_results/h_p_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr18<-as.data.frame(read.delim("./MoBa_results/h_p_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr19<-as.data.frame(read.delim("./MoBa_results/h_p_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr20<-as.data.frame(read.delim("./MoBa_results/h_p_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr21<-as.data.frame(read.delim("./MoBa_results/h_p_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
h_p_chr22<-as.data.frame(read.delim("./MoBa_results/h_p_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

harvest_dad<-rbind(h_p_chr1,h_p_chr2,h_p_chr3,h_p_chr4,h_p_chr5,h_p_chr6,h_p_chr7,h_p_chr8,h_p_chr9,h_p_chr10,h_p_chr11,
                   h_p_chr12,h_p_chr13,h_p_chr1,h_p_chr15,h_p_chr16,h_p_chr17,h_p_chr18,h_p_chr19,h_p_chr20,h_p_chr21,h_p_chr22)
harvest_dad$LOG_OR<-log(harvest_dad$OR)
write.table(harvest_dad,"./MoBa_results/harvest_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(harvest_dad,file="./MoBa_results/harvest_dad.RData")


### ROTTERDAM1 - FATHERS ###

r1_p_chr1<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr2<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr3<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr4<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr5<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr6<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr7<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr8<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr9<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr10<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr11<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr12<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr13<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr14<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr15<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr16<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr17<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr18<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr19<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr20<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr21<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
r1_p_chr22<-as.data.frame(read.delim("./MoBa_results/r1_p_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

rotterdam1_dad<-rbind(r1_p_chr1,r1_p_chr2,r1_p_chr3,r1_p_chr4,r1_p_chr5,r1_p_chr6,r1_p_chr7,r1_p_chr8,r1_p_chr9,r1_p_chr10,r1_p_chr11,
                      r1_p_chr12,r1_p_chr13,r1_p_chr1,r1_p_chr15,r1_p_chr16,r1_p_chr17,r1_p_chr18,r1_p_chr19,r1_p_chr20,r1_p_chr21,r1_p_chr22)
rotterdam1_dad$LOG_OR<-log(rotterdam1_dad$OR)
write.table(rotterdam1_dad,"./MoBa_results/rotterdam1_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(rotterdam1_dad,file="./MoBa_results/rotterdam1_dad.RData")


### ROTTERDAM2 - FATHERS ###

r2_p_chr1<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr2<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr3<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr4<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr5<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr6<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr7<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr8<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr9<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr10<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr11<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr12<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr13<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr14<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr15<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr16<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr17<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr18<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr19<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr20<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr21<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
r2_p_chr22<-as.data.frame(read.delim("./MoBa_results/r2_p_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

rotterdam2_dad<-rbind(r2_p_chr1,r2_p_chr2,r2_p_chr3,r2_p_chr4,r2_p_chr5,r2_p_chr6,r2_p_chr7,r2_p_chr8,r2_p_chr9,r2_p_chr10,r2_p_chr11,
                      r2_p_chr12,r2_p_chr13,r2_p_chr1,r2_p_chr15,r2_p_chr16,r2_p_chr17,r2_p_chr18,r2_p_chr19,r2_p_chr20,r2_p_chr21,r2_p_chr22)
rotterdam2_dad$LOG_OR<-log(rotterdam2_dad$OR)
write.table(rotterdam2_dad,"./MoBa_results/rotterdam2_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(rotterdam2_dad,file="./MoBa_results/rotterdam2_dad.RData")


### NORMENT_FEB18 - FATHERS ###

nf18_p_chr1<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr2<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr3<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr4<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr5<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr6<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr7<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr8<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr9<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr10<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr11<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr12<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr13<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr14<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr15<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr16<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr17<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr18<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr19<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr20<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr21<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
nf18_p_chr22<-as.data.frame(read.delim("./MoBa_results/nf18_p_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

norment_feb18_dad<-rbind(nf18_p_chr1,nf18_p_chr2,nf18_p_chr3,nf18_p_chr4,nf18_p_chr5,nf18_p_chr6,nf18_p_chr7,nf18_p_chr8,nf18_p_chr9,nf18_p_chr10,nf18_p_chr11,
                         nf18_p_chr12,nf18_p_chr13,nf18_p_chr1,nf18_p_chr15,nf18_p_chr16,nf18_p_chr17,nf18_p_chr18,nf18_p_chr19,nf18_p_chr20,nf18_p_chr21,nf18_p_chr22)
norment_feb18_dad$LOG_OR<-log(norment_feb18_dad$OR)
write.table(norment_feb18_dad,"./MoBa_results/norment_feb18_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(norment_feb18_dad,file="./MoBa_results/norment_feb18_dad.RData")


### NORMENT_MAY16 - FATHERS ###

nm16_p_chr1<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_1.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr2<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_2.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr3<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_3.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr4<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_4.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr5<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_5.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr6<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_6.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr7<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_7.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr8<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_8.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr9<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_9.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr10<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_10.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr11<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_11.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr12<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_12.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr13<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_13.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr14<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_14.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr15<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_15.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr16<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_16.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr17<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_17.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr18<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_18.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr19<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_19.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr20<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_20.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr21<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_21.infertility.glm.logistic",header=TRUE,sep="\t"))
nm16_p_chr22<-as.data.frame(read.delim("./MoBa_results/nm16_p_chr_22.infertility.glm.logistic",header=TRUE,sep="\t"))

norment_may16_dad<-rbind(nm16_p_chr1,nm16_p_chr2,nm16_p_chr3,nm16_p_chr4,nm16_p_chr5,nm16_p_chr6,nm16_p_chr7,nm16_p_chr8,nm16_p_chr9,nm16_p_chr10,nm16_p_chr11,
                         nm16_p_chr12,nm16_p_chr13,nm16_p_chr1,nm16_p_chr15,nm16_p_chr16,nm16_p_chr17,nm16_p_chr18,nm16_p_chr19,nm16_p_chr20,nm16_p_chr21,nm16_p_chr22)
norment_may16_dad$LOG_OR<-log(norment_may16_dad$OR)
write.table(norment_may16_dad,"./MoBa_results/norment_may16_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)
save(norment_may16_dad,file="./MoBa_results/norment_may16_dad.RData")


###########################################################################
### FORMATING PLINK2 GWAS OUTPUT TXT FILES: COMPATIBLE WITH GWAMA-Linux ###
###########################################################################

z<-qnorm(1-0.05/2)

vars01<-c("harvest_mom","rotterdam1_mom","rotterdam2_mom","norment_feb18_mom","norment_may16_mom",
          "harvest_dad","rotterdam1_dad","rotterdam2_dad","norment_feb18_dad","norment_may16_dad")

for(i in 1:length(vars01))
{
  name<-paste("./MoBa_results/",vars01[i],".txt",sep="")
  name2<-paste("./MoBa_results/",vars01[i],"_gwama.txt",sep="")
  dat<-as.data.frame(read.delim(name,header=TRUE,sep="\t"))
  dat$OR_SE<-dat$OR*dat$LOG.OR._SE
  dat$OR_95L<-dat$OR-(z*dat$OR_SE)
  dat$OR_95U<-dat$OR+(z*dat$OR_SE)
  dat<-rename.vars(dat,from=c("ID","REF","ALT","OBS_CT"),to=c("MARKERNAME","NEA","EA","N"))
  dat<-dat[,c("MARKERNAME","EA","NEA","OR","OR_95L","OR_95U","N")]
  write.table(dat,name2,sep="\t",row.names=FALSE, quote=FALSE)
}


########################################################################
### META-ANALYSIS OF GWAS RESULTS IN 5 STUDY BADGES (GWAMA IN Linux) ###
########################################################################


#################################################################################
### EXTRACTION OF SNPs SIGNIFICANTLY LINKED TO BMI IN GIANT + TO SUBFERTILITY ###
#################################################################################

setwd("N:/data/durable/Projects/Magnus_MR_BMI/gwas")

# EXPOSURE : BMI #

bmi<-read.csv2("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/gwas_bmi_cojo_locke2018.csv",
                   header=TRUE,sep=",",dec=".")
names(bmi)<-tolower(names(bmi))
bmi$Phenotype<-c("BMI")
bmi$units<-c("units")
bmi$id<-c("Yengo 2018")
bmi<-rename.vars(bmi,
                 from=c("snp","chr","tested_allele","freq_tested_allele_in_hrs","p","n"),
                 to=c("SNP","chrom","effect_allele","eaf","pval","samplesize"))
bmi<-bmi[,c("Phenotype","SNP","beta","se","eaf","effect_allele","other_allele","pval",
            "units","samplesize","id")]
write.table(bmi,"./twosampleMR/bmi.txt",sep="\t",row.names=FALSE, quote=FALSE)

bmi_harm<-read_exposure_data(filename = "./twosampleMR/bmi.txt",
                                 clump = FALSE,
                                 sep="\t",
                                 phenotype_col = "Phenotype",
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "se",
                                 eaf_col = "eaf",
                                 effect_allele_col = "effect_allele",
                                 other_allele_col = "other_allele",
                                 pval_col = "pval",
                                 units_col = "units",
                                 samplesize_col = "samplesize",
                                 id_col = "id",
                                 min_pval = 1e-200,
                                 log_pval = FALSE)
save(bmi_harm,file="./twosampleMR/bmi_harm.RData")


# OUTCOME : Subfertility #

# Women #

subf_mom<-as.data.frame(read.delim("./metaanalysis/gwama_mom.txt.out",header=TRUE,sep="\t"))
subf_mom<-rename.vars(subf_mom,
                      from=c("rs_number","reference_allele","OR","OR_se","p.value","n_samples"),
                      to=c("SNP","effect_allele","or","or_se","pval","samplesize"))
subf_mom<-subf_mom[,c("SNP","effect_allele","other_allele","eaf","or","or_se","pval","samplesize")]
subf_mom$beta<-log(subf_mom$or)
subf_mom$se<-subf_mom$or_se/subf_mom$or
subf_mom$Units<-c("Units")
subf_mom$Phenotype<-c("Subfertility")
write.table(subf_mom,"./twosampleMR/subf_mom.txt",sep="\t",row.names=FALSE, quote=FALSE)

subf_mom_harm <- read_outcome_data(
  snps = bmi_harm$SNP,
  filename = "./twosampleMR/subf_mom.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

save(subf_mom_harm,file="./twosampleMR/subf_mom_harm.RData")

dat<-harmonise_data(bmi_harm, subf_mom_harm, action = 2)
save(dat,file="./twosampleMR/bmi_subf_mom.RData")


# Men #

subf_dad<-as.data.frame(read.delim("./metaanalysis/gwama_dad.txt.out",header=TRUE,sep="\t"))
subf_dad<-rename.vars(subf_dad,
                      from=c("rs_number","reference_allele","OR","OR_se","p.value","n_samples"),
                      to=c("SNP","effect_allele","or","or_se","pval","samplesize"))
subf_dad<-subf_dad[,c("SNP","effect_allele","other_allele","eaf","or","or_se","pval","samplesize")]
subf_dad$beta<-log(subf_dad$or)
subf_dad$se<-subf_dad$or_se/subf_dad$or
subf_dad$Units<-c("Units")
subf_dad$Phenotype<-c("Subfertility")
write.table(subf_dad,"./twosampleMR/subf_dad.txt",sep="\t",row.names=FALSE, quote=FALSE)

subf_dad_harm <- read_outcome_data(
  snps = bmi_harm$SNP,
  filename = "./twosampleMR/subf_dad.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

save(subf_dad_harm,file="./twosampleMR/subf_dad_harm.RData")

dat<-harmonise_data(bmi_harm, subf_dad_harm, action = 2)
save(dat,file="./twosampleMR/bmi_subf_dad.RData")


### TWO SAMPLE MRs ###

setwd("N:/data/durable/Projects/Magnus_MR_BMI/gwas/twosampleMR")

# Women #

load("./bmi_subf_mom.RData")
mr_results<-mr(dat,
               method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
mr_heterogenq<-mr_heterogeneity(dat)
mr_q<-rbind(round(mr_heterogenq[2,6],2),round(mr_heterogenq[1,6],2),NA,NA)
mr_qpval<-rbind(pval_guapa(mr_heterogenq[2,8]),pval_guapa(mr_heterogenq[1,8]),NA,NA)
mr_pleiotropytest<-mr_pleiotropy_test(dat)
mr_pleiot<-rbind(NA,pval_guapa(mr_pleiotropytest[1,7]),NA,NA)

jpeg("./results/bmi_subf_mom_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/bmi_subf_mom_forestplot.jpg", width = 4000, height = 50000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/bmi_subf_mom_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()

mr_results$pval<-pval_guapa(mr_results$pval)
mr_results$beta<-paste(guapa(mr_results$b)," (",guapa(mr_results$se),")",sep="")
mr_results$q<-mr_q
mr_results$q_pval<-mr_qpval
mr_results$pleiotropy<-mr_pleiot
mr_results<-mr_results[,c("method","beta","pval","q","q_pval","pleiotropy")]
write.table(mr_results,file="./results/bmi_subf_mom.csv",sep=";",row.names=FALSE)

# Men #

load("./bmi_subf_dad.RData")
mr_results<-mr(dat,
               method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
mr_heterogenq<-mr_heterogeneity(dat)
mr_q<-rbind(round(mr_heterogenq[2,6],2),round(mr_heterogenq[1,6],2),NA,NA)
mr_qpval<-rbind(pval_guapa(mr_heterogenq[2,8]),pval_guapa(mr_heterogenq[1,8]),NA,NA)
mr_pleiotropytest<-mr_pleiotropy_test(dat)
mr_pleiot<-rbind(NA,pval_guapa(mr_pleiotropytest[1,7]),NA,NA)

jpeg("./results/bmi_subf_dad_scatterplot.jpg", width = 9000, height = 9000, res=1200)
mr_scatter_plot(mr_results, dat)
dev.off()

jpeg("./results/bmi_subf_dad_forestplot.jpg", width = 4000, height = 50000, res=600)
mr_forest_plot(mr_singlesnp(dat))
dev.off()

jpeg("./results/bmi_subf_dad_funnelplot.jpg", width = 6000, height = 6000, res=1200)
mr_funnel_plot(mr_singlesnp(dat))
dev.off()

mr_results$pval<-pval_guapa(mr_results$pval)
mr_results$beta<-paste(guapa(mr_results$b)," (",guapa(mr_results$se),")",sep="")
mr_results$q<-mr_q
mr_results$q_pval<-mr_qpval
mr_results$pleiotropy<-mr_pleiot
mr_results<-mr_results[,c("method","beta","pval","q","q_pval","pleiotropy")]
write.table(mr_results,file="./results/bmi_subf_dad.csv",sep=";",row.names=FALSE)

