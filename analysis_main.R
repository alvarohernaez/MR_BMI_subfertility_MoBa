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

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

# Fractional polynomial models (Stephen Burgess) #
source("N:/data/durable/Projects/Magnus_MR_BMI/R/Old/nlme_summ_aes.R")


#####################
### MAIN ANALYSES ###
#####################

setwd("N:/data/durable/Projects/Magnus_MR_BMI/R")

dir.create("./Outputs")
dir.create("./Outputs/descriptive")
dir.create("./Outputs/results")
dir.create("./Outputs/nlmr")
dir.create("./Outputs/nlmr/source")

setwd("N:/data/durable/Projects/Magnus_MR_BMI/R/Outputs")


### POPULATION DESCRIPTION ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")
dat$bmi_mom_cat<-with(dat,ifelse(bmi_mom<20,1,
                                 ifelse(bmi_mom<25,2,
                                        ifelse(bmi_mom<30,3,
                                               ifelse(bmi_mom>=30,4,NA)))))

datx<-subset2(dat,"dat$anthrop_na_mom==0 & !is.na(dat$bmi_grs_mom)")

xxx<-datx[,c("agedelivery_mom","eduyears_mom","smoking_mom","parity","bmi_mom","bmi_mom_cat","subfertility_12plus")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.
                               -subfertility_12plus,
                               xxx, method=c("bmi_mom"=2,"smoking_mom"=3,"bmi_mom_cat"=3,"parity"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=NA)

comp<-NULL
comp<-createTable(compareGroups(subfertility_12plus~.
                                -sel,
                                xxx, method=c("bmi_mom"=2,"smoking_mom"=3,"bmi_mom_cat"=3,"parity"=3)),
                  show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=NA)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1],comp$descr))
colnames(tab1)<-c("Mothers-All","Mothers-Non-subfertile","Mothers-Subfertile","Mothers-P-value","Mothers-N")

dat$bmi_dad_cat<-with(dat,ifelse(bmi_dad<20,1,
                                 ifelse(bmi_dad<25,2,
                                        ifelse(bmi_dad<30,3,
                                               ifelse(bmi_dad>=30,4,NA)))))
datx<-subset2(dat,"dat$anthrop_na_dad==0 & !is.na(dat$bmi_grs_dad)")

xxx<-datx[,c("agedelivery_dad","eduyears_dad","smoking_dad","parity","bmi_dad","bmi_dad_cat","subfertility_12plus")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.
                               -subfertility_12plus,
                               xxx, method=c("bmi_dad"=2,"smoking_dad"=3,"bmi_dad_cat"=3,"parity"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=NA)

comp<-NULL
comp<-createTable(compareGroups(subfertility_12plus~.
                                -sel,
                                xxx, method=c("bmi_dad"=2,"smoking_dad"=3,"bmi_dad_cat"=3,"parity"=3)),
                  show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=NA)

tab2<-NULL
tab2<-as.data.frame(cbind(all$descr[,1],comp$descr))
colnames(tab2)<-c("Fathers-All","Fathers-Non-subfertile","Fathers-Subfertile","Fathers-P-value","Fathers-N")

tab<-NULL
tab<-cbind(tab1,tab2)
rownames(tab)<-c("Age at delivery (years)","Education years",
                 "Smoking: 0","Smoking: 1","Smoking: 2","Smoking: 3",
                 "Parity: 0","Parity: 1","Parity: 2","Parity: 3",
                 "BMI GRS (units)","Low weight","Normal weight","Overweight","Obesity")

write.table(tab,file="./descriptive/descriptive_genotyped.csv",sep=";")


### ASSOCIATION BETWEEN BMI/HEIGHT/EDUYEARS AND RISK SCORE ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

# vars04: check the number of SNPs used for each GRS (BMI n=896)

vars01<-c("bmi_mom","bmi_dad")
vars02<-c("bmi_grs_mom","bmi_grs_dad")
vars03<-c("m_id_2374","f_id_2374")
vars04<-c("896","896")
vars05<-c("anthrop_na_mom","anthrop_na_dad")

tab<-NULL

for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars05[i]]==0 & !is.na(dat[,vars02[i]])")
  mod01<-lm_robust(datx[,vars01[i]]~datx[,vars02[i]], 
                   data=datx, clusters=datx[,vars03[i]], se_type="stata")
  coef<-paste(guapa(summary(mod01)$coefficients[2,1])," [",
              guapa(summary(mod01)$coefficients[2,5]),"; ",
              guapa(summary(mod01)$coefficients[2,6]),"]",sep="")
  fstat<-round(mod01$fstatistic[1],0)
  r2<-paste("'",guapa(mod01$adj.r.squared*100),sep="")
  
  datx<-subset2(datx,"!is.na(datx[,vars01[i]]) & !is.na(datx[,vars02[i]])")
  n_obs<-length(which(!is.na(datx[,vars03[i]])))
  n_snps<-vars04[i]
  mean_grs<-paste(guapa(mean(datx[,vars02[i]]))," (",guapa(sd(datx[,vars02[i]])),")",sep="")
  tab<-rbind(tab,cbind(n_obs,n_snps,mean_grs,coef,fstat,r2))
}

rownames(tab)<-vars01
write.table(tab,file="./grs_anthrop_assoc/assoc_grs_bmi.csv",sep=";")


### ASSOCIATION BETWEEN BMI/HEIGHT/EDUYEARS AND RISK SCORE: SMOOTHED SPLINES ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

vars01<-c("bmi_mom","bmi_dad")
vars02<-c("BMI (kg/m2), mothers","BMI (kg/m2), fathers")
vars03<-c("bmi_grs_mom","bmi_grs_dad")
vars04<-c("BMI GRS, mothers","BMI GRS, fathers")
vars05<-c("anthrop_na_mom","anthrop_na_dad")

for(i in 1:length(vars01))
  
{
  aaa<-dat[,vars03[i]]
  bbb<-dat[,vars01[i]]
  dat2<-subset2(dat,"!is.na(aaa) & !is.na(bbb) & dat[,vars05[i]]==0")
  aaa<-dat2[,vars03[i]]
  model<-glm(dat2[,vars01[i]]~bs(aaa),data=dat2)
  mod01<-glm(dat2[,vars01[i]]~aaa,data=dat2,family=gaussian)
  p_lin<-guapa(summary(mod01)$coefficients[2,4])
  p_nonlin<-guapa(lrtest(model,mod01)[2,5])
  
  z<-qnorm(1-0.05/2)
  res<-termplot(model,term='bs(aaa)',rug=FALSE,se=TRUE,plot=FALSE)
  res<-res$aaa
  corr<-res$y[1]
  ajus<-mean(dat2[,vars01[i]][which(dat2[,vars03[i]]==res$x[1])],na.rm=TRUE)
  res$y<-res$y-corr+ajus
  ci<-cbind(res$y,res$y-z*res$se,res$y+z*res$se)
  
  alc_n<-paste("./grs_anthrop_assoc/spline_",vars01[i],".jpg",sep="")
  xlabel<-vars04[i]
  ylabel<-vars02[i]
  
  jpeg(filename=alc_n,width=3000,height=3000,res=600,pointsize=11.5)
  
  par(las=1,cex=1,mar=c(5,5,2,2))
  matplot(res$x, ci, lty=c(1,0,0),lwd=2,type="l",col="black",xlab=xlabel,ylab=ylabel)
  polygon(c(res$x,rev(res$x)),c(ci[,2],rev(ci[,3])),col=grey(0.8),border="white")
  matpoints(res$x,ci,lty=c(1,0,0),lwd=2,type="l",col="black")
  abline(h=0, lty=3, col="black")
  abline(v=0, lty=3, col="black")
  leg<-paste("Linear component: P-value = ",p_lin,"\nNon-linear component: P-value = ",p_nonlin,sep="")
  legend("bottomleft",legend=leg, bty="n", cex=0.7, inset=c(0.025,0.025))
  
  dev.off()
  
}


### ASSOCIATION OF GRS WITH COVARIATES ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

datx<-subset2(dat,"dat$anthrop_na_mom==0 & !is.na(dat$bmi_grs_mom)")
datx$bmi_grs_mom_q<-as.numeric(ntile(datx$bmi_grs_mom, 4))
datx$height_grs_mom_q<-as.numeric(ntile(datx$height_grs_mom, 4))
datx$smoking_mom2<-with(datx,ifelse(smoking_mom==0,0,1))
datx$parity2<-with(datx,ifelse(parity==0,0,1))

xxx<-datx[,c("agedelivery_mom","eduyears_mom","smoking_mom2","parity2","bmi_grs_mom_q")]

bmi_grs_mom<-NULL
bmi_grs_mom<-createTable(compareGroups(bmi_grs_mom_q~.
                                       -height_grs_mom_q,
                                       xxx, method=c("smoking_mom2"=3,"parity2"=3)),
                         show.n=FALSE, show.p.overall=TRUE, show.p.trend=TRUE, hide.no=NA)
export2csv(bmi_grs_mom,file="./descriptive/descr_bmi_grs_mom.csv",sep=";")

datx<-subset2(dat,"dat$anthrop_na_dad==0 & !is.na(dat$bmi_grs_dad)")
datx$bmi_grs_dad_q<-as.numeric(ntile(datx$bmi_grs_dad, 4))
datx$smoking_dad2<-with(datx,ifelse(smoking_dad==0,0,1))
datx$parity2<-with(datx,ifelse(parity==0,0,1))

xxx<-datx[,c("agedelivery_dad","eduyears_dad","smoking_dad2","parity2","bmi_grs_dad_q")]

bmi_grs_dad<-NULL
bmi_grs_dad<-createTable(compareGroups(bmi_grs_dad_q~.
                                       -height_grs_dad_q,
                                       xxx, method=c("smoking_dad2"=3,"parity2"=3)),
                         show.n=FALSE, show.p.overall=TRUE, show.p.trend=TRUE, hide.no=NA)
export2csv(bmi_grs_dad,file="./descriptive/descr_bmi_grs_dad.csv",sep=";")


### LOGISTIC REGRESSION / MENDELIAN RANDOMIZATION: LINEAR ASSOCIATIONS ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

dat$var_null<-1

vars01<-c("gen_pred_bmi_mom","gen_pred_bmi_dad")
vars02<-c("bmi_mom","bmi_dad")
vars05<-c("agedelivery_mom","agedelivery_dad")
vars06<-c("eduyears_mom","eduyears_dad")
vars07<-c("smoking_mom","smoking_dad")
vars08<-c("pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374")
vars19<-c("anthrop_na_mom","anthrop_na_dad")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars01[i]])")
  mod01<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars02[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod01)[2,1])
  se<-as.numeric(summary(mod01)[2,2])
  coef01<-risk_se_ic_guapa(estimate,se)
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  
  mod02<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars02[i]]
                               +datx[,vars05[i]]+datx[,vars06[i]]+datx[,vars07[i]]+parity+datx[,vars20[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod02)[2,1])
  se<-as.numeric(summary(mod02)[2,2])
  coef02<-risk_se_ic_guapa(estimate,se)
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  
  mod03<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod03)[2,1])
  se<-as.numeric(summary(mod03)[2,2])
  coef03<-risk_se_ic_guapa(estimate,se)
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~datx[,vars01[i]]
                               +datx[,vars08[i]]+datx[,vars09[i]]+datx[,vars10[i]]+datx[,vars11[i]]
                               +datx[,vars12[i]]+datx[,vars13[i]]+datx[,vars14[i]]
                               +datx[,vars15[i]]+datx[,vars16[i]]+datx[,vars17[i]],
                               data=datx, cluster=datx[,vars18[i]], family="binomial")
  estimate<-as.numeric(summary(mod04)[2,1])
  se<-as.numeric(summary(mod04)[2,2])
  coef04<-risk_se_ic_guapa(estimate,se)
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  tab<-rbind(tab,cbind(coef01,pval01,coef02,pval02,coef03,pval03,coef04,pval04))
}

colnames(tab)<-c("Det-OR (raw)","Det-pval (raw)","Det-OR (adj.)","Det-pval (adj.)",
                 "GRS-OR (raw)","GRS-pval (raw)","GRS-OR (adj.)","GRS-pval (adj.)")
rownames(tab)<-c("BMI, mothers","BMI, fathers","Height, mothers","Height, fathers")
write.table(tab,file="./results/subfertility_linear.csv",sep=";")


### LOGISTIC REGRESSION: NON-LINEAR ASSOCIATIONS (SMOOTHED SPLINES) ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")
dat$var_null<-1

# EXECUTE MANUALLY: i=1, then the loop, repeat for i=2/3/4

vars01<-c("bmi_mom","bmi_dad")
vars02<-c("bmi_grs_mom_z","bmi_grs_dad_z")
vars04<-c("agedelivery_mom","agedelivery_dad")
vars05<-c("eduyears_mom","eduyears_dad")
vars06<-c("smoking_mom","smoking_dad")
vars07<-c("pc1_mom","pc1_dad")
vars08<-c("pc2_mom","pc2_dad")
vars09<-c("pc3_mom","pc3_dad")
vars10<-c("pc4_mom","pc4_dad")
vars11<-c("pc5_mom","pc5_dad")
vars12<-c("pc6_mom","pc6_dad")
vars13<-c("pc7_mom","pc7_dad")
vars14<-c("pc8_mom","pc8_dad")
vars15<-c("pc9_mom","pc9_dad")
vars16<-c("pc10_mom","pc10_dad")
vars17<-c("Body mass index (kg/m2)","Body mass index (kg/m2)")
vars19<-c("anthrop_na_mom","anthrop_na_dad")
vars20<-c(25,25)
vars21<-c(17,17)
vars22<-c(41,41)
vars23<-c(0.5,0.5)
vars24<-c(2.5,2.5)


for(i in 1:length(vars01))
  
{
  varstot<-c(vars01[i],vars02[i],vars04[i],vars05[i],vars06[i],vars19[i],"parity","subfertility_12plus",vars25[i])
  dat2<-na.omit(dat[,varstot])
  dat2<-subset2(dat2,"dat2[,vars19[i]]==0")
  aaa<-dat2[,vars01[i]]
  
  mod01<-glm(formula=as.factor(subfertility_12plus)~bs(aaa,degree=3),
             data=dat2, family="binomial")
  mod02<-glm(formula=as.factor(subfertility_12plus)~aaa,
             data=dat2, family="binomial")
  
  ptemp<-termplot(mod01,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,vars20[i])
  center<-with(temp, y[x==value])
  z<-qnorm(1-0.05/2)
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  ci<-exp(ytemp-center)
  name<-paste("./results/spline_",vars01[i],"_raw.jpg",sep="")
  labely<-c("Subfertility (odds ratio, unadjusted)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  
  figure<-ggplot(plot.data, aes(x=x, ymax=5))
  figure<-figure + geom_hline(aes(yintercept=1), colour="black", linetype=2) + 
    geom_line(aes(y=yest), color="black") + 
    geom_line(aes(y=lci), color="grey") + 
    geom_line(aes(y=uci), color="grey") + 
    scale_x_continuous(limits = c(vars21[i],vars22[i])) +
    scale_y_continuous(limits = c(vars23[i],vars24[i])) +
    theme_bw() +
    labs(x=vars17[i],y=labely) +
    theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20),
          axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + 
    geom_point(aes(x=vars20[i], y=1),data=plot.data, colour="black", size=3) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  jpeg(filename=name,width = 8000, height = 8000, res=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
  
  
  mod01<-glm(formula=as.factor(subfertility_12plus)~bs(aaa,degree=3)
             +dat2[,vars04[i]]+dat2[,vars05[i]]+dat2[,vars06[i]]+parity+dat2[,vars25[i]],
             data=dat2, family="binomial")
  mod02<-glm(formula=as.factor(subfertility_12plus)~aaa
             +dat2[,vars04[i]]+dat2[,vars05[i]]+dat2[,vars06[i]]+parity+dat2[,vars25[i]],
             data=dat2, family="binomial")
  p_lin<-guapa(summary(mod02)$coefficients[2,4])
  p_nonlin<-guapa(lrtest(mod01,mod02)[2,5])
  p_lin
  p_nonlin
  
  ptemp<-termplot(mod01,term=1,se=TRUE,plot=FALSE)
  temp<-ptemp$aaa
  value<-closest(temp$x,vars20[i])
  #min_val<-temp$x[which(temp$y==min(temp$y,na.rm=TRUE))]
  center<-with(temp, y[x==value])
  z<-qnorm(1-0.05/2)
  ytemp<-temp$y+outer(temp$se,c(0,-z,z),'*')
  min_val<-guapa(temp$x[which(temp$y==min(temp$y,na.rm=TRUE))])
  ci<-exp(ytemp-center)
  name<-paste("./results/spline_",vars01[i],"_adj.jpg",sep="")
  labely<-c("Subfertility (odds ratio, adjusted)")
  
  plot.data<-as.data.frame(cbind(temp$x,ci))
  colnames(plot.data)<-c("x","yest","lci","uci")
  
  figure<-ggplot(plot.data, aes(x=x, ymax=5))
  figure<-figure + geom_hline(aes(yintercept=1), colour="black", linetype=2) + 
    geom_line(aes(y=yest), color="black") + 
    geom_line(aes(y=lci), color="grey") + 
    geom_line(aes(y=uci), color="grey") + 
    scale_x_continuous(limits = c(vars21[i],vars22[i])) +
    scale_y_continuous(limits = c(vars23[i],vars24[i])) +
    theme_bw() +
    labs(x=vars17[i],y=labely) +
    theme(axis.title.x = element_text(vjust=0.5, size=20), axis.title.y = element_text(vjust=0.5, size=20),
          axis.text.x=element_text(size=18), axis.text.y=element_text(size=18)) + 
    geom_point(aes(x=vars20[i], y=1),data=plot.data, colour="black", size=3) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  jpeg(filename=name,width = 8000, height = 8000, res=1200)
  par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
  figure
  dev.off()
}


### LINEAR LOGISTIC REGRESSION IN STRATA ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

vars01<-c("bmi_mom","bmi_dad")
vars02<-c("bmi_grs_mom","bmi_grs_dad")
vars05<-c("agedelivery_mom","agedelivery_dad")
vars06<-c("eduyears_mom","eduyears_dad")
vars07<-c("smoking_mom","smoking_dad")
vars18<-c("m_id_2374","f_id_2374")
vars19<-c("anthrop_na_mom","anthrop_na_dad")

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars02[i]])")
  dat03<-subset2(datx,"datx[,vars01[i]]<20")
  dat04<-subset2(datx,"datx[,vars01[i]]>=20 & datx[,vars01[i]]<25")
  dat05<-subset2(datx,"datx[,vars01[i]]>=25 & datx[,vars01[i]]<30")
  dat06<-subset2(datx,"datx[,vars01[i]]>=30")

  mod03<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat03[,vars01[i]]
                               +dat03[,vars05[i]]+dat03[,vars06[i]]+dat03[,vars07[i]]+parity,
                               data=dat03, cluster=dat03[,vars18[i]], family="binomial")
  coef03<-risk_se_ic_guapa(summary(mod03)[2,1],summary(mod03)[2,2])
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat04[,vars01[i]]
                               +dat04[,vars05[i]]+dat04[,vars06[i]]+dat04[,vars07[i]]+parity,
                               data=dat04, cluster=dat04[,vars18[i]], family="binomial")
  coef04<-risk_se_ic_guapa(summary(mod04)[2,1],summary(mod04)[2,2])
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  mod05<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat05[,vars01[i]]
                               +dat05[,vars05[i]]+dat05[,vars06[i]]+dat05[,vars07[i]]+parity,
                               data=dat05, cluster=dat05[,vars18[i]], family="binomial")
  coef05<-risk_se_ic_guapa(summary(mod05)[2,1],summary(mod05)[2,2])
  pval05<-pval_guapa(as.numeric(summary(mod05)[2,4]))
  
  mod06<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat06[,vars01[i]]
                               +dat06[,vars05[i]]+dat06[,vars06[i]]+dat06[,vars07[i]]+parity,
                               data=dat06, cluster=dat06[,vars18[i]], family="binomial")
  coef06<-risk_se_ic_guapa(summary(mod06)[2,1],summary(mod06)[2,2])
  pval06<-pval_guapa(as.numeric(summary(mod06)[2,4]))
  
  tab<-cbind(tab,rbind(cbind(coef03,pval03),cbind(coef04,pval04),cbind(coef05,pval05),cbind(coef06,pval06)))
}

colnames(tab)<-c("Women-OR","Women-pval","Men-OR","Men-pval")
rownames(tab)<-c("<20.0","20.0-24.9","25.0-29.9",">=30.0")
write.table(tab,file="./results/subfertility_logreg_strata.csv",sep=";")


### NON-LINEAR MR ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

vars00<-c("environ_bmi_mom","environ_bmi_dad")
vars01<-c("environ_bmi_mom_q","environ_bmi_dad_q")
vars02<-c("bmi_mom","bmi_dad")
vars03<-c("bmi_grs_mom_z","bmi_grs_dad_z")
vars04<-c("environ_bmi_mom","environ_bmi_dad")
vars08<-c("pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374")
vars19<-c("anthrop_na_mom","anthrop_na_dad")
vars20<-c(100,100,1,1)

for(i in 1:length(vars01))
  
{
  series<-1:vars20[i]
  tab<-NULL
  tab2<-NULL
  
  for(j in 1:length(series))
    
  {
    seq<-as.numeric(series[j])
    dat2<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars03[i]])")
    dat2[,vars01[i]]<-as.numeric(ntile(dat2[,vars00[i]], vars20[i]))
    dat2<-subset2(dat2,"dat2[,vars01[i]]==seq")
    mean<-mean(dat2[,vars02[i]],na.rm=TRUE)
    
    mod01<-lm_robust(dat2[,vars02[i]]~dat2[,vars03[i]]
                     +dat2[,vars08[i]]+dat2[,vars09[i]]+dat2[,vars10[i]]+dat2[,vars11[i]]+dat2[,vars12[i]]
                     +dat2[,vars13[i]]+dat2[,vars14[i]]+dat2[,vars15[i]]+dat2[,vars16[i]]+dat2[,vars17[i]],
                     data=dat2, clusters=dat2[,vars18[i]], se_type="stata")
    beta_exp<-as.numeric(summary(mod01)$coefficients[2,1])
    se_exp<-as.numeric(summary(mod01)$coefficients[2,2])
    mod02<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat2[,vars03[i]]
                                 +dat2[,vars08[i]]+dat2[,vars09[i]]+dat2[,vars10[i]]+dat2[,vars11[i]]+dat2[,vars12[i]]
                                 +dat2[,vars13[i]]+dat2[,vars14[i]]+dat2[,vars15[i]]+dat2[,vars16[i]]+dat2[,vars17[i]],
                                 data=dat2, cluster=dat2[,vars18[i]], family="binomial")
    beta_outc<-as.numeric(summary(mod02)[2,1])
    se_outc<-as.numeric(summary(mod02)[2,2])
    tab<-rbind(tab,cbind(mean,beta_exp,se_exp,beta_outc,se_outc))
    
    mod03<-lm_robust(dat2[,vars02[i]]~dat2[,vars03[i]],
                     data=dat2, clusters=dat2[,vars18[i]], se_type="stata")
    beta_exp<-as.numeric(summary(mod03)$coefficients[2,1])
    se_exp<-as.numeric(summary(mod03)$coefficients[2,2])
    mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat2[,vars03[i]],
                                 data=dat2, cluster=dat2[,vars18[i]], family="binomial")
    beta_outc<-as.numeric(summary(mod04)[2,1])
    se_outc<-as.numeric(summary(mod04)[2,2])
    tab2<-rbind(tab2,cbind(mean,beta_exp,se_exp,beta_outc,se_outc))
    
    
  }
  
  namefile<-paste("./nlmr/source/",vars04[i],".csv",sep="")
  write.table(tab,file=namefile,sep=";")
  
  namefile2<-paste("./nlmr/source/",vars04[i],"_raw.csv",sep="")
  write.table(tab2,file=namefile2,sep=";")
  
}


### NON-LINEAR MR GRAPHS ###

data  <-as.data.frame(read.csv("./nlmr/source/environ_bmi_mom.csv",header=TRUE,sep=";",dec="."))
data<-data[order(data$mean),]

by    <-data$beta_outc  # genetic associations with outcome per quantile
byse  <-data$se_outc    # standard errors of genetic associations with outcome
bx    <-data$beta_exp   # genetic associations with exposure per quantile
bxse  <-data$se_exp     # standard errors of genetic associations with exposure
xmean <-data$mean       # mean BMI in each quantile

floorn   <-floor(min(data$mean,na.rm=TRUE))
ceilingn <-ceiling(max(data$mean,na.rm=TRUE))
xmean1   <-xmean-floorn
refx     <-25-floorn

results<-frac_poly_summ_mr(by, bx, byse, bxse, xmean1, 
                           pd=0.5, ref=refx, d=2, ylim_lower=0.8, ylim_upper=8, offset=floorn, 
                           xlim_upper=ceilingn, fig=TRUE, 
                           pref_x="Body mass index (kg/m2)", 
                           pref_y="Subfertility (odds ratio, adjusted)", 
                           breaks=c(1,2,4,8))
tab_res<-results$figure$data
summary(tab_res$x)

vars01<-c(18,20,22.5,25,27.5,30,32.5,35,37.5,40)

tab<-NULL
for(i in 1:length(vars01))
  
{
  line<-tab_res[which(tab_res$x==closest(tab_res$x,vars01[i])),]
  coef<-ic_guapa(guapa(line$yest),guapa(line$lci),guapa(line$uci))
  zval<-log(line$yest)/line$yse
  pval<-pval_guapa(exp(-0.717*zval-(0.416*(zval^2))))
  value<-paste(guapa(closest(tab_res$x,vars01[i]))," kg/m2",sep="")
  tab<-rbind(tab,cbind(value,coef,pval))
}

write.table(tab,file="./nlmr/coefs_nlmr_bmi_mom.csv",sep=";")

jpeg("./nlmr/nlmr_bmi_mom.jpg", width = 8000, height = 8000, res=1200)
par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
results
dev.off()

results$p_tests
guapa(tab_res[which(tab_res$yest==min(tab_res$yest,na.rm=TRUE)),]$x)


data  <-as.data.frame(read.csv("./nlmr/source/environ_bmi_dad.csv",header=TRUE,sep=";",dec="."))

by    <-data$beta_outc  # genetic associations with outcome per quantile
byse  <-data$se_outc    # standard errors of genetic associations with outcome
bx    <-data$beta_exp   # genetic associations with exposure per quantile
bxse  <-data$se_exp     # standard errors of genetic associations with exposure
xmean <-data$mean       # mean BMI in each quantile

floorn   <-floor(min(data$mean,na.rm=TRUE))
ceilingn <-ceiling(max(data$mean,na.rm=TRUE))
xmean1   <-xmean-floorn
refx     <-25-floorn

results<-frac_poly_summ_mr(by, bx, byse, bxse, xmean1, 
                           pd=0.5, ref=refx, d=2, ylim_lower=0.8, ylim_upper=16, offset=floorn, 
                           xlim_upper=ceilingn, fig=TRUE, 
                           pref_x="Body mass index (kg/m2)", 
                           pref_y="Subfertility (odds ratio, adjusted)", 
                           breaks=c(1,2,4,8,16))
tab_res<-results$figure$data
summary(tab_res$x)

vars01<-c(18,20,22.5,25,27.5,30,32.5,35,37.5,40)

tab<-NULL
for(i in 1:length(vars01))
  
{
  line<-tab_res[which(tab_res$x==closest(tab_res$x,vars01[i])),]
  coef<-ic_guapa(guapa(line$yest),guapa(line$lci),guapa(line$uci))
  zval<-log(line$yest)/line$yse
  pval<-pval_guapa(exp(-0.717*zval-(0.416*(zval^2))))
  value<-paste(guapa(closest(tab_res$x,vars01[i]))," kg/m2",sep="")
  tab<-rbind(tab,cbind(value,coef,pval))
}

write.table(tab,file="./nlmr/coefs_nlmr_bmi_dad.csv",sep=";")

jpeg("./nlmr/nlmr_bmi_dad.jpg", width = 8000, height = 8000, res=1200)
par(las=1,cex=1.2,mar=c(6,6,2,0),bty="n",lheight=0.9)
results
dev.off()

results$p_tests
guapa(tab_res[which(tab_res$yest==min(tab_res$yest,na.rm=TRUE)),]$x)


### LINEAR MR IN STRATA ACCORDING TO RESIDUAL BMI ###

load("N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")

# WOMEN q200: <20.0 (1-21), 20.0-24.9 (22-140), 25.0-29.9 (141-183), >=30.0-34.9 (184-200) - THRESHOLDS: 22-141-184 #
# MEN q200: <20.0 (1-2), 20.0-24.9 (3-88), 25.0-29.9 (89-180), >=30.0-34.9 (181-200) - THRESHOLDS:  3-89-181 #

vars00<-c("environ_bmi_mom","environ_bmi_dad")
vars01<-c("environ_bmi_mom_q","environ_bmi_dad_q")
vars03<-c("gen_pred_bmi_mom","gen_pred_bmi_dad")
vars08<-c("pc1_mom","pc1_dad")
vars09<-c("pc2_mom","pc2_dad")
vars10<-c("pc3_mom","pc3_dad")
vars11<-c("pc4_mom","pc4_dad")
vars12<-c("pc5_mom","pc5_dad")
vars13<-c("pc6_mom","pc6_dad")
vars14<-c("pc7_mom","pc7_dad")
vars15<-c("pc8_mom","pc8_dad")
vars16<-c("pc9_mom","pc9_dad")
vars17<-c("pc10_mom","pc10_dad")
vars18<-c("m_id_2374","f_id_2374")
vars19<-c("anthrop_na_mom","anthrop_na_dad")
vars20<-c(22,3)
vars21<-c(141,89)
vars22<-c(184,181)

tab<-NULL
for(i in 1:length(vars01))
  
{
  datx<-subset2(dat,"dat[,vars19[i]]==0 & !is.na(dat[,vars03[i]])")
  datx[,vars01[i]]<-as.numeric(ntile(datx[,vars00[i]], 200))
  
  dat01<-subset2(datx,"datx[,vars01[i]]<vars20[i]")
  dat02<-subset2(datx,"datx[,vars01[i]]>=vars20[i] & datx[,vars01[i]]<vars21[i]")
  dat03<-subset2(datx,"datx[,vars01[i]]>=vars21[i] & datx[,vars01[i]]<vars22[i]")
  dat04<-subset2(datx,"datx[,vars01[i]]>=vars22[i]")
  
  mod01<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat01[,vars03[i]]
                               +dat01[,vars08[i]]+dat01[,vars09[i]]+dat01[,vars10[i]]+dat01[,vars11[i]]+dat01[,vars12[i]]
                               +dat01[,vars13[i]]+dat01[,vars14[i]]+dat01[,vars15[i]]+dat01[,vars16[i]]+dat01[,vars17[i]],
                               data=dat01, cluster=dat01[,vars18[i]], family="binomial")
  coef01<-risk_se_ic_guapa(summary(mod01)[2,1],summary(mod01)[2,2])
  pval01<-pval_guapa(as.numeric(summary(mod01)[2,4]))
  
  mod02<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat02[,vars03[i]]
                               +dat02[,vars08[i]]+dat02[,vars09[i]]+dat02[,vars10[i]]+dat02[,vars11[i]]+dat02[,vars12[i]]
                               +dat02[,vars13[i]]+dat02[,vars14[i]]+dat02[,vars15[i]]+dat02[,vars16[i]]+dat02[,vars17[i]],
                               data=dat02, cluster=dat02[,vars18[i]], family="binomial")
  coef02<-risk_se_ic_guapa(summary(mod02)[2,1],summary(mod02)[2,2])
  pval02<-pval_guapa(as.numeric(summary(mod02)[2,4]))
  
  mod03<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat03[,vars03[i]]
                               +dat03[,vars08[i]]+dat03[,vars09[i]]+dat03[,vars10[i]]+dat03[,vars11[i]]+dat03[,vars12[i]]
                               +dat03[,vars13[i]]+dat03[,vars14[i]]+dat03[,vars15[i]]+dat03[,vars16[i]]+dat03[,vars17[i]],
                               data=dat03, cluster=dat03[,vars18[i]], family="binomial")
  coef03<-risk_se_ic_guapa(summary(mod03)[2,1],summary(mod03)[2,2])
  pval03<-pval_guapa(as.numeric(summary(mod03)[2,4]))
  
  mod04<-miceadds::glm.cluster(formula=as.factor(subfertility_12plus)~dat04[,vars03[i]]
                               +dat04[,vars08[i]]+dat04[,vars09[i]]+dat04[,vars10[i]]+dat04[,vars11[i]]+dat04[,vars12[i]]
                               +dat04[,vars13[i]]+dat04[,vars14[i]]+dat04[,vars15[i]]+dat04[,vars16[i]]+dat04[,vars17[i]],
                               data=dat04, cluster=dat04[,vars18[i]], family="binomial")
  coef04<-risk_se_ic_guapa(summary(mod04)[2,1],summary(mod04)[2,2])
  pval04<-pval_guapa(as.numeric(summary(mod04)[2,4]))
  
  tab<-cbind(tab,rbind(cbind(coef01,pval01),cbind(coef02,pval02),cbind(coef03,pval03),cbind(coef04,pval04)))
}

colnames(tab)<-c("Women-OR","Women-pval","Men-OR","Men-pval")
rownames(tab)<-c("<20.0","20.0-24.9","25.0-29.9",">=30.0")
write.table(tab,file="./results/subfertility_mr_strata.csv",sep=";")


