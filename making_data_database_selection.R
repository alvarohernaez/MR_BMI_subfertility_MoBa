rm(list=ls())

# Load required libraries #

# Shortcut format functions #

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

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }


# Fractional polynomial models (Prof. Stephen Burgess)
source("N:/data/durable/Projects/Magnus_MR_BMI/R/Old/nlme_summ_aes.R")

setwd("N:/data/durable/Projects/Magnus_MR_BMI/R")


######################################
### GENERATION OF WORKING DATABASE ###
######################################

### MERGING MOM/DAD GRSs ###

load("./bmi_grs.RData")
bmi_grs<-dat
load("./eduyears_grs.RData")
eduyears_grs<-dat
load("./smkinit_grs.RData")
smkinit_grs<-dat
load("./MoBa_raw.RData")
dat$sentrixid_mom<-with(dat,ifelse(sentrixid_mom=="",NA,sentrixid_mom))
dat$sentrixid_dad<-with(dat,ifelse(sentrixid_dad=="",NA,sentrixid_dad))

bmi_grs<-rename.vars(bmi_grs,from=c("id","bmi_grs"),to=c("sentrixid_mom","bmi_grs_mom"))
dat<-merge2(dat,bmi_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
bmi_grs<-rename.vars(bmi_grs,from=c("sentrixid_mom","bmi_grs_mom"),to=c("sentrixid_dad","bmi_grs_dad"))
dat<-merge2(dat,bmi_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

eduyears_grs<-rename.vars(eduyears_grs,from=c("id","eduyears_grs"),to=c("sentrixid_mom","eduyears_grs_mom"))
dat<-merge2(dat,eduyears_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
eduyears_grs<-rename.vars(eduyears_grs,from=c("sentrixid_mom","eduyears_grs_mom"),to=c("sentrixid_dad","eduyears_grs_dad"))
dat<-merge2(dat,eduyears_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

smkinit_grs<-rename.vars(smkinit_grs,from=c("id","smkinit_grs"),to=c("sentrixid_mom","smkinit_grs_mom"))
dat<-merge2(dat,smkinit_grs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
smkinit_grs<-rename.vars(smkinit_grs,from=c("sentrixid_mom","smkinit_grs_mom"),to=c("sentrixid_dad","smkinit_grs_dad"))
dat<-merge2(dat,smkinit_grs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

bmi_grs<-NULL
eduyears_grs<-NULL
smkinit_grs<-NULL

length(which(!is.na(dat$bmi_grs_mom)))
length(which(!is.na(dat$eduyears_grs_mom)))
length(which(!is.na(dat$smkinit_grs_mom)))
length(which(!is.na(dat$bmi_grs_dad)))
length(which(!is.na(dat$eduyears_grs_dad)))
length(which(!is.na(dat$smkinit_grs_dad)))


### MERGING QC INDICATORS ###

qc_indic<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/Stata/Marker_QC_Ind_Alex_group.txt",
                                   header=TRUE,sep=","))
qc_indic$V2<-NULL
qc_indic<-rename.vars(qc_indic,from=c("V1"),to=c("sentrixid_mom"))
qc_indic$qc_gen_mom<-1
dat<-merge2(dat,qc_indic,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
qc_indic<-rename.vars(qc_indic,from=c("sentrixid_mom","qc_gen_mom"),to=c("sentrixid_dad","qc_gen_dad"))
dat<-merge2(dat,qc_indic,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)
dat$qc_gen_mom<-with(dat,ifelse(is.na(qc_gen_mom),0,qc_gen_mom))
dat$qc_gen_dad<-with(dat,ifelse(is.na(qc_gen_dad),0,qc_gen_dad))

moms<-as.data.frame(read_dta("N:/data/durable/RAW/MoBaGenomics/_key/2019_10_01_MoBaGenetics_mother_2374.dta"))
names(moms)<-tolower(names(moms))
moms<-rename.vars(moms,from=c("sentrixid_m","batch_m"),to=c("sentrixid_mom","batch_mom"))
moms<-moms[,c("sentrixid_mom","batch_mom")]
dat<-merge2(dat,moms,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)
dads<-as.data.frame(read_dta("N:/data/durable/RAW/MoBaGenomics/_key/2019_10_01_MoBaGenetics_father_2374.dta"))
names(dads)<-tolower(names(dads))
dads<-rename.vars(dads,from=c("sentrixid_f","batch_f"),to=c("sentrixid_dad","batch_dad"))
dads<-dads[,c("sentrixid_dad","batch_dad")]
dat<-merge2(dat,dads,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

length(which(!is.na(dat$batch_mom)))
length(which(!is.na(dat$batch_dad)))


### BMI AND HEIGHT MEASURED VALUES ###

# height_mom == aa87

dat$height_mom<-dat$aa87
dat$weight_mom<-dat$aa85
dat$bmi_mom<-dat$weight_mom/((dat$height_mom/100)^2)

dat$height_dad<-dat$aa88
dat$weight_dad<-dat$aa89
dat$bmi_dad<-dat$weight_dad/((dat$height_dad/100)^2)


### HIGHEST EDUCATIONAL LEVEL, COMPLETED OR ONGOING ###

# Mothers #

dat$aa1125<-with(dat,ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125>aa1124),aa1125,
                            ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125<aa1124),NA,aa1125)))
dat$aa1128x<-with(dat,ifelse(aa1128==1 | aa1129==1,1,0))
dat$aa1128x<-with(dat,ifelse(is.na(aa1128x),0,aa1128x))

dat$edu_mom<-with(dat,ifelse(aa1124==1,1,
                             ifelse(aa1124==2,2,
                                    ifelse(aa1124==3,3,
                                           ifelse(aa1124==4,4,
                                                  ifelse(aa1124==5,5,
                                                         ifelse(aa1124==6,6,NA)))))))
dat$edu_mom<-with(dat,ifelse(is.na(aa1125),edu_mom,
                             ifelse(aa1125==1,1,
                                    ifelse(aa1125==2,2,
                                           ifelse(aa1125==3,3,
                                                  ifelse(aa1125==4,4,
                                                         ifelse(aa1125==5,5,
                                                                ifelse(aa1125==6,6,NA))))))))
dat$edu_mom<-with(dat,ifelse(is.na(edu_mom),0,edu_mom))

dat$edu_mom<-with(dat,ifelse(edu_mom==0 & aa1128x==1,3,
                             ifelse((edu_mom>0 & edu_mom<3) & aa1128x==1,3,
                                    ifelse(edu_mom>=3 & aa1128x==1,edu_mom,
                                           ifelse(edu_mom==0 & aa1128x==0,NA,
                                                  ifelse((edu_mom>0 & edu_mom<3) & aa1128x==0,edu_mom,
                                                         ifelse(edu_mom>=3 & aa1128x==0,edu_mom,NA)))))))

dat$eduyears_mom<-with(dat,ifelse(edu_mom==1,10,
                                  ifelse(edu_mom==2,10,
                                         ifelse(edu_mom==3,13,
                                                ifelse(edu_mom==4,13,
                                                       ifelse(edu_mom==5,19,
                                                              ifelse(edu_mom==6,20,
                                                                     ifelse(is.na(edu_mom),NA,NA))))))))

# Fathers #

dat$aa1127<-with(dat,ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127>aa1126),aa1127,
                            ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127<aa1126),NA,aa1127)))
dat$aa1130x<-with(dat,ifelse(aa1130==1 | aa1131==1,1,0))
dat$aa1130x<-with(dat,ifelse(is.na(aa1130x),0,aa1130x))

dat$edu_dad<-with(dat,ifelse(aa1126==1,1,
                             ifelse(aa1126==2,2,
                                    ifelse(aa1126==3,3,
                                           ifelse(aa1126==4,4,
                                                  ifelse(aa1126==5,5,
                                                         ifelse(aa1126==6,6,NA)))))))
dat$edu_dad<-with(dat,ifelse(is.na(aa1127),edu_dad,
                             ifelse(aa1127==1,1,
                                    ifelse(aa1127==2,2,
                                           ifelse(aa1127==3,3,
                                                  ifelse(aa1127==4,4,
                                                         ifelse(aa1127==5,5,
                                                                ifelse(aa1127==6,6,NA))))))))
dat$edu_dad<-with(dat,ifelse(is.na(edu_dad),0,edu_dad))

dat$edu_dad<-with(dat,ifelse(edu_dad==0 & aa1130x==1,3,
                             ifelse((edu_dad>0 & edu_dad<3) & aa1130x==1,3,
                                    ifelse(edu_dad>=3 & aa1130x==1,edu_dad,
                                           ifelse(edu_dad==0 & aa1130x==0,NA,
                                                  ifelse((edu_dad>0 & edu_dad<3) & aa1130x==0,edu_dad,
                                                         ifelse(edu_dad>=3 & aa1130x==0,edu_dad,NA)))))))

dat$eduyears_dad<-with(dat,ifelse(edu_dad==1,10,
                                  ifelse(edu_dad==2,10,
                                         ifelse(edu_dad==3,13,
                                                ifelse(edu_dad==4,13,
                                                       ifelse(edu_dad==5,19,
                                                              ifelse(edu_dad==6,20,
                                                                     ifelse(is.na(edu_dad),NA,NA))))))))

### DEFINITION OF EVER SMOKERS (smkinit) ###

# Mothers #

dat$smkinit_mom<-with(dat,ifelse(smoking_mom==0,0,
                                 ifelse(smoking_mom>0,1,NA)))

dat$smkinit2_mom<-with(dat,ifelse(aa1355==2 | aa1356>1 | aa1357!=0 | aa1358!=0 | aa1359>1 | aa1360!=0 | aa1361!=0 | 
                                    !is.na(aa1362) | aa1363==2 | !is.na(aa1364) | !is.na(aa1365),1,0))
dat$smkinit2_mom<-with(dat,ifelse(is.na(smkinit2_mom),0,smkinit2_mom))
dat$smkinit2_mom<-with(dat,ifelse(is.na(aa1355) & is.na(aa1356) & is.na(aa1357) & is.na(aa1358) & is.na(aa1359) & is.na(aa1360) & is.na(aa1361) & 
                                    is.na(aa1362) & is.na(aa1363) & is.na(aa1364) & is.na(aa1365),NA,smkinit2_mom))
dat$smkinit2_mom<-with(dat,ifelse(is.na(smkinit2_mom) & smoking_mom==0,0,
                                  ifelse(is.na(smkinit2_mom) & smoking_mom>0,1,
                                         ifelse(is.na(smkinit2_mom) & is.na(smoking_mom),NA,smkinit2_mom))))

# Fathers #

dat$smkinit_dad<-with(dat,ifelse(smoking_dad==0,0,
                                 ifelse(smoking_dad>0,1,NA)))

dat$smkinit2_dad<-with(dat,ifelse(aa1353==2 | aa1354==2 | ff214==2 | ff215>1 | ff216!=0 | ff217!=0 | ff218>1 | ff219!=0 | ff220!=0,1,0))
dat$smkinit2_dad<-with(dat,ifelse(is.na(smkinit2_dad),0,smkinit2_dad))
dat$smkinit2_dad<-with(dat,ifelse(is.na(aa1353) & is.na(aa1354) & is.na(ff214) & is.na(ff215) & is.na(ff216) & is.na(ff217) & is.na(ff218) & is.na(ff219) & is.na(ff220),NA,smkinit2_dad))
dat$smkinit2_dad<-with(dat,ifelse(is.na(smkinit2_dad) & smoking_dad==0,0,
                                  ifelse(is.na(smkinit2_dad) & smoking_dad>0,1,
                                         ifelse(is.na(smkinit2_dad) & is.na(smoking_dad),NA,smkinit2_dad))))


### DEFINITION OF PARITY (number of previous deliveries) ###

dat$parity<-with(dat,ifelse(paritet_5==0,0,
                            ifelse(paritet_5==1,1,
                                   ifelse(paritet_5==2,2,
                                          ifelse(paritet_5==3,3,
                                                 ifelse(paritet_5==4,3,NA))))))

length(which(!is.na(dat$parity)))


### MERGING PRINCIPAL COMPONENTS ###

pcs<-as.data.frame(read.csv("N:/data/durable/Projects/Magnus_MR_BMI/dat/chr_tmp/combined_filtered_indep_indep.csv",
                            header=TRUE,sep=",",dec=".")[,3:13])
names(pcs)<-c("sentrixid_mom","pc1_mom","pc2_mom","pc3_mom","pc4_mom","pc5_mom",
              "pc6_mom","pc7_mom","pc8_mom","pc9_mom","pc10_mom")
pcs$sentrixid_mom<-as.character(pcs$sentrixid_mom)
dat<-merge2(dat,pcs,by.id=c("sentrixid_mom"),all.x=TRUE,sort=FALSE)

pcs<-rename.vars(pcs,
                 from=c("sentrixid_mom","pc1_mom","pc2_mom","pc3_mom","pc4_mom","pc5_mom",
                        "pc6_mom","pc7_mom","pc8_mom","pc9_mom","pc10_mom"),
                 to=c("sentrixid_dad","pc1_dad","pc2_dad","pc3_dad","pc4_dad","pc5_dad",
                      "pc6_dad","pc7_dad","pc8_dad","pc9_dad","pc10_dad"))
pcs$sentrixid_dad<-as.character(pcs$sentrixid_dad)
dat<-merge2(dat,pcs,by.id=c("sentrixid_dad"),all.x=TRUE,sort=FALSE)

length(which(!is.na(dat$pc1_mom)))
length(which(!is.na(dat$pc2_mom)))

length(which(!is.na(dat$pc1_dad)))
length(which(!is.na(dat$pc2_dad)))


### DEFINITION OF SUBFERTILITY ###

dat$exclude<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))
dat$art_nona<-with(dat,ifelse(!is.na(art) & art>0,1,0))

# Threshold for subfertility: >12 (n=9732), >=12 (n=11850)
dat$subfertility_12plus<-with(dat,ifelse(aa48<12,0,
                                         ifelse(aa48>=12,1,NA)))
dat$subfertility_12plus<-with(dat,ifelse(art_nona==0,subfertility_12plus,
                                         ifelse(art_nona==1,1,NA)))
dat$subfertility_12plus<-with(dat,ifelse(exclude==0 & is.na(subfertility_12plus),0,
                                         ifelse(exclude==0 & subfertility_12plus==0,0,
                                                ifelse(exclude==0 & subfertility_12plus==1,1,
                                                       ifelse(exclude==1 & is.na(subfertility_12plus),NA,
                                                              ifelse(exclude==1 & subfertility_12plus==0,NA,
                                                                     ifelse(exclude==1 & subfertility_12plus==1,NA,NA)))))))

dat$preg_plan<-with(dat,ifelse(is.na(aa46),9,
                               ifelse(aa46==0,9,
                                      ifelse(aa46==1,0,
                                             ifelse(aa46==2,1,9)))))
attributes(dat$preg_plan)$value.label<-c("0=No","1=Yes","9=NA")
attributes(dat$preg_plan)$label<-c("Planned pregnancy")


### FINAL FORMAT OF DATA ###

dat<-rename.vars(dat,
                 from=c("mors_alder","fars_alder"),
                 to=c("agedelivery_mom","agedelivery_dad"))

dat$exclude<-with(dat,ifelse(is.na(subfertility),1,
                             ifelse(flerfodsel==1,1,0)))
dat$exclude<-with(dat,ifelse(is.na(exclude),0,exclude))
dat$flerfodsel<-with(dat,ifelse(is.na(flerfodsel),0,flerfodsel))

dat$exclude_batch_mom<-with(dat,ifelse(batch_mom=="TED",1,
                                       ifelse(batch_mom=="NORMENT-JAN15",1,
                                              ifelse(batch_mom=="NORMENT-JUN15",1,0))))
dat$exclude_batch_mom<-with(dat,ifelse(is.na(exclude_batch_mom),0,exclude_batch_mom))

dat$exclude_batch_dad<-with(dat,ifelse(batch_dad=="TED",1,
                                       ifelse(batch_dad=="NORMENT-JAN15",1,
                                              ifelse(batch_dad=="NORMENT-JUN15",1,0))))
dat$exclude_batch_dad<-with(dat,ifelse(is.na(exclude_batch_dad),0,exclude_batch_dad))

dat$exclude_skjema_tbl1<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))

dat$exclude_subfertility_na<-with(dat,ifelse(is.na(subfertility),1,0))

dat$exclude_flerfodsel<-with(dat,ifelse(flerfodsel==1,1,0))


nam_ok<-c("sentrixid_mom","sentrixid_dad","m_id_2374","f_id_2374","preg_id_2374","flerfodsel","preg_plan","art_nona",
          "bmi_mom","bmi_dad",
          "agedelivery_mom","agedelivery_dad","smkinit_mom","smkinit_dad","smkinit2_mom","smkinit2_dad",
          "eduyears_mom","eduyears_dad",
          "parity","subfertility","subfertility_12plus","batch_mom","batch_dad",
          "exclude","exclude_batch_mom","exclude_batch_dad","exclude_skjema_tbl1","exclude_flerfodsel","exclude_subfertility_na",
          "bmi_grs_mom","bmi_grs_dad","eduyears_grs_mom","eduyears_grs_dad","smkinit_grs_mom","smkinit_grs_dad",
          "pc1_mom","pc2_mom","pc3_mom","pc4_mom","pc5_mom","pc6_mom","pc7_mom","pc8_mom","pc9_mom","pc10_mom",
          "pc1_dad","pc2_dad","pc3_dad","pc4_dad","pc5_dad","pc6_dad","pc7_dad","pc8_dad","pc9_dad","pc10_dad")
dat<-dat[,nam_ok]


######################################
### DEFINITION OF STUDY FLOW CHART ###
######################################

# Total registers
length(which(!is.na(dat$preg_id_2374)))

# Singleton pregnancies (flerfodsel==0)
dat<-subset2(dat,"dat$exclude_flerfodsel==0")
dim(dat)[1]

# Exclusion by version_skjema_tbl1
table(dat$exclude_skjema_tbl1)
table(dat$exclude_skjema_tbl1,by=dat$subfertility_12plus)
dat<-subset2(dat,"dat$exclude_skjema_tbl1==0")
dim(dat)[1]

# Exclusion of mothers/fathers without any anthropometrical measurement
dat$anthrop_na_mom<-with(dat,ifelse(is.na(bmi_mom),1,0))
table(dat$anthrop_na_mom)

dat$anthrop_na_dad<-with(dat,ifelse(is.na(bmi_dad),1,0))
table(dat$anthrop_na_dad)


# Outlier exclusion

plot(compareGroups(~bmi_mom,data=dat))
sort(dat$bmi_mom,decreasing=FALSE)
sort(dat$bmi_mom,decreasing=TRUE)
dat$bmi_mom<-with(dat,ifelse(bmi_mom<9 | bmi_mom>90,NA,bmi_mom))

jpeg(filename="./Outputs/descriptive/bmi_mom.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_mom,data=dat))
dev.off()


plot(compareGroups(~bmi_dad,data=dat))
sort(dat$bmi_dad,decreasing=FALSE)
sort(dat$bmi_dad,decreasing=TRUE)
dat$bmi_dad<-with(dat,ifelse(bmi_dad<9 | bmi_dad>90,NA,bmi_dad))

jpeg(filename="./Outputs/descriptive/bmi_dad.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_dad,data=dat))
dev.off()


# Mothers/fathers without any anthropometrical measurement after outlier exclusion
dat$anthrop_na_mom<-with(dat,ifelse(is.na(bmi_mom),1,0))
table(dat$anthrop_na_mom)
dat$anthrop_na_dad<-with(dat,ifelse(is.na(bmi_dad),1,0))
table(dat$anthrop_na_dad)

# Unique IDs of mothers/fathers with any anthropometrical measurement
length(unique(dat[which(dat$anthrop_na_mom==0),c("m_id_2374")]))
length(unique(dat[which(dat$anthrop_na_dad==0),c("f_id_2374")]))

# IDs of mothers/fathers with genotype data
length(dat[which((!is.na(dat$height_grs_mom) | !is.na(dat$bmi_grs_mom)) & dat$anthrop_na_mom==0 & !is.na(dat$sentrixid_mom)),c("sentrixid_mom")])
length(dat[which((!is.na(dat$height_grs_dad) | !is.na(dat$bmi_grs_dad)) & dat$anthrop_na_dad==0 & !is.na(dat$sentrixid_dad)),c("sentrixid_dad")])

# Exclusion by batch_m / batch_d or lack of 
table(dat$exclude_batch_mom)
table(dat$exclude_batch_dad)
dat$exclude_genotype_mom<-with(dat,ifelse(exclude_batch_mom==1 | is.na(dat$sentrixid_mom),1,0))
dat$exclude_genotype_dad<-with(dat,ifelse(exclude_batch_dad==1 | is.na(dat$sentrixid_dad),1,0))

# Genotyped mothers and fathers with anthropometric measurements & quality OK

length(dat[which((!is.na(dat$height_grs_mom) | !is.na(dat$bmi_grs_mom)) & dat$anthrop_na_mom==0 & dat$exclude_genotype_mom==0),c("sentrixid_mom")])
length(dat[which((!is.na(dat$height_grs_dad) | !is.na(dat$bmi_grs_dad)) & dat$anthrop_na_dad==0 & dat$exclude_genotype_dad==0),c("sentrixid_dad")])

length(unique(dat[which((!is.na(dat$height_grs_mom) | !is.na(dat$bmi_grs_mom)) & dat$anthrop_na_mom==0 & dat$exclude_genotype_mom==0),c("sentrixid_mom")]))
length(unique(dat[which((!is.na(dat$height_grs_dad) | !is.na(dat$bmi_grs_dad)) & dat$anthrop_na_dad==0 & dat$exclude_genotype_dad==0),c("sentrixid_dad")]))

table(dat$exclude_subfertility_na)

# Mothers/fathers without educational data
dat$edu_na_mom<-with(dat,ifelse(is.na(edu_mom),1,0))
table(dat$edu_na_mom)
dat$edu_na_dad<-with(dat,ifelse(is.na(edu_dad),1,0))
table(dat$edu_na_dad)

# Mothers/fathers without smkinit data
dat$smkinit_na_mom<-with(dat,ifelse(is.na(smkinit_mom),1,0))
table(dat$smkinit_na_mom)
dat$smkinit_na_dad<-with(dat,ifelse(is.na(smkinit_dad),1,0))
table(dat$smkinit_na_dad)


# Calculation of genetically-determined and residual BMI #

vars01<-c("bmi_mom","bmi_dad")
vars02<-c("bmi_grs_mom","bmi_grs_dad")
vars03<-c("bmi_mom_z","bmi_dad_z")
vars04<-c("bmi_grs_mom_z","bmi_grs_dad_z")
vars05<-c("m_id_2374","f_id_2374")
vars06<-c("gen_pred_bmi_mom","gen_pred_bmi_dad")
vars07<-c("environ_bmi_mom","environ_bmi_dad")
vars08<-c("anthrop_na_mom","anthrop_na_dad")
vars09<-c("exclude_genotype_mom","exclude_genotype_dad")

for(i in 1:length(vars01))
  
{
  dat$exclude_var<-with(dat,ifelse(is.na(dat[,vars02[i]]) | dat[,vars08[i]]==1 | dat[,vars09[i]]==1,1,0))
  dat[,vars02[i]]<-with(dat,ifelse(exclude_var==0,dat[,vars02[i]],
                                   ifelse(exclude_var==1,NA,dat[,vars02[i]])))
  dat[,vars03[i]]<-as.numeric(with(dat,scale(dat[,vars01[i]])))
  dat[,vars04[i]]<-as.numeric(with(dat,scale(dat[,vars02[i]])))
  datx<-subset2(dat,"dat$exclude_var==0")
  mod01<-lm_robust(datx[,vars01[i]]~datx[,vars04[i]], data=datx, clusters=datx[,vars05[i]], se_type="stata")
  interc<-as.numeric(summary(mod01)$coefficients[1,1])
  slope<-as.numeric(summary(mod01)$coefficients[2,1])
  dat[,vars06[i]]<-(dat[,vars04[i]]*slope)+interc
  dat[,vars06[i]]<-with(dat,ifelse(exclude_var==0,dat[,vars06[i]],
                                   ifelse(exclude_var==1,NA,dat[,vars06[i]])))
  dat[,vars07[i]]<-dat[,vars01[i]]-dat[,vars06[i]]
  dat[,vars07[i]]<-with(dat,ifelse(exclude_var==0,dat[,vars07[i]],
                                   ifelse(exclude_var==1,NA,dat[,vars07[i]])))
  dat$exclude_var<-NULL
}

datx<-subset2(dat,"!is.na(dat$bmi_grs_mom) & dat$anthrop_na_mom==0 & dat$exclude_genotype_mom==0")
jpeg(filename="./Outputs/descriptive/bmi_grs_mom.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_grs_mom,data=datx))
dev.off()

jpeg(filename="./Outputs/descriptive/gen_pred_bmi_mom.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~gen_pred_bmi_mom,data=datx))
dev.off()

jpeg(filename="./Outputs/descriptive/environ_bmi_mom.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~environ_bmi_mom,data=datx))
dev.off()


datx<-subset2(dat,"!is.na(dat$bmi_grs_dad) & dat$anthrop_na_dad==0 & dat$exclude_genotype_dad==0")
jpeg(filename="./Outputs/descriptive/bmi_grs_dad.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~bmi_grs_dad,data=datx))
dev.off()

jpeg(filename="./Outputs/descriptive/gen_pred_bmi_dad.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~gen_pred_bmi_dad,data=datx))
dev.off()

jpeg(filename="./Outputs/descriptive/environ_bmi_dad.jpg",
     width=9000,height=9000,res=1200,pointsize=11.5)
par(las=1,cex=1,mar=c(5,5,2,2))
plot(compareGroups(~environ_bmi_dad,data=datx))
dev.off()


vars06<-c("gen_pred_bmi_mom","gen_pred_bmi_dad")
vars07<-c("environ_bmi_mom","environ_bmi_dad")
vars08<-c("gen_pred_bmi_mom_z","gen_pred_bmi_dad_z")
vars09<-c("environ_bmi_mom_z","environ_bmi_dad_z")

for(i in 1:length(vars06))
  
{
  dat[,vars08[i]]<-as.numeric(with(dat,scale(dat[,vars06[i]])))
  dat[,vars09[i]]<-as.numeric(with(dat,scale(dat[,vars07[i]])))
}

attributes(dat$agedelivery_mom)$label<-c("agedelivery_mom")
attributes(dat$agedelivery_dad)$label<-c("agedelivery_dad")
attributes(dat$bmi_mom)$label<-c("bmi_mom")
attributes(dat$bmi_dad)$label<-c("bmi_dad")
attributes(dat$height_mom)$label<-c("height_mom")
attributes(dat$height_dad)$label<-c("height_dad")


save(dat,file="N:/data/durable/Projects/Magnus_MR_BMI/R/MoBa_bmi_height.RData")


