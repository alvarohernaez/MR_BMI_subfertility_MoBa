rm(list=ls())

# Load required libraries #

setwd("N:/data/durable/Projects/Magnus_MR_BMI/R")


#############################
### GENERATION OF BMI GRS ###
#############################

# MoBa SELECTED SNPs ON BMI #

chr01<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr1.raw",header=TRUE,sep="\t"))
chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr4.raw",header=TRUE,sep="\t"))
chr05<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr5.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr9.raw",header=TRUE,sep="\t"))
chr10<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr10.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr11.raw",header=TRUE,sep="\t"))
chr12<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr12.raw",header=TRUE,sep="\t"))
chr13<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr13.raw",header=TRUE,sep="\t"))
chr14<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr14.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr15.raw",header=TRUE,sep="\t"))
chr16<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr16.raw",header=TRUE,sep="\t"))
chr17<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr17.raw",header=TRUE,sep="\t"))
chr18<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr18.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr20.raw",header=TRUE,sep="\t"))
chr21<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr21.raw",header=TRUE,sep="\t"))
chr22<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/chr22.raw",header=TRUE,sep="\t"))

chr01<-chr01[,c(1,7:dim(chr01)[2])]
chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr05<-chr05[,c(1,7:dim(chr05)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr10<-chr10[,c(1,7:dim(chr10)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr12<-chr12[,c(1,7:dim(chr12)[2])]
chr13<-chr13[,c(1,7:dim(chr13)[2])]
chr14<-chr14[,c(1,7:dim(chr14)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr16<-chr16[,c(1,7:dim(chr16)[2])]
chr17<-chr17[,c(1,7:dim(chr17)[2])]
chr18<-chr18[,c(1,7:dim(chr18)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr21<-chr21[,c(1,7:dim(chr21)[2])]
chr22<-chr22[,c(1,7:dim(chr22)[2])]

dat<-merge2(chr01,chr02,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr05,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr10,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr12,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr13,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr14,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr16,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr17,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr18,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr21,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr22,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO BMI, META-ANALYSIS #

bmi_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_bmi/gwas_bmi_cojo_locke2018.csv",header=TRUE,sep=",",dec=".")

names(bmi_ma)<-tolower(names(bmi_ma))
bmi_ma<-rename.vars(bmi_ma,
                    from=c("snp","tested_allele","freq_tested_allele_in_hrs","p"),
                    to=c("rsid","tested_allele_ma","maf_ma","pval"))
bmi_ma<-bmi_ma[,c("rsid","tested_allele_ma","maf_ma","beta","se","pval")]
bmi_ma$tested_allele_ma<-tolower(bmi_ma$tested_allele_ma)


### CHECK THAT THE 896 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(bmi_ma$rsid,decreasing=FALSE))

dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))

length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND bmi_ma MATCH


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
bmi_ma<-merge2(bmi_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
bmi_ma<-na.omit(bmi_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

bmi_ma$same_coding<-with(bmi_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                       ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

bmi_ma$flip<-with(bmi_ma,ifelse(beta>=0 & same_coding==1,0,
                                ifelse(beta<0 & same_coding==0,0,
                                       ifelse(beta>=0 & same_coding==0,1,
                                              ifelse(beta<0 & same_coding==1,1,NA)))))

bmi_ma$maf<-with(bmi_ma,ifelse(same_coding==1,maf_ma,
                               ifelse(same_coding==0,1-maf_ma,NA)))
bmi_ma$coef<-with(bmi_ma,ifelse(beta>0,beta,
                                ifelse(beta<0,-beta,NA)))


### CALCULATION OF BMI GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-bmi_ma[bmi_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-bmi_ma[bmi_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$bmi_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","bmi_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_BMI/R/bmi_grs.RData")


##################################
### GENERATION OF EDUYEARS GRS ###
##################################

# MoBa SELECTED SNPs ON EDUYEARS #

chr01<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr1.raw",header=TRUE,sep="\t"))
chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr4.raw",header=TRUE,sep="\t"))
chr05<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr5.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr9.raw",header=TRUE,sep="\t"))
chr10<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr10.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr11.raw",header=TRUE,sep="\t"))
chr12<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr12.raw",header=TRUE,sep="\t"))
chr13<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr13.raw",header=TRUE,sep="\t"))
chr14<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr14.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr15.raw",header=TRUE,sep="\t"))
chr16<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr16.raw",header=TRUE,sep="\t"))
chr17<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr17.raw",header=TRUE,sep="\t"))
chr18<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr18.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr20.raw",header=TRUE,sep="\t"))
chr21<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr21.raw",header=TRUE,sep="\t"))
chr22<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/chr22.raw",header=TRUE,sep="\t"))

chr01<-chr01[,c(1,7:dim(chr01)[2])]
chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr05<-chr05[,c(1,7:dim(chr05)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr10<-chr10[,c(1,7:dim(chr10)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr12<-chr12[,c(1,7:dim(chr12)[2])]
chr13<-chr13[,c(1,7:dim(chr13)[2])]
chr14<-chr14[,c(1,7:dim(chr14)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr16<-chr16[,c(1,7:dim(chr16)[2])]
chr17<-chr17[,c(1,7:dim(chr17)[2])]
chr18<-chr18[,c(1,7:dim(chr18)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr21<-chr21[,c(1,7:dim(chr21)[2])]
chr22<-chr22[,c(1,7:dim(chr22)[2])]

dat<-merge2(chr01,chr02,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr05,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr10,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr12,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr13,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr14,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr16,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr17,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr18,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr21,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr22,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO EDUYEARS, META-ANALYSIS LEE JJ ET AL, NAT GENET, 2018 #

eduyears_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_eduyears/gwas_eduyears_lee2018.csv",header=TRUE,sep=",",dec=".")

names(eduyears_ma)<-tolower(names(eduyears_ma))
eduyears_ma<-rename.vars(eduyears_ma,
                         from=c("snp","allele.1","frequency.allele.1","effect.size","p.value"),
                         to=c("rsid","tested_allele_ma","maf_ma","beta","pval"))
eduyears_ma$tested_allele_ma<-tolower(eduyears_ma$tested_allele_ma)


### CHECK THAT THE 1159 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(eduyears_ma$rsid,decreasing=FALSE))

dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))

length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND eduyears_ma MATCH


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
eduyears_ma<-merge2(eduyears_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
eduyears_ma<-na.omit(eduyears_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

eduyears_ma$same_coding<-with(eduyears_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                                 ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

eduyears_ma$flip<-with(eduyears_ma,ifelse(beta>=0 & same_coding==1,0,
                                          ifelse(beta<0 & same_coding==0,0,
                                                 ifelse(beta>=0 & same_coding==0,1,
                                                        ifelse(beta<0 & same_coding==1,1,NA)))))

eduyears_ma$maf<-with(eduyears_ma,ifelse(same_coding==1,maf_ma,
                                         ifelse(same_coding==0,1-maf_ma,NA)))
eduyears_ma$coef<-with(eduyears_ma,ifelse(beta>0,beta,
                                          ifelse(beta<0,-beta,NA)))


### CALCULATION OF EDUYEARS GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-eduyears_ma[eduyears_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-eduyears_ma[eduyears_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$eduyears_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","eduyears_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_BMI/R/eduyears_grs.RData")


#################################
### GENERATION OF SMKINIT GRS ###
#################################

# MoBa SELECTED SNPs ON SmkInit #

chr01<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr1.raw",header=TRUE,sep="\t"))
chr02<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr2.raw",header=TRUE,sep="\t"))
chr03<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr3.raw",header=TRUE,sep="\t"))
chr04<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr4.raw",header=TRUE,sep="\t"))
chr05<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr5.raw",header=TRUE,sep="\t"))
chr06<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr6.raw",header=TRUE,sep="\t"))
chr07<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr7.raw",header=TRUE,sep="\t"))
chr08<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr8.raw",header=TRUE,sep="\t"))
chr09<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr9.raw",header=TRUE,sep="\t"))
chr10<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr10.raw",header=TRUE,sep="\t"))
chr11<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr11.raw",header=TRUE,sep="\t"))
chr12<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr12.raw",header=TRUE,sep="\t"))
chr13<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr13.raw",header=TRUE,sep="\t"))
chr14<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr14.raw",header=TRUE,sep="\t"))
chr15<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr15.raw",header=TRUE,sep="\t"))
chr16<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr16.raw",header=TRUE,sep="\t"))
chr17<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr17.raw",header=TRUE,sep="\t"))
chr18<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr18.raw",header=TRUE,sep="\t"))
chr19<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr19.raw",header=TRUE,sep="\t"))
chr20<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr20.raw",header=TRUE,sep="\t"))
chr21<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr21.raw",header=TRUE,sep="\t"))
chr22<-as.data.frame(read.delim("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/chr22.raw",header=TRUE,sep="\t"))

chr01<-chr01[,c(1,7:dim(chr01)[2])]
chr02<-chr02[,c(1,7:dim(chr02)[2])]
chr03<-chr03[,c(1,7:dim(chr03)[2])]
chr04<-chr04[,c(1,7:dim(chr04)[2])]
chr05<-chr05[,c(1,7:dim(chr05)[2])]
chr06<-chr06[,c(1,7:dim(chr06)[2])]
chr07<-chr07[,c(1,7:dim(chr07)[2])]
chr08<-chr08[,c(1,7:dim(chr08)[2])]
chr09<-chr09[,c(1,7:dim(chr09)[2])]
chr10<-chr10[,c(1,7:dim(chr10)[2])]
chr11<-chr11[,c(1,7:dim(chr11)[2])]
chr12<-chr12[,c(1,7:dim(chr12)[2])]
chr13<-chr13[,c(1,7:dim(chr13)[2])]
chr14<-chr14[,c(1,7:dim(chr14)[2])]
chr15<-chr15[,c(1,7:dim(chr15)[2])]
chr16<-chr16[,c(1,7:dim(chr16)[2])]
chr17<-chr17[,c(1,7:dim(chr17)[2])]
chr18<-chr18[,c(1,7:dim(chr18)[2])]
chr19<-chr19[,c(1,7:dim(chr19)[2])]
chr20<-chr20[,c(1,7:dim(chr20)[2])]
chr21<-chr21[,c(1,7:dim(chr21)[2])]
chr22<-chr22[,c(1,7:dim(chr22)[2])]

dat<-merge2(chr01,chr02,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr03,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr04,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr05,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr06,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr07,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr08,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr09,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr10,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr11,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr12,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr13,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr14,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr15,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr16,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr17,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr18,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr19,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr20,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr21,by.id=c("FID"),all.x=TRUE,sort=FALSE)
dat<-merge2(dat,chr22,by.id=c("FID"),all.x=TRUE,sort=FALSE)


# SNPs RELATED TO SmkInit, META-ANALYSIS LIU ET AL, NAT GENET, 2018 #

smkinit_ma<-read.csv2("N:/data/durable/Projects/Magnus_MR_BMI/R/gwas_smoking/gwas_smkinit_liu2019.csv",header=TRUE,sep=",",dec=".")

names(smkinit_ma)<-tolower(names(smkinit_ma))
smkinit_ma<-rename.vars(smkinit_ma,
                        from=c("tested.allele","tested.allele.freq","pvalue"),
                        to=c("tested_allele_ma","maf_ma","beta","pval"))
smkinit_ma$tested_allele_ma<-tolower(smkinit_ma$tested_allele_ma)


### CHECK THAT THE 355 SNPs IN DAT BELONG TO THE META-ANALYSIS RESULTS ###

bbb<-as.character(sort(smkinit_ma$rsid,decreasing=FALSE))

dat_inv<-header.true(data.frame(t(dat)))
ccc<-as.character(sort(sub("\\_.*", "", rownames(dat_inv)),decreasing=FALSE))

length(which(ccc%in%bbb)) ### VARIABLES IN dat_inv AND smkinit_ma MATCH


# EXTRACT THE INFORMATION ON THE ALLELE TESTED IN MOBA #

rsid<-as.character(sub("\\_.*", "", rownames(dat_inv)))
tested_allele_moba<-sub('.*_', '', rownames(dat_inv))
moba_info<-as.data.frame(cbind(rsid,tested_allele_moba))
moba_info$tested_allele_moba<-tolower(moba_info$tested_allele_moba)
smkinit_ma<-merge2(smkinit_ma,moba_info,by.id=c("rsid"),all.x=TRUE,sort=FALSE)
smkinit_ma<-na.omit(smkinit_ma)


# SNP FLIPPING IF THE ALLELE ORDER IN MoBa AND THE META-ANALYSIS DIFFER #

smkinit_ma$same_coding<-with(smkinit_ma,ifelse(tested_allele_moba==tested_allele_ma,1,
                                               ifelse(tested_allele_moba!=tested_allele_ma,0,NA)))

smkinit_ma$flip<-with(smkinit_ma,ifelse(beta>=0 & same_coding==1,0,
                                        ifelse(beta<0 & same_coding==0,0,
                                               ifelse(beta>=0 & same_coding==0,1,
                                                      ifelse(beta<0 & same_coding==1,1,NA)))))

smkinit_ma$maf<-with(smkinit_ma,ifelse(same_coding==1,maf_ma,
                                       ifelse(same_coding==0,1-maf_ma,NA)))
smkinit_ma$coef<-with(smkinit_ma,ifelse(beta>0,beta,
                                        ifelse(beta<0,-beta,NA)))


### CALCULATION OF smkinit GRS ###

n_snps<-dim(dat)[2]-1
colnames(dat)<-sub("\\_.*", "", colnames(dat))
dat<-rename.vars(dat,from=c("FID"),to=c("id"))

vars01<-colnames(dat)[2:dim(dat)[2]]
vars02<-paste(vars01,"_flip",sep="")
vars03<-paste(vars01,"_coef",sep="")
vars04<-paste(vars01,"_weight",sep="")

for(i in 1:length(vars01))
  
{
  dat[,vars02[i]]<-smkinit_ma[smkinit_ma$rsid==vars01[i],"flip"]
  dat[,vars03[i]]<-smkinit_ma[smkinit_ma$rsid==vars01[i],"coef"]
  dat[,vars04[i]]<-with(dat,ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==0,0*dat[,vars03[i]],
                                   ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                          ifelse(dat[,vars02[i]]==0 & dat[,vars01[i]]==2,2*dat[,vars03[i]],
                                                 ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==0,2*dat[,vars03[i]],
                                                        ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==1,1*dat[,vars03[i]],
                                                               ifelse(dat[,vars02[i]]==1 & dat[,vars01[i]]==2,0*dat[,vars03[i]],NA)))))))
}

dat_coef<-dat[,vars03]
dat$coef_sum<-(dat_coef %>%
                 replace(is.na(.), 0) %>%
                 mutate(coef_sum = rowSums(.)))$coef_sum[1]

vars05<-append("id",vars04,after=1)
dat_weight<-dat[,vars05]
dat_weight$weight_sum<-(dat_weight[,2:dim(dat_weight)[2]] %>%
                          replace(is.na(.), 0) %>%
                          mutate(weight_sum = rowSums(.)))$weight_sum
dat_weight<-dat_weight[,c("id","weight_sum")]
dat<-merge(dat,dat_weight,by="id")
dat$smkinit_grs<-(dat$weight_sum/dat$coef_sum)*n_snps
dat<-dat[,c("id","smkinit_grs")]

save(dat,file="N:/data/durable/Projects/Magnus_MR_BMI/R/smkinit_grs.RData")

