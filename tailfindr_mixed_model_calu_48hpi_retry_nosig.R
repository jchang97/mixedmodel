library(lme4)
library(dplyr)
library(stringr)
library(ggplot2)
library(gtools)
library(purrr)
library(tidyr)
#install.packages("ggpubr")
library(ggpubr)
#install.packages("writexl")
library(writexl)
library(base)
library(lmerTest)
library(xlsx)
library(lattice)
setwd("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/mixed_model/calu_48hpi/retry/retry_2_nosig")

##FUNCTIONS


lmAnalysis<-function(x, lmer=T){
  x$replicate = factor(x$replicate)
  if(lmer){
    m = try(lmer(tail_length ~ condition + (1 | replicate), data=x))
  }else{
    m = try(lm(x$tail_length ~ x$condition))
  }

  if(!inherits(m,"try-error")){
    return(summary(m))
  }else{
    return(NULL)
  }
  
}



lmtestAnalysis<-function(x, lmer=T){
  x$replicate = factor(x$replicate)
  if(lmer){
    o = try(lmerTest::lmer(tail_length ~ condition + (1 | replicate), data=x))
  }else{
    o = try(lmerTest::lm(x$tail_length ~ x$condition))
  }
  
  if(!inherits(o,"try-error")){
    return(summary(o))
  }else{
    return(NULL)
  }
  
}

#import datasets
calu_c1 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_48hpi/c1/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_c1.tsv")
calu_c2 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_48hpi/c2/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_c2.tsv")
calu_c3 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_48hpi/c3/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_c3.tsv")
calu_i1 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_48hpi/i1/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_i1.tsv")
calu_i2 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_48hpi/i2/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_i2.tsv")
calu_i3 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_48hpi/i3/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_i3.tsv")


#remove duplicates
calu_c1 <- as_tibble(calu_c1)
calu_c2 <- as_tibble(calu_c2)
calu_c3 <- as_tibble(calu_c3)
calu_i1 <- as_tibble(calu_i1)
calu_i2 <- as_tibble(calu_i2)
calu_i3 <- as_tibble(calu_i3)


calu_c1_nd <- calu_c1 %>% distinct()
calu_c2_nd <- calu_c2 %>% distinct()
calu_c3_nd <- calu_c3 %>% distinct()
calu_i1_nd <- calu_i1 %>% distinct()
calu_i2_nd <- calu_i2 %>% distinct()
calu_i3_nd <- calu_i3 %>% distinct()




#add condition and replicate info
calu_c1_nd$replicate <- "c1"
calu_c2_nd$replicate <- "c2"
calu_c3_nd$replicate <- "c3"
calu_i1_nd$replicate <- "i1"
calu_i2_nd$replicate <- "i2"
calu_i3_nd$replicate <- "i3"


calu_c1_nd$condition <- "control"
calu_c2_nd$condition <- "control"
calu_c3_nd$condition <- "control"
calu_i1_nd$condition <- "infected"
calu_i2_nd$condition <- "infected"
calu_i3_nd$condition <- "infected"




#merge all datasets by row
calu_all <- rbind(calu_c1_nd, calu_c2_nd, calu_c3_nd, calu_i1_nd, calu_i2_nd, calu_i3_nd)


#subset data by tail_is_valid = TRUE
subset <- calu_all[calu_all$tail_is_valid == TRUE, ]


#subset data by tail_is_valid = TRUE and "ENSG"=TRUE
ens_subset <- subset %>% 
  filter(str_detect(ORFs, "ENS"))

#subset by polya or polyt

subset_polya <- ens_subset[ens_subset$read_type == "polyA",]
subset_polyt <- ens_subset[ens_subset$read_type == "polyT",]


#remove NA's
subset_polya_na <- subset_polya[!is.na(subset_polya$tail_length), ]
subset_polyt_na <- subset_polyt[!is.na(subset_polyt$tail_length), ]


####remove entries with just one value

group_subset_polya <- split(subset_polya_na, subset_polya_na$ORFs)
group_subset_polyt <- split(subset_polyt_na, subset_polyt_na$ORFs)


group_subset_polya_r <- group_subset_polya[lapply(group_subset_polya,nrow)>1]
group_subset_polyt_r <- group_subset_polyt[lapply(group_subset_polyt,nrow)>1]

names_polya <- names(group_subset_polya_r)
names_polyt <- names(group_subset_polyt_r)

subset_polya_final <- subset_polya_na[subset_polya_na$ORFs %in% names_polya,]
subset_polyt_final <- subset_polyt_na[subset_polyt_na$ORFs %in% names_polyt,]

#logtransform
subset_polyt_final$tail_length=log(subset_polyt_final$tail_length)

boxplot(subset_polyt_final$tail_length)
qqmath(subset_polyt_final$tail_length, id = 0.05)
plot(subset_polyt_final$tail_length, type = c("p", "smooth"))
qnorm(seq(0.01,0.99,0.01))
quantile(rnorm(200),probs = seq(0.01,0.99,0.01))


#subset by MT or non MT
subset_mt <- subset_polyt_final[grep("MT",subset_polyt_final$chrom),]
subset_nonmt <- subset_polyt_final[subset_polyt_final$chrom != "MT",]

#remove NA's
subset_na_mt <- subset_mt[!is.na(subset_mt$tail_length), ]
subset_na_nonmt <- subset_nonmt[!is.na(subset_nonmt$tail_length), ]

#get list of unique ENSG's
ensg_mt <- subset_na_mt$ORFs
ensg_nonmt <- subset_na_nonmt$ORFs


n_ensg_mt <- unique(unlist(strsplit(ensg_mt, " ")))
n_ensg_nonmt <- unique(unlist(strsplit(ensg_nonmt, " ")))


#split by ensg
split_mt = list()
for (i in n_ensg_mt) {
  split_mt[[i]] <- subset_na_mt[which(subset_na_mt$"ORFs"==i),] 
}


split_nonmt = list()
for (i in n_ensg_nonmt) {
  split_nonmt[[i]] <- subset_na_nonmt[which(subset_na_nonmt$"ORFs"==i),] 
}



####remove entries with <6 values

split_r_mt <- split_mt[lapply(split_mt,nrow)>5]
split_r_nonmt <- split_nonmt[lapply(split_nonmt,nrow)>5]



##Mixed model (lmerTest) p-value

mt2 = lapply(split_r_mt, lmtestAnalysis)
mt_coeff2=lapply(mt2,function(x) x$coefficients[2,5])
mt2_f <- unlist(mt_coeff2)


nonmt2 = lapply(split_r_nonmt, lmtestAnalysis)
nonmt_coeff2=lapply(nonmt2,function(x) x$coefficients[2,5])
nonmt2_f <- unlist(nonmt_coeff2)


##Mixed model (lmerTest) intercept (infected)

#mt2 = lapply(split_r_mt, lmtestAnalysis)
mt_coeff3=lapply(mt2,function(x) x$coefficients[2,1])
mt3_f <- unlist(mt_coeff3)


#nonmt2 = lapply(split_r_nonmt, lmtestAnalysis)
nonmt_coeff3=lapply(nonmt2,function(x) x$coefficients[2,1])
nonmt3_f <- unlist(nonmt_coeff3)



###Mixed model (lmerTest) intercept (control)

#mt_coeff4=lapply(mt2,function(x) x$coefficients[1,1])
#mt4_f <- unlist(mt_coeff4)



#nonmt_coeff4=lapply(nonmt2,function(x) x$coefficients[1,1])
#nonmt4_f <- unlist(nonmt_coeff4)


###Mixed model (lmerTest) SE (infected)

#mt_coeff5=lapply(mt2,function(x) x$coefficients[2,2])
#mt5_f <- unlist(mt_coeff5)



#nonmt_coeff5=lapply(nonmt2,function(x) x$coefficients[2,2])
#nonmt5_f <- unlist(nonmt_coeff5)

###Mixed model (lmerTest) SE (control)

#mt_coeff6=lapply(mt2,function(x) x$coefficients[1,2])
#mt6_f <- unlist(mt_coeff6)



#nonmt_coeff6=lapply(nonmt2,function(x) x$coefficients[1,2])
#nonmt6_f <- unlist(nonmt_coeff6)

##padj

mt2_p <- p.adjust(mt2_f, method="BH", n=length(mt2_f))
nonmt2_p <- p.adjust(nonmt2_f, method="BH", n=length(nonmt2_f))

#mt2_pf <- subset(mt2_p, mt2_p < 0.05)
#nonmt2_pf <- subset(nonmt2_p, nonmt2_p < 0.05)

#write.xlsx(mt2_pf, 'calu_48hpi_mt.xlsx')
#write.xlsx(nonmt2_pf, 'calu_48hpi_nonmt.xlsx')

write.xlsx(mt2_p, 'calu_48hpi_mt_log_nofilter.xlsx')
write.xlsx(nonmt2_p,'calu_48hpi_nonmt_log_nofilter.xlsx')

#Intercept (infected)

write.xlsx(mt3_f, 'calu_48hpi_mt_int.xlsx')
write.xlsx(nonmt3_f, 'calu_48hpi_nonmt_int.xlsx')


#Intercept (control)
#write.xlsx(mt4_f, 'calu_48hpi_mt_int_control.xlsx')
#write.xlsx(nonmt4_f, 'calu_48hpi_nonmt_int_control.xlsx')


#SE (infected)
#write.xlsx(mt5_f, 'calu_48hpi_mt_SE_infected.xlsx')
#write.xlsx(nonmt5_f, 'calu_48hpi_nonmt_SE_infected.xlsx')

#SE (control)
#write.xlsx(mt6_f, 'calu_48hpi_mt_SE_control.xlsx')
#write.xlsx(nonmt6_f, 'calu_48hpi_nonmt_SE_control.xlsx')

