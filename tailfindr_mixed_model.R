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

setwd("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi")

#import datasets
calu_c1 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi/c1/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_c1.tsv")
calu_c2 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi/c2/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_c2.tsv")
calu_c3 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi/c3/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_c3.tsv")
calu_i1 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi/i1/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_i1.tsv")
calu_i2 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi/i2/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_i2.tsv")
calu_i3 <- read.delim("Z:/User/Jessie/Projects/Corona_in_vitro_time_course/Analysis/host_analysis/tailfindr/calu_24hpi/i3/calu_dcDNA_readToCluster_cut_2.readPolyaMerged_i3.tsv")


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

#median by ORFs per replicate

group_subset_polya_m <- subset_polya_final %>% 
  group_by(ORFs,replicate,condition) %>% 
  summarise(median=median(tail_length))


group_subset_polyt_m <- subset_polyt_final %>% 
  group_by(ORFs,replicate,condition) %>% 
  summarise(median=median(tail_length))


#plot
boxplot(group_subset_polyt_m$median)

# ---------------------------------------------------

#mixed model (includes info about replicate effect)
gpa_mixed = lmer(median ~ condition + (1 | replicate), data = group_subset_polyt_m)
summary(gpa_mixed)
confint(gpa_mixed)

#convert t value to p value
df = 1
p <- 2*pt(q=3.033, df=df, lower.tail=FALSE)

#------------------------------------------------------
#logtransform
group_subset_polyt_m$log_median=log(group_subset_polyt_m$median)

#mixed model - logtransform (includes info about replicate effect)
gpa_mixed_log = lmer(log_median ~ condition + (1 | replicate), data = group_subset_polyt_m)
summary(gpa_mixed_log)
confint(gpa_mixed_log)


#plot
boxplot(group_subset_polyt_m$log_median)

#convert t value to p value
df = 1
p <- 2*pt(q=4.094, df=df, lower.tail=FALSE)

#---------------------------------------


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


##Mixed model
split_r_mt_un <- for (i in n_ensg_mt) {
  unlist(split_r_mt[[i]]$tail_length)
}
split_r_mt <- unlist(split_r_mt$)
mt = list()
mm <- for (i in n_ensg_mt) {
  mt[[i]] <- lmer(split_r_mt[[i]]$tail_length ~ split_r_mt[[i]]$condition + (1 | split_r_mt[[i]]$replicate))
}

split_r_mt$ENSG00000198938$tail_length

mm <- subset_na_mt %>% 
  group_by(ORFs) %>% 
  lmer(tail_length ~ condition + (1 | replicate), data = split_na_mt)

mm <- subset_na_mt %>%
  nest(-ORFs) %>%
  mutate(results = map(ORFs, ~ lmer(tail_length ~ condition + (1 | replicate), data = .)))


set.seed(1234)
nested_fit <- split_r_mt %>%
  mutate(fit   = map(lmer(tail_length ~ condition + (1 | replicate), data = split_r_mt))
         


splmt_m_mt <- split_r_mt %>%
  lmer(tail_length ~ [[i]]$condition + (1 | [[i]]$replicate), data = split_r_mt)) %>%
  summarise(total_ss =
split_m_mt <- sapply(n_ensg_mt, function(i) 
     mixed[[i]] <- lmer([[i]]$tail_length ~ [[i]]$condition + (1 | [[i]]$replicate), data = split_r_mt))
    



split_m_nonmt <- sapply(n_ensg_nonmt, function(i) 
  tryCatch(
    wilcox.test(split_c_r_nonmt[[i]]$tail_length, 
                split_i_r_nonmt[[i]]$tail_length, 
                na.action(na.omit)$p.value),
    warning = function(w) return(NA),
    error = function (e) return(NA)
    
  )
)


sink("lm_caco_24hpi.txt")
print(summary(lm(cars$speed ~ cars$dist)))
sink()  # returns output to the console






