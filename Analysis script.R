library(tidyverse)
library(performance)
library(lme4)
library(mgcv)
library(flextable)
library(psych)
library(ggpubr)
library(psych)
library(broom)


#### Prepare data ####

## read data

bl <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.BL.csv")


fu <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.FUP2.csv")


# remove cancer cases, breast & gastric cancer, duplicates

bl2 <- bl %>% 
  filter(!study %in% c("BC","GAC")) %>% 
  filter(!case == 1) %>% 
  filter(!duplicated(h2_id))


# merge baseline and FU

full <- inner_join(bl2, fu, by = "h2_id", suffix = c("_bl","_fu"))



# load corrected age data

load("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.Age.rda")


## Add corrected age variable into merged file 

LP.Age <- LP.Age %>% 
  filter(timepoint_cde == "Baseline") %>% 
  select(h2_id, atnd_age,age_fup2) %>% 
  rename(age_bl_correct = atnd_age, 
         age_fu_correct = age_fup2)

full <- left_join(full, LP.Age, by = "h2_id")


## time to follow-up

full <- full %>% mutate(time_fu = age_fu_correct - age_bl_correct)

hist(full$time_fu, breaks = 50) # minimum of 9 years 


### Check distributions of biomarkers to ensure no errors


# vector of all biomarker names

marker_names <- as.character(
  expression(neopt_d_bl,crp_g_bl,cnct_g_bl,s100at_g_bl,saat_g_bl,il6_msd_bl,
             il8_msd_bl,il10_msd_bl,ifng_msd_bl,tnfa_msd_bl,trp_d_bl,kyn_d_bl,
             hk_d_bl,ka_d_bl,xa_d_bl,aa_d_bl,haa_d_bl,pic_d_bl,qa_d_bl,KTr_bl,
             PAr_bl,HKXAr_bl,neopt_d_fu,crp_g_fu,cnct_g_fu,s100at_g_fu,
             saat_g_fu,il6_msd_fu,il8_msd_fu,il10_msd_fu,ifng_msd_fu,
             tnfa_msd_fu,trp_d_fu,kyn_d_fu,hk_d_fu,ka_d_fu,xa_d_fu,aa_d_fu,
             haa_d_fu,pic_d_fu,qa_d_fu,KTr_fu,PAr_fu,HKXAr_fu, cysta_d_bl, cysta_d_fu))


# histograms

#multi.hist(full[,marker_names])


## Replacing zero values 

# qqplot

#ggplot(full, aes(sample = log(crp_g_bl + 0.03))) +
  stat_qq(distribution = qnorm) + stat_qq_line(distribution = qnorm)


# replace zero's with 0.5 * limit of detection 

full$crp_g_bl2 = ifelse(full$crp_g_bl < 0.001, 0.06 / 2, full$crp_g_bl)

full$crp_g_fu2 = ifelse(full$crp_g_fu < 0.001, 0.06 / 2, full$crp_g_fu)

full$saat_g_bl2 <- ifelse(full$saat_g_bl < 0.001, 0.04 / 2, full$saat_g_bl)

full$il10_msd_bl2 <- ifelse(full$il10_msd_bl < 0.001, 0.045, full$il10_msd_bl)

full$il10_msd_fu2 <- ifelse(full$il10_msd_fu < 0.001, 0.045, full$il10_msd_fu)

full$ifng_msd_bl2 <- ifelse(full$ifng_msd_bl < 0.001, 0.396/2, full$ifng_msd_bl)

full$il6_msd_fu2 <- ifelse(full$il6_msd_fu < 0.001, 0.09/2, full$il6_msd_fu)


# update marker names

marker_names <- as.character(
  expression(neopt_d_bl,crp_g_bl2,cnct_g_bl,s100at_g_bl,saat_g_bl2,il6_msd_bl,
             il8_msd_bl,il10_msd_bl2,ifng_msd_bl2,tnfa_msd_bl,trp_d_bl,kyn_d_bl,
             hk_d_bl,ka_d_bl,xa_d_bl,aa_d_bl,haa_d_bl,pic_d_bl,qa_d_bl,KTr_bl,
             PAr_bl,HKXAr_bl,neopt_d_fu,crp_g_fu2,cnct_g_fu,s100at_g_fu,
             saat_g_fu,il6_msd_fu2,il8_msd_fu,il10_msd_fu2,ifng_msd_fu,
             tnfa_msd_fu,trp_d_fu,kyn_d_fu,hk_d_fu,ka_d_fu,xa_d_fu,aa_d_fu,
             haa_d_fu,pic_d_fu,qa_d_fu,KTr_fu,PAr_fu,HKXAr_fu, cysta_d_bl, cysta_d_fu))


### Log transform and winsorise biomarker variables 

# create winsorise function

winsorise <- function(x){
  
  maxin <- max(x[scale(x) < 3],na.rm=T)
  
  minin <- min(x[scale(x) > -3],na.rm=T)
  
  
  x <- ifelse(x > maxin & !is.na(x), maxin, x)
  
  x <- ifelse(x < minin & !is.na(x), minin, x)
  
  x
  
}


# winsorise biomarkers and AA

full <- full %>% 
  mutate(across(marker_names, ~ log2(.x), .names = "logW_{.col}")) %>% 
  mutate(across(starts_with("logW"), winsorise)) %>% 
  mutate(across(contains("ageacc"), winsorise, .names = "W_{.col}"))


#### Table 2 ####

## Baseline and follow-up means ##

# median and IQR at baseline

bl_res <- full %>% 
  select(marker_names) %>% 
  summarise(across(contains("_bl"), 
                   list(medianBL = ~ median(.,na.rm=T), IQRBL = ~ IQR(.,na.rm=T)))) %>% 
  pivot_longer(cols = everything(), names_sep = "_(?!.*_)", 
               names_to = c("variable",".value")) %>% 
  mutate(variable = sub("_[^_]+$", "", .$variable))


# median and IQR at fu

fu_res <- full %>% 
  select(bioms) %>% 
  summarise(across(contains("_fu"), 
                   list(medianFU = ~ median(.,na.rm=T), IQRFU = ~ IQR(.,na.rm=T)))) %>% 
  pivot_longer(cols = everything(), names_sep = "_(?!.*_)", 
               names_to = c("variable",".value")) %>% 
  mutate(variable = sub("_[^_]+$", "", .$variable))


# merge bl and follow-up res

res <- full_join(bl_res, fu_res, by = "variable")


# write output to table

res$BL <- sprintf("%2.2f (%2.2f)", res$medianBL, res$IQRBL)

res$FU <- sprintf("%2.2f (%2.2f)", res$medianFU, res$IQRFU)

res <- res %>% select(variable,BL,FU)


table1 <- flextable(res)

print(table1, preview = "docx")



### Change over time ###

# CRP

full$change_CRP = full$logW_crp_g_fu2 - full$logW_crp_g_bl2

m_crp <- lm(change_CRP ~ time_fu, data = full)

tidy(m_crp, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)


# Neopterin 

full$change_neopt = full$logW_neopt_d_fu - full$logW_neopt_d_bl

m_neo <- lm(change_neopt ~ time_fu, data = full)

tidy(m_neo, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

# cystatin

full$change_cysta = full$logW_cysta_d_fu - full$logW_cysta_d_bl

m_cys <- lm(change_cysta ~ time_fu, data = full)

tidy(m_cys, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)
  

# serum amyloid A

full$change_saat = full$logW_saat_g_fu - full$logW_saat_g_bl2

m_saat <- lm(change_saat ~ time_fu, data = full)

tidy(m_saat, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)


# interleukin 6

full$change_il6 = full$logW_il6_msd_fu2 - full$logW_il6_msd_bl

m_il6 <- lm(change_il6 ~ time_fu, data = full)

tidy(m_il6, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)


#### Table 3 ####

# confounders

full$sex_cde_bl <- as.factor(full$sex_cde_bl)

levels(full$sex_cde_bl) <- c("Male","Female")

summary(as.factor(full$cigst_cde))

summary(as.factor(full$educlvl_ord))

### CRP models

## Pheno

# model 1

m1 <- lm(W_AgeAccelPheno_fu ~ scale(logW_crp_g_fu2) + cob_cde_bl + sex_cde_bl + 
           bmi_rrto + cigst_cde + educlvl_ord + age_fu_correct, data = full)

tidy(m1, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

# model 2 (sensitivity analysis)

m2 <- lm(W_AgeAccelPheno_fu ~ scale(logW_crp_g_fu2) + age_fu_correct + 
           cob_cde_bl + sex_cde_bl + bmi_rrto + seifa_10_fu + 
           as.factor(phys.act) + tot_alclw_mrt + cigst_cde + tot_alclw_mrt + 
           logpack + educlvl_ord, data = full)

tidy(m2, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

## Grim

# model 1

m1 <- lm(W_AgeAccelGrim_fu ~ scale(logW_crp_g_fu2) + age_fu_correct + cob_cde_bl + 
           sex_cde_bl + bmi_rrto + cigst_cde + educlvl_ord, data = full)

tidy(m1, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

# model 2 (sensitivity analysis)

m2 <- lm(W_AgeAccelGrim_fu ~ scale(logW_crp_g_fu2) + age_fu_correct + 
           cob_cde_bl + sex_cde_bl + bmi_rrto + seifa_10_fu + 
           as.factor(phys.act) + tot_alclw_mrt + cigst_cde + tot_alclw_mrt + 
           logpack + educlvl_ord, data = full)

tidy(m2, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)


## Dunedin

m1 <- lm(scale(DunedinPoAm) ~ scale(logW_crp_g_fu2) + cob_cde_bl + sex_cde_bl + 
           bmi_rrto + cigst_cde + educlvl_ord + age_fu_correct, data = full)

tidy(m1, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)



## Interleukin-6 

# model 1

m1_il6 <- lm(W_AgeAccelGrim_bl ~ scale(logW_il6_msd_bl) + cob_cde_bl + sex_cde_bl + 
           bmi_rrto + cigst_cde + educlvl_ord, data = full)

tidy(m1_il6, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

# model 2 (sensitivity analysis)

m2_il6 <- lm(W_AgeAccelGrim_bl ~ scale(logW_il6_msd_bl) + age_bl_correct + 
           cob_cde_bl + sex_cde_bl + base_bmi_rrto + seifa_10_bl + 
           as.factor(phys.act) + tot_alclw_mrt + cigst_cde + tot_alclw_mrt + 
           logpack + educlvl_ord, data = full)

tidy(m2_il6, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)









