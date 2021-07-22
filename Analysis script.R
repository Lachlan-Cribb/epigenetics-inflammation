library(tidyverse)
library(performance)
library(lme4)
library(mgcv)
library(flextable)
library(psych)
library(ggpubr)
library(psych)
library(broom)
library(car)


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


#### Table 1 demographics ####

psych::describe(full$age_fu_correct)

full %>% group_by(sex_cde_bl) %>% tally() %>% mutate(percent = n/sum(n) * 100)

psych::describe(full$bmi_rrto)

full %>% group_by(cob_cde_bl) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full %>% group_by(sample_type) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full %>% group_by(educlvl_ord) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full %>% group_by(cigst_cde) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full <- full %>% mutate(packyears = if_else(logpack > 0, exp(logpack), 0))


#### process biomarkers ####

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
  mutate(across(contains("ageacc"), winsorise, .names = "W_{.col}")) %>% 
  mutate(across(c("AA.Zhang","AA.DunedinPoAm"), winsorise, .names = "W_{.col}"))

## calculate inflammaging signature



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

### Change over time) ###

# CRP

full$change_CRP = full$logW_crp_g_fu2 - full$logW_crp_g_bl2

m_crp <- lm(change_CRP ~ time_fu, data = full)

tidy(m_crp, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

crPlots(m_crp)

# Neopterin 

full$change_neopt = full$logW_neopt_d_fu - full$logW_neopt_d_bl

m_neo <- lm(change_neopt ~ time_fu, data = full)

tidy(m_neo, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

crPlots(m_neo)

# cystatin

full$change_cysta = full$logW_cysta_d_fu - full$logW_cysta_d_bl

m_cys <- lm(change_cysta ~ time_fu, data = full)

tidy(m_cys, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

crPlots(m_cys)

# serum amyloid A

full$change_saat = full$logW_saat_g_fu - full$logW_saat_g_bl2

m_saat <- lm(change_saat ~ time_fu, data = full)

tidy(m_saat, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

crPlots(m_saat)

# interleukin 6

full$change_il6 = full$logW_il6_msd_fu2 - full$logW_il6_msd_bl

m_il6 <- lm(change_il6 ~ time_fu, data = full)

tidy(m_il6, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

crPlots(m_il6)

#### Figure 1 ####

library(knitr)

cor_dat <- full %>% 
  select(starts_with("logW")) %>% 
  select(ends_with("fu"), ends_with("fu2"))


cors_spearman <- function(df) { 
  M <- Hmisc::rcorr(as.matrix(df), type = "spearman")
  Mdf <- map(M, ~data.frame(.x))
}

plot_data <- cors_spearman(cor_dat) %>%
  map(~rownames_to_column(.x, var="measure1")) %>%
  map(~pivot_longer(.x, -measure1, "measure2")) %>%
  bind_rows(.id = "id") %>%
  pivot_wider(names_from = id, values_from = value) %>% 
  mutate(across(c("measure1","measure2"), str_replace, "logW_", "")) %>% 
  mutate(measure1 = sub("_.*", "", plot_data$measure1)) %>% 
  mutate(measure2 = sub("_.*", "", plot_data$measure2)) %>% 
  mutate(across(contains("measure"), str_to_upper))

plot_data %>%
  ggplot(aes(measure1, measure2, fill=r, label=round(r,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text(size =2.5) +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90))

ggsave("correlation heatmap.png", device = "png", dpi = 450)

#### figure 2 ####

aa_dat <- full %>% 
  select(starts_with("W_")) %>% 
  select(W_AA.Zhang, W_AA.DunedinPoAm, W_AgeAccelGrim_bl, W_AgeAccelGrim_fu, 
         W_AgeAccelPheno_bl, W_AgeAccelPheno_fu)

cors <- function(df) { 
  M <- Hmisc::rcorr(as.matrix(df))
  Mdf <- map(M, ~data.frame(.x))
}

aa_data <- cors(aa_dat) %>%
  map(~rownames_to_column(.x, var="measure1")) %>%
  map(~pivot_longer(.x, -measure1, "measure2")) %>%
  bind_rows(.id = "id") %>%
  pivot_wider(names_from = id, values_from = value) %>% 
  mutate(across(contains("measure"), recode, 
                W_AA.Zhang = "Zhang AA", 
                W_AA.DunedinPoAm = "DunedinPoAm AA",
                W_AgeAccelGrim_bl = "BL GrimAge AA",
                W_AgeAccelGrim_fu = "FU GrimAge AA",
                W_AgeAccelPheno_bl = "BL PhenoAge AA",
                W_AgeAccelPheno_fu = "FU PhenoAge AA"))

aa_data %>%
  ggplot(aes(measure1, measure2, fill=r, label=round(r,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text(size =3) +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90))

ggsave("AA heatmap.png", device = "png", dpi = 300)




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

crPlots(m1,~ scale(logW_crp_g_fu2))


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

crPlots(m1,~ scale(logW_crp_g_fu2))

# model 2 (sensitivity analysis)

m2 <- lm(W_AgeAccelGrim_fu ~ scale(logW_crp_g_fu2) + age_fu_correct + 
           cob_cde_bl + sex_cde_bl + bmi_rrto + seifa_10_fu + 
           as.factor(phys.act) + tot_alclw_mrt + cigst_cde + tot_alclw_mrt + 
           logpack + educlvl_ord, data = full)

tidy(m2, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)


## Dunedin
# model 1

m1 <- lm(scale(W_AA.DunedinPoAm) ~ scale(logW_crp_g_fu2) + cob_cde_bl + sex_cde_bl + 
           bmi_rrto + cigst_cde + educlvl_ord + age_fu_correct, data = full)

tidy(m1, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

crPlots(m1,~ scale(logW_crp_g_fu2))

# model 2 

m2 <- lm(scale(W_AA.DunedinPoAm) ~ scale(logW_crp_g_fu2) + age_fu_correct + 
           cob_cde_bl + sex_cde_bl + bmi_rrto + seifa_10_fu + 
           as.factor(phys.act) + tot_alclw_mrt + cigst_cde + tot_alclw_mrt + 
           logpack + educlvl_ord, data = full)

tidy(m2, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

## Zhang
# model 1

m1 <- lm(W_AA.Zhang ~ scale(logW_crp_g_fu2) + cob_cde_bl + sex_cde_bl + 
           bmi_rrto + cigst_cde + educlvl_ord + age_fu_correct, data = full)

tidy(m1, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

crPlots(m1,~ scale(logW_crp_g_fu2))

# model 2 

m2 <- lm(W_AA.Zhang ~ scale(logW_crp_g_fu2) + age_fu_correct + 
           cob_cde_bl + sex_cde_bl + bmi_rrto + seifa_10_fu + 
           as.factor(phys.act) + tot_alclw_mrt + cigst_cde + tot_alclw_mrt + 
           logpack + educlvl_ord, data = full)

tidy(m2, conf.int = T) %>% select(term, estimate, conf.low, conf.high, p.value)

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



### R squared analysis 

trim <- full %>% select(starts_with("logW"), W_AgeAccelGrim_fu, age_fu_correct, cob_cde_bl, sex_cde_bl, bmi_rrto, educlvl_ord) %>% 
  select(ends_with("fu"), ends_with("fu2"), W_AgeAccelGrim_fu, age_fu_correct, cob_cde_bl, sex_cde_bl, bmi_rrto, educlvl_ord)

m1 <- lm(W_AgeAccelGrim_fu ~ age_fu_correct + cob_cde_bl + sex_cde_bl + bmi_rrto + educlvl_ord, data = trim)
m2 <- lm(W_AgeAccelGrim_fu ~ ., data = trim)

summary(m1)

summary(m2)
