# Load the required packages

library(tidyverse)
library(flextable)
library(psych)
library(ggpubr)
library(psych)
library(broom)
library(ggthemes)
library(knitr)
library(car)
library(glmnet)

#### Prepare data ####

### read data

## baseline

bl_new <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\LP.BL.csv")

bl_old <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.BL.csv")

aa_bl <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\AA.BL.csv")

bl_imp <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\LP.BL_imp.csv")

# Merge baseline data files 

bl <- left_join(bl_new, bl_old, by = "h2_id", suffix = c("","_old"))

bl <- bl %>% select(-ends_with("_old")) # remove redundant variables

bl2 <- left_join(bl, aa_bl, by = "h2_id", suffix = c("_old", ""))

bl2 <- bl2 %>% select(-ends_with("_old")) # remove redundant variables

# add in imputed data

bl3 <- left_join(bl2, bl_imp, by = "h2_id", suffix = c("","_imp"))

# remove cancer cases, breast & gastric cancer, duplicates

bl3 <- bl3 %>% 
  filter(!study %in% c("BC","GAC")) %>% 
  filter(!case == 1) %>% 
  filter(!duplicated(h2_id))

# add _bl suffix to all variables (except participant ID)

bl3 <- bl3 %>% 
  rename_at(vars(-h2_id), function(x) paste0(x,"_bl"))


## Follow-up

fu_new <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\LP.FUP2.csv")

fu_old <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.FUP2.csv")

fu_imp <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\LP.FUP2_imp.csv")

# merge follow-up files

fu <- left_join(fu_new, fu_old, by = "h2_id", suffix = c("","_old"))

fu <- fu %>% select(-ends_with("_old")) # remove redundant variables

fu2 <- left_join(fu, fu_imp, by = "h2_id", suffix = c("", "_imp")) # add in imputed data

# add fu suffix to all variables (except participant ID)

fu2 <- fu2 %>% 
  rename_at(vars(-h2_id), function(x) paste0(x,"_fu"))

## Merge baseline and follow-up 

full <- inner_join(bl3, fu2, by = "h2_id")

# load corrected age data

load("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.Age.rda")

## Add corrected age variable into merged file 

LP.Age <- LP.Age %>% 
  filter(timepoint_cde == "Baseline") %>% 
  select(h2_id, atnd_age,age_fup2) %>% 
  rename(age_bl_correct = atnd_age, 
         age_fu_correct = age_fup2)

full <- left_join(full, LP.Age, by = "h2_id")

# remove rows without follow-up epigenetic age

full <- full %>% filter(!is.na(DNAmGrimAge_fu))

# Calculate time to follow-up

full <- full %>% mutate(time_fu = age_fu_correct - age_bl_correct)

hist(full$time_fu, breaks = 50) # minimum of 9 years 

# Calculate age acceleration variables

full <- full %>% 
  mutate(AgeAccelGrim_bl = residuals(lm(DNAmGrimAge_bl ~ age_bl_correct, data = full)),
         AgeAccelGrim_fu = residuals(lm(DNAmGrimAge_fu ~ age_fu_correct, data = full)),
         AgeAccelZhang_bl = residuals(lm(score.Zhang.cont_bl ~ age_bl_correct, data = full)),
         AgeAccelZhang_fu = residuals(lm(score.Zhang.cont_fu ~ age_fu_correct, data = full)),
         AgeAccelDunedin_bl = residuals(lm(DunedinPoAm_bl ~ age_bl_correct, data = full)),
         AgeAccelDunedin_fu = residuals(lm(DunedinPoAm_fu ~ age_fu_correct, data = full)),
         AgeAccelPheno_bl = residuals(lm(DNAmPhenoAge_bl ~ age_bl_correct, data = full)),
         AgeAccelPheno_fu = residuals(lm(DNAmPhenoAge_fu ~ age_fu_correct, data = full)))

#### Table 1: baseline demographics ####

psych::describe(full$age_bl_correct)

full %>% group_by(sex_cde_bl) %>% tally() %>% mutate(percent = n/sum(n) * 100)

psych::describe(full$bmi_rrto_bl)

psych::describe(full$time_fu)

full %>% group_by(cob_cde_bl) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full %>% group_by(sample_type_bl) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full %>% group_by(cigst_cde_bl) %>% tally() %>% mutate(percent = n/sum(n) * 100)

full <- full %>% mutate(packyears_bl = if_else(logpack_bl > 0, exp(logpack_bl) / 20, 0))

hist(full$packyears_bl, breaks = 50)

full %>% filter(packyears_bl > 0) %>% summarise(median = median(packyears_bl),
                                                IQR = IQR(packyears_bl))

full %>% summarise(median = median(seifa_10_bl,na.rm=T), 
                   IQR = IQR(seifa_10_bl,na.rm=T))

psych::describe(full$time_fu)

#### Transform and prepare biomarkers ####

# vector of all biomarker names

marker_names <- as.character(
  expression(neopt_d_bl,crp_g_bl,cnct_g_bl,s100at_g_bl,saat_g_bl,il6_msd_bl,
             il8_msd_bl,il10_msd_bl,ifng_msd_bl,tnfa_msd_bl,trp_d_bl,kyn_d_bl,
             hk_d_bl,ka_d_bl,xa_d_bl,aa_d_bl,haa_d_bl,pic_d_bl,qa_d_bl,KTr_bl,
             PAr_bl,HKXAr_bl,neopt_d_fu,crp_g_fu,cnct_g_fu,s100at_g_fu,
             saat_g_fu,il6_msd_fu,il8_msd_fu,il10_msd_fu,ifng_msd_fu,
             tnfa_msd_fu,trp_d_fu,kyn_d_fu,hk_d_fu,ka_d_fu,xa_d_fu,aa_d_fu,
             haa_d_fu,pic_d_fu,qa_d_fu,KTr_fu,PAr_fu,HKXAr_fu, cysta_d_bl, cysta_d_fu))

## count missing values from each marker 

# baseline

full %>% 
  select(marker_names) %>% select(ends_with("_bl")) %>% 
  summarise(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  summarise(median = median(value),
            min = min(value),
            max = max(value))

# follow-up

full %>% 
  select(marker_names) %>% select(ends_with("_fu")) %>% 
  summarise(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  summarise(median = median(value),
            min = min(value),
            max = max(value))

## Count missing data in confounders

confounders <- c("cob_cde_bl", "cob_cde_fu", "sex_cde_bl","sex_cde_fu", 
                 "bmi_rrto_bl", "bmi_rrto_fu", "cigst_cde_bl","cigst_cde_fu",
                 "educlvl_ord_bl","educlvl_ord_fu","age_bl_correct", "age_fu_correct")

# baseline 

full %>% 
  select(confounders) %>% 
  summarise(across(everything(), ~ sum(is.na(.)))) %>% 
  pivot_longer(everything()) %>% 
  summarise(median = median(value),
            min = min(value),
            max = max(value))

### From here on, using *imputed* data 

marker_names <- str_replace(marker_names, "_bl", "_imp_bl")

marker_names <- str_replace(marker_names, "_fu", "_imp_fu")

## Count zero values (below LOD) for each biomarker variable

full %>% 
  select(marker_names) %>%  
  summarise(across(everything(), ~ sum(. == 0))) %>% 
  pivot_longer(everything()) %>% 
  summarise(median = median(value),
            min = min(value),
            max = max(value))


# Plot histograms to check for outliers/issues

full %>% select(all_of(marker_names)) %>% 
  select(neopt_d_imp_bl:xa_d_imp_bl) %>% 
  multi.hist(global = F)

full %>% select(all_of(marker_names)) %>% 
  select(aa_d_imp_bl:il10_msd_imp_fu) %>% 
  multi.hist(global = F)

full %>% select(all_of(marker_names)) %>% 
  select(ifng_msd_imp_fu:cysta_d_imp_fu) %>% 
  multi.hist(global = F)

## Several biomarkers have zero values indicating values below the LOD. 
## these are replaced with half the LOD for each assay

full$crp_g_imp_bl2 = ifelse(full$crp_g_imp_bl < 0.001, 0.06 / 2, full$crp_g_imp_bl)

full$crp_g_imp_fu2 = ifelse(full$crp_g_imp_fu < 0.001, 0.06 / 2, full$crp_g_imp_fu)

full$saat_g_imp_bl2 <- ifelse(full$saat_g_imp_bl < 0.001, 0.04 / 2, full$saat_g_imp_bl)

full$il10_msd_imp_bl2 <- ifelse(full$il10_msd_imp_bl < 0.001, 0.09 / 2, full$il10_msd_imp_bl)

full$il10_msd_imp_fu2 <- ifelse(full$il10_msd_imp_fu < 0.001, 0.09 / 2, full$il10_msd_imp_fu)

full$ifng_msd_imp_bl2 <- ifelse(full$ifng_msd_imp_bl < 0.001, 0.396/2, full$ifng_msd_imp_bl)

full$il6_msd_imp_fu2 <- ifelse(full$il6_msd_imp_fu < 0.001, 0.09/2, full$il6_msd_imp_fu)

### Transform and winsorise biomarkers
# create function to winsorise variables at 3 standard deviations. 
# variables beyond the 3 SD threshold are replaced with the largest 
# absolute value within this threshold. 

winsorise <- function(x){
  maxin <- max(x[scale(x) < 3],na.rm=T)
  minin <- min(x[scale(x) > -3],na.rm=T)
  x <- ifelse(x > maxin & !is.na(x), maxin, x)
  x <- ifelse(x < minin & !is.na(x), minin, x)
  x
}

# update marker names 

marker_names <- as.character(
  expression(neopt_d_bl,crp_g_bl2,cnct_g_bl,s100at_g_bl,saat_g_bl2,il6_msd_bl,
             il8_msd_bl,il10_msd_bl2,ifng_msd_bl2,tnfa_msd_bl,trp_d_bl,kyn_d_bl,
             hk_d_bl,ka_d_bl,xa_d_bl,aa_d_bl,haa_d_bl,pic_d_bl,qa_d_bl,KTr_bl,
             PAr_bl,HKXAr_bl,neopt_d_fu,crp_g_fu2,cnct_g_fu,s100at_g_fu,
             saat_g_fu,il6_msd_fu2,il8_msd_fu,il10_msd_fu2,ifng_msd_fu,
             tnfa_msd_fu,trp_d_fu,kyn_d_fu,hk_d_fu,ka_d_fu,xa_d_fu,aa_d_fu,
             haa_d_fu,pic_d_fu,qa_d_fu,KTr_fu,PAr_fu,HKXAr_fu, cysta_d_bl, cysta_d_fu))

marker_names <- str_replace(marker_names, "_bl", "_imp_bl")

marker_names <- str_replace(marker_names, "_fu", "_imp_fu")

# log transform, winsorise, and scale biomarkers to z scores

full <- full %>% 
  mutate(across(marker_names, ~ log(.x), .names = "wlog_{.col}")) %>% 
  mutate(across(starts_with("wlog_"), winsorise)) %>% 
  mutate(across(starts_with("wlog_"), scale, .names = "z{.col}"))

# winsorise epigenetic age variables

epi_vars <- full %>% 
  select(starts_with("AgeAccel")) %>% 
  select(-contains("Residual")) %>% 
  names()

full <- full %>% 
  mutate(across(all_of(epi_vars), winsorise,
                .names = "w{.col}"))

# standardize epigenetic ageing measures 

full <- full %>% 
  mutate(across(c("wAgeAccelZhang_fu", "wAgeAccelDunedin_fu", "wAgeAccelGrim_fu",
                  "wAgeAccelPheno_fu", "wAgeAccelZhang_bl", "wAgeAccelDunedin_bl",
                  "wAgeAccelGrim_bl", "wAgeAccelPheno_bl"), scale,
                .names = "z{.col}"))

# histograms of transformed markers

full %>% select(wlog_neopt_d_imp_bl:wlog_xa_d_imp_bl) %>% 
  multi.hist(global = F)

full %>% select(wlog_aa_d_imp_bl:wlog_il10_msd_imp_fu2) %>% 
  multi.hist(global = F)

full %>% select(wlog_ifng_msd_imp_fu:wlog_cysta_d_imp_fu) %>% 
  multi.hist(global = F)

full %>% 
  select(wAgeAccelGrim_bl, wAgeAccelGrim_fu,
         wAgeAccelPheno_bl, wAgeAccelPheno_fu,
         wAgeAccelDunedin_bl, wAgeAccelDunedin_fu,
         wAgeAccelZhang_bl, wAgeAccelZhang_fu) %>% 
  multi.hist(global = F)

### calculate inflammaging signature from baseline and follow-up inflam data
## Pierre's LogZ (natural log) variables are used instead of the newly created
## variables. This is to be sure the variables match those from the
## 2021 Journals of Gerontology: Series A paper

# using the mean age at each visit as the intercept 

full <- full %>% 
  mutate(across(starts_with("LogZ_"), winsorise, .names = "W{.col}")) %>% 
  rowwise() %>% 
  mutate(inflamm_sig_bl = 57.5 - 1.33*WLogZ_nam_d_imp_bl - 0.22*WLogZ_trp_d_imp_bl + 0.28*WLogZ_aa_d_imp_bl + 
           0.90*WLogZ_qa_d_imp_bl + 0.76*WLogZ_neopt_d_imp_bl + 0.65*WLogZ_cysta_d_imp_bl + 
           0.09*WLogZ_il6_msd_imp_bl + 0.19*WLogZ_il8_msd_imp_bl + 0.63*WLogZ_PAr_imp_bl + 0.76*WLogZ_HKAAr_imp_bl) %>% 
  mutate(inflamm_sig_fu = 68.9 - 1.33*WLogZ_nam_d_imp_fu - 0.22*WLogZ_trp_d_imp_fu + 0.28*WLogZ_aa_d_imp_fu + 
           0.90*WLogZ_qa_d_imp_fu + 0.76*WLogZ_neopt_d_imp_fu + 0.65*WLogZ_cysta_d_imp_fu + 
           0.09*WLogZ_il6_msd_imp_fu + 0.19*WLogZ_il8_msd_imp_fu + 0.63*WLogZ_PAr_imp_fu + 0.76*WLogZ_HKAAr_imp_fu) %>%
  ungroup()

full %>% select(contains("inflamm")) %>% multi.hist()

# winsorise and scale to Z scores inflammaging variables

full <- full %>% 
  mutate(across(contains("inflamm"), winsorise, .names = "w{.col}")) %>% 
  mutate(across(contains("winflamm"), scale, .names = "z{.col}"))


#### Table 2 ####

## Baseline and follow-up biomarker summary ##

# median and IQR at baseline

bl_res <- full %>% 
  select(marker_names) %>% 
  summarise(across(contains("_bl"), 
                   list(medianBL = ~ median(.,na.rm=T), IQRBL = ~ IQR(.,na.rm=T)))) %>% 
  pivot_longer(cols = everything(), names_sep = "_(?!.*_)", 
               names_to = c("variable",".value")) %>% 
  mutate(variable = sub("_[^_]+$", "", .$variable))

bl_res

# median and IQR at fu

fu_res <- full %>% 
  select(marker_names) %>% 
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

#print(table1, preview = "docx")

### Correlation with age at follow-up ###

# Names of all the biomarker variables at FU

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  select(-zwlog_cysta_d_imp_fu) %>% 
  names()

# Spearman's correlation between each biomarker and age at FU 

cors <- map(xs,
     ~ cor.test(full[[.x]], full$age_fu_correct, method = "spearman")) %>% 
  map_df(tidy) %>% 
  cbind(xs, .)
  
knitr::kable(cors, digits = 2)


#### Regression models ####

# sex as a factor

full$sex_cde_bl <- as.factor(full$sex_cde_bl)

full$sex_cde_fu <- as.factor(full$sex_cde_fu)

#### Figure 1 - Cross sectional association at BL #### 

# create vector including all baseline biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()


# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_bl", "zwAgeAccelGrim_bl", "zwAgeAccelZhang_bl", "zwAgeAccelDunedin_bl")

# regression for each combination of biomarker and AgeAccel variable 


bl_cs_results <- 
  
  # every combination of biomarker and AgeAccel
  
  crossing(Var1 = xs, Var2 = ys) %>%
  
  # perform linear regression adjusted for age and sex and save results
  
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_bl + cob_cde_bl + "),  
         models = map(frm, 
                       ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  
  # Drop unneeded parameters and rename variables for plotting 
  
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(!term %in% c("(Intercept)", "sex_cde_blMale", "cob_cde_blNorthern Europe", "cob_cde_blSouthern Europe")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_bl = "AgeAccelDunedin",
                               zwAgeAccelZhang_bl = "AgeAccelZhang",
                               zwAgeAccelGrim_bl = "AgeAccelGrim",
                               zwAgeAccelPheno_bl = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

bl_cs_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("Baseline forest.png", device = "png",
       height = 7, width = 10)



#### Figure 2 - Cross sectional association at FU #### 

# create vector including all FU biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  select(-zwlog_cysta_d_imp_fu) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 


fu_cs_results <- crossing(Var1 = xs, Var2 = ys) %>%
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_fu + cob_cde_fu + "),  
         models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(!term %in% c("(Intercept)", "sex_cde_fuMale", "cob_cde_fuNorthern Europe", "cob_cde_fuSouthern Europe")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

fu_cs_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("follow-up forest.png", device = "png",
       height = 7, width = 10)


#### Figure 3 - Prospective associations #### 

# create vector including all BL biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 

prosp_results <- 
  
  # all combinations of biomarker and AgeAccel vars
  
  crossing(Var1 = xs, Var2 = ys) %>%
  
  # create formula for each regression model 
  
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_fu + cob_cde_fu + ")) %>% 
  
  # Add baseline AgeAccel variable to formula 
  
  mutate(Var3 = sub("_fu", "_bl", Var2)) %>% 
  mutate(frm = str_c(frm, Var3, sep = " + ")) %>% 
  select(-Var3) %>% 
  
  # Regressions and save results 
  
  mutate(models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  
  # format results for plotting 
  
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(str_detect(term, "_bl"), !term %in% c("zwAgeAccelDunedin_bl", "zwAgeAccelZhang_bl", "zwAgeAccelGrim_bl", "zwAgeAccelPheno_bl")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

prosp_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("Prospective forest.png", device = "png",
       height = 7, width = 10)


#### Figure 4 - Change over time ####

# create vector including all BL biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 

change_results <- 
  
  # create formulas with confounder variables 
  
  crossing(Var1 = xs, Var2 = ys) %>%
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_fu + cob_cde_fu + ")) %>% 
  
  # Add baseline AgeAccel and baseline biomarker variables to each model
  
  mutate(Var3 = sub("_fu", "_bl", Var2)) %>% 
  mutate(Var4 = sub("_bl", "_fu", Var1)) %>% 
  
  # Correct variable names
  
  mutate(Var4 = dplyr::recode(Var4,
    zwlog_ifng_msd_imp_fu2 = 'zwlog_ifng_msd_imp_fu',
    zwlog_saat_g_imp_fu2 = 'zwlog_saat_g_imp_fu',
    zwlog_il6_msd_imp_fu = 'zwlog_il6_msd_imp_fu2'
  )) %>% 

  # add confounders to regression
  
  mutate(frm = str_c(frm, Var3, sep = " + ")) %>% 
  mutate(frm = str_c(frm, Var4, sep = " + ")) %>% 
  select(-Var3, -Var4) %>% 
  
  # regressions 
  
  mutate(models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-Var1, -frm) %>% 
  
  # prepare data for plotting 
  
  rename('clock' = Var2) %>% 
  filter(str_detect(term, "imp_fu") | str_detect(term, "sig_fu")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

change_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))

ggsave("Longitudinal forest.png", device = "png",
       height = 7, width = 10)


#### Comparing estimates ####

## Baseline ##

# strongest results 

bl_cs_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value)

# proportion with FDR p < 0.05


bl_cs_results %>% 
  mutate(fdr_p = p.adjust(bl_cs_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# proportion with estimate > 0 

bl_cs_results %>% 
  mutate(fdr_p = p.adjust(bl_cs_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  mutate(positive = if_else(estimate > 0, "Yes", "No")) %>% 
  group_by(sig, positive) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

## Follow-up ##

# strongest results 

fu_cs_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value)

# proportion with FDR p < 0.05

fu_cs_results %>% 
  mutate(fdr_p = p.adjust(fu_cs_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# proportion with estimate > 0 

fu_cs_results %>% 
  mutate(fdr_p = p.adjust(fu_cs_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  mutate(positive = if_else(estimate > 0, "Yes", "No")) %>% 
  group_by(sig, positive) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# Result with estimate < 0

fu_cs_results %>% 
  mutate(fdr_p = p.adjust(fu_cs_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  filter(sig == "Yes" & estimate < 0)

## comparing baseline and follow-up 

## comparison ##

estimates <- tibble(
  term = bl_cs_results$term2,
  clock = bl_cs_results$clock,
  bl = bl_cs_results$estimate,
  fu = fu_cs_results$estimate,
  prosp = prosp_results$estimate,
  change = change_results$estimate
)

# largest differences from baseline to FU

estimates %>% 
  mutate(change_est = fu - bl) %>% 
  arrange(desc(abs(change_est)))

# Baseline and follow-up correlation

cor(x = bl_cs_results$estimate, y = fu_cs_results$estimate)

## propsective ##

# strongest results 

prosp_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value)

# proportion with p < Bonferroni threshold (0.0021)

prosp_results %>% 
  mutate(fdr_p = p.adjust(prosp_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)


## change model ##

# strongest results 

change_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value)

# proportion with p < Bonferroni threshold (0.0021)

change_results %>% 
  mutate(fdr_p = p.adjust(change_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# proportion with estimate > 0 

change_results %>% 
  mutate(fdr_p = p.adjust(change_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  mutate(positive = if_else(estimate > 0, "Yes", "No")) %>% 
  group_by(sig, positive) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# negative estimate

change_results %>% 
  mutate(fdr_p = p.adjust(change_results$p.value, method = "fdr")) %>% 
  mutate(sig = if_else(fdr_p < 0.05, "Yes", "No")) %>% 
  filter(sig == "Yes" & estimate < 0)




#### Variance explained (ridge regression) ####

### Baseline ###

## Grim

# predictors

X <- as.matrix(full %>% 
  select(starts_with("zwlog_"), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")))

# outcome

grim <- as.matrix(full$zwAgeAccelGrim_bl)

# range of lambda values

lambda_seq <- 10^seq(2, -2, by = -.05)

# find best lambda

grim_ridge <- cv.glmnet(x = X, y = grim, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(grim_ridge)

best_lambda <- grim_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = grim, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 


## Pheno ##

# outcome

pheno <- as.matrix(full$zwAgeAccelPheno_bl)

# find best lambda

pheno_ridge <- cv.glmnet(x = X, y = pheno, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(pheno_ridge)

best_lambda <- pheno_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = pheno, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 


## Zhang ##

# outcome

zhang <- as.matrix(full$zwAgeAccelZhang_bl)

# find best lambda

zhang_ridge <- cv.glmnet(x = X, y = zhang, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(zhang_ridge)

best_lambda <- zhang_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = zhang, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 

## Dunedin ##

# outcome

dunedin <- as.matrix(full$zwAgeAccelDunedin_bl)

# find best lambda

dunedin_ridge <- cv.glmnet(x = X, y = dunedin, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(dunedin_ridge)

best_lambda <- dunedin_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = dunedin, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 



### Follow-up ###

## Grim

# predictors

X <- as.matrix(full %>% 
                 select(starts_with("zwlog_"), zwinflamm_sig_fu) %>% 
                 select(ends_with("fu"), ends_with("fu2")))

# outcome

grim <- as.matrix(full$zwAgeAccelGrim_fu)

# find best lambda

grim_ridge <- cv.glmnet(x = X, y = grim, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(grim_ridge)

best_lambda <- grim_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = grim, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 


## Pheno ##

# outcome

pheno <- as.matrix(full$zwAgeAccelPheno_fu)

# find best lambda

pheno_ridge <- cv.glmnet(x = X, y = pheno, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(pheno_ridge)

best_lambda <- pheno_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = pheno, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 


## Zhang ##

# outcome

zhang <- as.matrix(full$zwAgeAccelZhang_fu)

# find best lambda

zhang_ridge <- cv.glmnet(x = X, y = zhang, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(zhang_ridge)

best_lambda <- zhang_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = zhang, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 

## Dunedin ##

# outcome

dunedin <- as.matrix(full$zwAgeAccelDunedin_fu)

# find best lambda

dunedin_ridge <- cv.glmnet(x = X, y = dunedin, family = "gaussian", lambda = lambda_seq, alpha = 0)

plot(dunedin_ridge)

best_lambda <- dunedin_ridge$lambda.min

# best model

fit_best <- glmnet(x = X, y = dunedin, family = "gaussian", lambda = best_lambda, alpha = 0)

coef(fit_best)

fit_best$dev.ratio # r squared 



#### Supplementary Figure 1 -correlation matrix ####
## Correlation matrix of inflammation/TK biomarkers and age accel variables
## at follow-up 

aa_vars <- full %>% 
  select(wAgeAccelZhang_fu, wAgeAccelDunedin_fu, wAgeAccelGrim_fu,
         wAgeAccelPheno_fu, wAgeAccelZhang_bl, wAgeAccelDunedin_fu,
         wAgeAccelGrim_bl, wAgeAccelPheno_bl) %>% 
  names()

cor_dat <- full %>% 
  select(starts_with("wlog_"), all_of(aa_vars), winflamm_sig_fu) %>% 
  select(ends_with("_fu"), ends_with("fu2")) %>% 
  select(-wlog_cysta_d_imp_fu)

cors_spearman <- function(df) { 
  M <- Hmisc::rcorr(as.matrix(df), type = "spearman")
  Mdf <- map(M, ~data.frame(.x))
}

plot_data <- cors_spearman(cor_dat) %>%
  map(~rownames_to_column(.x, var="measure1")) %>%
  map(~pivot_longer(.x, -measure1, "measure2")) %>%
  bind_rows(.id = "id") %>%
  pivot_wider(names_from = id, values_from = value) %>% 
  mutate(across(c("measure1","measure2"), str_replace, "wlog_", "")) %>% 
  mutate(measure1 = sub("_.*", "", .$measure1)) %>% 
  mutate(measure2 = sub("_.*", "", .$measure2)) %>% 
  mutate(across(contains("measure"), str_to_upper)) %>% 
  mutate(across(contains("measure"), ~ dplyr::recode(.,
                                                     WAGEACCELZHANG = "AgeAccelZhang",
                                                     WAGEACCELDUNEDIN = "AgeAccelDunedin",
                                                     WAGEACCELGRIM = "AgeAccelGrim",
                                                     WAGEACCELPHENO = "AgeAccelPheno",
                                                     WINFLAMM = "Inflammaging signature",
                                                     NEOPT = "Neopterin",
                                                     CNCT = "Cystatin C",
                                                     S100AT = "Calprotectin",
                                                     SAAT = "Serum amyloid A",
                                                     IL8 = "Interleukin-8",
                                                     IL10 = "Interleukin-10",
                                                     IFNG = "Interferon-g",
                                                     TNFA = "TNF-a",
                                                     TRP = "Tryptophan",
                                                     KYN = "Kynurenine",
                                                     HK = "3-Hydroxykynurenine",
                                                     KA = "Kynurenic acid",
                                                     XA = "Xanthurenic acid",
                                                     AA = "Anthranilic acid",
                                                     HAA = "3-Hydroxyanthranilic acid",
                                                     PIC = "Piconilic acid",
                                                     QA = "Quinolinic acid",
                                                     KTR = "KTr",
                                                     PAR = "PAr",
                                                     HKXAR = "HK:XA",
                                                     CRP = "C-reactive protein",
                                                     IL6 = "Interleukin-6")))

order <-  c("AgeAccelZhang", "AgeAccelDunedin", "AgeAccelGrim",
            "AgeAccelPheno", "Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

plot_data %>%
  ggplot(aes(factor(measure1, levels = order), factor(measure2, levels = rev(order)), 
             fill=r, label=round(r,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nCorrelation") + 
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text(size =2.5) +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 10),
        axis.text.y = element_text(size = 10))

ggsave("correlation heatmap.png", device = "png", dpi = 450,
       width = 12, height = 8)

#### Supplementary Figure 2 - age correlations ####

# create text to annotate plot 

ann_text <- data.frame(
  label = c("r = 0.22", "r = 0.85", "r = 0.74", "r = 0.41", "r = 0.50"),
  variable = c("DunedinPoAm", "GrimAge","PhenoAge","Zhang", "Inflammaging"),
  x = rep(-2, 5),
  y = rep(88, 5)
)

full %>% 
  mutate(across(c("wscore.Zhang.cont_fu", "wDunedinPoAm_fu", 
                  "wDNAmGrimAge_fu", "wDNAmPhenoAge_fu", "winflamm_sig_fu"), scale)) %>% 
  pivot_longer(cols = c("wscore.Zhang.cont_fu", "wDunedinPoAm_fu", 
                        "wDNAmGrimAge_fu", "wDNAmPhenoAge_fu", "winflamm_sig_fu"),
               names_to = "variable", values_to = "ageaccel") %>%
  mutate(variable = dplyr::recode(variable,
                                  wscore.Zhang.cont_fu = "Zhang",
                                  wDunedinPoAm_fu = "DunedinPoAm",
                                  wDNAmGrimAge_fu = "GrimAge",
                                  wDNAmPhenoAge_fu = "PhenoAge",
                                  winflamm_sig_fu = "Inflammaging")) %>%
  ggplot(aes(x = ageaccel, y = age_fu_correct)) +
  geom_point() + geom_smooth(method = "lm", colour = "darkblue") + 
  facet_wrap(~ factor(variable,
                      levels = c("Zhang", "DunedinPoAm",
                                 "GrimAge", "PhenoAge",
                                 "Inflammaging"))) +
  labs(x = "Epigenetic ageing/inflammaging biomarker (scaled)", y = "Chronological age") +
  theme_stata() +
  geom_text(
    data    = ann_text,
    mapping = aes(x = x, y = y, label = label)
  )

ggsave("clock age correlation.png", device = "png",
       width = 7, height = 5)

bl_cs_results %>% 
  filter(term == "CRP")
#### Supplementary Figures 3-7 - Sensitivity confounders ####

full$cigst_cde_imp_bl <- as.factor(full$cigst_cde_imp_bl)

levels(full$cigst_cde_imp_bl) <- c("Never", "Current", "Former")

full$cigst_cde_imp_fu <- as.factor(full$cigst_cde_imp_fu)

levels(full$cigst_cde_imp_fu) <- c("Never", "Current", "Former")

### Figure S1 - Cross sectional association at BL ###

# create vector including all baseline biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()


# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_bl", "zwAgeAccelGrim_bl", "zwAgeAccelZhang_bl", "zwAgeAccelDunedin_bl")

# regression for each combination of biomarker and AgeAccel variable 

bl_sf3 <- 
  
  # every combination of biomarker and AgeAccel
  
  crossing(Var1 = xs, Var2 = ys) %>%
  
  # perform linear regression adjusted for age and sex and save results
  
  mutate(frm = str_c(Var2, Var1, sep = " ~ bmi_rrto_imp_bl + cigst_cde_imp_bl + seifa_10_imp_bl + sex_cde_bl + cob_cde_bl + "),  
         models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  
  # Drop unneeded parameters and rename variables for plotting 
  
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(!term %in% c("(Intercept)", "sex_cde_blMale", "cob_cde_blNorthern Europe", "cob_cde_blSouthern Europe",
                      "bmi_rrto_imp_bl", "cigst_cde_imp_blCurrent", "cigst_cde_imp_blFormer", "seifa_10_imp_bl")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_bl = "AgeAccelDunedin",
                               zwAgeAccelZhang_bl = "AgeAccelZhang",
                               zwAgeAccelGrim_bl = "AgeAccelGrim",
                               zwAgeAccelPheno_bl = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

bl_sf3 %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("BL supplementary figure 3.png", device = "png",
       height = 7, width = 10)

## Compare to unadjusted results

bl_sf3 %>% 
  mutate(unadj_est = bl_cs_results$estimate) %>% 
  mutate(dif = estimate - unadj_est) %>% 
  group_by(clock) %>% 
  summarise(mean = mean(dif))

### Supplementary figure 4 - Cross sectional association at FU ###

# create vector including all FU biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  select(-zwlog_cysta_d_imp_fu) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 


fu_cs_results <- crossing(Var1 = xs, Var2 = ys) %>%
  mutate(frm = str_c(Var2, Var1, sep = " ~ bmi_rrto_imp_fu + cigst_cde_imp_fu + seifa_10_imp_fu + sex_cde_fu + cob_cde_fu + "),  
         models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(!term %in% c("(Intercept)", "sex_cde_fuMale", "cob_cde_fuNorthern Europe", "cob_cde_fuSouthern Europe")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))

# strongest results 

fu_cs_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value)

# proportion with p < Bonferroni threshold (0.0021)

fu_cs_results %>% 
  mutate(sig = if_else(p.value < 0.0021, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# proportion with estimate > 0 

fu_cs_results %>% 
  mutate(sig = if_else(p.value < 0.0021, "Yes", "No")) %>% 
  mutate(positive = if_else(estimate > 0, "Yes", "No")) %>% 
  group_by(sig, positive) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

fu_cs_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("follow-up forest.png", device = "png",
       height = 7, width = 10)


#### Figure 3 - Prospective associations #### 

# create vector including all BL biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 

prosp_results <- 
  
  # all combinations of biomarker and AgeAccel vars
  
  crossing(Var1 = xs, Var2 = ys) %>%
  
  # create formula for each regression model 
  
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_fu + cob_cde_fu + ")) %>% 
  
  # Add baseline AgeAccel variable to formula 
  
  mutate(Var3 = sub("_fu", "_bl", Var2)) %>% 
  mutate(frm = str_c(frm, Var3, sep = " + ")) %>% 
  select(-Var3) %>% 
  
  # Regressions and save results 
  
  mutate(models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  
  # format results for plotting 
  
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(str_detect(term, "_bl"), !term %in% c("zwAgeAccelDunedin_bl", "zwAgeAccelZhang_bl", "zwAgeAccelGrim_bl", "zwAgeAccelPheno_bl")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))

# strongest results 

prosp_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value)

# proportion with p < Bonferroni threshold (0.0021)

prosp_results %>% 
  mutate(sig = if_else(p.value < 0.0021, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

prosp_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("Prospective forest.png", device = "png",
       height = 7, width = 10)


#### Figure 4 - Change over time ####

# create vector including all BL biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 

change_results <- 
  
  # create formulas with confounder variables 
  
  crossing(Var1 = xs, Var2 = ys) %>%
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_fu + cob_cde_fu + ")) %>% 
  
  # Add baseline AgeAccel and baseline biomarker variables to each model
  
  mutate(Var3 = sub("_fu", "_bl", Var2)) %>% 
  mutate(Var4 = sub("_bl", "_fu", Var1)) %>% 
  
  # Correct variable names
  
  mutate(Var4 = dplyr::recode(Var4,
                              zwlog_ifng_msd_imp_fu2 = 'zwlog_ifng_msd_imp_fu',
                              zwlog_saat_g_imp_fu2 = 'zwlog_saat_g_imp_fu',
                              zwlog_il6_msd_imp_fu = 'zwlog_il6_msd_imp_fu2'
  )) %>% 
  
  # add confounders to regression
  
  mutate(frm = str_c(frm, Var3, sep = " + ")) %>% 
  mutate(frm = str_c(frm, Var4, sep = " + ")) %>% 
  select(-Var3, -Var4) %>% 
  
  # regressions 
  
  mutate(models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-Var1, -frm) %>% 
  
  # prepare data for plotting 
  
  rename('clock' = Var2) %>% 
  filter(str_detect(term, "imp_fu") | str_detect(term, "sig_fu")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWINFLAMM = "Inflammaging signature",
                               NEOPT = "Neopterin",
                               CNCT = "Cystatin C",
                               S100AT = "Calprotectin",
                               SAAT = "Serum amyloid A",
                               IL8 = "Interleukin-8",
                               IL10 = "Interleukin-10",
                               IFNG = "Interferon-g",
                               TNFA = "TNF-a",
                               TRP = "Tryptophan",
                               KYN = "Kynurenine",
                               HK = "3-Hydroxykynurenine",
                               KA = "Kynurenic acid",
                               XA = "Xanthurenic acid",
                               AA = "Anthranilic acid",
                               HAA = "3-Hydroxyanthranilic acid",
                               PIC = "Piconilic acid",
                               QA = "Quinolinic acid",
                               KTR = "KTr",
                               PAR = "PAr",
                               HKXAR = "HK:XA",
                               CRP = "C-reactive protein",
                               IL6 = "Interleukin-6"))

# strongest results 

change_results %>% 
  select(clock, term2, estimate, p.value) %>% 
  arrange(p.value) %>% 
  print(n = 50)

# proportion with p < Bonferroni threshold (0.0021)

change_results %>% 
  mutate(sig = if_else(p.value < 0.0021, "Yes", "No")) %>% 
  group_by(sig) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

# proportion with estimate > 0 

change_results %>% 
  mutate(sig = if_else(p.value < 0.0021, "Yes", "No")) %>% 
  mutate(positive = if_else(estimate > 0, "Yes", "No")) %>% 
  group_by(sig, positive) %>% 
  tally() %>% 
  mutate(percent = n/sum(n) * 100)

## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

change_results %>% 
  ggplot(aes(x = estimate, y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 4) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))

ggsave("Longitudinal forest.png", device = "png",
       height = 7, width = 10)
