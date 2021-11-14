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
library(brms)

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

# Add in follow-up smoking variable

BLFUP2 <- read.csv("S:/MNHS-SCS-Medicine/PM-Users/Lachlan/data/BLFUP2.csv")

smokdat <- BLFUP2 %>% 
  select(h2_id, smokingstop)

full <- full %>% 
  left_join(smokdat, by = "h2_id")

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

full$sex_cde_fu <- factor(full$sex_cde_fu, levels = c("Male", "Female"))

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
  mutate(across(marker_names, ~ log2(.x), .names = "wlog_{.col}")) %>% 
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


## residuals of inflamm_sig

full$AAinflammaging_bl <- 
  residuals(lm(inflamm_sig_bl ~ age_bl_correct, data = full))

full$AAinflammaging_fu <- 
  residuals(lm(inflamm_sig_fu ~ age_fu_correct, data = full))

full %>% select(contains("AAinflamm")) %>% multi.hist()

# winsorise and scale to Z scores inflammaging variables

full <- full %>% 
  mutate(across(contains("AAinflamm"), winsorise, .names = "w{.col}")) %>% 
  mutate(across(contains("wAAinflamm"), scale, .names = "z{.col}"))


### Fold change in inflammatory biomarkers from baseline to follow-up ###

full <- full %>% 
  mutate(change_neopt = wlog_neopt_d_imp_fu - wlog_neopt_d_imp_bl,
         change_crp = wlog_crp_g_imp_fu2 - wlog_crp_g_imp_bl2,
         change_cnct = wlog_cnct_g_imp_fu - wlog_cnct_g_imp_bl,
         change_s100at = wlog_s100at_g_imp_fu - wlog_s100at_g_imp_bl,
         change_saat = wlog_saat_g_imp_fu - wlog_saat_g_imp_bl2,
         change_il6 = wlog_il6_msd_imp_fu2 - wlog_il6_msd_imp_bl,
         change_il8 = wlog_il8_msd_imp_fu - wlog_il8_msd_imp_bl,
         change_il10 = wlog_il10_msd_imp_fu2 - wlog_il10_msd_imp_bl2,
         change_ifng = wlog_ifng_msd_imp_fu - wlog_ifng_msd_imp_bl2,
         change_tnfa = wlog_tnfa_msd_imp_fu - wlog_tnfa_msd_imp_bl,
         change_trp = wlog_trp_d_imp_fu - wlog_trp_d_imp_bl,
         change_kyn = wlog_kyn_d_imp_fu - wlog_kyn_d_imp_bl,
         change_hk = wlog_hk_d_imp_fu - wlog_hk_d_imp_bl,
         change_ka = wlog_ka_d_imp_fu - wlog_ka_d_imp_bl,
         change_xa = wlog_xa_d_imp_fu - wlog_xa_d_imp_bl,
         change_aa = wlog_aa_d_imp_fu - wlog_aa_d_imp_bl,
         change_haa = wlog_haa_d_imp_fu - wlog_haa_d_imp_bl,
         change_pic = wlog_pic_d_imp_fu - wlog_pic_d_imp_bl,
         change_qa = wlog_qa_d_imp_fu - wlog_qa_d_imp_bl,
         change_KTr = wlog_KTr_imp_fu - wlog_KTr_imp_bl,
         change_PAr = wlog_PAr_imp_fu - wlog_PAr_imp_bl,
         change_HKXAr = wlog_HKXAr_imp_fu - wlog_HKXAr_imp_bl,
         change_inflam = scale(zwAAinflammaging_fu - zwAAinflammaging_bl))


## winsorise changes cores

full <- full %>% 
  mutate(across(starts_with("change_"), ~ winsorise(.)))

full %>% 
  select(starts_with("change_")) %>% 
  multi.hist(global = F)


## Change in AgeAccel

full <- full %>% 
  mutate(change_grim = wAgeAccelGrim_fu - wAgeAccelGrim_bl,
         change_pheno = wAgeAccelPheno_fu - wAgeAccelPheno_bl,
         change_zhang = wAgeAccelZhang_fu - wAgeAccelZhang_bl,
         change_dunedin = wAgeAccelDunedin_fu - wAgeAccelDunedin_bl)


## standardize change scores

full <- full %>% 
  mutate(across(starts_with("change_"), scale, .names = "z{.col}"))


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

full$sex_cde_bl <- factor(full$sex_cde_bl)

#### Figure 1 - Cross sectional association at BL #### 

# create vector including all baseline biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwAAinflammaging_bl) %>% 
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
                               ZWAAINFLAMMAGING = "AAinflammaging",
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

order <-  c("AAinflammaging", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")


bl_cs_results %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 1.5) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Coefficient, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("Baseline forest v2.png", device = "png",
       height = 5, width = 9)


#### Figure 2 - Cross sectional association at FU #### 

# create vector including all FU biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwAAinflammaging_fu) %>% 
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
  filter(!term %in% c("(Intercept)", "sex_cde_fuFemale", "cob_cde_fuNorthern Europe", "cob_cde_fuSouthern Europe")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWAAINFLAMMAGING = "AAinflammaging",
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

order <-  c("AAinflammaging", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

fu_cs_results %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 1.5) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))

ggsave("follow-up forest v2.png", device = "png",
       height = 5, width = 9)

#### Figure 3 - Prospective associations #### 

# create vector including all BL biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwAAinflammaging_bl) %>% 
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
                               ZWAAINFLAMMAGING = "AAinflammaging",
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

order <-  c("AAinflammaging", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")


prosp_results %>% 

  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 1.5) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))

ggsave("Prospective forest v2.png", device = "png",
       height = 5, width = 9)


#### Figure 4 - Change over time ####

# create vector including all BL biomarker variables names

bl_markers <- paste("zwlog_", marker_names, sep = "")

# drop IL8 & calprotectin due to almost perfect correlation with change 

xs <- full %>% 
  select(all_of(bl_markers), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl, -zwlog_il8_msd_imp_bl, -zwlog_s100at_g_imp_bl) %>% 
  names()

xs

# vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 


change_results <- 
  
  # create formulas with confounder variables 
  
  crossing(Var1 = xs, Var2 = ys) %>%
  mutate(Var3 = str_to_lower(Var2)) %>% 
  mutate(Var3 = sub("zwageaccel", "zchange_", Var3)) %>% 
  mutate(Var3 = sub("_fu", "", Var3)) %>% 
  mutate(Var4 = sub("zwlog_", "zchange", Var1)) %>% 
  mutate(Var4 = sub("_.*", "", .$Var4)) %>% 
  mutate(Var4 = sub("change", "change_", Var4)) %>% 
  mutate(Var4 = sub("zw", "zchange_", Var4)) %>% 
  mutate(Var4 = sub("change_inflamm", "change_inflam", Var4)) %>% 
  mutate(frm = str_c(Var3, Var1, sep = " ~ sex_cde_fu + cob_cde_fu + ")) %>% 
  mutate(Var2 = sub("_fu", "_bl", Var2)) %>% 
  mutate(frm = str_c(frm, Var2, sep = " + ")) %>% 
  mutate(frm = str_c(frm, Var4, sep = " + ")) %>% 
  select(-Var1, -Var2) %>% 
  
  
  # regressions 
  
  mutate(models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-frm) %>% 
  
  # prepare data for plotting 
  
  rename('clock' = Var3) %>% 
  filter(str_detect(term, "change_")) %>% 
  mutate(term = sub("change_", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zchange_dunedin = "AgeAccelDunedin",
                               zchange_zhang = "AgeAccelZhang",
                               zchange_grim = "AgeAccelGrim",
                               zchange_pheno = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZINFLAM = "Inflammaging signature",
                               ZNEOPT = "Neopterin",
                               ZCNCT = "Cystatin C",
                               ZS100AT = "Calprotectin",
                               ZSAAT = "Serum amyloid A",
                               ZIL8 = "Interleukin-8",
                               ZIL10 = "Interleukin-10",
                               ZIFNG = "Interferon-g",
                               ZTNFA = "TNF-a",
                               ZTRP = "Tryptophan",
                               ZKYN = "Kynurenine",
                               ZHK = "3-Hydroxykynurenine",
                               ZKA = "Kynurenic acid",
                               ZXA = "Xanthurenic acid",
                               ZAA = "Anthranilic acid",
                               ZHAA = "3-Hydroxyanthranilic acid",
                               ZPIC = "Piconilic acid",
                               ZQA = "Quinolinic acid",
                               ZKTR = "KTr",
                               ZPAR = "PAr",
                               ZHKXAR = "HK:XA",
                               ZCRP = "C-reactive protein",
                               ZIL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

change_results %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 2) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))

ggsave("Longitudinal forest v2.png", device = "png",
       height = 5, width = 9)

#### Comparing estimates ####

## Baseline ##

# strongest results 

bl_cs_results %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high) %>% 
  arrange(p.value) %>% 
  print(n = 50)

# inflammaging 

bl_cs_results %>% 
  filter(term2 == "AAinflammaging") %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high)


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
  select(clock, term2, estimate, p.value, conf.low, conf.high) %>% 
  arrange(p.value) %>% 
  print(n = 100)

# inflammaging 

fu_cs_results %>% 
  filter(term2 == "AAinflammaging") %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high)

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
  filter(estimate < 0)
  

## comparing baseline and follow-up 

## comparison ##

estimates <- tibble(
  term = bl_cs_results$term2,
  clock = bl_cs_results$clock,
  bl = bl_cs_results$estimate,
  fu = fu_cs_results$estimate,
)

# largest differences from baseline to FU

estimates %>% 
  mutate(change_est = fu - bl) %>% 
  arrange(desc(abs(change_est)))

# Baseline and follow-up correlation

cor(x = bl_cs_results$estimate, y = fu_cs_results$estimate)

t.test(x = bl_cs_results$estimate, y = fu_cs_results$estimate, paired = T)

## propsective ##

# strongest results 

prosp_results %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high) %>% 
  arrange(p.value)

# inflammaging 

prosp_results %>% 
  filter(term2 == "AAinflammaging") %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high) %>% 
  arrange(p.value)


## change model ##

# strongest results 

change_results %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high) %>% 
  arrange(p.value) %>% 
  print(n = 100)

# inflammaging

change_results %>% 
  filter(term2 == "Inflammaging signature") %>% 
  select(clock, term2, estimate, p.value, conf.low, conf.high) %>% 
  arrange(p.value)

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


#### Supplementary figure 4 - correlation of coefficients ####

fu_cs_cor <- fu_cs_results %>% 
  pivot_wider(id_cols = term, names_from = clock, values_from = estimate) %>% 
  select(-term)

cors_spearman <- function(df) { 
  M <- Hmisc::rcorr(as.matrix(df), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))
}


plot_data <- cors_spearman(fu_cs_cor) %>%
  map(~rownames_to_column(.x, var="measure1")) %>%
  map(~pivot_longer(.x, -measure1, "measure2")) %>%
  bind_rows(.id = "id") %>%
  pivot_wider(names_from = id, values_from = value)


plot_data %>%
  ggplot(aes(factor(measure1), factor(measure2), 
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

ggsave("coefficient correlation.png", device = "png", dpi = 450,
       width = 6, height = 4)




#### Variance explained (Bayesian horseshoe regression) ####

## Grim

# formula

form1 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelGrim_bl ~", .)

# fit model

bl_grim <- 
  brm(form1,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "grim_bl_hs",
      backend = "cmdstanr")

# model results

summary(bl_grim)

bl_grim_R2 <- loo_R2(bl_grim) # LOOCV R squared 

bl_grim_R2 <- bl_grim_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "GrimAge", visit = "BL")

## Pheno

# formula

form2 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelPheno_bl ~", .)

# fit model

bl_pheno <- 
  brm(form2,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "pheno_bl_hs",
      backend = "cmdstanr")

# model results

summary(bl_pheno)

bl_pheno_R2 <- loo_R2(bl_pheno) # LOOCV R squared 

bl_pheno_R2 <- bl_pheno_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "PhenoAge", visit = "BL")

## Zhang

# formula

form3 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelZhang_bl ~", .)

# fit model

bl_zhang <- 
  brm(form3,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "zhang_bl_hs",
      backend = "cmdstanr")

# model results

summary(bl_zhang)

bl_zhang_R2 <- loo_R2(bl_zhang) # LOOCV R squared 

bl_zhang_R2 <- bl_zhang_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "Zhang", visit = "BL")

## Dunedin

# formula

form4 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelDunedin_bl ~", .)

# fit model

bl_dunedin <- 
  brm(form4,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "dunedin_bl_hs",
      backend = "cmdstanr")

# model results

summary(bl_dunedin)

bl_dunedin_R2 <- loo_R2(bl_dunedin) # LOOCV R squared 

bl_dunedin_R2 <- bl_dunedin_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "Dunedin", visit = "BL")


#### follow-up variance explained ####

## Grim

# formula

form1 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelGrim_fu ~", .)

# fit model

fu_grim <- 
  brm(form1,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.999, max_treedepth = 15),
      file = "grim_fu_hs",
      backend = "cmdstanr")

# model results

summary(fu_grim)

fu_grim_R2 <- loo_R2(fu_grim) # LOOCV R squared 

fu_grim_R2 <- fu_grim_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "GrimAge", visit = "FU")

## Pheno

# formula

form2 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelPheno_fu ~", .)

# fit model

fu_pheno <- 
  brm(form2,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "pheno_fu_hs",
      backend = "cmdstanr")

# model results

summary(fu_pheno)

fu_pheno_R2 <- loo_R2(fu_pheno) # LOOCV R squared 

fu_pheno_R2 <- fu_pheno_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "PhenoAge", visit = "FU")


## Zhang

# formula

form3 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelZhang_fu ~", .)

# fit model

fu_zhang <- 
  brm(form3,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "zhang_fu_hs",
      backend = "cmdstanr")

# model results

summary(fu_zhang)

fu_zhang_R2 <- loo_R2(fu_zhang) # LOOCV R squared 

fu_zhang_R2 <- fu_zhang_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "Zhang", visit = "FU")

## Dunedin

# formula

form4 <- full %>% 
  select(starts_with("zwlog_"), zwAAinflammaging_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  names() %>% 
  paste(., collapse = " + ") %>% 
  paste("zwAgeAccelDunedin_fu ~", .)

# fit model

fu_dunedin <- 
  brm(form4,
      data = full,
      prior = c(prior(horseshoe(df = 3, par_ratio = 0.2)),
                prior(normal(0,0.5), class = "Intercept")),
      iter = 4000,
      warmup = 1000,
      chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      file = "dunedin_fu_hs",
      backend = "cmdstanr")

# model results

summary(fu_dunedin)

fu_dunedin_R2 <- loo_R2(fu_dunedin) # LOOCV R squared 

fu_dunedin_R2 <- fu_dunedin_R2 %>% 
  as_tibble() %>% 
  mutate(outcome = "Dunedin", visit = "FU")


#### R squared plot ####

R2_dat <- bind_rows(
  bl_grim_R2,
  bl_pheno_R2,
  bl_zhang_R2,
  bl_dunedin_R2,
  fu_grim_R2,
  fu_pheno_R2,
  fu_zhang_R2,
  fu_dunedin_R2
)


R2_dat %>% 
  mutate(outcome = fct_recode(outcome,
    DunedinAgeAccel = 'Dunedin',
    GrimAgeAccel = 'GrimAge',
    PhenoAgeAccel = 'PhenoAge',
    ZhangAgeAccel = 'Zhang')) %>% 
  mutate(Q2.5 = if_else(Q2.5 < 0, 0, Q2.5)) %>% 
  mutate(visit = if_else(visit == "BL", "Baseline", "Follow-up")) %>% 
  ggplot(aes(x = outcome, y = Estimate, fill = visit)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.75) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5),
                position = position_dodge(0.75),
                width = 0.2,
                alpha = 0.6) +
  ggthemes::theme_stata() +
  ylim(c(0,0.15)) +
  labs(x = "\nOutcome variable", y = "Bayes R-squared, 95% CI")

ggsave("Bayes R squared.png", device = "png",
       width = 6, height = 6)


#### Supplementary Figure 1 -correlation matrix ####
## Correlation matrix of inflammation/TK biomarkers and age accel variables
## at follow-up 

aa_vars <- full %>% 
  select(wAgeAccelZhang_fu, wAgeAccelDunedin_fu, wAgeAccelGrim_fu,
         wAgeAccelPheno_fu, wAgeAccelDunedin_fu) %>% 
  names()

cor_dat <- full %>% 
  select(starts_with("wlog_"), all_of(aa_vars), zwAAinflammaging_fu) %>% 
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
                                                     ZWAAINFLAMMAGING = "AAinflammaging",
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
            "AgeAccelPheno", "AAinflammaging", "Neopterin", "C-reactive protein",
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
       width = 14, height = 8)

#### Supplementary Figure 2 - age correlations ####

# create text to annotate plot 

ann_text <- data.frame(
  label = c("r = 0.22", "r = 0.85", "r = 0.74", "r = 0.41", "r = 0.50"),
  variable = c("DunedinPoAm", "GrimAge","PhenoAge","Zhang", "Inflammaging signature"),
  x = rep(-2, 5),
  y = rep(88, 5)
)

full %>% 
  mutate(across(c("score.Zhang.cont_fu", "DunedinPoAm_fu", 
                "DNAmGrimAge_fu", "DNAmPhenoAge_fu", "inflamm_sig_fu"), 
                winsorise, .names = "w{.col}")) %>% 
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
                                  winflamm_sig_fu = "Inflammaging signature")) %>%
  ggplot(aes(x = ageaccel, y = age_fu_correct)) +
  geom_point() + geom_smooth(method = "lm", colour = "darkblue") + 
  facet_wrap(~ factor(variable,
                      levels = c("Zhang", "DunedinPoAm",
                                 "GrimAge", "PhenoAge",
                                 "Inflammaging signature"))) +
  labs(x = "Epigenetic ageing/inflammaging biomarker (z score)", y = "Chronological age") +
  theme_stata() +
  geom_text(
    data    = ann_text,
    mapping = aes(x = x, y = y, label = label)
  )

ggsave("clock age correlation.png", device = "png",
       width = 7, height = 6)

## grim age only

full <- full %>% 
  mutate(across(c("score.Zhang.cont_fu", "DunedinPoAm_fu", 
                  "DNAmGrimAge_fu", "DNAmPhenoAge_fu", "inflamm_sig_fu"), 
                winsorise, .names = "w{.col}"))

ann_text <- data.frame(
  label = "r = 0.85",
  variable = "GrimAge",
  x = 52,
  y = 82
)

full %>% 
  ggplot(aes(x = wDNAmGrimAge_fu, y = age_fu_correct)) +
  geom_point() + geom_smooth(method = "lm", colour = "darkblue") +
  labs(x = "GrimAge", y = "Chronological age") +
  theme_stata() +
  geom_text(
    data    = ann_text,
    mapping = aes(x = x, y = y, label = label),
    size = 6
  )

ggsave("grim age correlation.png", device = "png",
       width = 4, height = 5)


## GrimAge residuals


m1 <- lm(wDNAmGrimAge_fu ~ age_fu_correct, data = full)

full %>% 
  mutate(resids = residuals(m1)) %>% 
  mutate(fitted = fitted(m1)) %>% 
  select(h2_id, resids, wDNAmGrimAge_fu, age_fu_correct, fitted) %>% 
  print(n = 30)

full %>% 
  ggplot(aes(x = age_fu_correct, y = wDNAmGrimAge_fu)) +
  geom_point() +
  geom_point(data = full %>% filter(h2_id == 11059),
             pch = 21,
             size = 5,
             colour = "red",
             stroke = 2) + 
  geom_point(data = full %>% filter(h2_id == 11402),
             pch = 21,
             size = 5,
             stroke = 2,
             colour = "green") +
  geom_smooth(method = "lm", colour = "darkblue") +
  labs(x = "Chronological age", y = "GrimAge") +
  theme_stata() +
  geom_segment(aes(x = 66.1, y= 70.3, yend = 62.0, xend = 66.1),
               colour = "red", size = 1) +
  geom_segment(aes(x = 70.3, y= 57.7, yend = 65.6, xend = 70.3),
               colour = "green", size = 1)

ggsave("grim age residual.png", device = "png",
       width = 5, height = 5)

#### Supplementary Figures 6-9 - Sensitivity confounders ####

hist(full$seifa_10_imp_bl)

hist(full$seifa_10_f2_imp_fu)

hist(full$bmi_rrto_imp_bl)

hist(full$bmi_rrto_f2_imp_fu)

full$cigst_cde_bl <- factor(full$cigst_cde_bl)

full$smokingstop <- factor(full$smokingstop)

hist(full$alc_mrat_imp_bl)

hist(full$alc_mrat_f2_imp_fu)

#### Figure S5 - baseline cross sectional #### 

# create vector including all baseline biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwAAinflammaging_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_bl", "zwAgeAccelGrim_bl", "zwAgeAccelZhang_bl", "zwAgeAccelDunedin_bl")

# regression for each combination of biomarker and AgeAccel variable 


bl_cs_results_s5 <- 
  
  # every combination of biomarker and AgeAccel
  
  crossing(Var1 = xs, Var2 = ys) %>%
  
  # perform linear regression adjusted for age and sex and save results
  
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_bl + cob_cde_bl + seifa_10_imp_bl + bmi_rrto_imp_bl + cigst_cde_bl + alc_mrat_imp_bl + "),  
         models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  
  # Drop unneeded parameters and rename variables for plotting 
  
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(!term %in% c("(Intercept)", "sex_cde_blMale", "cob_cde_blNorthern Europe", "cob_cde_blSouthern Europe", "seifa_10_imp_bl",
                      "bmi_rrto_imp_bl", "cigst_cde_bl1", "cigst_cde_bl2", "alc_mrat_imp_bl")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_bl = "AgeAccelDunedin",
                               zwAgeAccelZhang_bl = "AgeAccelZhang",
                               zwAgeAccelGrim_bl = "AgeAccelGrim",
                               zwAgeAccelPheno_bl = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWAAINFLAMMAGING = "AAinflammaging",
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

order <-  c("AAinflammaging", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")


bl_cs_results_s5 %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 1.5) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Coefficient, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("Baseline sensitivty fig S5.png", device = "png",
       height = 5, width = 9)

#### Figure S6 - sensitivity cross sectional at FU #### 

# create vector including all FU biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwAAinflammaging_fu) %>% 
  select(ends_with("fu"), ends_with("fu2")) %>% 
  select(-zwlog_cysta_d_imp_fu) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 

fu_cs_results_s6 <- crossing(Var1 = xs, Var2 = ys) %>%
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_bl + cob_cde_bl + seifa_10_f2_imp_fu + bmi_rrto_f2_imp_fu + smokingstop + alc_mrat_f2_imp_fu + "),  
         models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-Var1, -frm) %>% 
  rename('clock' = Var2) %>% 
  filter(!term %in% c("(Intercept)", "sex_cde_blMale", "cob_cde_blNorthern Europe", "cob_cde_blSouthern Europe","seifa_10_f2_imp_fu", "bmi_rrto_f2_imp_fu", "smokingstopStill/started smoking", "smokingstopStopped smoking", "alc_mrat_f2_imp_fu")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWAAINFLAMMAGING = "AAinflammaging",
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

order <-  c("AAinflammaging", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

fu_cs_results_s6 %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 1.5) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))


ggsave("follow-up sensitivity fig s6.png", device = "png",
       height = 5, width = 9)

#### Figure S7 - sensitivity prospective #### 

# create vector including all BL biomarker variables names

xs <- paste("zwlog_", marker_names, sep = "")

xs <- full %>% 
  select(all_of(xs), zwAAinflammaging_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl) %>% 
  names()

# create vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 

prosp_results_s7 <- 
  
  # all combinations of biomarker and AgeAccel vars
  
  crossing(Var1 = xs, Var2 = ys) %>%
  
  # create formula for each regression model 
  
  mutate(frm = str_c(Var2, Var1, sep = " ~ sex_cde_bl + cob_cde_bl + seifa_10_imp_bl + bmi_rrto_imp_bl + cigst_cde_bl + alc_mrat_imp_bl + ")) %>% 
  
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
  filter(str_detect(term, "_bl"), !term %in% c("(Intercept)", "sex_cde_blMale", "cob_cde_blNorthern Europe", "cob_cde_blSouthern Europe", "seifa_10_imp_bl","bmi_rrto_imp_bl", 
                                               "cigst_cde_bl1", "cigst_cde_bl2", "alc_mrat_imp_bl", "zwAgeAccelDunedin_bl", "zwAgeAccelGrim_bl",
                                               "zwAgeAccelPheno_bl", "zwAgeAccelZhang_bl")) %>% 
  mutate(term = str_replace(term, "zwlog_", "")) %>% 
  mutate(term = sub("_.*", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zwAgeAccelDunedin_fu = "AgeAccelDunedin",
                               zwAgeAccelZhang_fu = "AgeAccelZhang",
                               zwAgeAccelGrim_fu = "AgeAccelGrim",
                               zwAgeAccelPheno_fu = "AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZWAAINFLAMMAGING = "AAinflammaging",
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

order <-  c("AAinflammaging", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

prosp_results_s7 %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low,
                      xmax = conf.high),
                  fatten = 1.5) +
  facet_wrap(~ clock, nrow = 1) +
  theme_stata() +
  theme(axis.text.y = element_text(angle = 0)) +
  labs(x = "Estimate, 95% CI", y = "") +
  xlim(c(-0.4, 0.4)) +
  theme(panel.spacing.x = unit(1, "lines"))

ggsave("Prospective sensitivity fig s7.png", device = "png",
       height = 5, width = 9)


#### Figure S8 - Change over time sensitivity ####

# create vector including all BL biomarker variables names

bl_markers <- paste("zwlog_", marker_names, sep = "")

# drop IL8 & calprotectin due to almost perfect correlation with change 

xs <- full %>% 
  select(all_of(bl_markers), zwinflamm_sig_bl) %>% 
  select(ends_with("bl"), ends_with("bl2")) %>% 
  select(-zwlog_cysta_d_imp_bl, -zwlog_il8_msd_imp_bl, -zwlog_s100at_g_imp_bl) %>% 
  names()

xs

# vector of all epigenetic ageing variables 

ys <- c("zwAgeAccelPheno_fu", "zwAgeAccelGrim_fu", "zwAgeAccelZhang_fu", "zwAgeAccelDunedin_fu")

# regression for each combination of biomarker and epigenetic ageing variable 


change_results <- 
  
  # create formulas with confounder variables 
  
  crossing(Var1 = xs, Var2 = ys) %>%
  mutate(Var3 = str_to_lower(Var2)) %>% 
  mutate(Var3 = sub("zwageaccel", "zchange_", Var3)) %>% 
  mutate(Var3 = sub("_fu", "", Var3)) %>% 
  mutate(Var4 = sub("zwlog_", "zchange", Var1)) %>% 
  mutate(Var4 = sub("_.*", "", .$Var4)) %>% 
  mutate(Var4 = sub("change", "change_", Var4)) %>% 
  mutate(Var4 = sub("zw", "zchange_", Var4)) %>% 
  mutate(Var4 = sub("change_inflamm", "change_inflam", Var4)) %>% 
  mutate(frm = str_c(Var3, Var1, sep = " ~ sex_cde_bl + cob_cde_bl + seifa_10_imp_bl + bmi_rrto_imp_bl + cigst_cde_bl + alc_mrat_imp_bl + ")) %>% 
  mutate(Var2 = sub("_fu", "_bl", Var2)) %>% 
  mutate(frm = str_c(frm, Var2, sep = " + ")) %>% 
  mutate(frm = str_c(frm, Var4, sep = " + ")) %>% 
  select(-Var1, -Var2) %>% 
  select(frm)
  
  
  # regressions 
  
  mutate(models = map(frm, 
                      ~tidy(lm(as.formula(.x), data=full), conf.int=T))) %>% 
  unnest(cols = c(models)) %>% 
  select(-frm) %>% 
  
  # prepare data for plotting 
  
  rename('clock' = Var3) %>% 
  filter(str_detect(term, "change_")) %>% 
  mutate(term = sub("change_", "", .$term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  mutate(clock = dplyr::recode(clock,
                               zchange_dunedin = "Change AgeAccelDunedin",
                               zchange_zhang = "Change AgeAccelZhang",
                               zchange_grim = "Change AgeAccelGrim",
                               zchange_pheno = "Change AgeAccelPheno")) %>% 
  mutate(term2 = dplyr::recode(term,
                               ZINFLAM = "Inflammaging signature",
                               ZNEOPT = "Neopterin",
                               ZCNCT = "Cystatin C",
                               ZS100AT = "Calprotectin",
                               ZSAAT = "Serum amyloid A",
                               ZIL8 = "Interleukin-8",
                               ZIL10 = "Interleukin-10",
                               ZIFNG = "Interferon-g",
                               ZTNFA = "TNF-a",
                               ZTRP = "Tryptophan",
                               ZKYN = "Kynurenine",
                               ZHK = "3-Hydroxykynurenine",
                               ZKA = "Kynurenic acid",
                               ZXA = "Xanthurenic acid",
                               ZAA = "Anthranilic acid",
                               ZHAA = "3-Hydroxyanthranilic acid",
                               ZPIC = "Piconilic acid",
                               ZQA = "Quinolinic acid",
                               ZKTR = "KTr",
                               ZPAR = "PAr",
                               ZHKXAR = "HK:XA",
                               ZCRP = "C-reactive protein",
                               ZIL6 = "Interleukin-6"))


## plot ## 

order <-  c("Inflammaging signature", "Neopterin", "C-reactive protein",
            "Calprotectin", "Cystatin C", "Serum amyloid A",
            "Interferon-g","TNF-a", "Interleukin-6", 
            "Interleukin-8", "Interleukin-10","Kynurenine",
            "Tryptophan", "3-Hydroxykynurenine","Kynurenic acid",
            "Xanthurenic acid", "Anthranilic acid", "3-Hydroxyanthranilic acid",
            "Piconilic acid", "Quinolinic acid", "KTr", "PAr", "HK:XA")

change_results %>% 
  ggplot(aes(x = estimate, 
             y = factor(term2, levels = rev(order)),
             colour = I(ifelse(abs(estimate) > 0.1, "mediumblue", "grey20")))) +
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

ggsave("Longitudinal forest v2.png", device = "png",
       height = 5, width = 9)


#### Comparison between sensitivity and unadjusted results ####

sens_res <- bl_cs_results %>% 
  rename(estimate_bl = estimate) %>% 
  mutate(estimate_bl_s = bl_cs_results_s5$estimate,
         estimate_fu = fu_cs_results$estimate,
         estimate_fu_s = fu_cs_results_s6$estimate,
         estimate_prosp = prosp_results$estimate,
         estimate_prosp_s = prosp_results_s7$estimate)

# baseline #

sens_res %>% 
  group_by(term2) %>% 
  summarise(mean_bl = mean(estimate_bl),
            mean_bl_s = mean(estimate_bl_s),
            dif = mean(estimate_bl_s - estimate_bl),
            sd_dif = sd(estimate_bl_s - estimate_bl)) %>% 
  mutate(percent_change = ((mean_bl_s - mean_bl) / mean_bl) * 100) %>% 
  arrange(desc(abs(dif)))

# follow-up #

sens_res %>% 
  group_by(term2) %>% 
  summarise(mean_fu = mean(estimate_fu),
            mean_fu_s = mean(estimate_fu_s),
            dif = mean(estimate_fu_s - estimate_fu),
            sd_dif = sd(estimate_fu_s - estimate_fu)) %>% 
  mutate(percent_change = ((mean_fu_s - mean_fu) / mean_fu) * 100) %>% 
  arrange(desc(abs(dif)))

# Prospective #

sens_res %>%
  group_by(term2) %>% 
  summarise(mean_prosp = mean(estimate_prosp),
            mean_prosp_s = mean(estimate_prosp_s),
            dif = mean(estimate_prosp_s - estimate_prosp),
            sd_dif = sd(estimate_prosp_s - estimate_prosp)) %>% 
  mutate(percent_change = ((mean_prosp_s - mean_prosp) / mean_prosp) * 100) %>% 
  arrange(desc(abs(dif)))

#### Survival ####


library(survival)

# survival function 

survive <- Surv(time = full$age_fu_correct, time2 = full$AGE_DTH_fu, event = full$death_COD_imp_fu)


# AgeAccel model


m1 <- coxph(survive ~ zwAgeAccelGrim_fu + zwAgeAccelPheno_fu + zwAgeAccelZhang_fu + zwAgeAccelDunedin_fu + 
              age_fu_correct + sex_cde_fu, 
               x = T, y = T, data = full)


# AgeAccel and inflammaging model 

m2 <- coxph(survive ~ zwAgeAccelGrim_fu + zwAgeAccelPheno_fu + zwAgeAccelZhang_fu + zwAgeAccelDunedin_fu + 
              age_fu_correct + sex_cde_fu + zwinflamm_sig_fu, 
            x = T, y = T, data = full)


anova(m1, m2)


# R squared comparison 

m1 <- rms::cph(survive ~ zwAgeAccelGrim_fu + zwAgeAccelPheno_fu + zwAgeAccelZhang_fu + zwAgeAccelDunedin_fu + 
              age_fu_correct + sex_cde_fu, 
            x = T, y = T, data = full)

m2 <- rms::cph(survive ~ zwAgeAccelGrim_fu + zwAgeAccelPheno_fu + zwAgeAccelZhang_fu + zwAgeAccelDunedin_fu + 
              age_fu_correct + sex_cde_fu + zwinflamm_sig_fu, 
            x = T, y = T, data = full)


print(m1)

print(m2)

full$bmi_rrto_f2_imp_fu

summary(full$bmi_rrto_fu)

full %>% 
  summarise(mean = mean(bmi_rrto_fu),
            sd = sd(bmi_rrto_fu))

summary(factor(full$smokingstop_fu))
