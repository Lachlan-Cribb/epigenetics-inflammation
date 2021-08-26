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

# add bl_suffix to all variables (except participant ID)

bl3 <- bl3 %>% 
  rename_at(vars(-h2_id), function(x) paste0(x,"_bl"))


## Follow-up

fu_new <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\LP.FUP2.csv")

fu_old <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.FUP2.csv")

fu_imp <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\new\\LP.FUP2_imp.csv")

# merge follow-up files

fu <- left_join(fu_new, fu_old, by = "h2_id", suffix = c("","_old"))

fu <- fu %>% select(-ends_with("_old"))

fu2 <- left_join(fu, fu_imp, by = "h2_id", suffix = c("", "_imp"))

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

full <- full %>% filter(!is.na(AgeAccelGrim_fu))

# Calculate time to follow-up

full <- full %>% mutate(time_fu = age_fu_correct - age_bl_correct)

hist(full$time_fu, breaks = 50) # minimum of 9 years 


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
  select(contains("grim"), contains("pheno"),
         contains("dunedin"), contains("zhang")) %>% 
  names()

full <- full %>% 
  mutate(across(all_of(epi_vars), winsorise,
                .names = "w{.col}"))

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
         wAA.DunedinPoAm_bl, wAA.DunedinPoAm_fu,
         wAA.Zhang.cont_bl, wAA.Zhang.cont_fu) %>% 
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

# winsorise inflammaging variables

full <- full %>% 
  mutate(across(contains("inflamm"), winsorise, .names = "w{.col}"))

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


### Change over time for table 2 ###
# We regress change in biomarker on
# follow-up time

# mean center time to follow-up

full$ctime_fu <- full$time_fu - mean(full$time_fu)

# CRP

ggplot(full, aes(x = age_fu_correct, y = wlog_crp_g_imp_fu2)) +
  geom_jitter() + geom_smooth()

m_crp <- lm(wlog_crp_g_imp_fu2 ~ 1 + ctime_fu + wlog_crp_g_imp_bl2, data = full)

sjPlot::plot_model(m_crp, type = "pred", terms = "ctime_fu")
cor(full$age_fu_correct, full$wlog_crp_g_imp_fu2)

# Neopterin 

ggplot(full, aes(x = age_fu_correct, y = wlog_neopt_d_imp_fu)) +
  geom_jitter() + geom_smooth()

exp(mean(full$wlog_neopt_d_imp_fu - full$wlog_neopt_d_imp_bl))

m_neo <- lm(wlog_neopt_d_imp_fu ~ ctime_fu + wlog_neopt_d_imp_bl, data = full)

tidy(m_neo, conf.int = T) %>% 
  mutate(percent = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100) %>% 
  select(term, percent, conf.low, conf.high)

cor(full$age_fu_correct, full$wlog_neopt_d_imp_fu)
# cystatin

full$change_cnct = full$wlog_cnct_g_imp_fu - full$wlog_cnct_g_imp_bl

m_cys <- lm(change_cnct ~ time_fu + wlog_cnct_g_imp_bl, data = full)

tidy(m_cys, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

# serum amyloid A

full$change_saat = full$wlog_saat_g_imp_fu - full$wlog_saat_g_imp_bl2

m_saat <- lm(wlog_saat_g_imp_fu ~ time_fu + wlog_saat_g_imp_bl2, data = full)

tidy(m_saat, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

# interleukin 6

full$change_il6 = full$wlog_il6_msd_imp_fu2 - full$wlog_il6_msd_imp_bl

m_il6 <- lm(wlog_il6_msd_imp_fu2 ~ time_fu + wlog_il6_msd_imp_bl, data = full)

tidy(m_il6, conf.int = T) %>% mutate(estimate = (exp(estimate)-1) * 100) %>% 
  mutate(conf.low = (exp(conf.low)-1) * 100, 
         conf.high = (exp(conf.high)-1) * 100)

ggplot(full, aes(x = age_bl_correct, y = wlog_crp_g_imp_bl2)) +
  geom_point() + geom_smooth()

#### Figure 1 ####
## Correlation matrix of inflammation/TK biomarkers and age accel variables
## at follow-up 

aa_vars <- full %>% 
  select(wAA.Zhang.cont_fu, wAA.DunedinPoAm_fu, wAgeAccelGrim_fu,
         wAgeAccelPheno_fu, wAA.Zhang.cont_bl, wAA.DunedinPoAm_fu,
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
                           WAA.ZHANG.CONT = "AgeAccelZhang",
                           WAA.DUNEDINPOAM = "AgeAccelDunedin",
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
       width = 12, height = 10)

#### Figure 3 ####

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


#### Table 3 ####

## confounders

full$sex_cde_bl <- as.factor(full$sex_cde_bl)

levels(full$sex_cde_bl) <- c("Male","Female")

summary(as.factor(full$cigst_cde_bl))

summary(as.factor(full$educlvl_ord_bl))

## standardize epigenetic ageing measures 

full <- full %>% 
  mutate(across(c("wAA.Zhang.cont_fu", "wAA.DunedinPoAm_fu", "wAgeAccelGrim_fu",
                  "wAgeAccelPheno_fu", "wAA.Zhang.cont_bl", "wAA.DunedinPoAm_fu",
                  "wAgeAccelGrim_bl", "wAgeAccelPheno_bl"), scale,
                .names = "z{.col}"))

## CRP models

# Pheno

m1 <- lm(zwAgeAccelPheno_fu ~ zwlog_crp_g_imp_fu2 + cob_cde_fu + sex_cde_fu + 
         bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m1, terms = "zwlog_crp_g_imp_fu2")

tidy(m1, conf.int = T) %>% 
  mutate(clock = "pheno") %>% 
  filter(str_detect(term, "zwlog")) -> m1_res

# grim

m2 <- lm(zwAgeAccelGrim_fu ~ zwlog_crp_g_imp_fu2 + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m2, terms = "zwlog_crp_g_imp_fu2")

tidy(m2, conf.int = T) %>% 
  mutate(clock = "grim") %>% 
  filter(str_detect(term, "zwlog")) -> m2_res

# Zhang

m3 <- lm(zwAA.Zhang.cont_fu ~ zwlog_crp_g_imp_fu2 + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m3, terms = "zwlog_crp_g_imp_fu2")

tidy(m3, conf.int = T) %>% 
  mutate(clock = "zhang") %>% 
  filter(str_detect(term, "zwlog")) -> m3_res

# Dunedin 

m4 <- lm(zwAA.DunedinPoAm_fu ~ zwlog_crp_g_imp_fu2 + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m4, terms = "zwlog_crp_g_imp_fu2")

tidy(m4, conf.int = T) %>% 
  mutate(clock = "dunedin") %>% 
  filter(str_detect(term, "zwlog")) -> m4_res

## Neopterin 

# Pheno

m5 <- lm(zwAgeAccelPheno_fu ~ zwlog_neopt_d_imp_fu + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m5, terms = "zwlog_neopt_d_imp_fu")

tidy(m5, conf.int = T) %>% 
  mutate(clock = "pheno") %>% 
  filter(str_detect(term, "zwlog")) -> m5_res

# grim

m6 <- lm(zwAgeAccelGrim_fu ~ zwlog_neopt_d_imp_fu + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m6, terms = "zwlog_neopt_d_imp_fu")

tidy(m6, conf.int = T) %>% 
  mutate(clock = "grim") %>% 
  filter(str_detect(term, "zwlog")) -> m6_res

# Zhang

m7 <- lm(zwAA.Zhang.cont_fu ~ zwlog_neopt_d_imp_fu + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m7, terms = "zwlog_neopt_d_imp_fu")

tidy(m7, conf.int = T) %>% 
  mutate(clock = "zhang") %>% 
  filter(str_detect(term, "zwlog")) -> m7_res

# Dunedin 

m8 <- lm(zwAA.DunedinPoAm_fu ~ zwlog_neopt_d_imp_fu + cob_cde_fu + sex_cde_fu + 
           bmi_rrto_fu + cigst_cde_fu + educlvl_ord_fu + age_fu_correct, 
         data = full)

crPlots(m8, terms = "zwlog_neopt_d_imp_fu")

tidy(m8, conf.int = T) %>% 
  mutate(clock = "dunedin") %>% 
  filter(str_detect(term, "zwlog")) -> m8_res






#### Forest plot ####

plot_data <- rbind(crp_res, cys_res, il6_res, saat_res)

plot_data <- plot_data %>% mutate(term = gsub("\\(|\\)","",term)) %>% 
  mutate(term = sub("scalelogW_", "", term)) %>% 
  mutate(term = gsub("_.*$", "", term)) %>% 
  mutate(term = str_to_upper(term)) %>% 
  rename(lower = conf.low, upper = conf.high)

plot_data %>% 
  ggplot(aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = lower,
                 xmax = upper)) +
  bayesplot_theme_set() +
  labs(x = "Estimate with 95% CI", y = "Inflammatory/kynurenine biomarker")


#### overall R squared ####

trim <- full %>% select(starts_with("logW"), W_AgeAccelGrim_fu, age_fu_correct, cob_cde_bl, sex_cde_bl, bmi_rrto, educlvl_ord) %>% 
  select(ends_with("fu"), ends_with("fu2"), W_AgeAccelGrim_fu, age_fu_correct, cob_cde_bl, sex_cde_bl, bmi_rrto, educlvl_ord)

m1 <- lm(W_AgeAccelGrim_fu ~ age_fu_correct + cob_cde_bl + sex_cde_bl + bmi_rrto + educlvl_ord, data = trim)
m2 <- lm(W_AgeAccelGrim_fu ~ ., data = trim)

summary(m1)

summary(m2)
