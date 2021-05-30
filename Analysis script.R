library(tidyverse)
library(blandr)
library(performance)
library(lme4)
library(mgcv)
library(flextable)
library(psych)
library(ggpubr)


#### Prepare data ####

## read data

bl <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.BL.csv")

fu <- read.csv("S:\\MNHS-SCS-Medicine\\PM-Users\\Lachlan\\data\\LP.FUP2.csv")


## merge baseline and FU


bl2 <- bl %>% 
  filter(!study %in% c("BC","GAC")) %>% 
  filter(!case == 1) %>% 
  filter(!duplicated(h2_id))


full <- inner_join(bl2, fu, by = "h2_id", suffix = c("_bl","_fu"))



## time to follow-up

full <- full %>% mutate(time_fu = fup_age - age)


# concerning cases 

full %>% 
  filter(time_fu < 9) %>% 
  select(age, fup_age, h2_id)



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
  mutate(across(bioms, ~ log(.x + 0.1), .names = "logW_{.col}")) %>% 
  mutate(across(starts_with("logW"), winsorise)) %>% 
  mutate(across(contains("ageacc"), winsorise, .names = "W_{.col}"))




### Table 2. Change over time


# vector of all biomarker names

marker_names <- as.character(
  expression(neopt_d_bl,crp_g_bl,cnct_g_bl,s100at_g_bl,saat_g_bl,il6_msd_bl,
             il8_msd_bl,il10_msd_bl,ifng_msd_bl,tnfa_msd_bl,trp_d_bl,kyn_d_bl,
             hk_d_bl,ka_d_bl,xa_d_bl,aa_d_bl,haa_d_bl,pic_d_bl,qa_d_bl,KTr_bl,
             PAr_bl,HKXAr_bl,neopt_d_fu,crp_g_fu,cnct_g_fu,s100at_g_fu,
             saat_g_fu,il6_msd_fu,il8_msd_fu,il10_msd_fu,ifng_msd_fu,
             tnfa_msd_fu,trp_d_fu,kyn_d_fu,hk_d_fu,ka_d_fu,xa_d_fu,aa_d_fu,
             haa_d_fu,pic_d_fu,qa_d_fu,KTr_fu,PAr_fu,HKXAr_fu))


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



## Change scores as function of follow-up time





## CRP LOW VALUES

full %>% 
  filter(crp_g_bl > -1) %>% select(crp_g_bl) %>% 
  arrange(crp_g_bl)







#### spaghetti plots ####


## PhenoAge

hist(full$W_AgeAccelPheno_fu, breaks = 20)


# plot


plot_pheno <- full %>% 
  pivot_longer(cols = c("W_AgeAccelPheno_bl", "W_AgeAccelPheno_fu"), 
               names_to = "visit", values_to = "pheno", 
               names_prefix = "W_AgeAccelPheno_") %>%
  mutate(visit = if_else(visit == "bl", 0, 1)) %>% 
  ggplot(aes(x = as.numeric(visit), y = pheno, group = h2_id)) +
  geom_line() + scale_x_continuous(breaks = c(0,1),
                                   labels = c("BL","FU")) +
  geom_smooth(method = "lm", group = 1, colour = "red") +
  xlab("Visit") + ylim(c(-25,25)) + 
  ylab("PhenoAge AA") + theme_bw()

plot_pheno



## GrimAge

hist(full$W_AgeAccelGrim_fu, breaks = 20)


#plot 

plot_grim <- full %>% 
  pivot_longer(cols = c("W_AgeAccelGrim_bl", "W_AgeAccelGrim_fu"), 
               names_to = "visit", values_to = "grim", 
               names_prefix = "W_AgeAccelGrim_") %>%
  mutate(visit = if_else(visit == "bl", 0, 1)) %>% 
  ggplot(aes(x = as.numeric(visit), y = grim, group = h2_id)) +
  geom_line() + scale_x_continuous(breaks = c(0,1),
                                   labels = c("BL","FU")) +
  geom_smooth(method = "lm", group = 1, colour = "red") +
  xlab("Visit") + ylab("GrimAge AA") + theme_bw() +
  ylim(c(-15,15))


plot_grim


fig1 <- ggarrange(plot_pheno, plot_grim, ncol = 2)

fig1

sum(is.na(full$crp_g_bl))

## CRP

hist(full$crp_g_fu, breaks = 50)


# plot


plot_crp <- full %>% 
  pivot_longer(cols = c("logW_crp_g_bl", "logW_crp_g_fu"), 
               names_to = "visit", values_to = "crp", 
               names_prefix = "logW_crp_g_") %>%
  mutate(visit = if_else(visit == "bl", 0, 1)) %>% 
  ggplot(aes(x = visit, y = crp, group = h2_id)) +
  geom_line() + scale_x_continuous(breaks = c(0,1),
                                   labels = c("BL","FU")) +
  geom_smooth(method = "lm", group = 1, colour = "red") +
  xlab("Visit") + ylab("log(hsCRP)") + theme_bw() +
  ylim(c(-3,4))

plot_crp



## neopterin


hist(full$neopt_d_bl, breaks = 50)


# plot


plot_neopt <- full %>% 
  pivot_longer(cols = c("logW_neopt_d_bl", "logW_neopt_d_fu"), 
               names_to = "visit", values_to = "Neopterin", 
               names_prefix = "logW_neopt_d_") %>%
  mutate(visit = if_else(visit == "bl", 0, 1)) %>% 
  ggplot(aes(x = visit, y = Neopterin, group = h2_id)) +
  geom_line() + scale_x_continuous(breaks = c(0,1),
                                   labels = c("BL","FU")) +
  geom_smooth(method = "lm", group = 1, colour = "red") +
  xlab("Visit") + ylab("log(Neopterin)") + theme_bw() +
  ylim(c(1, 4))

plot_neopt


## Cystatin C


hist(full$cnct_g_fu, breaks = 50)


#plot

plot_cnct <- full %>% 
  pivot_longer(cols = c("logW_cnct_g_bl", "logW_cnct_g_fu"), 
               names_to = "visit", values_to = "CNCT", 
               names_prefix = "logW_cnct_g_") %>%
  mutate(visit = if_else(visit == "bl", 0, 1)) %>% 
  ggplot(aes(x = visit, y = CNCT, group = h2_id)) +
  geom_line() + scale_x_continuous(breaks = c(0,1),
                                   labels = c("BL","FU")) +
  geom_smooth(method = "lm", group = 1, colour = "red") +
  xlab("Visit") + ylab("log(Cystatin C)") + theme_bw() +
  ylim(c(-0.75, 1))

plot_cnct



## calprotectin

hist(full$s100at_g_fu, breaks = 50)


plot_s100 <- full %>% 
  pivot_longer(cols = c("logW_s100at_g_bl", "logW_s100at_g_fu"), 
               names_to = "visit", values_to = "s100", 
               names_prefix = "logW_s100at_g_") %>%
  mutate(visit = if_else(visit == "bl", 0, 1)) %>% 
  ggplot(aes(x = visit, y = s100, group = h2_id)) +
  geom_line() + scale_x_continuous(breaks = c(0,1),
                                   labels = c("BL","FU")) +
  geom_smooth(method = "lm", group = 1, colour = "red") +
  xlab("Visit") + ylab("log(Calprotectin)") + theme_bw() +
  ylim(c(-1.5, 3))

plot_s100


# make figure

fig2 <- ggarrange(plot_crp, plot_neopt, plot_cnct, plot_s100, ncol = 2,
                  nrow = 2)

ggsave("spaghetti 2.png", device = "png")


### ISSUE WITH CALPROTECTIN? ###

full %>% 
  pivot_longer(cols = c("s100at_g_bl", "s100at_g_fu"), 
               names_to = "visit", values_to = "s100", 
               names_prefix = "s100at_g_") %>%
  filter(s100 < 10) %>% 
  ggplot(aes(x = s100)) + geom_histogram() + 
  facet_wrap(~ visit) + xlab("Calprotectin") theme_bw()

ggsave("calprotectin histogram.png", device = "png")


# log

full %>% 
  pivot_longer(cols = c("logW_s100at_g_bl", "logW_s100at_g_fu"), 
               names_to = "visit", values_to = "s100", 
               names_prefix = "logW_s100at_g_") %>%
  filter(s100 < 10) %>% 
  ggplot(aes(x = s100)) + geom_histogram() + 
  facet_wrap(~ visit) + xlab("log(calprotectin)") + theme_bw()

ggsave("log calprotectin histogram.png", device = "png")