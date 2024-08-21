rm(list=ls())

#The data is preprocessed:
#the fiations less than 80 ms are removed, 
#filler and practice trials are removed, 
#failed calibration trials removed
#split into regions
#merged with frequency data
shhh <- suppressPackageStartupMessages # It's a library, so shhh!

shhh(library( mgcv ))
shhh(library(dplyr))
shhh(library(ggplot2))
shhh(library(lme4))
shhh(library(tidymv))
shhh(library(gamlss))
shhh(library(gsubfn))
shhh(library(lmerTest))
shhh(library(tidyverse))
shhh(library(boot))
shhh(library(rsample))
shhh(library(plotrix))
shhh(library(ggrepel))
shhh(library(mgcv))
shhh(library(brms))
shhh(library(bayesplot))
shhh(library(tidyr))
shhh(library(car))
shhh(library(HDInterval))
shhh(library(gridExtra))
shhh(library(posterior))

shhh(library(coda))
shhh(library(cmdstanr))
shhh(library(rstan))
shhh(library(rstantools))

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_bw())
options(digits=4)
options(scipen=999)
set.seed(444)

library(readxl)
library(plyr)
################################################################################
#IMPORT DATA BY REGIONS#
################################################################################

R2 <- read_excel("R2.xlsx")
R3 <- read_excel("R3.xlsx")
R4 <- read_excel("R4.xlsx")
R5 <- read_excel("R5.xlsx")

#################Descriptive stats###############################################

# REGION 2###

FF <- ddply(R2, c("gender_match"), summarise,
            N    = length(IA_FIRST_FIXATION_DURATION),
            mean = mean(IA_FIRST_FIXATION_DURATION, na.rm = T),
            sd   = sd(IA_FIRST_FIXATION_DURATION, na.rm = T))
FF$DV<-"FF"


GD <- ddply(R2, c("gender_match"), summarise,
            N    = length(IA_FIRST_RUN_DWELL_TIME),
            mean = mean(IA_FIRST_RUN_DWELL_TIME, na.rm = T),
            sd   = sd(IA_FIRST_RUN_DWELL_TIME, na.rm = T))
GD$DV<-"GD"



SF <- ddply(R2, c("gender_match"), summarise,
            N    = length(SFD),
            mean = mean(SFD, na.rm = T),
            sd   = sd(SFD, na.rm = T))

SF$DV<-"SF"


TT <- ddply(R2, c("gender_match"), summarise,
            N    = length(IA_DWELL_TIME),
            mean = mean(IA_DWELL_TIME, na.rm = T),
            sd   = sd(IA_DWELL_TIME, na.rm = T))

TT$DV<-"TT"


#descr table for durations
rbind(FF, SF, GD, TT)->descr_dur_R2
descr_dur_R2$SE<-descr_dur_R2$sd/sqrt(descr_dur_R2$N)
#create region column
descr_dur_R2$Region<-"R2"



ROf <- ddply(R2, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_OUT_FULL),
             mean = mean(IA_REGRESSION_OUT_FULL, na.rm = T),
             sd   = sd(IA_REGRESSION_OUT_FULL, na.rm = T))

ROf$DV<-"ROf"


Rin <- ddply(R2, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_IN),
             mean = mean(IA_REGRESSION_IN, na.rm = T),
             sd   = sd(IA_REGRESSION_IN, na.rm = T))

Rin$DV<-"Rin"

View(R2)

#descr table for durations
rbind(skip, RO, ROf, Rin)->descr_prob_R2
descr_prob_R2$SE<-descr_prob_R2$sd/sqrt(descr_prob_R2$N)
descr_prob_R2$Region<-"R2"

############MODELING#############################################################


#Frequency was taken from http://stimul.cognitivestudies.ru/ru_stimul/ for lemmas(for R2s and adverbs, 
#conjuctions, particles adj in Sg Nom Masc) and wordforms *everything else)
#and from https://skell.sketchengine.eu/#home?lang=ru for phrases


###########predictors in full model:############################################
gender_match*type+ target_gender+
  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX+
  (1 + gender_match*type| DATA_FILE) + (1 | item)

#################################################################################

library(lme4)

#FFD
ffd.mod <- lmer(log(IA_FIRST_FIXATION_DURATION) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 + gender_match| DATA_FILE) + (1 | item),
                data = R2, 
                control = lmerControl(optimizer = "bobyqa"))

summary(ffd.mod)

#SFD
sfd.mod <- lmer(log(SFD) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 + gender_match| DATA_FILE) + (1 | item),
                data = R2, 
                control = lmerControl(optimizer = "bobyqa"))

#GD
gd.mod <- lmer(log(IA_FIRST_RUN_DWELL_TIME) ~ gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 | DATA_FILE) + (1 | item),
               data = R2, 
               control = lmerControl(optimizer = "bobyqa"))


summary(gd.mod)

#TT
R2$IA_DWELL_TIME<-as.numeric(R2$IA_DWELL_TIME)
R2$lg_frequency<-as.numeric(R2$lg_frequency)
data.tt <- R2[R2$IA_DWELL_TIME != 0,]
tt.mod <- lmer(log(IA_DWELL_TIME) ~gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 + gender_match| DATA_FILE) + (1 | item),
               data = data.tt, 
               control = lmerControl(optimizer = "bobyqa"))

summary(tt.mod)

#Rout

ROf.mod <- glmer(IA_REGRESSION_OUT_FULL ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 + gender_match| DATA_FILE) + (1 | item),
                 data = R2, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))

#Rin
Rin.mod <- glmer(IA_REGRESSION_IN ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 + gender_match| DATA_FILE) + (1 | item),
                 data = R2, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))

summary(Rin.mod)

#CREATE A TABLE WITH ALL MODELS

library(sjPlot)

tab_model(ffd.mod, sfd.mod, gd.mod, tt.mod,ROf.mod, Rin.mod,
          string.se = "SE",digits = 3, show.ci = FALSE, show.icc = FALSE, show.est = TRUE, p.adjust = "bonferroni",
          show.se = TRUE, show.p = TRUE, dv.labels = c("Log FFD", "Log SFD", "Log GD", "Log TT", "ROF","Rin"))

#SOME GRAPHING AND POSTHOCS##

library(interactions)
cat_plot(tt.mod, pred = gender_match,  y.label = "TT") 
cat_plot(Rin.mod, pred = gender_match,  y.label = "Rin") 
cat_plot(tt.mod, pred = type,  y.label = "TT") 
cat_plot(Rin.mod, pred = type,  y.label = "Rin") 

library(emmeans)
emmeans(tt.mod, pairwise ~ type, type="response", adjust ="bonferroni") #EMMs
# emmeans
# type          response   SE   df lower.CL upper.CL
# stim_adj           439 29.2 78.1      384      501
# stim_pred_adj      512 35.7 84.0      446      588
# stim_verb          396 26.7 76.8      347      453

Results are averaged over the levels of: gender_match, target_gender 
Degrees-of-freedom method: kenward-roger 
Confidence level used: 0.95 
Intervals are back-transformed from the log scale 

# $contrasts #ODDS SCALE
# contrast                  ratio     SE   df null t.ratio p.value
# stim_adj / stim_pred_adj  0.857 0.0633 51.8    1  -2.091  0.1244
# stim_adj / stim_verb      1.108 0.0761 42.0    1   1.490  0.4313
# stim_pred_adj / stim_verb 1.293 0.0981 48.3    1   3.386  0.0042

# $contrasts logit
# contrast                  estimate     SE   df t.ratio p.value
# stim_adj - stim_pred_adj    -0.155 0.0739 51.8  -2.091  0.1244
# stim_adj - stim_verb         0.102 0.0687 42.0   1.490  0.4313
# stim_pred_adj - stim_verb    0.257 0.0759 48.3   3.386  0.0042

emmeans(Rin.mod, pairwise ~ type, type="response", adjust ="bonferroni") #EMMs
# 
# $emmeans
# type           prob     SE  df asymp.LCL asymp.UCL
# stim_adj      0.350 0.0392 Inf     0.278     0.430
# stim_pred_adj 0.473 0.0455 Inf     0.385     0.562
# stim_verb     0.281 0.0355 Inf     0.217     0.355

Results are averaged over the levels of: gender_match, target_gender 
Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 
# 
# $contrasts
# contrast                  odds.ratio    SE  df null z.ratio p.value
# stim_adj / stim_pred_adj       0.602 0.127 Inf    1  -2.415  0.0472
# stim_adj / stim_verb           1.381 0.262 Inf    1   1.704  0.2653
# stim_pred_adj / stim_verb      2.295 0.492 Inf    1   3.877  0.0003

# $contrasts logit
# contrast                  estimate    SE  df z.ratio p.value
# stim_adj - stim_pred_adj    -0.508 0.210 Inf  -2.415  0.0472
# stim_adj - stim_verb         0.323 0.190 Inf   1.704  0.2653
# stim_pred_adj - stim_verb    0.831 0.214 Inf   3.877  0.0003


####################################################################################
#REGION 3 #

#DESCRIPTIVE STATS FOR REGION 3 ####
FF <- ddply(R3, c("gender_match"), summarise,
            N    = length(IA_FIRST_FIXATION_DURATION),
            mean = mean(IA_FIRST_FIXATION_DURATION, na.rm = T),
            sd   = sd(IA_FIRST_FIXATION_DURATION, na.rm = T))
FF$DV<-"FF"


GD <- ddply(R3, c("gender_match"), summarise,
            N    = length(IA_FIRST_RUN_DWELL_TIME),
            mean = mean(IA_FIRST_RUN_DWELL_TIME, na.rm = T),
            sd   = sd(IA_FIRST_RUN_DWELL_TIME, na.rm = T))
GD$DV<-"GD"



SF <- ddply(R3, c("gender_match"), summarise,
            N    = length(SFD),
            mean = mean(SFD, na.rm = T),
            sd   = sd(SFD, na.rm = T))

SF$DV<-"SF"


TT <- ddply(R3, c("gender_match"), summarise,
            N    = length(IA_DWELL_TIME),
            mean = mean(IA_DWELL_TIME, na.rm = T),
            sd   = sd(IA_DWELL_TIME, na.rm = T))

TT$DV<-"TT"


#descr table for durations
rbind(FF, SF, GD, TT)->descr_dur_R3
descr_dur_R3$SE<-descr_dur_R3$sd/sqrt(descr_dur_R3$N)
descr_dur_R3$Region<-"R3"


ROf <- ddply(R3, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_OUT_FULL),
             mean = mean(IA_REGRESSION_OUT_FULL, na.rm = T),
             sd   = sd(IA_REGRESSION_OUT_FULL, na.rm = T))

ROf$DV<-"ROf"

Rin <- ddply(R3, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_IN),
             mean = mean(IA_REGRESSION_IN, na.rm = T),
             sd   = sd(IA_REGRESSION_IN, na.rm = T))

Rin$DV<-"Rin"

#descr table for probabilities
rbind(skip, RO, ROf, Rin)->descr_prob_R3
descr_prob_R3$SE<-descr_prob_R3$sd/sqrt(descr_prob_R3$N)
descr_prob_R3$Region<-"R3"

#MODELING REGION 3##############################################################


ffd.mod <- lmer(log(IA_FIRST_FIXATION_DURATION) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 | DATA_FILE) + (1 | item),
                data = R3, 
                control = lmerControl(optimizer = "bobyqa"))

summary(ffd.mod)



sfd.mod <- lmer(log(SFD) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 | DATA_FILE) + (1 | item),
                data = R3, 
                control = lmerControl(optimizer = "bobyqa"))


gd.mod <- lmer(log(IA_FIRST_RUN_DWELL_TIME) ~ gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 + gender_match| DATA_FILE) + (1 | item),
               data = R3, 
               control = lmerControl(optimizer = "bobyqa"))


summary(gd.mod)


R3$IA_DWELL_TIME<-as.numeric(R3$IA_DWELL_TIME)
R3$lg_frequency<-as.numeric(R3$lg_frequency)
data.tt <- R3[R3$IA_DWELL_TIME != 0,]

tt.mod <- lmer(log(IA_DWELL_TIME) ~gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 + gender_match| DATA_FILE) + (1 | item),
               data = data.tt, 
               control = lmerControl(optimizer = "bobyqa"))

summary(tt.mod)




Rof.mod <- glmer(IA_REGRESSION_OUT_FULL ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 + gender_match| DATA_FILE) + (1 | item),
                 data = R3, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))

Rin.mod <- glmer(IA_REGRESSION_IN ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 | DATA_FILE) + (1 | item),
                 data = R3, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))
summary(Rin.mod)

library(sjPlot)

tab_model(ffd.mod, sfd.mod, gd.mod, tt.mod,Rof.mod, Rin.mod,
          string.se = "SE",digits = 3, show.ci = FALSE, show.icc = FALSE, show.est = TRUE, p.adjust = "bonferroni",
          show.se = TRUE, show.p = TRUE, dv.labels = c("Log FFD", "Log SFD", "Log GD", "Log TT", "ROF","Rin"))

#SOME GRAPHING AND POSTHOCS Region 3#

library(interactions)

#Diff between match/mismatch, but only at pred_adj level


mytheme <- theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    # legend.key = element_blank(),
    text = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 14),
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 14),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "gray96")
  )

cat_plot(sfd.mod, pred = gender_match, modx = type, y.label = "SFD", interval = TRUE,  
         point.shape = TRUE, 
         legend.main= "Agreement type", 
         modx.labels	= c("Long adjective", "Short adjective", "Verb"))->c

c+ theme(legend.position = "right") +
  mytheme +
  scale_y_continuous("SFD (Log)", limits = c(5, 6)) +
  scale_x_discrete("Grammaticality", labels = c("Match", "Mismatch")) +
  guides(fill = FALSE)->c


tiff("region2_SFD_interaction.tiff", units="in", width=7, height=5, res=300)
c
dev.off()


library(emmeans)
emmeans(sfd.mod, pairwise ~ type|gender_match, type="response", adjust ="bonferroni") #EMMs
$emmeans
# gender_match = Match:
#   type          response    SE    df lower.CL upper.CL
# stim_adj           233 10.15  75.5      214      254
# stim_pred_adj      235 10.50  85.8      215      257
# stim_verb          220  9.30  69.9      202      239
# 
# gender_match = Mis:
#   type          response    SE    df lower.CL upper.CL
# stim_adj           242 11.00  89.1      221      265
# stim_pred_adj      280 13.41 111.3      254      308
# stim_verb          221  9.82  84.8      203      242
# 
# Results are averaged over the levels of: target_gender 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# gender_match = Match:
#   contrast                  ratio     SE    df null t.ratio p.value
# stim_adj / stim_pred_adj  0.991 0.0412  82.4    1  -0.216  1.0000
# stim_adj / stim_verb      1.059 0.0426  62.0    1   1.430  0.4734
# stim_pred_adj / stim_verb 1.069 0.0417  84.8    1   1.705  0.2757
# 
# gender_match = Mis:
#   contrast                  ratio     SE    df null t.ratio p.value
# stim_adj / stim_pred_adj  0.864 0.0407 123.0    1  -3.114  0.0069
# stim_adj / stim_verb      1.092 0.0484  89.3    1   1.993  0.1478
# stim_pred_adj / stim_verb 1.265 0.0560 136.4    1   5.305  <.0001

# $contrasts logit
# gender_match = Match:
#   contrast                  estimate     SE    df t.ratio p.value
# stim_adj - stim_pred_adj  -0.00898 0.0416  82.4  -0.216  1.0000
# stim_adj - stim_verb       0.05749 0.0402  62.0   1.430  0.4734
# stim_pred_adj - stim_verb  0.06647 0.0390  84.8   1.705  0.2757
# 
# gender_match = Mis:
#   contrast                  estimate     SE    df t.ratio p.value
# stim_adj - stim_pred_adj  -0.14659 0.0471 123.0  -3.114  0.0069
# stim_adj - stim_verb       0.08839 0.0443  89.3   1.993  0.1478
# stim_pred_adj - stim_verb  0.23497 0.0443 136.4   5.305  <.0001
# 
# Results are averaged over the levels of: target_gender 
Degrees-of-freedom method: kenward-roger 
P value adjustment: bonferroni method for 3 tests 
Tests are performed on the log scale 


cat_plot(tt.mod, pred = gender_match,  y.label = "TT") 
cat_plot(RO.mod, pred = gender_match,  y.label = "RO") 
cat_plot(Rin.mod, pred = gender_match, y.label = "Rin") 


cat_plot(tt.mod, pred = type,  y.label = "TT") 
cat_plot(RO.mod, pred = type,  y.label = "RO") 

library(emmeans)
emmeans(tt.mod, pairwise ~ type, type="response", adjust ="bonferroni") #EMMs
# $emmeans
# type          response   SE   df lower.CL upper.CL
# stim_adj           353 21.4 71.8      313      398
# stim_pred_adj      455 26.8 69.0      404      512
# stim_verb          349 20.8 69.4      310      393
# 
# Results are averaged over the levels of: gender_match, target_gender 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Intervals are back-transformed from the log scale 
# 
# $contrasts
# contrast                  ratio     SE   df null t.ratio p.value
# stim_adj / stim_pred_adj  0.776 0.0432 44.9    1  -4.557  0.0001
# stim_adj / stim_verb      1.011 0.0586 42.5    1   0.189  1.0000
# stim_pred_adj / stim_verb 1.303 0.0674 42.3    1   5.115  <.0001

# $contrastslogit
# contrast                  estimate     SE   df t.ratio p.value
# stim_adj - stim_pred_adj    -0.254 0.0556 44.9  -4.557  0.0001
# stim_adj - stim_verb         0.011 0.0580 42.5   0.189  1.0000
# stim_pred_adj - stim_verb    0.265 0.0517 42.3   5.115  <.0001


emmeans(Rof.mod, pairwise ~ type, type="response", adjust ="bonferroni") #EMMs
# $emmeans
# type           prob     SE  df asymp.LCL asymp.UCL
# stim_adj      0.261 0.0344 Inf     0.199     0.333
# stim_pred_adj 0.417 0.0401 Inf     0.341     0.497
# stim_verb     0.245 0.0325 Inf     0.187     0.314
# 
# Results are averaged over the levels of: gender_match, target_gender 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 
# 
# $contrasts
# contrast                  odds.ratio    SE  df null z.ratio p.value
# stim_adj / stim_pred_adj       0.493 0.102 Inf    1  -3.412  0.0019
# stim_adj / stim_verb           1.085 0.239 Inf    1   0.373  1.0000
# stim_pred_adj / stim_verb      2.201 0.426 Inf    1   4.077  0.0001
# 
# Results are averaged over the levels of: gender_match, target_gender 
# P value adjustment: bonferroni method for 3 tests 
# Tests are performed on the log odds ratio scale 

# $contrasts logit
# contrast                  estimate    SE  df z.ratio p.value
# stim_adj - stim_pred_adj   -0.7069 0.207 Inf  -3.412  0.0019
# stim_adj - stim_verb        0.0819 0.220 Inf   0.373  1.0000
# stim_pred_adj - stim_verb   0.7888 0.193 Inf   4.077  0.0001
# 

#################################################################################################
#REGION 4######################

#################Descriptive stats###############################################

# REGION 4###

FF <- ddply(R4, c("gender_match"), summarise,
            N    = length(IA_FIRST_FIXATION_DURATION),
            mean = mean(IA_FIRST_FIXATION_DURATION, na.rm = T),
            sd   = sd(IA_FIRST_FIXATION_DURATION, na.rm = T))
FF$DV<-"FF"


GD <- ddply(R4, c("gender_match"), summarise,
            N    = length(IA_FIRST_RUN_DWELL_TIME),
            mean = mean(IA_FIRST_RUN_DWELL_TIME, na.rm = T),
            sd   = sd(IA_FIRST_RUN_DWELL_TIME, na.rm = T))
GD$DV<-"GD"



SF <- ddply(R4, c("gender_match"), summarise,
            N    = length(SFD),
            mean = mean(SFD, na.rm = T),
            sd   = sd(SFD, na.rm = T))

SF$DV<-"SF"


TT <- ddply(R4, c("gender_match"), summarise,
            N    = length(IA_DWELL_TIME),
            mean = mean(IA_DWELL_TIME, na.rm = T),
            sd   = sd(IA_DWELL_TIME, na.rm = T))

TT$DV<-"TT"


#descr table for durations
rbind(FF, SF, GD, TT)->descr_dur_R4
descr_dur_R4$SE<-descr_dur_R4$sd/sqrt(descr_dur_R4$N)
#create region column
descr_dur_R4$Region<-"R4"


ROf <- ddply(R4, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_OUT_FULL),
             mean = mean(IA_REGRESSION_OUT_FULL, na.rm = T),
             sd   = sd(IA_REGRESSION_OUT_FULL, na.rm = T))

ROf$DV<-"ROf"

Rin <- ddply(R4, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_IN),
             mean = mean(IA_REGRESSION_IN, na.rm = T),
             sd   = sd(IA_REGRESSION_IN, na.rm = T))

Rin$DV<-"Rin"


#descr table for probabilities
rbind(skip, RO, ROf, Rin)->descr_prob_R4
descr_prob_R4$SE<-descr_prob_R4$sd/sqrt(descr_prob_R4$N)
descr_prob_R4$Region<-"R4"



#MODELING REGION 4#################################################################


ffd.mod <- lmer(log(IA_FIRST_FIXATION_DURATION) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 + gender_match| DATA_FILE) + (1 | item),
                data = R4, 
                control = lmerControl(optimizer = "bobyqa"))

summary(ffd.mod)





sfd.mod <- lmer(log(SFD) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 | DATA_FILE) + (1 | item),
                data = R4, 
                control = lmerControl(optimizer = "bobyqa"))

plot(fitted(sfd.mod), residuals(sfd.mod))





gd.mod <- lmer(log(IA_FIRST_RUN_DWELL_TIME) ~ gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 + gender_match| DATA_FILE) + (1 | item),
               data = R4, 
               control = lmerControl(optimizer = "bobyqa"))


summary(gd.mod)


R4$IA_DWELL_TIME<-as.numeric(R4$IA_DWELL_TIME)

R4$lg_frequency<-as.numeric(R4$lg_frequency)
data.tt <- R4[R4$IA_DWELL_TIME != 0,]

tt.mod <- lmer(log(IA_DWELL_TIME) ~gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 + gender_match| DATA_FILE) + (1 | item),
               data = data.tt, 
               control = lmerControl(optimizer = "bobyqa"))

summary(tt.mod)


summary(skip.mod)



Rof.mod <- glmer(IA_REGRESSION_OUT_FULL ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 | DATA_FILE) + (1 | item),
                 data = R4, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))

Rin.mod <- glmer(IA_REGRESSION_IN ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 | DATA_FILE) + (1 | item),
                 data = R4, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))
summary(Rin.mod)

library(sjPlot)

tab_model(ffd.mod, sfd.mod, gd.mod, tt.mod,Rof.mod, Rin.mod,
          string.se = "SE",digits = 3, show.ci = FALSE, show.icc = FALSE, 
          show.est = TRUE, p.adjust = "bonferroni",
          show.se = TRUE, show.p = TRUE,
          dv.labels = c("Log FFD", "Log SFD", "Log GD", "Log TT", "ROF","Rin"))

##POTHOCS AND GRAPHING###################################################

library(interactions)
#Diff between match/mismatch, but only at pred_adj level in match compared to long adjetcives
cat_plot(ffd.mod, pred = gender_match, modx=type, y.label = "FFD") 
library(emmeans)
emmeans(ffd.mod, pairwise ~ type|gender_match, type="response", adjust ="bonferroni") #EMMs

# gender_match = Match:
#   type          response   SE   df lower.CL upper.CL
# stim_adj           215 7.70 83.3      200      231
# stim_pred_adj      236 8.08 77.9      220      252
# stim_verb          224 7.60 76.0      209      240
# 
# gender_match = Mis:
#   type          response   SE   df lower.CL upper.CL
# stim_adj           218 7.40 87.1      204      233
# stim_pred_adj      217 7.00 82.6      203      231
# stim_verb          225 7.20 81.1      211      240

Results are averaged over the levels of: target_gender 
Degrees-of-freedom method: kenward-roger 
Confidence level used: 0.95 
Intervals are back-transformed from the log scale 

# $contrasts
# gender_match = Match:
#   contrast                  ratio     SE   df null t.ratio p.value
# stim_adj / stim_pred_adj  0.912 0.0335 75.6    1  -2.492  0.0447
# stim_adj / stim_verb      0.961 0.0345 74.8    1  -1.103  0.8213
# stim_pred_adj / stim_verb 1.053 0.0336 89.4    1   1.633  0.3177
# 
# gender_match = Mis:
#   contrast                  ratio     SE   df null t.ratio p.value
# stim_adj / stim_pred_adj  1.006 0.0369 75.0    1   0.171  1.0000
# stim_adj / stim_verb      0.970 0.0348 74.8    1  -0.848  1.0000
# stim_pred_adj / stim_verb 0.964 0.0307 89.0    1  -1.153  0.7554

# $contrasts logit
# gender_match = Match:
#   contrast                  estimate     SE   df t.ratio p.value
# stim_adj - stim_pred_adj  -0.09157 0.0367 75.6  -2.492  0.0447
# stim_adj - stim_verb      -0.03953 0.0359 74.8  -1.103  0.8213
# stim_pred_adj - stim_verb  0.05204 0.0319 89.4   1.633  0.3177
# 
# gender_match = Mis:
#   contrast                  estimate     SE   df t.ratio p.value
# stim_adj - stim_pred_adj   0.00626 0.0367 75.0   0.171  1.0000
# stim_adj - stim_verb      -0.03043 0.0359 74.8  -0.848  1.0000
# stim_pred_adj - stim_verb -0.03669 0.0318 89.0  -1.153  0.7554




cat_plot(tt.mod, pred = gender_match, y.label = "TT") 
cat_plot(RO.mod, pred = gender_match, y.label = "Ro")
cat_plot(Rof.mod, pred = gender_match, y.label = "Rof") 



#################################################################################################
#REGION 5######################

#################Descriptive stats###############################################

# REGION 5###

FF <- ddply(R5, c("gender_match"), summarise,
            N    = length(IA_FIRST_FIXATION_DURATION),
            mean = mean(IA_FIRST_FIXATION_DURATION, na.rm = T),
            sd   = sd(IA_FIRST_FIXATION_DURATION, na.rm = T))
FF$DV<-"FF"


GD <- ddply(R5, c("gender_match"), summarise,
            N    = length(IA_FIRST_RUN_DWELL_TIME),
            mean = mean(IA_FIRST_RUN_DWELL_TIME, na.rm = T),
            sd   = sd(IA_FIRST_RUN_DWELL_TIME, na.rm = T))
GD$DV<-"GD"



SF <- ddply(R5, c("gender_match"), summarise,
            N    = length(SFD),
            mean = mean(SFD, na.rm = T),
            sd   = sd(SFD, na.rm = T))

SF$DV<-"SF"


TT <- ddply(R5, c("gender_match"), summarise,
            N    = length(IA_DWELL_TIME),
            mean = mean(IA_DWELL_TIME, na.rm = T),
            sd   = sd(IA_DWELL_TIME, na.rm = T))

TT$DV<-"TT"


#descr table for durations
rbind(FF, SF, GD, TT)->descr_dur_R5
descr_dur_R5$SE<-descr_dur_R5$sd/sqrt(descr_dur_R5$N)
#create region column
descr_dur_R5$Region<-"R5"



ROf <- ddply(R5, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_OUT_FULL),
             mean = mean(IA_REGRESSION_OUT_FULL, na.rm = T),
             sd   = sd(IA_REGRESSION_OUT_FULL, na.rm = T))

ROf$DV<-"ROf"

Rin <- ddply(R5, c("gender_match"), summarise,
             N    = length(IA_REGRESSION_IN),
             mean = mean(IA_REGRESSION_IN, na.rm = T),
             sd   = sd(IA_REGRESSION_IN, na.rm = T))

Rin$DV<-"Rin"



#descr table for probs
rbind(skip, RO, ROf, Rin)->descr_prob_R5
descr_prob_R5$SE<-descr_prob_R5$sd/sqrt(descr_prob_R5$N)
descr_prob_R5$Region<-"R5"



#MODELING REGION 5####################################################################

ffd.mod <- lmer(log(IA_FIRST_FIXATION_DURATION) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 | DATA_FILE) + (1 | item),
                data = R5, 
                control = lmerControl(optimizer = "bobyqa"))

summary(ffd.mod)





sfd.mod <- lmer(log(SFD) ~ gender_match*type+ target_gender+
                  lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                  (1 | DATA_FILE) + (1 | item),
                data = R5, 
                control = lmerControl(optimizer = "bobyqa"))

plot(fitted(sfd.mod), residuals(sfd.mod))





gd.mod <- lmer(log(IA_FIRST_RUN_DWELL_TIME) ~ gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 | DATA_FILE) + (1 | item),
               data = R5, 
               control = lmerControl(optimizer = "bobyqa"))


summary(gd.mod)


R5$IA_DWELL_TIME<-as.numeric(R5$IA_DWELL_TIME)

R5$lg_frequency<-as.numeric(R5$lg_frequency)
data.tt <- R5[R5$IA_DWELL_TIME != 0,]

tt.mod <- lmer(log(IA_DWELL_TIME) ~gender_match*type+ target_gender+
                 lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                 (1 | DATA_FILE) + (1 | item),
               data = data.tt, 
               control = lmerControl(optimizer = "bobyqa"))

summary(tt.mod)



summary(skip.mod)



Rof.mod <- glmer(IA_REGRESSION_OUT_FULL ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 + gender_match| DATA_FILE) + (1 | item),
                 data = R5, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))

Rin.mod <- glmer(IA_REGRESSION_IN ~gender_match*type+ target_gender+
                   lg_frequency+len.scaled+prev.len.scaled+next.len.scaled+ TRIAL_INDEX_sc+
                   (1 | DATA_FILE) + (1 | item),
                 data = R5, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))
summary(Rin.mod)

library(sjPlot)


tab_model(ffd.mod, sfd.mod, gd.mod, tt.mod,Rof.mod, Rin.mod,
          string.se = "SE",digits = 3, show.ci = FALSE, show.icc = FALSE, show.est = TRUE, p.adjust = "bonferroni",
          show.se = TRUE, show.p = TRUE, dv.labels = c("Log FFD", "Log SFD", "Log GD", "Log TT",  "ROF","Rin"))

####POSTHOCS AND GRAPHING##################

library(interactions)
cat_plot(Rof.mod, pred = type, y.label = "Rof") 
library(emmeans)
emmeans(Rof.mod, pairwise ~ type, type="response", adjust ="bonferroni") #EMMs

# emmeans
# type           prob     SE  df asymp.LCL asymp.UCL
# stim_adj      0.162 0.0247 Inf     0.119     0.216
# stim_pred_adj 0.274 0.0310 Inf     0.217     0.338
# stim_verb     0.268 0.0327 Inf     0.209     0.336
# 
# Results are averaged over the levels of: gender_match, target_gender 
# Confidence level used: 0.95 
# Intervals are back-transformed from the logit scale 
# 
# $contrasts
# contrast                  odds.ratio    SE  df null z.ratio p.value
# stim_adj / stim_pred_adj       0.512 0.117 Inf    1  -2.927  0.0103
# stim_adj / stim_verb           0.528 0.132 Inf    1  -2.549  0.0325
# stim_pred_adj / stim_verb      1.031 0.216 Inf    1   0.145  1.0000

# $contrasts logit
# contrast                  estimate    SE  df z.ratio p.value
# stim_adj - stim_pred_adj   -0.6696 0.229 Inf  -2.927  0.0103
# stim_adj - stim_verb       -0.6392 0.251 Inf  -2.549  0.0325
# stim_pred_adj - stim_verb   0.0304 0.210 Inf   0.145  1.0000




#################################################################################################
#GRAPHING DESCRIPTIVE STATS######################
rbind(descr_dur_R2, descr_dur_R3, descr_dur_R4, descr_dur_R5)->dur_R2345

library(ggplot2)
#descr chart


mytheme <- theme_bw() + theme(axis.title = element_text(size = 20),
                              strip.text = element_text(size = 20),
                              legend.key = element_blank(),
                              text = element_text(size = 20),
                              legend.text = element_text(size = 20), 
                              axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 20),
                              plot.title = element_text(size=20, hjust = 0.5, vjust = 0.5),
                              axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 20),
                              strip.background = element_rect(fill = 'gray96'))


# Reorder factor levels
library(plyr)
dur_R2345$DV <- factor(dur_R2345$DV, levels=c("FF", "SF", "GD", "TT") )


tiff("descr_dur_R2345.tiff", units="in", width=20, height=20, res=300)
ggplot(data = dur_R2345) +mytheme+
  aes(x = gender_match, y = mean,  fill = DV) + geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2, show.legend=FALSE,position=position_dodge(.9)) + 
  theme(legend.title = element_blank())+
  ggtitle("")+
  xlab("") + 
  ylab("Duration (ms)")+
  scale_fill_manual(values=c('lightgrey','darkgrey','black', 'grey'))->l
l + facet_grid(. ~ Region)->l
dev.off()


rbind(descr_prob_R2, descr_prob_R3, descr_prob_R4, descr_prob_R5)->descr_prob2345



tiff("descr_prob_R2345.tiff", units="in", width=20, height=7, res=300)
ggplot(data = descr_prob2345) +mytheme+
  aes(x = gender_match, y = mean,  fill = DV) + geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2, show.legend=FALSE,position=position_dodge(.9)) + 
  theme(legend.title = element_blank())+
  ggtitle("")+
  xlab("") + 
  ylab("Proportion")+
  scale_fill_manual(values=c('lightgrey','darkgrey','black', 'grey'))->k
k + facet_grid(. ~ Region)->
  k
dev.off()



library(ggpubr)
fg<-ggarrange(l, k,
              labels = c("A", "B"),  font.label = list(size = 16),label.x = 0.00, common.legend = FALSE, legend = "right",
              ncol = 1, nrow = 2)

fg

tiff("Figure_descriptives.tiff", units="in", width=20, height=20, res=300)
fg
dev.off()




