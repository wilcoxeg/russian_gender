## ==============================================
## this code process and analyzes the SPR data
## from the study on Russian internal vs external agreement
## last updated 9/26/23
## ==============================================

library(tidyr)
library(ibextor)
library(readr)
library(ggplot2)
library(magrittr)
library(lme4)
library(lmerTest)

# ibextor is a package written by Anton Malko. 
# If not downloaded, do the following:
# install.packages("pkgload")
# devtools::install_github("antonmalko/ibextor")

# ---- Load data ----

results_file <- "results_13Apr2020.txt" # indicate USER's path / working dir.

sent_dat <- get_results_ds(results_file, encoding = "UTF-8")
quest_dat <- get_results_q(results_file, encoding = "UTF-8")
subj_info <- get_subj_info(results_file, encoding = "UTF-8")

# ---- Subject info ----

# Add subj number - it does not get added due to a bug in ibextor package,
# but people are collected from the results file in order of their appearance,
# so we can just add a sequence of numbers
subj_info$subj <- 1:nrow(subj_info)

# Some people wrote "18 years" (in Russian) in the "age" field; clean it up
subj_info$age <- gsub("[^[:digit:]]","", subj_info$age)

subj_info$age <- as.numeric(subj_info$age)
subj_info$sex <- factor(subj_info$sex)

summary(subj_info)

subj_info %>%
  dplyr::arrange(age) %T>% # Tee operator from magrittr
  print() %>%
  dplyr::pull(age) %>%
  boxplot()
# 1 subjects is older than 40, exclude them later

quest_acc <- quest_dat %>%
  dplyr::group_by(subj) %>%
  dplyr::summarize(n_quest = n(),
                   n_correct = sum(is_correct),
                   prop_correct = n_correct/n_quest)
summary(quest_acc)
# All subjects have at least 82% of correct answers, no reasons for
# exclusion here



# ---- Assign experimental factors ----

stim_dat <- sent_dat[sent_dat$type %in% c("a","b","c","d"),]
fillers_dat <- sent_dat[sent_dat$type == "filler",]

# Encode agreeing element type
stim_dat$agr_elem_type <- ifelse(stim_dat$item %in% 1:16, "adj",
                                 ifelse(stim_dat$item %in% 17:32, "verb", "pred_adj"))

# Encode gender match
stim_dat$gender_match <- ifelse(stim_dat$type %in% c("a","c"), "match", "mismatch")
stim_dat$gram <- ifelse(stim_dat$gender_match == "match", "gram", "ungram")

# Add dummy columns to fillers so that we can use them later during residualisation
fillers_dat$agr_elem_type <- NA
fillers_dat$gender_match <- NA

# Encode fillers grammaticality
fillers_dat$gram <- ifelse(fillers_dat$item %in% 101:145, "gram", "ungram")

# ---- Check data integrity ----

# Check that we have the expected number of sentences per subject ----

# Only experimental - expect 48 per subject (thus, 48 should be the only number in the output
# of the code below)
stim_dat %>%
  dplyr::group_by(subj) %>%
  dplyr::filter(region == 1, type %in% c("a","b","c","d")) %>%
  dplyr::count() %>%
  dplyr::pull(n) %>%
  unique()

# Experimental stimuli: by condition -- should be 12: four per each PoS (verb, adj, predicative adj)
stim_dat %>%
  dplyr::group_by(type,subj) %>%
  dplyr::filter(region == 1, type %in% c("a","b","c","d")) %>%
  dplyr::count() %>%
  dplyr::pull(n) %>%
  unique()

# Only fillers - expect 90 per subject # ZZF note: I am getting 0 instead of 90
stim_dat %>%
  dplyr::group_by(subj) %>%
  dplyr::filter(region == 1, type == "filler") %>%
  dplyr::count() %>%
  dplyr::pull(n) %>%
  unique()

# Check that the experimental factors were assigned correctly ----

# check that verbs and adjs do only have these words in pre-critical region.
# Replace the filtering condition with "verb", "adj" or "pred_adj" to check all
# three groups
stim_dat %>%
  dplyr::filter(agr_elem_type == "pred_adj", region == 2) %>%
  dplyr::pull(word) %>%
  unique()

# Check grammaticality (gram vs  ungram). Replace the filtering clause with "gram" or "ungram"
# to check both conditions
stim_dat %>%
  dplyr::filter(gram == "gram", region %in% c(2,3)) %>%
  dplyr::select(subj, item, region, word) %>%
  tidyr::pivot_wider(id_cols = c(subj,item), names_from = region, names_prefix = "reg_", values_from = "word") %>%
  dplyr::mutate(comb = paste(reg_2, reg_3)) %>%
  dplyr::pull(comb) %>%
  unique()

# ---- Check for overly long reading times ----
summary(stim_dat$rt)
small_rts <- stim_dat[stim_dat$rt < 50,]
extreme_rts <- stim_dat[stim_dat$rt > 1000,]

tmp <- boxplot(stim_dat$rt)
# There are reaction times from 0 ms (just one) up to 37 seconds

# Let's check pre-critical, critical and spillover
extreme_rts_rois <- stim_dat[stim_dat$rt > 1000 & stim_dat$region %in% 2:4,]

# ==== Data cleanup ====

# exclude subjects older than the majority of the group

bad_subjects <- c(55) # age > 40 
stim_dat_clean <- stim_dat[!stim_dat$subj %in% bad_subjects,]
fillers_dat_clean <- fillers_dat[!fillers_dat$subj %in% bad_subjects,]
# original analysis exluded outliers; we will not (Matuschek et al. 2017)
    # exclude RTs > 5000 ms and shorter than 50 ms
    # stim_dat_clean <- stim_dat_clean[stim_dat_clean$rt > 50 & stim_dat_clean$rt <= 5000,]
    # fillers_dat_clean <- fillers_dat_clean[fillers_dat_clean$rt > 50 & fillers_dat_clean$rt <= 5000,]

# ==== Data analysis ====

# Statistical analyses
#
# Long-transform RTs and residualise based on length and by-subject random intercept
stim_and_fillers_clean <- rbind(stim_dat_clean, fillers_dat_clean)
stim_and_fillers_clean %<>%
  dplyr::mutate(word_length = nchar(word),
                log_rt = log(rt))

# because we didn't remove outliers, we need to remove -Inf values
stim_and_fillers_clean <- stim_and_fillers_clean %>% filter(rt>0)
aux_model_for_resid <- lmer(log_rt ~ word_length + (1|subj), data = stim_and_fillers_clean)
stim_and_fillers_clean$resid_log_rt <- resid(aux_model_for_resid)

# Notice that we are replacing an existing object. But in fact, it
# will be the same as the old one with one additional column for residualized log_rts
stim_dat_clean <- stim_and_fillers_clean[stim_and_fillers_clean$type %in% c("a","b","c","d"),]

# Set contrast - we'll use sum coding
stim_dat_clean$agr_elem_type <- factor(stim_dat_clean$agr_elem_type)
contrasts(stim_dat_clean$agr_elem_type) <- contr.sum(3)
    #exploratory: other coding schemes
    #contrasts(stim_dat_clean$agr_elem_type) <- contr.helmert(3) # added by ZZF
    #stim_dat_clean$agr_elem_type <- factor(stim_dat_clean$agr_elem_type, levels=c("pred_adj","verb","adj"))
    #stim_dat_clean$agr_elem_type <- factor(stim_dat_clean$agr_elem_type, levels=c("adj","verb","pred_adj"))
    #contrasts(stim_dat_clean$agr_elem_type) <- contr.helmert(3)

stim_dat_clean$gram <- factor(stim_dat_clean$gram)
contrasts(stim_dat_clean$gram) <- contr.sum(2)


# ---- Models of interest -- on residualized data ----
# models <- list() # Anton's code

# Order of model simplification:
# 1) interaction by items
# 2) interaction by subjects
# 3) agr_elem_type by item
# 4) agr_elem_type by subj
# 5) gram by item
# 6) gram by subj

# pre-critical region

prec.mod1 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram*agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,]) # singular fit
prec.mod2 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,]) # singular fit
prec.mod3 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,]) # singular fit
prec.mod4 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,]) # singular fit
prec.mod5 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1+gram| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,]) # singular fit
prec.mod6 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,])
prec.mod7 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1| subj) + (1| item),
                  data = stim_dat_clean[stim_dat_clean$region == 2,])

anova(prec.mod6, prec.mod7) # the random slope for gram grouped by
# subject does not improve the model. Use prec.mod7.
summary(prec.mod7)  

    # model validation
    qqnorm(resid(prec.mod7))
    qqline(resid(prec.mod7))
    plot(resid(prec.mod7), predict(prec.mod7))

# critical region
crit.mod1 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram*agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) # singular fit
crit.mod2 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) # singular fit
crit.mod3 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) # singular fit
crit.mod4 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) # singular fit
crit.mod5 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1+gram| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) # singular fit
crit.mod6 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) # singular fit
crit.mod7 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 | subj) + (1| item),
                  data = stim_dat_clean[stim_dat_clean$region == 3,]) 

summary(crit.mod7)
    # model validation
    qqnorm(resid(crit.mod7))
    qqline(resid(crit.mod7))
    plot(resid(crit.mod7), predict(crit.mod7))
    # two overly influential data points, potentially
    lsmeans(crit.mod7, ~gram)
    lsmeans(crit.mod7, ~agr_elem_type)
    
# spillover region 1  
spill1.mod1 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram*agr_elem_type| item),
                  data = stim_dat_clean[stim_dat_clean$region == 4,]) # singular fit
spill1.mod2 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                    data = stim_dat_clean[stim_dat_clean$region == 4,]) # singular fit
spill1.mod3 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                    data = stim_dat_clean[stim_dat_clean$region == 4,]) # singular fit
spill1.mod4 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram| item),
                    data = stim_dat_clean[stim_dat_clean$region == 4,]) # singular fit
spill1.mod5 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1+gram| item),
                    data = stim_dat_clean[stim_dat_clean$region == 4,]) # singular fit
spill1.mod6 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1| item),
                    data = stim_dat_clean[stim_dat_clean$region == 4,]) # singular fit
spill1.mod7 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 | subj) + (1| item),
                    data = stim_dat_clean[stim_dat_clean$region == 4,]) 

summary(spill1.mod7)
    # model validation
    qqnorm(resid(spill1.mod7))
    qqline(resid(spill1.mod7))
    plot(resid(spill1.mod7), predict(spill1.mod7))
    # two overly influential data points, potentially
    lsmeans(spill1.mod7, ~gram)
    lsmeans(spill1.mod7, ~agr_elem_type)

# spillover region 2
spill2.mod1 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram*agr_elem_type| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) # singular fit
spill2.mod2 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram*agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) # singular fit
spill2.mod3 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram+agr_elem_type| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) # singular fit
spill2.mod4 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram+agr_elem_type| subj) + (1+gram| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) # singular fit
spill2.mod5 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1+gram| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) # singular fit
spill2.mod6 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 + gram| subj) + (1| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) # singular fit
spill2.mod7 <- lmer(resid_log_rt ~ gram*agr_elem_type + (1 | subj) + (1| item),
                    data = stim_dat_clean[stim_dat_clean$region == 5,]) 

summary(spill2.mod7)
lsmeans(spill2.mod7, ~agr_elem_type)

# ---- Adjectives vs verbs ----

stim_dat_clean$synt_type <- factor(ifelse(stim_dat_clean$agr_elem_type == "adj", "mod", "pred")) # modifier vs predicates
stim_dat_clean$pos <- factor(ifelse(stim_dat_clean$agr_elem_type == "verb", "verb", "adj")) # adjectives vs verbs
stim_dat_clean$pred_adj <- factor(ifelse(stim_dat_clean$agr_elem_type == "pred_adj", "pa", "non-pa")) # pred. adj vs mod.adj&verb

contrasts(stim_dat_clean$synt_type) <- contr.sum(2)
contrasts(stim_dat_clean$pos) <- contr.sum(2)
contrasts(stim_dat_clean$pred_adj) <- contr.sum(2)

# precritical region
pos.prec.mod1 <- lmer(resid_log_rt ~ gram*pos + (1 | subj) + (1| item),
                      data = stim_dat_clean[stim_dat_clean$region == 2,])
pos.prec.mod2 <- lmer(resid_log_rt ~ gram*pos + (1 +gram| subj) + (1| item),
                      data = stim_dat_clean[stim_dat_clean$region == 2,])
anova(pos.prec.mod1, pos.prec.mod2) #random slope for gram does not improve model,
                                    #use poc.prec.mod1
summary(pos.prec.mod1)

# critical region
pos.crit.mod1 <- lmer(resid_log_rt ~ gram*pos + (1 | subj) + (1| item),
                      data = stim_dat_clean[stim_dat_clean$region == 3,])
    # any more complexity in random effect structure -> singular fit
summary(pos.crit.mod1)

# spillover 1
pos.spill1.mod1 <- lmer(resid_log_rt ~ gram*pos + (1 | subj) + (1| item),
                      data = stim_dat_clean[stim_dat_clean$region == 4,])
    # any more complexity in random effect structure -> singular fit
summary(pos.spill1.mod1)

# spillover 2
pos.spill2.mod1 <- lmer(resid_log_rt ~ gram*pos + (1 | subj) + (1| item),
                        data = stim_dat_clean[stim_dat_clean$region == 5,])
    # any more complexity in random effect structure -> singular fit
summary(pos.spill2.mod1)

# ---- Modifiers vs predicates ----

# precritical region
synt.prec.mod1 <- lmer(resid_log_rt ~ gram*synt_type + (1 | subj) + (1| item),
                      data = stim_dat_clean[stim_dat_clean$region == 2,])
synt.prec.mod2 <- lmer(resid_log_rt ~ gram*synt_type + (1+gram | subj) + (1| item),
                       data = stim_dat_clean[stim_dat_clean$region == 2,])
anova(synt.prec.mod2, synt.prec.mod1) #random slope for gram does not improve model,
                                      #use synt.prec.mod1
summary(synt.prec.mod1)

# critical region
synt.crit.mod1 <- lmer(resid_log_rt ~ gram*synt_type + (1 | subj) + (1| item),
                       data = stim_dat_clean[stim_dat_clean$region == 3,])
    # any more complexity in random effect structure -> singular fit
summary(synt.crit.mod1)

# spillover 1 
synt.spill1.mod1 <- lmer(resid_log_rt ~ gram*synt_type + (1 | subj) + (1| item),
                       data = stim_dat_clean[stim_dat_clean$region == 4,])
    # any more complexity in random effect structure -> singular fit
summary(synt.spill1.mod1)

# spillover 2 
synt.spill2.mod1 <- lmer(resid_log_rt ~ gram*synt_type + (1 | subj) + (1| item),
                         data = stim_dat_clean[stim_dat_clean$region == 5,])
# any more complexity in random effect structure -> singular fit
summary(synt.spill2.mod1)
lsmeans(synt.spill2.mod1, ~synt_type)


# ==== Visualising data ====

se <- function(x) sqrt(var(x)/length(x))

# ---- Sanity check: gram vs. ungramm fillers ----

fillers_dat_clean %>%
  dplyr::group_by(region, gram) %>%
  dplyr::summarize(mean_rt = mean(rt)) %>%
  ggplot(aes(x = region, y = mean_rt, color = gram, group = gram)) +
  geom_point() +
  geom_line() +
  ggtitle("Fillers: gram vs ungram") +
  ggsave(filename = "fillers_gram_ungram.png", plots$fillers) # ZZF note: object "plots" not found

# ---- Experimental conditions: gram vs. ungram

# Raw RTs
stim_dat_clean %>%
  dplyr::group_by(region, gram) %>%
  dplyr::summarize(mean_rt = mean(rt), std_err = se(rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = gram, group = gram)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: gram vs ungram (Raw RTs)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25) 

#Residualised log RTs
stim_dat_clean %>%
  dplyr::group_by(region, gram) %>%
  dplyr::summarize(mean_rt = mean(resid_log_rt), std_err = se(resid_log_rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, linetype = gram, group = gram)) +
  geom_point() + geom_line() +
  ylab("Residualized log(RT)") + xlab("Region") + labs(linetype="Grammaticality") +
  ggtitle("Reading times for grammatical vs ungrammatical") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25) 

#Residualised log RTs with error bars
stim_dat_clean %>%
  dplyr::group_by(region, gram) %>%
  dplyr::summarize(mean_rt = mean(resid_log_rt), std_err = se(resid_log_rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = gram, group = gram)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: gram vs ungram (Residualised log RTs)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25) 

# log RTs
stim_dat_clean %>%
  dplyr::group_by(region, gram) %>%
  dplyr::summarize(mean_rt = mean(log_rt), std_err = se(log_rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = gram, group = gram)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: gram vs ungram (Log RTs)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25) 


# ---- Experimental conditions: agreeing element type ----

# Raw RTs
stim_dat_clean %>%
  dplyr::group_by(region, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(rt), std_err = se(rt))%>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = agr_elem_type)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: agreeing element type (Raw RTs") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25)

# Residualised log RTs
stim_dat_clean %>%
  dplyr::group_by(region, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(resid_log_rt), std_err = se(resid_log_rt))%>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = agr_elem_type)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: agreeing element type (Residualised log RTs") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25)


# Residualised log RTs with error bars
stim_dat_clean %>%
  dplyr::group_by(region, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(resid_log_rt), std_err = se(resid_log_rt))%>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = agr_elem_type)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: agreeing element type (Residualised log RTs") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25) 

# Log RTs
stim_dat_clean %>%
  dplyr::group_by(region, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(log_rt), std_err = se(log_rt))%>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = agr_elem_type)) +
  geom_point() +
  geom_line() +
  ggtitle("Stimuli: agreeing element type (Log RTs") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25)

# ---- Experimental conditions: gram x type ----

# Raw RTs
stim_dat_clean %>%
  dplyr::group_by(region, gram, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(rt), std_err = se(rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = interaction(gram, agr_elem_type))) +
  geom_point() +
  geom_line(aes(linetype = gram)) +
  ggtitle("Stimuli: experimental manipulation (Raw RTs)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25)

# Residualised log RTs
stim_dat_clean %>%
  dplyr::group_by(region, gram, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(resid_log_rt), std_err = se(resid_log_rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = interaction(gram, agr_elem_type))) +
  geom_point() +
  geom_line(aes(linetype = gram)) +
  ggtitle("Stimuli: experimental manipulation (Residualised log RTs)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25)


# Residualised log RTs with error bars
stim_dat_clean %>%
  dplyr::group_by(region, gram, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(resid_log_rt), std_err = se(resid_log_rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = interaction(gram, agr_elem_type))) +
  geom_point() + 
  geom_line(aes(linetype = gram)) +
  ylab("Residualized log(RT)") + xlab("Region") + labs(linetype="Grammaticality", color = "Agreeing element") + 
  ggtitle("Reading times: grammaticality x agreeing element") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25) +
  scale_color_discrete(labels = c("long adj.", "short adj.", "verb"))  +
  scale_linetype_discrete(labels = c("grammatical", "ungrammatical"))

stim_dat_clean %>%
  dplyr::group_by(region, gram, agr_elem_type) %>%
  dplyr::summarize(mean_rt = mean(log_rt), std_err = se(log_rt)) %>%
  dplyr::mutate(ymin = mean_rt - std_err, ymax = mean_rt + std_err) %>%
  ggplot(aes(x = region, y = mean_rt, color = agr_elem_type, group = interaction(gram, agr_elem_type))) +
  geom_point() +
  geom_line(aes(linetype = gram)) +
  ggtitle("Stimuli: experimental manipulation (Log RTs)") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width= 0.25)

# ==== EXPLORATORY ====

# Try taking frequency into account on the pre-critical and spillover.
freq_data <- readxl::read_excel("C:/Users/zuzan/OneDrive/Documents/Stimuli.xlsx", sheet = "freq_to_use", trim_ws = TRUE) %>%
  dplyr::filter(!is.na(freq))

stim_dat_clean_freq <- stim_dat_clean[stim_dat_clean$region %in% 2:3,] %>%
  dplyr::left_join(freq_data)

aux_model_for_resid_freq <- lmer(log(rt) ~ word_length + freq + (1|subj), data = stim_dat_clean_freq)
stim_dat_clean_freq$resid_log_rt_freq <- resid(aux_model_for_resid_freq)

# Pre-crit
models$freq$precrit <- lmer(resid_log_rt_freq ~ gram*agr_elem_type + (1 | item),
                            data = stim_dat_clean_freq[stim_dat_clean_freq$region == 2,])

models$freq$crit <- lmer(resid_log_rt_freq ~ gram*agr_elem_type + (1 | item),
                         data = stim_dat_clean_freq[stim_dat_clean_freq$region == 3,])

lm(resid_log_rt_freq ~ gram*pos,
   data = stim_dat_clean_freq[stim_dat_clean_freq$region == 3,]) %>% summary()

# --- compare verbs and predicative adjectives separately ---
# quick check, so need to vary region manually; output goes directly to "summary()"

stim_dat_clean_no_adj <- stim_dat_clean[stim_dat_clean$agr_elem_type != "adj",]
stim_dat_clean_no_adj <- droplevels(stim_dat_clean_no_adj)
contrasts(stim_dat_clean_no_adj$gram) <- contr.sum(2)
contrasts(stim_dat_clean_no_adj$agr_elem_type) <- contr.sum(2)

lmer(resid_log_rt ~ gram*agr_elem_type + (1 | subj) + (1 | item),
     data = stim_dat_clean_no_adj[stim_dat_clean_no_adj$region == 4,]) %>% summary()

# --- compare verbs and short adjectives separately ---

stim_dat_clean_no_pred_adj <- stim_dat_clean[stim_dat_clean$agr_elem_type != "pred_adj",]
stim_dat_clean_no_pred_adj <- droplevels(stim_dat_clean_no_pred_adj)
contrasts(stim_dat_clean_no_pred_adj$gram) <- contr.sum(2)
contrasts(stim_dat_clean_no_pred_adj$agr_elem_type) <- contr.sum(2)

lmer(resid_log_rt ~ gram*agr_elem_type + (1 | subj) + (1 | item),
     data = stim_dat_clean_no_pred_adj[stim_dat_clean_no_pred_adj$region == 4,]) %>% summary()

# --- compare full adjectives and short adjectives separately ---

stim_dat_clean_no_verb <- stim_dat_clean[stim_dat_clean$agr_elem_type != "verb",]
stim_dat_clean_no_verb <- droplevels(stim_dat_clean_no_verb)
contrasts(stim_dat_clean_no_verb$gram) <- contr.sum(2)
contrasts(stim_dat_clean_no_verb$agr_elem_type) <- contr.sum(2)

lmer(resid_log_rt ~ gram*agr_elem_type + (1 | subj) + (1 | item),
     data = stim_dat_clean_no_verb[stim_dat_clean_no_verb$region == 4,]) %>% summary()


