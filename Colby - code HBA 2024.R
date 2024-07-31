###########
## packages
###########

library(brms)
library(dplyr)
library(lme4)

##########################
## data read in and set up
##########################

setwd("/Users/acolby/Desktop/HBA 2024")
data = read.csv("data_mar_10.csv")

data$OT_SG = as.numeric(data$OT_SG)
data$log_OT = log(data$OT_SG) ## log OT 
data$age = as.numeric(data$age)
data$afb_observed = as.numeric(data$afb_observed)
data$n_births_observed = as.numeric(data$n_births_observed)

data$PID = as.factor(data$PID)
data$sex = as.factor(data$sex)
data$context = as.factor(data$context)
data$extraction_batch = as.factor(data$extraction_batch) 
data$assay_batch = as.factor(data$assay_batch)
data$assay_lot = as.factor(data$assay_lot)
data$BF = as.factor(data$BF)

#######################
## repeatability models
#######################

data$age.z=scale(data$age)

## age scaled - both sexes

rep.1 <- bf(log_OT | cens(censored) ~ 1 + s(age.z, by = sex) + BF + context + (1|PID) + (1|extraction_batch))

rep.model.1 <- brm(rep.1, data = data, 
                   prior = c(
                     prior(normal(0,2), class = "Intercept"),
                     prior(normal(0,1), class = "b"),
                     prior(exponential(2), class = "sd")),
                   chains=3, cores=3, warmup=500, iter=3000,
                   control = list(adapt_delta = 0.999))

saveRDS(rep.model.1, "repeatability.model.RDS")

##########
## results
##########

summary(rep.model.1)
conditional_effects(rep.model.1)
plot(rep.model.1)

###########
## variance
###########

## posterior samples
post.rep.model.1 <- as.data.frame(as.matrix(rep.model.1))
## PID
sd_PID_1 <- post.rep.model.1$sd_PID__Intercept
variance_PID_1 = sd_PID_1^2
## residual
sd_residual_1 <- post.rep.model.1$sigma
variance_residual_1 = sd_residual_1^2
## extraction batch
sd_exba_1 = post.rep.model.1$sd_extraction_batch__Intercept
variance_extraction_1 = sd_exba_1^2
## total variance
variance_total_1 <- variance_PID_1 + variance_residual_1 + variance_extraction_1 
## repeatability
repeatability_1 <- variance_PID_1 / variance_total_1
median(repeatability_1) ## median = 0.3448974
mean(repeatability_1) ## mean = 0.0.3440946

################
## HPD Intervals
################

## confidence interval
confidence_level <- 0.95
# sort the posterior variances
sorted_variances_1 <- sort(repeatability_1)
# calculate the cumulative probability
cumulative_probs_1 <- cumsum(sorted_variances_1) / sum(sorted_variances_1)
# find the indices corresponding to the lower and upper bounds
lower_bound_index_1 <- which(cumulative_probs_1 >= (1 - confidence_level)/2)[1]
upper_bound_index_1 <- which(cumulative_probs_1 >= 1 - (1 - confidence_level)/2)[1]
# extract the lower and upper bounds
hpd_lower_1 <- sorted_variances_1[lower_bound_index_1]
hpd_upper_1 <- sorted_variances_1[upper_bound_index_1]
# display the result
cat("HPD Interval:", hpd_lower_1, "to", hpd_upper_1, "\n") ## HPD Interval: 0.2263581 to 0.4902073 

############################################
## estimates and standard deviations for PID
############################################

randomeffects = ranef(rep.model.1)
PID = as.data.frame(randomeffects$PID)
PID$PID <- rownames(PID)
PID <- PID[, c("PID", names(PID)[-ncol(PID)])]
colnames(PID) <- c("PID", "meanOT", "sdOT", "2.5Q", "97.5Q")

new_df <- merge(PID, data, by = "PID", all = TRUE)
new_df$meanOT = as.numeric(new_df$meanOT)
new_df$sdOT = as.numeric(new_df$sdOT)
no.dup.df =  new_df %>%
  distinct(PID, .keep_all = TRUE)

##################
## fertility model
##################

fert.data <- no.dup.df %>%
  filter(!is.na(n_births_observed))

## males and females

fert.mod.1 = bf(scale(n_births_observed) ~ 1 + s(scale(age)) + me(meanOT, sdOT) + sex)

fert.model.1 <- brm(fert.mod.1, data = fert.data, 
                    prior = c(
                      prior(normal(0,2), class = "Intercept"),
                      prior(normal(0,1), class = "b")),
                    chains=3, cores=3, warmup=500, iter=3000,
                    control = list(adapt_delta = 0.999))
saveRDS(fert.model.1, "fertility.model.RDS")

summary(fert.model.1)
conditional_effects(fert.model.1)
plot(fert.model.1)
posterior_samples_fert = as.data.frame(as.matrix(fert.model.1))
sum_positive <- sum(posterior_samples_fert$bsp_memeanOTsdOT > 0)
sum_positive/7500 ## 0.3242667

############
## afr model 
############

afr.data <- no.dup.dfa %>%
  filter(!is.na(afb_observed))

## males and females

afr.mod.1 = bf(scale(afb_observed) ~ 1 + me(meanOT, sdOT) + sex)

afr.model.1 <- brm(afr.mod.1, data = afr.data, 
                   prior = c(
                     prior(normal(0,2), class = "Intercept"),
                     prior(normal(0,1), class = "b")),
                   chains=3, cores=3, warmup=500, iter=3500,
                   control = list(adapt_delta = 0.999))
saveRDS(afr.model.1, "afr.model.RDS")

summary(afr.model.1)
conditional_effects(afr.model.1)
plot(afr.model.1)
posterior_samples_afr = as.data.frame(as.matrix(afr.model.1))
sum_negative <- sum(posterior_samples_afr$bsp_memeanOTsdOT < 0)
sum_negative/9000 ## 0.06711111

#####################
## extraversion model 
#####################

extra.data <- no.dup.df %>%
  filter(!is.na(extraversionM))

extra.mod.1 = bf(scale(extraversionM) ~ 1 + me(meanOT, sdOT))

extra.model.1 <- brm(extra.mod.1, data = extra.data, 
                   prior = c(
                     prior(normal(0,2), class = "Intercept"),
                     prior(normal(0,1), class = "b")),
                   chains=3, cores=3, warmup=500, iter=3000,
                   control = list(adapt_delta = 0.999))
saveRDS(extra.model.1, "extraversion.model.RDS")

summary(extra.model.1)
conditional_effects(extra.model.1)
plot(extra.model.1)
posterior_samples_extra = as.data.frame(as.matrix(extra.model.1))
sum_positive_extra <- sum(posterior_samples_extra$bsp_memeanOTsdOT > 0)
sum_positive_extra/7500 ## 0.4868


#################
## openness model 
#################

open.data <- no.dup.df %>%
  filter(!is.na(opennessM))

open.mod.1 = bf(scale(opennessM) ~ 1 + me(meanOT, sdOT))

open.model.1 <- brm(open.mod.1, data = open.data, 
                     prior = c(
                       prior(normal(0,2), class = "Intercept"),
                       prior(normal(0,1), class = "b")),
                     chains=3, cores=3, warmup=500, iter=3000,
                     control = list(adapt_delta = 0.999))
saveRDS(open.model.1, "openness.model.RDS")

summary(open.model.1)
conditional_effects(open.model.1)
plot(open.model.1)
posterior_samples_open = as.data.frame(as.matrix(open.model.1))
sum_positive_open <- sum(posterior_samples_open$bsp_memeanOTsdOT > 0)
sum_positive_open/7500 ## 0.5270667
