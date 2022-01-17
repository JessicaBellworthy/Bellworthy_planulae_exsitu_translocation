# Settlement and Mortality Ex Situ 2020 translocation expt.

library(dplyr)
library(ggplot2)
library(cowplot)


setwd("G:/.shortcut-targets-by-id/1m-AP2D-7C4qAOnoLoy3tuMO7AG8cQtdI/Eilat_NSF-BSF April 2021/R")
d = read.csv("Settlement&Survival2020_JB.csv")
str(d)

d$age_days = as.factor(d$age_days)
d$treatment = factor(d$treatment, levels = c("DD", "SD", "DS", "SS"))
d$well_no = as.factor(d$well_no)
d$swim = as.numeric(d$swim)
d$set = as.numeric(d$set)
d$dead = as.numeric(d$dead)


Sum_all <- d %>% 
  group_by(treatment, age_days) %>% 
  summarise_each(funs(median(., na.rm=TRUE), mean(., na.rm = TRUE), max(., na.rm = TRUE), sd(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), swim:pc_dead)

View(Sum_all)
write.csv(file="Set_Mort_Averages.csv", Sum_all)



# Subset df to have only wells with 5 or more planulae in
d_sub = subset(d, larvae_start_num >4)
View(d_sub)

Sum_sub <- d_sub %>% 
  group_by(treatment, age_days) %>% 
  summarise_each(funs(median(., na.rm=TRUE), mean(., na.rm = TRUE), max(., na.rm = TRUE), sd(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), swim:pc_dead)

View(Sum_sub)
write.csv(file="Set_Mort_Averages_5+.csv", Sum_sub)

mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black", size = 6), axis.title = element_text(size = 9))


# View settlement and mortality of all data faceted by age_days - all data, not subset
# Bar graphs from mean and se
set_all = ggplot(Sum_all, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  facet_wrap(facets = "age_days")+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()
set_all
ggsave("Settlement_ExSitu_StyloSpat2020_dayfacet.jpeg", plot = set_all, width = 15, height = 15,dpi=300, 
       units = "cm")


dead_all = ggplot(Sum_all, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  facet_wrap(facets = "age_days")+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()
dead_all
ggsave("Mortality_ExSitu_StyloSpat2020_dayfacet.jpeg", plot = dead_all, width = 15, height = 15,dpi=300, 
       units = "cm")



# View settlement and mortality of all data faceted by age_days but only for wells with 5 or more planulae - i.e. the subset data
set_all_sub = ggplot(Sum_sub, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  facet_wrap(facets = "age_days")+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()

set_all_sub
ggsave("Settlement_ExSitu_StyloSpat2020_dayfacet_+5.jpeg", plot = set_all_sub, width = 15, height = 15,dpi=300, 
       units = "cm")

dead_all_sub = ggplot(Sum_sub, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  facet_wrap(facets = "age_days")+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()
dead_all_sub
ggsave("Mortality_ExSitu_StyloSpat2020_dayfacet_+5.jpeg", plot = dead_all_sub, width = 15, height = 15,dpi=300, 
       units = "cm")




# Settlement and mortality at 8 days from 2020 data
d8= subset(Sum_all, age_days ==8)
View(d8)

set8 = ggplot(d8, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()

set8

dead8 = ggplot(d8, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()
dead8


set.plots8 <- plot_grid(set8, dead8, labels = c('A', 'B'),
                       label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Settlement_ExSitu_StyloSpat2020_8days.jpeg", plot = set.plots8, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_8days.pdf", plot = set.plots8, width = 12, height = 6,dpi=300, 
       units = "cm")

# Settlement and mortality at 8 days from 2020 data for wells with 5 or more planulae
d8_sub= subset(Sum_sub, age_days ==8)
View(d8_sub)

set8_sub = ggplot(d8_sub, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()

set8_sub

dead8_sub = ggplot(d8_sub, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()
dead8_sub


set.plots8_sub <- plot_grid(set8_sub, dead8_sub, labels = c('A', 'B'),
                        label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Settlement_ExSitu_StyloSpat2020_8days_+5.jpeg", plot = set.plots8_sub, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_8days_+5.pdf", plot = set.plots8_sub, width = 12, height = 6,dpi=300, 
       units = "cm")


# Settlement and mortality at 15 days from 2020 data (dont have all treatments for 14 days)
d15= subset(Sum_all, age_days ==15)
View(d15)

set15 = ggplot(d15, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()

set15

dead15 = ggplot(d15, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 85))+
  scale_fill_viridis_d()
dead15


set.plots15 <- plot_grid(set15, dead15, labels = c('A', 'B'),
                       label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Settlement_ExSitu_StyloSpat2020_15days.jpeg", plot = set.plots15, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_15days.pdf", plot = set.plots15, width = 12, height = 6,dpi=300, 
       units = "cm")


# Settlement and mortality at 20 days from 2020 data
d20= subset(Sum_all, age_days ==20)
View(d20)


p = ggplot(d20, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 80))+
  scale_fill_viridis_d()
p

p2 = ggplot(d20, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 80))+
  scale_fill_viridis_d()
p2


set.plots20 <- plot_grid(p, p2, labels = c('A', 'B'),
                       label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Settlement_ExSitu_StyloSpat2020_20days.jpeg", plot = set.plots20, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_20days.pdf", plot = set.plots20, width = 12, height = 6,dpi=300, 
       units = "cm")


# # Settlement and mortality at 20 days from 2020 data with wells only with 5 + planulae in
d20_sub= subset(Sum_sub, age_days ==20)
View(d20_sub)


p = ggplot(d20_sub, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_viridis_d()
p

p2 = ggplot(d20_sub, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_viridis_d()
p2


set.plots20_sub <- plot_grid(p, p2, labels = c('A', 'B'),
                         label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Settlement_ExSitu_StyloSpat2020_20days_+5.jpeg", plot = set.plots20_sub, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_20days_+5.pdf", plot = set.plots20_sub, width = 12, height = 6,dpi=300, 
       units = "cm")


### 2020 Set and mortality panel 
all.plots <- plot_grid(set.plots8, set.plots15, set.plots20,# labels = c('A', 'B', 'C'),
                       label_y = 0.985, label_size = 11,ncol = 1, align = "h", byrow = F)


ggsave("Settlement_ExSitu_StyloSpat2020_8,15,20days.jpeg", plot = all.plots, width = 12, height = 15,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_8,15,20days.pdf", plot = all.plots, width = 12, height = 15,dpi=300, 
       units = "cm")



# Find the age_days with the highest number of settled spat for each well_no
d_sorted = d[order(d$well_no, -d$set),]  # sort the df by set descending 
View(d_sorted)
ans <- d_sorted[!duplicated(d_sorted$well_no),] # keep only rows with max no. settled
View(ans)
str(ans)
ans$age_days = as.numeric(ans$age_days)
mean(ans$age_days) # 9.37
median(ans$age_days) # 9.0
aggregate(age_days ~ treatment, ans, mean) ##  DD 7.428571 SD 9.612903, DS 8.818182, SS 9.702703

# same as above but with the subset of data +5
d_sort_sub = d_sub[order(d_sub$well_no, -d_sub$set),]  # sort the df by set descending 
View(d_sort_sub)
ans_sub <- d_sort_sub[!duplicated(d_sort_sub$well_no),] # keep only rows with max no. settled
View(ans_sub)
ans_sub$age_days = as.numeric(ans_sub$age_days)
mean(ans_sub$age_days, na.rm = TRUE) # 9.08
median(ans_sub$age_days, na.rm = TRUE) # 8.5
aggregate(age_days ~ treatment, ans_sub, mean)  # DD 9.00, SD 9.40, DS 5.00, SS 9.25



# # Settlement and mortality at 9 days (time of max settlement) from 2020 data with wells only with 5 + planulae in
d9_sub= subset(Sum_sub, age_days ==9)
View(d9_sub)


p = ggplot(d9_sub, aes(x= treatment, y = pc_set_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_set_mean - pc_set_se, ymax = pc_set_mean + pc_set_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean settlement %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_viridis_d()
p

p2 = ggplot(d9_sub, aes(x= treatment, y = pc_dead_mean, fill = treatment))+
  geom_col()+
  geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mean mortality %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_viridis_d()
p2


set.plots9_sub <- plot_grid(p, p2, labels = c('A', 'B'),
                             label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Settlement_ExSitu_StyloSpat2020_9days_+5.jpeg", plot = set.plots9_sub, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Settlement_ExSitu_StyloSpat2020_9days_+5.pdf", plot = set.plots9_sub, width = 12, height = 6,dpi=300, 
       units = "cm")

# Already low n in DS and DD
# Decision to keep all data, not subset


#####################
# Kaplan Meier Curves

library(coxme)
library(purrr)
library(tidyr)
library(tidyverse)
library(survminer)
library(survival)
library(ggsci)
library(splitstackshape)
library(MuMIn)


################### Non filtered data, all data ############### 
k = read.csv("KM2020_new.csv")
#k = subset(k, age_days >4)
#k = subset(k, age_days <21)
View(k)
K = expandRows(k, "count")  #Expands row by count and Removes rows where count = 0
View(K)


### Fit KM curves
# Settlement
fitK = survfit(Surv(age_days, new_set_status)~ treatment, data = K)
summary(fitK)

names(fitK$strata) <- gsub("treatment=", "", names(fitK$strata))

setplot <- ggsurvplot(fitK, 
                      fun = "event",
                      palette = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
                      pval = TRUE, pval.method = TRUE, 
                      legend = c(0.85, 0.2), legend.title = "", xlab="Days", ylab = "Settlement probability", legend.labs = c("DD","DS", "SD", "SS"))

setplot

#Pairwise comparison between treatments, fixed effects only, p value adj. method BH
res <- pairwise_survdiff(Surv(age_days, new_set_status) ~ treatment,
                         data= K)
res   
#   DD   DS   SD 
#DS 0.88 -    -   
# SD 0.93 0.88 -   
# SS 0.88 0.88 0.32

d = subset(d, age_days == 20)
View(d)
set20 = lm(data = d, pc_set ~ treatment)
summary(set20)
anova(set20)

d$treatment = factor(d$treatment, levels = c("DD", "DS", "SD", "SS"))

p20 = ggplot(d, aes(x= treatment, y = pc_set, fill = treatment))+
  geom_boxplot(aes(fill = treatment), outlier.size = 0.2, fatten = 0.5, lwd = 0.2, alpha = 0.8)+
  geom_jitter(position = position_jitter(width = .25), size = 0.2)+
 # geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Settlement %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_viridis_d()
p20


# Mortality
fitM = survfit(Surv(age_days, new_dead_status)~ treatment, data = K)
summary(fitM)

names(fitM$strata) <- gsub("treatment=", "", names(fitM$strata))

mortplot <- ggsurvplot(fitM,
                       palette = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"),
                       pval = TRUE, pval.method = TRUE, 
                       legend = c(0.85, 0.8), legend.title = "", xlab="Days", ylab = "Survial Probability", legend.labs = c("DD","DS", "SD", "SS"))

mortplot

d = subset(d, age_days == 20)
View(d)
set20 = lm(data = d, pc_dead ~ treatment)
summary(set20)
anova(set20)

p20.d = ggplot(d, aes(x= treatment, y = pc_dead, fill = treatment))+
  geom_boxplot(aes(fill = treatment), outlier.size = 0.2, fatten = 0.5, lwd = 0.2, alpha = 0.8)+
  geom_jitter(position = position_jitter(width = .25), size = 0.2)+
  # geom_errorbar(aes(ymin = pc_dead_mean - pc_dead_se, ymax = pc_dead_mean + pc_dead_se), width=0)+
  mytheme +  guides(fill = FALSE)+
  labs(y = "Mortality %")+
  scale_y_continuous(limits = c(0, 100))+
  scale_fill_viridis_d()
p20.d



#Pairwise comparison between treatments, fixed effects only, p value adj. method BH
res <- pairwise_survdiff(Surv(age_days, new_dead_status) ~ treatment,
                         data= K)
res   
#   DD   DS   SD  
#DS 0.97 -    -   
#SD 0.97 0.97 -   
# SS 0.97 0.97 0.97


#### These plots to be used for inital paper submission and results

day.plots20 <- plot_grid(p20, p20.d, labels = c('A', 'B'),
                            label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Set.Dead_ExSitu_StyloSpat2020_20days_violin.jpeg", plot = day.plots20, width = 12, height = 6,dpi=300, 
       units = "cm")
ggsave("Set.Dead_ExSitu_StyloSpat2020_20days_violin.pdf", plot = day.plots20, width = 12, height = 6,dpi=300, 
       units = "cm")


### Composite figure 2
SM.plots <- plot_grid(setplot, p20,  mortplot,p20.d, labels = c('A', 'C', 'B', 'D'),
                         label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)
SM.plots
ggsave("Figure2_SettlementMortality.jpeg", plot = SM.plots, width = 12, height = 12,dpi=300, 
       units = "cm")
ggsave("Figure2_SettlementMortality.pdf", plot = SM.plots, width = 12, height = 12,dpi=300, 
       units = "cm")

############### KM with all data are not significant but try different models ######

# model selection for mortality 
# Starting with most complex

# Fixed treatment effect plus random effect of well plate and start date
mix_eff = coxme(Surv(age_days, status_dead) ~ treatment+(1|well_no)+(1|start_date),data=K)
summary(mix_eff)
#Random effects
#Group      Variable  Std Dev      Variance    
#well_no    Intercept 0.8282528564 0.6860027941
#start_date Intercept 0.0196052073 0.0003843642

# Fixed treatment effect plus random effect of larvae start number
mix_eff = coxme(Surv(age_days, new_dead_status) ~ treatment+(1|larvae_start_num),data=K)
summary(mix_eff)
#Random effects
#Group            Variable  Std Dev   Variance 
#larvae_start_num Intercept 0.3224603 0.1039806


# Fixed treatment effect plus random effect of well plate 
mix_eff1<- coxme(Surv(age_days, status_dead) ~ treatment+(1|well_no),data=K)
summary(mix_eff1)
# Random effects
#Group   Variable  Std Dev   Variance 
#well_no Intercept 0.8271708 0.6842115

# Fixed treatment effect only
fix_eff1<- coxph(Surv(age_days, status_dead) ~ treatment,data=K)
summary(fix_eff1)
#              coef exp(coef) se(coef)     z Pr(>|z|)  
# treatmentDS 0.3963    1.4863   0.2008 1.973   0.0485 *
# treatmentSD 0.2109    1.2348   0.1713 1.231   0.2183  
# treatmentSS 0.1029    1.1084   0.1713 0.601   0.5479 
ggforest(fix_eff1, data = K) # only available for fixed effect models


anova(mix_eff, mix_eff1, fix_eff1)
#loglik    Chisq Df P(>|Chi|)    
# 1 -23101                          
# 2 -23101   0.0252  1    0.8738    
# 3 -23332 462.8217  1    <2e-16 ***
AIC(mix_eff, mix_eff1, fix_eff1)
#             df      AIC
#mix_eff  67.65035 46070.24
#mix_eff1 67.62979 46070.28
#fix_eff1  3.00000 46669.79

# The mixed effect model with only well_no as random effect is the best model i.e. lowest AIC and the start_date effect was not significant

#Pairwise comparison between treatments, fixed effects only
res <- pairwise_survdiff(Surv(age_days, status_dead) ~ treatment,
                         data= K)
res #   DD    DS    SD   
# DS 0.064 -     -    
# SD 0.265 0.121 -    
# SS 0.494 0.028 0.028




# model selection for time to settlement
# Starting with most complex

# Fixed treatment effect plus random effect of well plate and start date
mix_eff = coxme(Surv(age_days, status_set) ~ treatment+(1|well_no)+(1|start_date),data=K)
summary(mix_eff)
#Random effects
#Group      Variable  Std Dev      Variance    
#well_no    Intercept 0.5176720 0.2679843
#start_date Intercept 0.1261400 0.0159113

# Fixed treatment effect plus random effect of well plate 
mix_eff1<- coxme(Surv(age_days, status_set) ~ treatment+(1|well_no),data=K)
summary(mix_eff1)
# Random effects
#Group   Variable  Std Dev   Variance 
#well_no Intercept 0.5250496 0.2756771

# Fixed treatment effect only
fix_eff1<- coxph(Surv(age_days, status_set) ~ treatment,data=K)
summary(fix_eff1)
#              coef exp(coef) se(coef)     z Pr(>|z|)  
# treatmentDS -0.08049   0.92266  0.16228 -0.496     0.62  
# treatmentSD -0.73249   0.48071  0.13173 -5.561 2.69e-08 ***
# ttreatmentSS -0.61736   0.53937  0.13088 -4.717 2.39e-06 ***
ggforest(fix_eff1, data = K) # only available for fixed effect models


anova(mix_eff, mix_eff1, fix_eff1)
#loglik    Chisq Df P(>|Chi|)    
# 1 -16498                          
# 2 -16498   0.4101  1    0.5219    
# 3 -16649 302.0328  1    <2e-16 ***
AIC(mix_eff, mix_eff1, fix_eff1)
#             df      AIC
#mix_eff  64.25684 32915.65
#mix_eff1 64.37569 32915.71
#fix_eff1  3.00000 33303.75

# The mixed effect model with only well_no as random effect is the best model i.e. lowest AIC and the start_date effect was not significant

#Pairwise comparison between treatments, fixed effects only
res <- pairwise_survdiff(Surv(age_days, status_set) ~ treatment,
                         data= K)
res #   DD    DS    SD   
# DS 0.686    -     -    
# SD 1.9e-07 2.7e-08 -  
# SS 1.6e-06 2.5e-07 0.017



## Next step to get median time to settlement for each treatment as a single value
## https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

#   https://github.com/ronenliberman/CC-disrupt-octocoral-reproductive-synchrony/blob/master/Survival%20stats%20for%20publication.R





### Kaplan Meier curves with subset data, i.e. more than 5 planulae per well

k = read.csv("KM2020_new.csv")
k_sub = subset(k, larvae_start_num >4)
k_sub = subset(k_sub, age_days >4)
k_sub = subset(k_sub, age_days <21)
View(k_sub)
K_sub = expandRows(k_sub, "count")  #Expands row by count and Removes rows where count = 0

# Fit KM curves

# Settlement Probability
fitK = survfit(Surv(age_days, new_set_status)~ treatment, data = K_sub)
fitK
summary(fitK)

setplot <- ggsurvplot(fitK, 
                      fun = "event",
                      pval = TRUE, pval.method = TRUE, 
                      legend = "bottom", xlab="Days", ylab = "Settlement", legend.labs = c("DD","DS", "SD", "SS")+ 
                        scale_fill_viridis_d(direction = -1),
                      surv.median.line = "hv")
# scale_fill_manual(values = c("#FDE725FF","#2C728EFF")))

setplot

# Mortality
fitM = survfit(Surv(age_days, new_dead_status)~ treatment, data = K_sub)
fitM
summary(fitM)

mortplot <- ggsurvplot(fitM,
                       pval = TRUE, pval.method = TRUE, 
                       legend = "bottom", xlab="Days", ylab = "Survial Probability", legend.labs = c("DD","DS", "SD", "SS")+ 
                         scale_fill_viridis_d(direction = -1),
                       surv.median.line = "hv")

mortplot

fitM_sub = survfit(Surv(age_days, status_dead)~ treatment, data = K_sub)
summary(fitM_sub)

mortplot_sub <- ggsurvplot(fitM_sub, K_sub,
                           pval = TRUE, legend = "bottom",
                           xlab="Age days", ylab = "Survival Probability", legend.labs = c("DD","DS", "SD", "SS")+ scale_fill_viridis_d())

mortplot_sub




################### Data exploration ############### 
# subset data to include wells with only 4 - 11 planulae
k = read.csv("KM2020_new.csv")
k = subset(k, larvae_start_num >3.9)
k = subset(k, larvae_start_num <11.1)
View(k)
K = expandRows(k, "count")  #Expands row by count and Removes rows where count = 0
View(K)


### Fit KM curves
# Settlement
fitK = survfit(Surv(age_days, new_set_status)~ treatment, data = K)
summary(fitK)

names(fitK$strata) <- gsub("treatment=", "", names(fitK$strata))

setplot <- ggsurvplot(fitK, 
                      fun = "event",
                      pval = TRUE, pval.method = TRUE, 
                      legend = "bottom", legend.title = "", xlab="Days", ylab = "Settlement probability", legend.labs = c("DD","DS", "SD", "SS")+ 
                        scale_fill_viridis_d(direction = -1))

setplot #not sig. in word file

#Pairwise comparison between treatments, fixed effects only, p value adj. method BH
res <- pairwise_survdiff(Surv(age_days, new_set_status) ~ treatment,
                         data= K)
res      # not significant


# Mortality
fitM = survfit(Surv(age_days, new_dead_status)~ treatment, data = K)
summary(fitM)

names(fitM$strata) <- gsub("treatment=", "", names(fitM$strata))

mortplot <- ggsurvplot(fitM,
                       pval = TRUE, pval.method = TRUE, 
                       legend = "bottom", legend.title = "", xlab="Days", ylab = "Survial Probability", legend.labs = c("DD","DS", "SD", "SS")+ 
                         scale_fill_viridis_d(direction = -1))

mortplot #not sig. in word file


#Pairwise comparison between treatments, fixed effects only, p value adj. method BH
res <- pairwise_survdiff(Surv(age_days, new_dead_status) ~ treatment,
                         data= K)
res   # not significant



