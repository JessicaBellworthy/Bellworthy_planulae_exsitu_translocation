library(tidyverse)
library(lubridate)
library(ggplot2)  
library(phytotools)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
rm(list=ls())

#calculating RLC parameters using phytotools (pratt 1980)
setwd("G:/Shared drives/Jessica Bellworthy/R/RapidLightCurves/TranslocationStyloSpat")
d <- read.csv('Translocation_ETR_alldata.csv')
View(d)

rlc.data <- r[c(1,2,3)] # Subset data and ensure that column order is: PAR - ETR - ID/COLONY
rlc.data<- setNames(rlc.data, c("par","etr", "id")) # Set column names
rlc.data$id <- as.factor(rlc.data$id)
head(rlc.data)

rlc.data$etr <- na_if(rlc.data$etr, 0) # Change any zeros in etr to NA
View(rlc.data)

ncurves <- length(unique(rlc.data$id)) # number of unique ids in the data 
ids <- unique(rlc.data$id) # store the unique ids 


rlc.parameters <- data.frame(
  id = ids, 
  alpha = 0, 
  beta = 0, 
  ETRmax = 0, 
  Ek = 0, 
  ps = 0
)


for (i in 1:ncurves){
  
  temp.id = ids[i] # extract the id of the curve to be fitted
  
  print(paste("Now fitting curve ", as.character(temp.id))) # to keep track what's happening if the data has many curves
  
  temp.rlc.data <- rlc.data[rlc.data$id==temp.id,] # extract the the data of a single curve into a temporary variable
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") # for more options and explanation see package phytotools manual
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  # store the parameters
  rlc.parameters$id[i] <- temp.id
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g.Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237). 
  # Note that the equation depends on the model fitted, the code below applies only to the PGH model! 
  # Model equations are documented in the phytotools package code examples (and in the original papers): https://cran.r-project.org/web/packages/phytotools/phytotools.pdf
  
  ETRmax = ps.rlc*(alpha.rlc/(alpha.rlc + beta.rlc))*(beta.rlc/(alpha.rlc+beta.rlc))^(beta.rlc/alpha.rlc)
  Ek = ETRmax/alpha.rlc 
  
  # store the variables
  rlc.parameters$ETRmax[i] <- ETRmax
  rlc.parameters$Ek[i] <- Ek
  
  #plotting the curve and fitted model into a tiff file. By default the file name is the id of the curve. 
  tiff(file=paste0(temp.id, ".tiff"), compression="lzw")
  
  # plot the data, 
  plot(x=PAR, y=ETR, main=temp.id) 
  
  # plot the model fit
  with(fit, {
    P <- ps.rlc*(1-exp(-1*alpha.rlc*PAR/ps.rlc))*exp(-1*beta.rlc*PAR/ps.rlc) # the PGH model equation
    lines(PAR,P)
  }
  ) # end of with
  dev.off() #close the plotting devide. if this is not done, the next run of the loop will override the plot. 
  
}

# now the data frame rlc.parameters contains the fitted values for each curve. Tiff plots should be in current working directory. 
rlc.parameters

warnings()

write.csv(file="rlc.parameters.insitu_4months_StyloSpat_recalc.csv",rlc.parameters)


#### check statistical differences
library(lme4)
library(nlme)
library(lmerTest)
library(MASS)
library(car)
library(predictmeans)
library(ggpubr)

###manually copy the FVFM from the raw data to the rlc parameter data file and metadata

d <- read_csv('rlc.parameters.TranslocationStyloSpat_fvfm.csv')

d$age <- as.factor(d$age)
d$plate <- as.factor(d$plate)

gghistogram(d, x = "ETRmax", rug = TRUE, fill = "treatment", bins = 10)
etrMAX <- lmer(ETRmax~age*treatment + (1|plate), data=d)
summary(etrMAX)
anova(etrMAX) 
# No sig diff p = 0.5074

#normality assumptions - OK
shapiro.test(d$ETRmax)
leveneTest(ETRmax~depth,d=d)


###alpha
Alpha <- lm(alpha~depth, data=d2) # with the anomaly removed p = 0.7109
anova(Alpha)
# No sig diff p = 0.321

#summary(glht(Alpha, linfct=mcp(treatment="Tukey")))

#normality assumptions -
shapiro.test(d2$alpha) #OK - with anomaly removed
leveneTest(alpha~depth,d=d) #OK

###eK
head(d)
eK <- lm(Ek~depth, data=d)
summary(eK)
anova(eK)
# No sig diff p = 0.538
#summary(glht(eK, linfct=mcp(treatment="Tukey")))

#normality assumptions -ok
shapiro.test(d$Ek)
leveneTest(Ek~depth,d=d)

###FV/FM
head(d)
fv <- lm(fvfm~depth, data=d)
anova(fv)
# No sig diff p = 0.132
#summary(glht(fv, linfct=mcp(treatment="Tukey")))

#normality assumptions
shapiro.test(d$fvfm) ### not normally distributed, two low deep outliers
leveneTest(fvfm~depth,d=d)
qqnorm(d$fvfm)


###beta
head(d)
b <- lm(beta~depth, data=d)
anova(b)
# No sig diff p = 0.6228


#normality assumptions
shapiro.test(d$beta) ### not normally distributed
leveneTest(beta~depth,d=d)
qqnorm(d$fvfm)


## summary of results 
library(plotrix)

Sum_all <- d %>% 
  group_by(age, treatment) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), n = sum(!is.na(.)), max(., na.rm = TRUE), min(., na.rm = TRUE),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), alpha:fvfm)

View(Sum_all)
write.csv(file="rlc.summarytable.StyloSpat.csv", Sum_all)


#### Graphics in Box plots #########
r = d

r$treatment = factor(r$treatment, levels = c("SD", "SS"))
#ETRmax
p.etr = ggplot(r, aes(y = ETRmax, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
 # facet_wrap(facets = "age")+
  theme_classic()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values = c('#3CBB75FF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.etr

p.alpha= ggplot(r, aes(y = alpha, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
 # facet_wrap(facets = "age")+
  theme_classic()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values = c('#3CBB75FF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.alpha

p.beta= ggplot(r, aes(y = beta, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  #facet_wrap(facets = "age")+
  theme_classic()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values = c('#3CBB75FF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.beta

p.eK= ggplot(r, aes(y = Ek, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  #facet_wrap(facets = "age")+
  theme_classic()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values = c('#3CBB75FF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.eK

p.fvfm= ggplot(r, aes(y = fvfm, x = treatment)) +
  geom_boxplot(outlier.shape = NA, aes(fill= treatment, alpha = 0.8), fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.25)+
  #facet_wrap(facets = "age")+
  theme_classic()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values = c('#3CBB75FF', '#FDE725FF'))+
  theme(axis.title.x = element_blank(), axis.text = element_text(colour = "black"), strip.background = element_rect(size = 0.5))+
  guides(fill = FALSE, alpha = FALSE)
p.fvfm


library(cowplot)

rlc.plots <- plot_grid(p.etr, p.alpha, p.beta, p.eK, p.fvfm, labels = c('A', 'B', 'C', 'D', 'E'), label_x = 0.12,
                       label_y = 0.985, label_size = 14,ncol = 1, align = "v", byrow = F, hjust =2)

rlc.plots
ggsave("PAM_StyloSpat.jpeg", plot = rlc.plots, width = 6, height = 19,dpi=300, 
       units = "cm")

ggsave("PAM_StyloSpat.pdf", plot = rlc.plots, width = 6, height = 19,dpi=300, 
       units = "cm")

