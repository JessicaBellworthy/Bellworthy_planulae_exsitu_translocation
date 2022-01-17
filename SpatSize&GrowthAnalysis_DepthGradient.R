# Depth Gradient Translocation Expt., Spat Size Analysis
library(lmerTest)
library(lme4)
library(ggplot2)
library(dplyr)
library("FSA")
library(cowplot)
library(car)
library(ggpubr)
library(ggpmisc)


setwd("G:/.shortcut-targets-by-id/1m-AP2D-7C4qAOnoLoy3tuMO7AG8cQtdI/Eilat_NSF-BSF April 2021/R")

## Analysis of Rachel's Bino photos for spat size and polyp number
d = read.csv("SpatCalyxSize_2020.csv")
View(d)
str(d)

mytheme = theme_classic()+
  theme(axis.text = element_text(colour = "black", size = 6), axis.title = element_text(size = 9))

d$age = as.numeric(d$age)
d$treatment = factor(d$treatment, levels = c("DD", "DS", "SD", "SS"))
d$colony = as.factor(d$colony)
d$polyps = as.numeric(d$polyps)


# Spat diameter with time by treatment
growth = ggplot(d, aes(x = age, y = (length_um/1000)))+
  geom_point(aes(colour = treatment))+
  facet_wrap(~ treatment)+
  geom_smooth(method='lm', lwd = 0.5, aes(colour = treatment), formula = y ~ x)+
  scale_colour_viridis_d()+
  mytheme+  
  labs(y = "Diameter (mm)", x = "Age (days)")+
  scale_x_continuous(limits = c(0, 60, 5))+
  scale_y_continuous(limits = c(-0.5, 3.5, 1))+
  guides(colour = FALSE)+
  stat_poly_eq(label.y = 0.9, aes(label = paste(..eq.label..)), coef.digits = 4, formula = y ~ x, size = 3)
growth

ggsave("SpatGrowth_StyloSpat2020.jpeg", plot = growth, width = 13, height = 13,dpi=300, 
       units = "cm")
ggsave("SpatGrowth_StyloSpat2020.pdf", plot = growth, width = 13, height = 13,dpi=300, 
       units = "cm")

# Calyx Diameter with age by treatment
calyx = ggplot(d, aes(x = age, y = (calyx_um/1000)))+
  geom_point(aes(colour = treatment), shape = 1)+
  facet_wrap(~ treatment)+
  geom_smooth(method='lm', lwd = 0.5, aes(colour = treatment), formula = y ~ x)+
  scale_colour_viridis_d()+
  mytheme+  
  labs(y = "Calyx Diameter (mm)", x = "Age (days)")+
  scale_x_continuous(limits = c(0, 60, 5))+
  scale_y_continuous(limits = c(-0.5, 3.5, 1))+
  guides(colour = FALSE)+
  stat_poly_eq(label.y = 0.5, aes(label = paste(..eq.label..)), coef.digits = 3, formula = y ~ x, size = 3)
calyx

ggsave("CalyxDiameter_StyloSpat2020.jpeg", plot = calyx, width = 13, height = 13,dpi=300, 
       units = "cm")
ggsave("CalyxDiamter_StyloSpat2020.pdf", plot = calyx, width = 13, height = 13,dpi=300, 
       units = "cm")

# are the slopes significantly different from zero?
DD = subset(d, treatment == "DD")
DS = subset(d, treatment == "DS")
SS = subset(d, treatment == "SS")
SD = subset(d, treatment == "SD")
summary(lm(data = DD, calyx_um ~ age))
summary(lm(data = DS, calyx_um ~ age))
summary(lm(data = SS, calyx_um ~ age))
summary(lm(data = SD, calyx_um ~ age))

model1 = lm(data = d, calyx_um ~ treatment)
anova(model1)
dunnTest(data = d, calyx_um ~ treatment)

#Comparison          Z      P.unadj        P.adj
#1    DD - DS  0.8755594 3.812696e-01 0.762539
#2    DD - SD -3.3090257 9.362125e-04 0.002808637
#3    DS - SD -5.3256491 1.005932e-07 0.0000005029659
#4    DD - SS -3.6289572 2.845684e-04 0.001138274
#5    DS - SS -5.9347905 2.942210e-09 0.00000001765326
#6    SD - SS -0.2540216 7.994788e-01 0.7994788e-01

# are the intercepts significantly different from each other?



# number of polps with age
polyps = ggplot(d, aes(x = age, y = polyps))+
  geom_boxplot(aes(group = age, fill = treatment), outlier.size = 0.1, fatten = 0.5, alpha = 0.9, lwd = 0.2)+
  facet_wrap(~ treatment)+
#  geom_smooth(method='lm', lwd = 0.5, aes(colour = treatment), formula = y ~ x)+
  scale_fill_viridis_d()+
  mytheme+ theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA))+
  labs(y = "No. of Polyps", y.position = "right", x = "Age (days)")+
  scale_x_continuous(limits = c(0, 60, 5))+
  scale_y_continuous(limits = c(0, 16, 2), position = "right")+
  guides(colour = FALSE, fill = FALSE)
polyps

ggsave("SpatPolypNo_StyloSpat2020.png",bg = "transparent",  plot = polyps, width = 13, height = 13,dpi=300, 
       units = "cm")
ggsave("SpatPolypNo_StyloSpat2020.pdf", plot = polyps, width = 13, height = 13,dpi=300, 
       units = "cm")


# combine polyp and size graphs
comb = ggplot(d, aes(x = age))+
  geom_point(aes(colour = treatment, y = (length_um/1000)))+
 geom_smooth(d, aes(y = (length_um/1000), x = age, color = treatment), method='lm', lwd = 0.5, formula = y~x)+
  stat_poly_eq(label.y = 0.9, aes(label = paste(..eq.label..)), coef.digits = 4, formula = y ~ x, size = 3)+
  labs(y = "Diameter (mm)", x = "Age (days)")+
  scale_x_continuous(limits = c(0, 60, 5))
 # scale_y_continuous(name = "Diameter (mm)", sec_axis = sec_axis(trans=~.*10, name = "No. of Polyps"))+
  #geom_boxplot(aes(y = polyps, group = age, fill = treatment), outlier.size = 0.1, fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  #facet_wrap(~ treatment)+
  #scale_fill_viridis_d()+
  #mytheme+  
  #labs(y = "No. of Polyps", x = "Age (days)")+
  #guides(colour = FALSE, fill = FALSE)
comb

ggsave("SpatPolypNo_StyloSpat2020.jpeg", plot = polyps, width = 14, height = 13,dpi=300, 
       units = "cm")
ggsave("SpatPolypNo_StyloSpat2020.pdf", plot = polyps, width = 14, height = 13,dpi=300, 
       units = "cm")

