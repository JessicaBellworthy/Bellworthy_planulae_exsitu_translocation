library(lmerTest)
library(lme4)
library(ggplot2)
library(dplyr)
library("FSA")
library(cowplot)
library(rstatix)
library(MuMIn)
library(ggpubr)
library(magick)


####  Larval Volume   ######
setwd("G:/.shortcut-targets-by-id/1m-AP2D-7C4qAOnoLoy3tuMO7AG8cQtdI/Eilat_NSF-BSF April 2021/R")
vol = read.csv("LarvalVolume.csv")
View(vol)

sh = subset(vol, depth == "S" & vol_mm3 <0.7) #remove 1
de = subset(vol, depth == "D"& vol_mm3 <0.45)   # remove 2


shapiro.test(sh$vol_mm3)# OK, with >0.7 removed
hist(sh$vol_mm3)
shapiro.test(de$vol_mm3) # OK, with >0.45 removed
hist(de$vol_mm3)
levene_test(vol, vol_mm3 ~ depth)


vol$day = as.factor(vol$day)
vol$depth = factor(vol$depth, levels = c("S", "D"))
vol$vol_mm3 = as.numeric(vol$vol_mm3)
str(vol)

vol = rbind(de, sh)
View(vol)
vol$depth = factor(vol$depth, levels = c("S", "D"))

Sum_all <- vol %>% 
  group_by(depth) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), n = sum(!is.na(.)),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), vol_mm3)

View(Sum_all)
write.csv(file="LarvalVolumeAverages.csv", Sum_all)


volume = lmer(vol_mm3~depth + (1|day), data=vol)
anova(volume) # F = 60.52 DenDF = 19.444 p = 2.167e-07 ***
summary(volume)
r.squaredGLMM(volume) #          R2m       R2c
                      #    0.4265383 0.4655294

volume_lm = lm(vol_mm3~depth, data=vol)
anova(volume_lm) # F = 60.52 DenDF = 19.444 p = 2.167e-07 ***  p = 0.000000216
anova(volume, volume_lm)


mytheme = theme_classic()+
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour = "black", size = 8),
        axis.text.y = element_text(colour = "black", size = 6), axis.title = element_text(size = 9))



mynames = c("Shallow", "Deep")

stat.test.vol <- vol %>%
  t_test(vol_mm3 ~ depth) %>%
  add_significance() %>%
  add_xy_position

vol.plot = ggplot(vol, aes(y = vol_mm3, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.2)+
  labs(y= ~Larval ~Volume ~(mm^3))+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#FDE725FF","#2C728EFF"))+
 scale_y_continuous(limits = c(0, 0.66), breaks = seq(0,0.6,0.2))+
  stat_pvalue_manual(stat.test.vol, label = "p.signif", tip.length = 0.01, y.position = 0.65)+
  mytheme+
  guides(fill = FALSE)

vol.plot


ggsave("LarvalVolumeEllip_DepthGradient2021.jpeg", plot = vol.plot, width = 4.5, height = 6,dpi=300, 
       units = "cm")
ggsave("LarvalVolumeEllip_DepthGradient2021.pdf", plot = vol.plot, width = 4.5, height = 6,dpi=300, 
       units = "cm")


####### Larval Physiology #######
phys = read.csv("Physio_r.csv")

sh = subset(phys, depth == "S")
de = subset(phys, depth == "D"& Sample != "DPe")   # removes DPe, low protein and cell count

shapiro.test(sh$cell.sample)# OK
hist(sh$cell.sample)
shapiro.test(de$cell.sample) # OK
hist(de$cell.sample)
levene_test(phys, cell.sample ~ depth) # OK

phys = rbind(de, sh)
phys$depth = factor(phys$depth, levels = c("S", "D"))

str(phys)


Sum_all <- phys %>% 
  group_by(depth) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), n = sum(!is.na(.)),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), cell.sample: cell.protein)

View(Sum_all)
write.csv(file="PlanulaePhysioAverages.csv", Sum_all)

# cells / planula significantly different between depths
cell = lm(cell.sample/n_planulae~depth, data=phys)
anova(cell) # F = 8.2732 DF = 1 p = 0.02063 *

# protein/ planula not(!) significantly different between depths
pro = lm(ugProtein.sample/n_planulae~depth, data=phys)
anova(pro) # F = 5.1802 DF = 1 p = 0.0524

# cells / protein significantly different between depths
cell.p = lm(cell.protein/n_planulae~depth, data=phys)
anova(cell.p) # F = 20.25 DF = 1 p = 0.002002 **


# Cells/ larva
stat.test.cell <- phys %>%
  t_test(cell.sample ~ depth) %>%
  add_significance() %>%
  add_xy_position

cell.plot = ggplot(phys, aes(y = (cell.sample/n_planulae), x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
 geom_point(position = position_jitter(width = .25), size = 0.2)+
labs(y= ~symbiont ~cells ~larva^-1)+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#FDE725FF","#2C728EFF"))+
 scale_y_continuous(limits = c(0, 15000), breaks = seq(0,15000,5000))+
 stat_pvalue_manual(stat.test.cell, label = "p.signif", tip.length = 0.01, y.position = 15000)+
  mytheme+
  guides(fill = FALSE)

cell.plot

# Protein/ larva
stat.test.pro <- phys %>%
  t_test(ugProtein.sample ~ depth) %>%
  add_significance() %>%
  add_xy_position

pro.plot = ggplot(phys, aes(y = (ugProtein.sample/n_planulae), x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_point(position = position_jitter(width = .25), size = 0.2)+
  labs(y= ~µg ~protein ~larva^-1)+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#FDE725FF","#2C728EFF"))+
  scale_y_continuous(limits = c(0, 20), breaks = seq(0,20,5))+
  stat_pvalue_manual(stat.test.pro, label = "p.signif", tip.length = 0.01, y.position = 19)+
  mytheme+
  guides(fill = FALSE)

pro.plot

# cell/ protein
stat.test.cp <- phys %>%
  t_test(cell.protein ~ depth) %>%
  add_significance() %>%
  add_xy_position

cp.plot = ggplot(phys, aes(y = (cell.protein), x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_point(position = position_jitter(width = .25), size = 0.2)+
  labs(y= ~cells/~µg ~protein)+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#FDE725FF","#2C728EFF"))+
  scale_y_continuous(limits = c(0, 1200), breaks = seq(0,1200,400))+
  stat_pvalue_manual(stat.test.cp, label = "p.signif", tip.length = 0.01, y.position = 1200)+
  mytheme+
  guides(fill = FALSE)

cp.plot


phys.plots <- plot_grid(cell.plot, pro.plot, cp.plot, labels = c('A', 'B', 'C'),
                       label_y = 0.985, label_size = 11,ncol = 3, align = "h", byrow = F)

ggsave("LarvalPhysio_DepthGradient2021.jpeg", plot = phys.plots, width = 11, height = 4.5,dpi=300, 
       units = "cm")
ggsave("LarvalPhysio_DepthGradient2021.pdf", plot = phys.plots, width = 11, height = 4.5,dpi=300, 
       units = "cm")





###### Spat Size ######
size = read.csv("SpatSize.csv")
View(size)


Sum_size <-  size%>% 
  group_by(treatment, in_ex, age_days) %>% 
  summarise_each(funs(mean(., na.rm=TRUE), sd = sd, n = sum(!is.na(.)),
                      se = sd(., na.rm=TRUE)/sqrt(sum(!is.na(.)))), diameter)
write.csv(file="SpatSizeAverages.csv", Sum_size)

# Ex Situ Experiment
ex_size = subset(size, in_ex == "ex")
ex_size$treatment = factor(ex_size$treatment, levels = c("DD", "SD","DS", "SS"))
ex_size$age_days = factor(ex_size$age_days)
str(ex_size)

SD.8 = subset(ex_size, treatment == "SD" & age_days == "8")
SS.8 = subset(ex_size, treatment == "SS" & age_days == "8")
DS.8 = subset(ex_size, treatment == "DS" & age_days == "8")
DD.8 = subset(ex_size, treatment == "DD" & age_days == "8")
day8_ex = subset(ex_size, age_days == "8")

SD.14 = subset(ex_size, treatment == "SD" & age_days == "14")
SS.14 = subset(ex_size, treatment == "SS" & age_days == "14")
DS.14 = subset(ex_size, treatment == "DS" & age_days == "14")
DD.14 = subset(ex_size, treatment == "DD" & age_days == "14")
day14_ex = subset(ex_size, age_days == "14")


shapiro.test(SS.8$diameter)
hist(SS.8$diameter)
shapiro.test(SD.8$diameter) # not normal
hist(SD.8$diameter)
shapiro.test(DS.8$diameter)
hist(DS.8$diameter)
shapiro.test(DD.8$diameter)  
hist(DD.8$diameter)

shapiro.test(SS.14$diameter)
hist(SS.14$diameter)
shapiro.test(SD.14$diameter) # not normal
hist(SD.14$diameter)
shapiro.test(DS.14$diameter) # sample size too small
hist(DS.14$diameter)
shapiro.test(DD.14$diameter) 
hist(DD.14$diameter)


exsize8_KS = kruskal.test(diameter~treatment, data=day8_ex)
exsize8_KS   # Kruskal-Wallis chi-squared = 10.804, df = 3, p-value = 0.01283
dunnTest(diameter~treatment, data=day8_ex, method = "bonferroni")  # only DD vs SD different

exsize14_KS = kruskal.test(diameter~treatment, data=day14_ex)
exsize14_KS   # Kruskal-Wallis chi-squared = 28.069, df = 3, p-value = 3.512e-06
dunnTest(diameter~treatment, data=day14_ex, method = "bonferroni")  # DD vs SD, DS vs SD, DD vs SS, 


ex.size.p = ggplot(ex_size, aes(y = diameter, x = age_days, fill = treatment)) +
  geom_boxplot(outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  facet_wrap(facets = "treatment", nrow = 1)+
  geom_point(position=position_jitterdodge(), size = 0.5)+
  labs(y= ~Spat ~Diameter ~(mm), x = ~Days ~Post ~Settlement)+
  scale_fill_viridis_d()+
  # scale_y_continuous(limits = c(0.5, 0.7), breaks = seq(0.5,0.7,0.1))+
  mytheme+
  guides(fill = FALSE)

ex.size.p

ggsave("SpatSizeByDay_ExSitu_StyloSpat2021.jpeg", plot = ex.size.p, width = 10, height = 6,dpi=300, 
       units = "cm")
ggsave("SpatSizeByDay_ExSitu_StyloSpat2021.pdf", plot = ex.size.p, width = 10, height = 6,dpi=300, 
       units = "cm")




### Respiration Planulae ######
r = read.csv("planulae_respiration_forR.csv")
str(r)
r$depth = as.factor(r$depth)
r.ind =subset(r, nmol.min.ind <0)
r.mm = subset(r, nmol.min.mm3 <0)


#by individual
stat.test.ind <- r.ind %>%
  t_test(nmol.min.ind ~ depth) %>%
  add_significance() %>%
  add_xy_position

ind.plot = ggplot(r.ind, aes(y = nmol.min.ind, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
 geom_jitter(position = position_jitter(width = .25), size = 0.2)+
labs(y= ~Respiration ~(nmol ~min^-1 ~larva^-1))+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#FDE725FF","#2C728EFF"))+
  scale_y_continuous(limits = c(-0.5, 0.1), breaks = seq(-0.5,0.1,0.2))+
  stat_pvalue_manual(stat.test.ind, label = "p.signif", tip.length = 0.01, y.position = 0.07)+
  mytheme+
  guides(fill = FALSE)

ind.plot


#by volume
stat.test.mm <- r.mm %>%
  t_test(nmol.min.mm3 ~ depth) %>%
  add_significance() %>%
  add_xy_position

mm.plot = ggplot(r.mm, aes(y = nmol.min.mm3, x = depth)) +
  geom_boxplot(aes(fill = depth), outlier.shape = "", fatten = 0.5, alpha = 0.7, lwd = 0.2)+
  geom_jitter(position = position_jitter(width = .25), size = 0.2)+
  labs(y= ~Larval ~Respiration ~(nmol ~min^1 ~mm^3))+
  scale_x_discrete(labels = mynames)+
  scale_fill_manual(values = c("#FDE725FF","#2C728EFF"))+
  #scale_y_continuous(limits = c(0, 0.66), breaks = seq(0,0.6,0.2))+
  stat_pvalue_manual(stat.test.mm, label = "p.signif", tip.length = 0.01, y.position = 0.2)+
  mytheme+
  guides(fill = FALSE)

mm.plot

r.plots <- plot_grid(ind.plot, mm.plot, 
                     #labels = c('A', 'B'),
                       label_y = 0.985, label_size = 11,ncol = 2, align = "h", byrow = F)

ggsave("Respiration_Planulae_2021.jpeg", plot = r.plots, width = 8, height = 6,dpi=300, 
       units = "cm")
ggsave("Respiration_Planulae_2021.pdf", plot = r.plots, width = 8, height = 6,dpi=300, 
       units = "cm")


### Combine plots for manuscript figure with images
planulae.jpg = ggdraw() + draw_image("Deep_Sha_Planulae.jpg")
nets.jpg = ggdraw() + draw_grob("MesophoticNets.jpg")

physio.plots <- plot_grid(vol.plot, pro.plot, cell.plot, ind.plot, 
                     labels = c('C', 'E', 'D', 'F'),
                     label_y = 0.985, label_size = 10, ncol = 2, align = "h", byrow = F)

ggsave("Figure2_PlanulaePhysiology.jpeg", plot = physio.plots, width = 11, height = 12,dpi=300, 
       units = "cm")
ggsave("Figure2_PlanulaePhysiology.pdf", plot = physio.plots, width = 11, height = 12,dpi=300, 
       units = "cm")

