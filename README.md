
## Newly identified Lentil (Lens culinaris Medik.) genotypes for Enhanced Adaptation to the Mediterranean Environment 
In this readme are reported all the analysis steps present in the paper of Rocchetti et al.,2025

> open libraries and set directory
```
setwd("/lustre/rocchettil/Paper_lentil/Manoscript_Lentil_editing")
datalentil46<- read.csv("lentil46whole.csv", header = TRUE)
install.packages("colorspace")
install.packages("metan")
library(metan)
```

> BLUP estimation for the whole set of 46 genotypes in è different environment

Considering the linear model:

 (1) $P = G + E + GEI$
 
the phenotypic value of a trait of any individual in a given environment can be written as the sum of its genetic effect G, the environmental effect E, the genotype by environment interaction (GEI) and e as the random residual effect within each environment following a normal distribution N (0, σ2).

Considering G, E and GEI as random factors we estimated the respective variance component using residual maximum likelihood (RELM).
From the equation (1) by modelling the genotypic effect (G) and GEI as random and (E) as fixed, the Best Linear Unbiased Prediction for genotypes (BLUPg) (Piepho et al., 1998; Piepho et al., 1994) or the genotypic value was estimated for the entire set of 46 genotypes.  
```
#Yield
mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = YLD, random = ("gen"), verbose = TRUE)#random = "all" means G;E;GEI are random if random = "gen" means G and GEI random, E and rep[E]fixed 
get_model_data(mixed_mod, "genpar")
plot(mixed_mod, type = "re")

get.variance <- get_model_data(mixed_mod,what = "vcomp")
print.table(get.variance)
print.table(mixed_mod)
library(ggplot2)
b <- plot_blup(mixed_mod, which = "gen", plot_theme = theme_metan_minimal())
b
```
For a more personalzed output using a ggplot
```

  blup_yld<-data.frame(mixed_mod$YLD$BLUPgen)
yld_average<- mean(datalentil46$YLD, na.rm = TRUE)
above<- data.frame(gen = blup_yld$GEN[which(blup_yld$LL>yld_average)])
below<- data.frame(gen= blup_yld$GEN[which(blup_yld$UL<yld_average)])


blup_yld$significance <- "Average"
blup_yld$significance[blup_yld$GEN%in%above$gen] <- "Above"
blup_yld$significance[blup_yld$GEN%in%below$gen] <- "Below"

library(ggplot2)
library(viridis)
library(ggpubr)
Yield<- ggplot(blup_yld, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 480)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("Yield (g/plot)") + ylab("genotypes")

Yield

#Canopy height
mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = CH, random = ("gen"), verbose = TRUE)
blup_ch<-data.frame(mixed_mod$CH$BLUPgen)
ch_average<- mean(datalentil46$CH, na.rm = TRUE)
above<- data.frame(gen = blup_ch$GEN[which(blup_ch$LL>ch_average)])
below<- data.frame(gen= blup_ch$GEN[which(blup_ch$UL<ch_average)])


blup_ch$significance <- "Average"
blup_ch$significance[blup_ch$GEN%in%above$gen] <- "Above"
blup_ch$significance[blup_ch$GEN%in%below$gen] <- "Below"
CH<- ggplot(blup_ch, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 23.8)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("Canopy Height (cm)") + ylab("genotypes")
CH

#First Flower

mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = First_flower, random = ("gen"), verbose = TRUE)
blup_f<-data.frame(mixed_mod$First_flower$BLUPgen)
f_average<- mean(datalentil46$First_flower, na.rm = TRUE)
above<- data.frame(gen = blup_f$GEN[which(blup_f$LL>f_average)])
below<- data.frame(gen= blup_f$GEN[which(blup_f$UL<f_average)])


blup_f$significance <- "Average"
blup_f$significance[blup_f$GEN%in%above$gen] <- "Above"
blup_f$significance[blup_f$GEN%in%below$gen] <- "Below"
F<- ggplot(blup_f, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 118)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("Fisrt flower(DAS)") + ylab("genotypes")
F

ggsave("blup46.jpeg", plot = together, device = "jpeg", width = 400, height = 200, units = "mm", dpi = 1000)

```

![blup46](https://github.com/user-attachments/assets/e391481e-6473-4d72-acfc-51f29bf4c019)

Phenotypic data collected from the 16 genotypes in common among all trials were used to estimate trait correlation and dissect Genotypes by Environment Interaction (GEI).
In the upcoming code we are going to estimate BLUPs for each environment following the model  $Y = G + rep + e$

```
#plot blup all together
datalentil_16<- read.csv("lentil_16.csv", header = TRUE)
#Metaponto autumn 2019
mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = YLD)
blup_yld_met_autu_2019<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_met_autu_2019)[3]<- paste("yld_met_autu_2019")
mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_met_autu_2019<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_met_autu_2019)[3]<- paste("FirstF_met_autu_2019")
data_c<- data.frame(merge(blup_yld_met_autu_2019, blup_FirstF_met_autu_2019, by = "GEN"))
```

