
## Newl-identified Lentil (Lens culinaris Medik.) genotypes for Enhanced Adaptation to the Mediterranean Environment 
In this readme are reported all the analysis's steps present in the paper of Rocchetti et al.,2025

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
library(ggplot2)
library(viridis)
library(ggpubr)

print("Yield")
mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = YLD, random = ("gen"), verbose = TRUE)
  blup_yld<-data.frame(mixed_mod$YLD$BLUPgen)
yld_average<- mean(datalentil46$YLD, na.rm = TRUE)
above<- data.frame(gen = blup_yld$GEN[which(blup_yld$LL>yld_average)])
below<- data.frame(gen= blup_yld$GEN[which(blup_yld$UL<yld_average)])

blup_yld$significance <- "Average"
blup_yld$significance[blup_yld$GEN%in%above$gen] <- "Above"
blup_yld$significance[blup_yld$GEN%in%below$gen] <- "Below"

Yield<- ggplot(blup_yld, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 480)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("Yield (g/plot)") + ylab("genotypes")

Yield

print("Seed Weight")
mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = SW, random = ("gen"), verbose = TRUE)
blup_sw<-data.frame(mixed_mod$SW$BLUPgen)
sw_average<- mean(datalentil46$SW, na.rm = TRUE)
above<- data.frame(gen = blup_sw$GEN[which(blup_sw$LL>sw_average)])
below<- data.frame(gen= blup_sw$GEN[which(blup_sw$UL<sw_average)])


blup_sw$significance <- "Average"
blup_sw$significance[blup_sw$GEN%in%above$gen] <- "Above"
blup_sw$significance[blup_sw$GEN%in%below$gen] <- "Below"
SW<- ggplot(blup_sw, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 23)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("1000 Seed Weight(g) ") + ylab("genotypes")
SW

print("Canopy height")
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

print("First Flower")

mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = First_flower, random = ("gen"), verbose = TRUE)
blup_f<-data.frame(mixed_mod$First_flower$BLUPgen)
f_average<- mean(datalentil46$First_flower, na.rm = TRUE)
above<- data.frame(gen = blup_f$GEN[which(blup_f$LL>f_average)])
below<- data.frame(gen= blup_f$GEN[which(blup_f$UL<f_average)])


blup_f$significance <- "Average"
blup_f$significance[blup_f$GEN%in%above$gen] <- "Above"
blup_f$significance[blup_f$GEN%in%below$gen] <- "Below"
F1<- ggplot(blup_f, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 118)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("First flower(DAS)") + ylab("genotypes")
F1

print("First pod")

mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = First_pod, random = ("gen"), verbose = TRUE)
blup_fp<-data.frame(mixed_mod$First_pod$BLUPgen)
fp_average<- mean(datalentil46$First_pod, na.rm = TRUE)
above<- data.frame(gen = blup_f$GEN[which(blup_fp$LL>fp_average)])
below<- data.frame(gen= blup_f$GEN[which(blup_fp$UL<fp_average)])


blup_fp$significance <- "Average"
blup_fp$significance[blup_f$GEN%in%above$gen] <- "Above"
blup_fp$significance[blup_f$GEN%in%below$gen] <- "Below"
FP<- ggplot(blup_fp, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 130)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("First pod(DAS)") + ylab("genotypes")
FP


 print("Plant height")
 mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = PH, random = ("gen"), verbose = TRUE)
 blup_ph<-data.frame(mixed_mod$PH$BLUPgen)
 ph_average<- mean(datalentil46$PH, na.rm = TRUE)
 above<- data.frame(gen = blup_ph$GEN[which(blup_ph$LL>ph_average)])
 below<- data.frame(gen= blup_ph$GEN[which(blup_ph$UL<ph_average)])
 
 
 blup_ph$significance <- "Average"
 blup_ph$significance[blup_ph$GEN%in%above$gen] <- "Above"
 blup_ph$significance[blup_ph$GEN%in%below$gen] <- "Below"
 PH<- ggplot(blup_ph, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
     geom_point(aes(col = significance), size=4)+
     geom_vline(xintercept = 35.10)+
     scale_color_manual(values = c("Blue", "gray", "red"))+
     geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
     theme_classic2(base_size = 11)+
     xlab("Plant Height (cm)") + ylab("genotypes")
 PH


print("First pod height")
mixed_mod<-gamem_met(datalentil46, env = ENV, gen = GEN, rep = REP, resp = FPH, random = ("gen"), verbose = TRUE)
blup_fph<-data.frame(mixed_mod$FPH$BLUPgen)
fph_average<- mean(datalentil46$FPH, na.rm = TRUE)
above<- data.frame(gen = blup_fph$GEN[which(blup_fph$LL>fph_average)])
below<- data.frame(gen= blup_fph$GEN[which(blup_fph$UL<fph_average)])


blup_fph$significance <- "Average"
blup_fph$significance[blup_fph$GEN%in%above$gen] <- "Above"
blup_fph$significance[blup_fph$GEN%in%below$gen] <- "Below"
FPH<- ggplot(blup_fph, aes(x=Predicted, y=reorder(GEN, Predicted), group =significance)) + 
    geom_point(aes(col = significance), size=4)+
    geom_vline(xintercept = 8.10)+
    scale_color_manual(values = c("Blue", "gray", "red"))+
    geom_errorbar(aes(xmin=LL, xmax=UL), width=.1)+
    theme_classic2(base_size = 11)+
    xlab("First pod height (cm)") + ylab("genotypes")
FPH

together<- ggarrange(Yield, SW, F1, FP, CH, PH, FPH, nrow = 2, ncol = 4)

ggsave("blup46.jpeg", plot = together, device = "jpeg", width = 400, height = 300, units = "mm", dpi = 1000)

```

![blup46](https://github.com/user-attachments/assets/e6a4c269-46ae-4d5c-a850-6eb73f0226ac)


Phenotypic data collected from the 16 genotypes in common among all trials were used to estimate trait correlation and dissect Genotypes by Environment Interaction (GEI).
In the upcoming code we are going to estimate BLUPs for each combination of trait and environment following the model the mixed model  $Y = G + rep + e$ where _G_ is considered random and _rep_ fixed

```
#plot blup all together
datalentil_16<- read.csv("lentil_16.csv", header = TRUE)

############## Metaponto autumn 2019
mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = YLD)
blup_yld_MA2019<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_MA2019)[3]<- paste("yld_MA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_MA2019<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_MA2019)[3]<- paste("FirstF_MA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_MA2019<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_MA2019)[3]<- paste("FirstP_MA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = CH)
blup_CH_MA2019<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_MA2019)[3]<- paste("CH_MA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = PH)
blup_PH_MA2019<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_MA2019)[3]<- paste("PH_MA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2019"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_MA2019<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_MA2019)[3]<- paste("FPH_MA2019")



################# Metaponto autumn 2020
mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2020"), gen = GEN, rep= REP, resp = YLD)
blup_yld_MA2020<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_MA2020)[3]<- paste("yld_MA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2020"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_MA2020<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_MA2020)[3]<- paste("FirstF_MA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2020"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_MA2020<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_MA2020)[3]<- paste("FirstP_MA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2020"), gen = GEN, rep= REP, resp = CH)
blup_CH_MA2020<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_MA2020)[3]<- paste("CH_MA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2020"), gen = GEN, rep= REP, resp = PH)
blup_PH_MA2020<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_MA2020)[3]<- paste("PH_MA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Metaponto_autumn_2020"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_MA2020<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_MA2020)[3]<- paste("FPH_MA2020")

################# Osimo_autumn_2019
mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2019"), gen = GEN, rep= REP, resp = YLD)
blup_yld_OA2019<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_OA2019)[3]<- paste("yld_OA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2019"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_OA2019<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_OA2019)[3]<- paste("FirstF_OA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2019"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_OA2019<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_OA2019)[3]<- paste("FirstP_OA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2019"), gen = GEN, rep= REP, resp = CH)
blup_CH_OA2019<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_OA2019)[3]<- paste("CH_OA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2019"), gen = GEN, rep= REP, resp = PH)
blup_PH_OA2019<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_OA2019)[3]<- paste("PH_OA2019")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2019"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_OA2019<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_OA2019)[3]<- paste("FPH_OA2019")

################# Osimo_autumn_2020
mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2020"), gen = GEN, rep= REP, resp = YLD)
blup_yld_OA2020<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_OA2020)[3]<- paste("yld_OA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2020"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_OA2020<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_OA2020)[3]<- paste("FirstF_OA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2020"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_OA2020<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_OA2020)[3]<- paste("FirstP_OA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2020"), gen = GEN, rep= REP, resp = CH)
blup_CH_OA2020<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_OA2020)[3]<- paste("CH_OA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2020"), gen = GEN, rep= REP, resp = PH)
blup_PH_OA2020<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_OA2020)[3]<- paste("PH_OA2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2020"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_OA2020<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_OA2020)[3]<- paste("FPH_OA2020")


################# Osimo_autumn_2021
mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2021"), gen = GEN, rep= REP, resp = YLD)
blup_yld_OA2021<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_OA2021)[3]<- paste("yld_OA2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2021"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_OA2021<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_OA2021)[3]<- paste("FirstF_OA2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2021"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_OA2021<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_OA2021)[3]<- paste("FirstP_OA2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2021"), gen = GEN, rep= REP, resp = CH)
blup_CH_OA2021<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_OA2021)[3]<- paste("CH_OA2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2021"), gen = GEN, rep= REP, resp = PH)
blup_PH_OA2021<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_OA2021)[3]<- paste("PH_OA2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_autumn_2021"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_OA2021<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_OA2021)[3]<- paste("FPH_OA2021")


################# Osimo_spring_2020
mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = YLD)
blup_yld_OS2020<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_OS2020)[3]<- paste("yld_OS2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_OS2020<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_OS2020)[3]<- paste("FirstF_OS2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_OS2020<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_OS2020)[3]<- paste("FirstP_OS2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = CH)
blup_CH_OS2020<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_OS2020)[3]<- paste("CH_OS2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = PH)
blup_PH_OS2020<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_OS2020)[3]<- paste("PH_OS2020")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_OS2020<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_OS2020)[3]<- paste("FPH_OS2020")

################# Osimo_spring_2021
mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2020"), gen = GEN, rep= REP, resp = YLD)
blup_yld_OS2021<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld_OS2021)[3]<- paste("yld_OS2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2021"), gen = GEN, rep= REP, resp = FirstF)
blup_FirstF_OS2021<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_FirstF_OS2021)[3]<- paste("FirstF_OS2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2021"), gen = GEN, rep= REP, resp = FirstP)
blup_FirstP_OS2021<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_FirstP_OS2021)[3]<- paste("FirstP_OS2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2021"), gen = GEN, rep= REP, resp = CH)
blup_CH_OS2021<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_CH_OS2021)[3]<- paste("CH_OS2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2021"), gen = GEN, rep= REP, resp = PH)
blup_PH_OS2021<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_PH_OS2021)[3]<- paste("PH_OS2021")

mixed_mod<-gamem(subset(datalentil_16, ENV == "Osimo_spring_2021"), gen = GEN, rep= REP, resp = FPH)
blup_FPH_OS2021<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_FPH_OS2021)[3]<- paste("FPH_OS2021")


#put all data frames into list
df_list <- list(blup_yld_MA2019, blup_FirstF_MA2019, blup_FirstP_MA2019, blup_CH_MA2019, blup_PH_MA2019, blup_FPH_MA2019, blup_yld_MA2020, blup_FirstF_MA2020, blup_FirstP_MA2020, blup_CH_MA2020, blup_PH_MA2020, blup_FPH_MA2020, blup_yld_OA2019, blup_FirstF_OA2019, blup_FirstP_OA2019, blup_CH_OA2019, blup_PH_OA2019, blup_FPH_OA2019, blup_yld_OA2020, blup_FirstF_OA2020, blup_FirstP_OA2020, blup_CH_OA2020, blup_PH_OA2020, blup_FPH_OA2020, blup_yld_OA2021, blup_FirstF_OA2021, blup_FirstP_OA2021, blup_CH_OA2021, blup_PH_OA2021, blup_FPH_OA2021, blup_yld_OS2020, blup_FirstF_OS2020, blup_FirstP_OS2020, blup_CH_OS2020, blup_PH_OS2020, blup_FPH_OS2020, blup_yld_OS2021, blup_FirstF_OS2021, blup_FirstP_OS2021, blup_CH_OS2021, blup_PH_OS2021, blup_FPH_OS2021)      

#merge all data frames together

library(tidyverse)
library(dplyr)
df2 <- df_list %>% reduce(inner_join, by='GEN')

# R base - Select columns from list
df<-df2[,c("GEN","yld_MA2019", "FirstF_MA2019", "FirstP_MA2019", "CH_MA2019", "PH_MA2019", "FPH_MA2019", "yld_MA2020", "FirstF_MA2020", "FirstP_MA2020", "CH_MA2020", "PH_MA2020", "FPH_MA2020", "yld_OA2019", "FirstF_OA2019", "FirstP_OA2019", "CH_OA2019", "PH_OA2019", "FPH_OA2019", "yld_OA2020", "FirstF_OA2020", "FirstP_OA2020", "CH_OA2020", "PH_OA2020", "FPH_OA2020", "yld_OA2021", "FirstF_OA2021", "FirstP_OA2021", "CH_OA2021", "PH_OA2021", "FPH_OA2021", "yld_OS2020", "FirstF_OS2020", "FirstP_OS2020", "CH_OS2020", "PH_OS2020", "FPH_OS2020", "yld_OS2021", "FirstF_OS2021", "FirstP_OS2021", "CH_OS2021", "PH_OS2021", "FPH_OS2021")]


```
The derive dataset contains all the BLUP combination for the specific combination of trait and environment. 
The obtained dataset was used to conducte a PCA to visualize trait's genetic correlation and to have a general agronomic overview of the material

```
#add biological status
pivot_table <- datalentil_16 %>%
  select(GEN, Bio_stat) %>%
  group_by(GEN, Bio_stat) %>%
  summarise( .groups = 'drop') 

readyPC<-left_join(df, pivot_table, by='GEN')

library(FactoMineR)
library(factoextra)
rownames(readyPC) <- readyPC$GEN
res.pca<-PCA(readyPC[,2:43], scale.unit = TRUE, ncp = 5, graph = TRUE)
ind <- get_pca_ind(res.pca)
var <- get_pca_var(res.pca)
var
```
The results are printed out and plotted using ggplot2
```
library(ggplot2)
library(ggrepel)

# Create a data frame for PCA results
pca_data <- as.data.frame(ind$coord)
pca_data$Bio_stat <- readyPC$Bio_stat
pca_data$GEN <- readyPC$GEN  # Assuming you want to label individuals with GEN

# Plotting
#score plot
score<- ggplot(pca_data, aes(x = Dim.1, y = Dim.2)) +
 geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = pca_data, aes(x = Dim.1, y = Dim.2, color = Bio_stat), size = 3, shape = 16) +  # Individual points
geom_label_repel(data = pca_data, aes(x=Dim.1, y=Dim.2, label = GEN), size = 2.5, family = "Times",max.overlaps = Inf) +
  scale_color_manual(values = c("darkred", "darkorange", "darkgreen")) +  # Custom color palette
theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=15),axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))+
    xlab("PC1: 33.4%") + ylab("PC1: 23.6%")+
  labs(title = "PCA score plot", color = "Groups")

TAB_var <- as.data.frame(var$coord)
dd<-data.frame(trait = c("yield", "flowering", "flowering", "architecture", "architecture","architecture", "yield", "flowering", "flowering", "architecture", "architecture","architecture", "yield", "flowering", "flowering", "architecture", "architecture","architecture", "yield", "flowering", "flowering", "architecture", "architecture","architecture", "yield", "flowering", "flowering", "architecture", "architecture","architecture", "yield", "flowering", "flowering", "architecture", "architecture","architecture", "yield", "flowering", "flowering", "architecture", "architecture","architecture"))
TAB_var <- cbind(TAB_var, dd)
#loading
loading<- ggplot(TAB_var, aes(x = Dim.1, y = Dim.2)) +
 geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_segment(data = TAB_var, aes(xend=Dim.1, yend=Dim.2, x=0, y=0, color = trait), size=0.4, linetype=1, arrow=arrow(length = unit(0.02, "npc")))+ # Individual points
geom_label_repel(data = TAB_var, aes(x=Dim.1, y=Dim.2, label = rownames(TAB_var)), size = 2.5, family = "Times",max.overlaps = Inf) +
  scale_color_manual(values = c("blue3",  "darkorange", "green3")) +  # Custom color palette
theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=15),axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))+
    xlab("PC1: 33.4%") + ylab("PC1: 23.6%")+
  labs(title = "PCA loading plot")

PCAplot<-ggarrange(score, loading, nrow=1, ncol=2)

ggsave("PCAplot.jpeg", plot = PCAplot, device = "jpeg", width = 400, height = 180, units = "mm", dpi = 1000, bg = "white")
```

![PCAplot](https://github.com/user-attachments/assets/7d9c74af-c30e-411f-952a-8a7a98f8f21a)



From the same dataset of 16 genotypes in 7 environments we are going to estimate BLUPs for flowering architectural and production traits using the mixed model (1) $P = G + E + GEI$. 

Traits genetic correlation will be evaluated using pairwise-Pearson correlation.

```
#Yield
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = YLD, random = ("gen"), verbose = TRUE)
  blup_yld<-data.frame(mixed_mod$YLD$BLUPgen)
names(blup_yld)[5]<- paste("blup_YLD")

#FirstF
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = FirstF, random = ("gen"), verbose = TRUE)
  blup_f1<-data.frame(mixed_mod$FirstF$BLUPgen)
names(blup_f1)[5]<- paste("blup_FirstF")

#FirstP
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = FirstP, random = ("gen"), verbose = TRUE)
  blup_fp<-data.frame(mixed_mod$FirstP$BLUPgen)
names(blup_fp)[5]<- paste("blup_FirstP")

#CH
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = CH, random = ("gen"), verbose = TRUE)
  blup_ch<-data.frame(mixed_mod$CH$BLUPgen)
names(blup_ch)[5]<- paste("blup_CH")

#PH
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = PH, random = ("gen"), verbose = TRUE)
  blup_ph<-data.frame(mixed_mod$PH$BLUPgen)
names(blup_ph)[5]<- paste("blup_PH")

#FPH
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = FPH, random = ("gen"), verbose = TRUE)
  blup_fph<-data.frame(mixed_mod$FPH$BLUPgen)
names(blup_fph)[5]<- paste("blup_FPH")

#SW
mixed_mod<-gamem_met(datalentil_16, env = ENV, gen = GEN, rep = REP, resp = SW, random = ("gen"), verbose = TRUE)
  blup_sw<-data.frame(mixed_mod$SW$BLUPgen)
names( blup_sw)[5]<- paste("blup_SW")

#put all data frames into list
df_list <- list(blup_yld,blup_f1, blup_fp, blup_ch, blup_ph, blup_fph, blup_sw)  
#merge all data frames together

library(tidyverse)
library(dplyr)
df <- df_list %>% reduce(inner_join, by='GEN')

# R base - Select columns from list
df<-df[,c("GEN","blup_YLD", "blup_FirstF", "blup_FirstP","blup_CH", "blup_PH", "blup_FPH", "blup_SW")]
cc<-corr_plot(df)
ggsave("corrplot.jpeg", plot = cc, device = "jpeg", width = 250, height = 180, units = "mm", dpi = 1000, bg = "white")
```

![corrplot](https://github.com/user-attachments/assets/d71abbe2-8bb7-4421-a753-099a416765bf)


The Additive Main Effect and Multiplicative Interaction (AMMI) model (Gauch, 1988; van Eeuwijk, 1995) was used as a fixed model framework to dissect GEI variance in different Interaction Principal Components (IPCAs), using singular value decomposition (SVD) procedure.
In the next chuck of code we are going to apply AMMI using the package Metan;

```
#perform AMMI1 biplot
model <- performs_ammi(datalentil_16, ENV, GEN, rep, resp = YLD, Verbose =FALSE)
a<- plot_scores(model, type = 1,
                col.env = "blue",
                col.segm.env = "blue",
                col.gen = "black", plot_theme = theme_metan_minimal(), title = TRUE)
a

data_ammi<- as.data.frame(a[["data"]])
pivot_table <- datalentil_16 %>%
  select(GEN, Bio_stat) %>%
  group_by(GEN, Bio_stat) %>%
  summarise( .groups = 'drop')
names(pivot_table)[1] <-paste("Code")

data_ammi<-left_join(data_ammi, pivot_table, by='Code')

ammi1<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=528, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = subset(data_ammi,type == "GEN"), aes(x=Y, y=PC1, color=Bio_stat), size = 4.5) +
  scale_color_manual(values = c("darkred", "darkorange", "chartreuse3")) + 
  geom_segment(data = subset(data_ammi,type == "ENV"), aes(xend=Y, yend=PC1, x=528, y=0), size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc")))+
 geom_label_repel(data = subset(data_ammi,type == "GEN"), aes(x=Y, y=PC1, label = Code), size = 2.5, family = "Times", max.overlaps = Inf) +
  geom_label_repel(data = subset(data_ammi,type == "ENV"), aes(x=Y, y=PC1, label = Code), size = 2.5, family = "Times, max.overlaps = Inf") +
  xlab("Yield (g/plot)") + ylab("PC1: 45%") +
  guides(color=guide_legend(title="Biological status")) +
  theme_bw(base_size = 10, base_family = "Times") +
   theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=15),axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))+
 labs(title = "AMMI 1")
ammi1


#perform AMMI2 biplot
model <- performs_ammi(datalentil_16, ENV, GEN, rep, resp = YLD, Verbose =FALSE)
b<- plot_scores(model, type = 2,
                col.env = "blue",
                col.segm.env = "blue",
                col.gen = "black", plot_theme = theme_metan_minimal(), title = TRUE)
b

data_ammi2<- as.data.frame(b[["data"]])
pivot_table <- datalentil_16 %>%
  select(GEN, Bio_stat) %>%
  group_by(GEN, Bio_stat) %>%
  summarise( .groups = 'drop')
names(pivot_table)[1] <-paste("Code")

data_ammi2<-left_join(data_ammi2, pivot_table, by='Code')

ammi2<-ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = subset(data_ammi,type == "GEN"), aes(x=PC1, y=PC2, color=Bio_stat), size = 4.5) +
  scale_color_manual(values = c("darkred", "darkorange", "chartreuse3")) + 
  geom_segment(data = subset(data_ammi,type == "ENV"), aes(xend=PC1, yend=PC2, x=0, y=0), size=0.2, linetype=1, arrow=arrow(length = unit(0.02, "npc")))+
geom_label_repel(data = subset(data_ammi,type == "GEN"), aes(x=PC1, y=PC2, label = Code), size = 2.5, family = "Times" ,max.overlaps = Inf) +
  geom_label_repel(data = subset(data_ammi,type == "ENV"), aes(x=PC1, y=PC2, label = Code), size = 2.5, family = "Times", max.overlaps = Inf) +
  xlab("PC1: 45%") + ylab("PC2: 34.7%") +
  guides(color=guide_legend(title="Biological status")) +
  theme_bw(base_size = 10, base_family = "Times") +
   theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(1)), strip.text = element_text(size=15),axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))+
 labs(title = "AMMI 2")
ammi2

ammi<-ggarrange(ammi1, ammi2, nrow=1, ncol=2)
ggsave("ammi.jpeg", plot = ammi, device = "jpeg", width = 400, height = 180, units = "mm", dpi = 1000, bg = "white")
```

![ammi](https://github.com/user-attachments/assets/ed5609e5-2c97-4497-9fbc-6bef31171c9a)
