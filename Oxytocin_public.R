## Complementary Code for the Paper "Salivary Oxytocin and Amygdalar Alterations in Functional Neurological Disorders
# Author: Samantha Weber (samantha.weber@bli.uzh.ch)
# Current R version: 4.2.0

# This code works with example data but won't produce the same results as reported in the paper, as the data
# can only be shared on request. This code should give a braod overview on how all the analyses in the manuscript
# were conducted. 

# v1.0

# Load Libraries

library(ISwR)
library(MASS)
library(ggplot2)
library(Hmisc)
library(plotly)
library(readxl)
library(car)
library(lawstat)
library(ggpubr)
library(rstatix)
library(lme4)
library(plyr)
library(multcomp)
library(dplyr)
library(tidyverse)
library(dplyr)
library(SNPassoc)
library(ggpubr)
library(corrplot)
library(Hmisc)
library(mediation)


# Define the folder, in which you save your plots
setwd(getwd()) #set path for outputfiles 
#####---------------------------------------------------------------- ##########

# 1. Load Data -----------------------------------------------------------------

# Read in the data
df <- read_excel("ExampleData.xlsx")
df$ctq<-df$ctq_emoab+df$ctq_physab+df$ctq_sexab+df$ctq_emoneg+df$ctq_physneg


# Save covariates as factors
df$p_code<-as.character(df$p_code)
df$group<-as.factor(df$group)
df$gender<-as.factor(df$gender)
df$psychMed<-as.factor(df$psychMed)
df$MensCycle<-as.factor(df$MensCycle)
df$menopause<-as.factor(df$menopause)
df$contraception<-as.factor(df$contraception)
df$smoke<-as.factor(df$smoke)


# 2. Statistics OXYTOCIN -------------------------------------------------------

## Plot Data
OTplot<-ggplot(df, aes(x=group, y=Oxytocin, fill=group)) +
  geom_boxplot()+
  theme_classic()+  
  scale_fill_manual(values = c("#868686FF" ,"#2BB08E"))+ 
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.title.x = element_blank()) +
  ggtitle("Oxytocin pg/ml") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+  
  theme(axis.title.x=element_blank(),axis.text.y = element_text(size = 14, face = "bold"), 
        axis.text.x = element_text(size = 18, face = "bold"), plot.title = element_text(hjust = 0.5))+
  stat_compare_means(label = "p.signif", comparisons = list(c("HC", "FND")), label.x=1)+ 
  geom_jitter(shape=16, position=position_jitter(0.2))

OTplot

table(df$group)
table(df$group,df$gender)


# Without covariates
OTlm <- lm(formula = Oxytocin ~ group,
           data = df) 

summary(OTlm)
OTAov <- aov(OTlm)
summary(OTAov) # 

TukeyHSD(OTAov, which ='group') 


## ANOVA on the fitted model/data
# We run a repeated measures Anova on the fitted data. So first a glm, then an Anova on the glm data in order to have covariates included. 
# Using "flagHigh can reduce the number of factors, here we use now the date of menstrual cycle, contraception and menopause
OTlm <- lm(formula = Oxytocin ~ group  + gender + psychMed +bdi +stai1 + date_diff + contraception +menopause + age,
           data = df) 

summary(OTlm) 
OTAov <- aov(OTlm)
summary(OTAov) # 
TukeyHSD(OTAov, which ='group') 


# 3. Genetic Association study -------------------------------------------------

#3.1  Quality control of data --------------------------------------------------

#we normally check per row and per column for missing data. If > 95% is missing in a row/column, then the row/column is excluded. 
# as we have only 10 SNPs and very few subjects, we have 0 % data.

names(df) #df_all = selected data
df$cmiss<-apply(df[,3:12],1,function(x) length(which(is.na(x))))
df$pmiss<-df$cmiss*100/10
head(df) # check column pmiss and cmiss
# remove colums/rows with more than 90% missing.


## Quality control using Hardy-Weinberg Equilibrium (HWE) 
# FUNDAMENTAL: CREATE setupSNP:
# we first specify which columns contain the genotyping data, and then what alleles they are without separation

myData_S<-setupSNP(data=df,colSNPs=2:11,sep="")

plotMissing(myData_S) #visualizes combination of SNPs x Subjects who are missing --> P60 in our case. 
class(myData_S)

#Make individual dataframes
df_FND<-filter(df, group =="FND")
df_HC<-filter(df, group =="HC")


myData_FND<-setupSNP(data=df_FND, colSNPs=2:11,sep="")
myData_HC<-setupSNP(data=df_HC,colSNPs=2:11,sep="")


## 3.2. Hardy-Weinberg Equilibrium Test ----------------------------------------
# CALCULATE frequency of genotypes, alleles y do HWE test. 

# summary view of various aspects of SNP quality control

summary_S<-summary(myData_S)

# Separate for HC and FND
summary_FND<-summary(myData_FND)
summary_HC<-summary(myData_HC)


##3.2 Association Analysis for one SNP -----------------------------------------
#TPH2 rs4570625
#TPH1 rs1800532
#BDNF rs6265
#BDNF rs1491850
#DRD4 rs3758653
#OXTR rs2254298
#OXTR rs53576
#DRD2 rs1799732
#FKBP51 rs1360780
#FKBP52 rs3800373
#TPH1 rs1800532

myData_S$group2<-myData_S$group
myData_S$group<-revalue(myData_S$group, c("FND"="1","HC"="0"))
myData_S$group<-as.factor(myData_S$group)
#re-order factor levels for region
myData_S$group <- factor(myData_S$group, levels=c('0', '1'))
levels(myData_S$group)

# Check all genes at once
asoc_S1<-WGassociation(group~1+ age + gender + bdi + stai2 + ctq, data=myData_S, model=c("codominant","log-additive","dominant","recessive"),genotypingRate=80)
asoc_S1 #Interpretation: we see that rs53576 has a significant association whith group, when using the recessive model (homozygote for "bad" allele)
WGstats(asoc_S1)# to see the details of the regression
# only rs53576 looks interesting
WGstats(asoc_S1)
asoc_S1

# Logistic regression with functions from SNPassoc:
# Can group be explained by SNP?

# No covariates
logistic_S<-association(group~rs53576, data=myData_S)
logistic_S

#covariates
logistic_SC<-association(group~rs53576 + age + gender + ctq + bdi + stai2, data=myData_S) 
logistic_SC #significant association between group and genotype --> this SNP is a risk factor for FND

# Linear Regression with oxytocin

#Oxytocin
asoc_OT<-WGassociation(Oxytocin~1+ gender + bdi +psychMed + stai2 + date_diff + menopause + contraception + age + ctq,data=myData_S, 
                       model=c("codominant","dominant","recessive"),genotypingRate=80)

asoc_OT # NO associations between OT and genotype

#In FND only, with and without covariates
asoc_OT_FND<-WGassociation(Oxytocin~1,data=myData_FND, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_FND
asoc_OT_FNDC<-WGassociation(Oxytocin~1 + gender + bdi +psychMed + stai2 + date_diff + menopause + contraception + age + ctq, data=myData_FND, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_FNDC # Nada

# In HC only, with and without covariates
asoc_OT_HC<-WGassociation(Oxytocin~1,data=myData_HC, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_HC
asoc_OT_HCC<-WGassociation(Oxytocin~1 + gender + bdi + stai2 + date_diff + menopause + contraception + age + ctq, data=myData_HC, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_HCC #even more nada


# 4. Methylation data -----------------------------------------------------------
# Let's check methylation data

df$OXTR_CpG_sum<-(df$OXTR_CpG_2 + df$OXTR_CpG_5 + df$OXTR_CpG_6.7 + df$OXTR_CpG_8)/4

stat.test <- df %>%
  wilcox_test(OXTR_CpG_sum ~ group) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()
stat.test

## 4.1. Analyses regression with methylation -----------------------------------

# We run a repeated measures Anova on the fitted data. So first a glm, then an Anova on the glm data in order to have covariates included. 
# Plot data
dfOT<-df

OXTRplot<-ggplot(dfOT, aes(x=group, y=OXTR_CpG_sum, fill=group)) +
  geom_boxplot()+
  theme_classic()+  
  scale_fill_manual(values = c("#868686FF" ,"#2BB08E"))+ 
  theme(legend.position = "none", axis.title.y = element_blank(), 
        axis.title.x = element_blank()) +
  ggtitle("OXTR Methylation [%]") + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+  
  theme(axis.title.x=element_blank(),axis.text.y = element_text(size = 14, face = "bold"), 
        axis.text.x = element_text(size = 18, face = "bold"), plot.title = element_text(hjust = 0.5))+ 
  stat_compare_means(label = "p.signif", comparisons = list(c("HC", "FND")), label.x=1)+ 
  geom_jitter(shape=16, position=position_jitter(0.2))

OXTRplot

## Figure 2
library(ggtext) # Load the ggtext package
dfOT_FND <- filter(dfOT, group == "FND")
dfOT_HC <- filter(dfOT, group == "HC")

#Change to violin plot
plota <- ggplot(dfOT_FND, aes_string(x = "rs53576", y = "Oxytocin", fill = "rs53576")) +
  geom_violin(show.legend = FALSE, alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", show.legend = FALSE, outlier.shape = NA) +
  ylim(0, 25) +
  scale_fill_manual(values = c("#440154FF", "#440154FF", "#440154FF")) +
  labs(y = "Oxytocin [pg/ml]", x = "rs53576 genotype", title = "<b>A</b> Oxytocin FND") +
  theme_classic() +
  theme(plot.title = element_markdown())

plotb <- ggplot(dfOT_HC, aes_string(x = "rs53576", y = "Oxytocin", fill = "rs53576")) +
  geom_violin(show.legend = FALSE, alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", show.legend = FALSE, outlier.shape = NA) +
  ylim(0, 25) +
  scale_fill_manual(values = c("#21908CFF", "#21908CFF", "#21908CFF")) +
  labs(y = "Oxytocin [pg/ml]", x = "rs53576 genotype", title = "<b>B</b> Oxytocin Controls") +
  theme_classic() +
  theme(plot.title = element_markdown())

plotc <- ggplot(dfOT_FND, aes_string(x = "rs53576", y = "OXTR_CpG_sum", fill = "rs53576")) +
  geom_violin(show.legend = FALSE, alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", show.legend = FALSE, outlier.shape = NA) +
  ylim(0.2, 0.8) +
  scale_fill_manual(values = c("#440154FF", "#440154FF", "#440154FF")) +
  labs(y = "OXTR Methylation [%]", x = "rs53576 genotype", title = "<b>C</b> OXTR Methylation FND") +
  theme_classic() +
  theme(plot.title = element_markdown()) +
  stat_compare_means(label = "p.signif", comparisons = list(c("GG", "AA"), c("GA", "AA")), label.x = 1, label.y = c(0.6, 0.7))

plotd <- ggplot(dfOT_HC, aes_string(x = "rs53576", y = "OXTR_CpG_sum", fill = "rs53576")) +
  geom_violin(show.legend = FALSE, alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", show.legend = FALSE, outlier.shape = NA) +
  ylim(0.2, 0.8) +
  scale_fill_manual(values = c("#21908CFF", "#21908CFF", "#21908CFF")) +
  labs(y = "OXTR Methylation [%]", x = "rs53576 genotype", title = "<b>D</b> OXTR Methylation Controls") +
  theme_classic() +
  theme(plot.title = element_markdown())


Fig2<-ggarrange(plota, plotb, plotc, plotd)
Fig2

## Replication Apazoglou paper
#check mean 
mean_meth<-dfOT%>%
  group_by(group) %>%
  get_summary_stats(OXTR_CpG_sum, type = "mean_sd")
mean_meth
mean_meth<-dfOT%>%
  group_by(group) %>%
  get_summary_stats(OXTR_CpG_2:OXTR_CpG_sum, type = "mean_sd")
mean_meth
mean_meth<-test%>%
  group_by(group) %>%
  get_summary_stats(OXTR_CpG_2:OXTR_CpG_sum, type = "mean_sd")
mean_meth

# Method: Group differences were assessed with two-sided Student’s t-test

res_meth<-t.test(dfOT$OXTR_CpG_sum[dfOT$group == "FND"], dfOT$OXTR_CpG_sum[dfOT$group == "HC"])
res_meth # Not significant

# Of note: Apazoglou reports higher values in FND patients (68%) compared ot HC (62%). This is a difference of 6%. Everything below 10% anyway has no
# biological significance... 


## Let's make all this a bit more elegant... 
# First without covariates
OXTRlm <- lm(formula = Oxytocin ~ group*OXTR_CpG_sum,
             data = dfOT) 

OTAov <- aov(OXTRlm)
summary(OTAov) # significant
TukeyHSD(OTAov, which ='group') 

# with covariates
OXTRlm <- lm(formula = OXTR_CpG_sum ~ group + psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
             data = dfOT)

#Interaction with ctq
OXTRlm <- lm(formula = OXTR_CpG_sum ~ group*ctq + psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
             data = dfOT) 
 
summary(aov(OXTRlm))

#Interaction with sex
OXTRlm <- lm(formula = OXTR_CpG_sum ~ group*sex + psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
             data = dfOT) 

summary(aov(OXTRlm))


# Interaction with genotype
OXTRlm <- aov(lm(formula = OXTR_CpG_sum ~ group*rs53576+ Oxytocin +  psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
                 data = dfOT))
summary(OXTRlm)

OXTRlm <- aov(lm(formula = OXTR_CpG_sum ~ rs53576 + Oxytocin +  psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
                 data = dfOT_FND))
summary(OXTRlm)

OXTRlm <- aov(lm(formula = OXTR_CpG_sum ~ rs53576+ Oxytocin +bdi +stai1 + date_diff + contraception + menopause + age,
                 data = dfOT_HC))
summary(OXTRlm) 


OXTRlm <- aov(lm(formula = OXTR_CpG_sum ~ group*rs53576*ctq + Oxytocin+  psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
                 data = dfOT))
summary(OXTRlm) 


### 4.1.1. Explore Interaction -------------------------------------------------
# Set a color-blind friendly palette from viridis
library(viridis)
my_palette <- viridis_pal(option = "D")(3)
ggscatter(dfOT, 
          x = "Oxytocin", 
          y = "OXTR_CpG_sum",
          fill = "group", 
          color = "group",
          shape = "group",  # Use the 'group' variable for shapes
          add = "reg.line", 
          conf.int = TRUE, 
          cor.coef = FALSE, 
          cor.method = "spearman",
          palette = my_palette,  # Set the custom color palette
          xlab = "Oxytocin [pg/ml]", 
          ylab = "OXTR Methylation [%]") +
  # stat_cor(aes(color = group), method = "spearman")
  theme_minimal() +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Group"))


## IN FND ONLY
# With covariates
OT_FND<-filter(dfOT, group =="FND")
OXTRlm <- lm(formula = Oxytocin ~ OXTR_CpG_sum*rs53576 + ctq_emoneg + bdi + psychMed +stai1 + date_diff + contraception + menopause + age,
             data = OT_FND) 

summary(OXTRlm)  
OTextAov <- aov(OXTRlm)
summary(OTextAov) # significant

TukeyHSD(OTextAov, which ='group')


CTQ_high<-filter(OT_FND, ctq > 35)
dim(CTQ_high)
CTQ_low<-filter(OT_FND, ctq < 36)
dim(CTQ_low)

CTQ_all<-CTQ_high
CTQ_all$group_ctq<-1
CTQ_low$group_ctq<-0
CTQ_all<-rbind(CTQ_all,CTQ_low)
dim(CTQ_all)
CTQ_all$group<-as.factor(CTQ_all$group)

OXTRlm <- lm(formula = group_ctq ~ Oxytocin * OXTR_CpG_sum + psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
             data = CTQ_all) 
summary(OXTRlm)  
OTextAov <- aov(OXTRlm)
OXTRlm <- lm(formula = Oxytocin ~ group_ctq* OXTR_CpG_sum + psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
             data = CTQ_all)
summary(OXTRlm)  
OTextAov <- aov(OXTRlm)
OXTRlm <- lm(formula = Oxytocin ~ rs53576*ctq + psychMed +bdi +stai1 + date_diff + contraception + menopause + age,
             data = OT_FND) 
summary(OXTRlm)  
OTextAov <- aov(OXTRlm)
OXTRlm <- lm(formula = Oxytocin ~ rs53576*ctq +bdi +stai1 + date_diff + contraception + menopause + age,
             data = OT_HC)
summary(OXTRlm)  
OTextAov <- aov(OXTRlm)


## 4.2 OT x OXTR Meth x Genotype Assoc -----------------------------------------

# Linear Regression with oxytocin
#Oxytocin
asoc_OT<-WGassociation(Oxytocin~1+ psychMed +bdi +stai1 + date_diff + contraception + menopause + age + ctq + OXTR_CpG_sum,data=myData_S, 
                       model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT

#In FND only
myData_FND <- filter(myData_S, group == "FND")
asoc_OT_FND<-WGassociation(Oxytocin~1+ OXTR_CpG_sum, data=myData_FND, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_FND

asoc_OT_FNDC<-WGassociation(Oxytocin~1 +psychMed +bdi +stai1 + date_diff + contraception + menopause + age + ctq + OXTR_CpG_sum,data=myData_FND, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_FNDC

# In HC only
myData_HC <- filter(myData_S, group == "HC")
asoc_OT_HC<-WGassociation(Oxytocin~1+ OXTR_CpG_sum,data=myData_HC, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_HC

asoc_OT_HCC<-WGassociation(Oxytocin~1 +bdi +stai1 + date_diff + contraception + menopause + age + ctq+ OXTR_CpG_sum,data=myData_HC, model=c("codominant","dominant","recessive"),genotypingRate=80)
asoc_OT_HCC





## Supplementary Material Figures

a <- ggscatter(dfOT, x = "ctq_emoab", y = "OXTR_CpG_sum", fill = "group", color = "group",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = FALSE, cor.method = "spearman",
               xlab = "CTQ - Emotional Abuse", ylab = "OXTR Methylation [%]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
a
  
b<-ggscatter(dfOT, x = "ctq_emoab", y = "Oxytocin",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Abuse", ylab = "Oxytocin [ng/ml]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
b
c<-ggscatter(dfOT, x = "ctq_emoneg", y = "OXTR_CpG_sum",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Neglect", ylab = "OXTR Methylation [%]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
c

d<-ggscatter(dfOT, x = "ctq_emoneg", y = "Oxytocin",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Neglect", ylab = "Oxytocin [ng/ml]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
d

e<-ggscatter(dfOT, x = "ctq_physab", y = "OXTR_CpG_sum",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Abuse", ylab = "OXTR Methylation [%]")+ 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
e
f<-ggscatter(dfOT, x = "ctq_physab", y = "Oxytocin",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Abuse", ylab = "Oxytocin [ng/ml]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
f
g<-ggscatter(dfOT, x = "ctq_physneg", y = "OXTR_CpG_sum",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Neglect", ylab = "OXTR Methylation [%]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
g

h<-ggscatter(dfOT, x = "ctq_physneg", y = "Oxytocin",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Neglect", ylab = "Oxytocin [ng/ml]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
h
i<-ggscatter(dfOT, x = "ctq_sexab", y = "OXTR_CpG_sum",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Sexual Abuse", ylab = "OXTR Methylation [%]")+ 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
i
j<-ggscatter(dfOT, x = "ctq_sexab", y = "Oxytocin",fill = "group", color = "group",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Sexual Abuse", ylab = "Oxytocin [ng/ml]") + 
  stat_cor(aes(color = group), method = "spearman", 
           label.x = Inf, label.y = c(Inf, Inf), 
           hjust = 1.1, vjust = c(2, 4)) + 
  scale_fill_manual(values = c("#440154FF", "#21908CFF")) +
  scale_color_manual(values = c("#440154FF", "#21908CFF"))
j
ggarrange(a,b)
ggarrange(c,d)
ggarrange(e,f)
ggarrange(g,h)
ggarrange(i,j)

#Emotional Neglect
a <- ggscatter(dfOT_FND, x = "ctq_emoneg", y = "OXTR_CpG_sum", fill = "rs53576", color = "rs53576",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = FALSE, cor.method = "spearman",
               xlab = "CTQ - Emotional Neglect", ylab = "OXTR Methylation [%]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.59, 0.56, 0.53), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.6)  # Extend y-axis to 0.6

a


b<-ggscatter(dfOT_HC, x = "ctq_emoneg", y = "OXTR_CpG_sum",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Neglect", ylab = "OXTR Methylation [%]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.66, 0.63, 0.60), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.7)  # Extend y-axis to 0.6
b

c<-ggscatter(dfOT_FND, x = "ctq_emoneg", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Neglect", ylab = "Oxytocin [pg/ml]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6
c

d<-ggscatter(dfOT_HC, x = "ctq_emoneg", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Neglect", ylab = "Oxytocin [pg/ml]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6

d
ggarrange(a,b,c,d)



# Emotional Abuse

#Emotional Neglect
a <- ggscatter(dfOT_FND, x = "ctq_emoab", y = "OXTR_CpG_sum", fill = "rs53576", color = "rs53576",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = FALSE, cor.method = "spearman",
               xlab = "CTQ - Emotional Abuse", ylab = "OXTR Methylation [%]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.59, 0.56, 0.53), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.6)  # Extend y-axis to 0.6
a

b<-ggscatter(dfOT_HC, x = "ctq_emoab", y = "OXTR_CpG_sum",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Abuse", ylab = "OXTR Methylation [%]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.66, 0.63, 0.60), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.7)  # Extend y-axis to 0.6
b

c<-ggscatter(dfOT_FND, x = "ctq_emoab", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Abuse", ylab = "Oxytocin [pg/ml]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6
c

d<-ggscatter(dfOT_HC, x = "ctq_emoab", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Emotional Abuse", ylab = "Oxytocin [pg/ml]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6

d
ggarrange(a,b,c,d)


#Physical Neglect

#Emotional Neglect
a <- ggscatter(dfOT_FND, x = "ctq_physneg", y = "OXTR_CpG_sum", fill = "rs53576", color = "rs53576",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = FALSE, cor.method = "spearman",
               xlab = "CTQ - Physical Neglect", ylab = "OXTR Methylation [%]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.59, 0.56, 0.53), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.6)  # Extend y-axis to 0.6
a

b<-ggscatter(dfOT_HC, x = "ctq_physneg", y = "OXTR_CpG_sum",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Neglect", ylab = "OXTR Methylation [%]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.66, 0.63, 0.60), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.7)  # Extend y-axis to 0.6
b

c<-ggscatter(dfOT_FND, x = "ctq_physneg", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Neglect", ylab = "Oxytocin [pg/ml]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6
c

d<-ggscatter(dfOT_HC, x = "ctq_physneg", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Neglect", ylab = "Oxytocin [pg/ml]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6

d
ggarrange(a,b,c,d)

#Physical abuse
a <- ggscatter(dfOT_FND, x = "ctq_physab", y = "OXTR_CpG_sum", fill = "rs53576", color = "rs53576",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = FALSE, cor.method = "spearman",
               xlab = "CTQ - Physical Abuse", ylab = "OXTR Methylation [%]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.59, 0.56, 0.53), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.6)  # Extend y-axis to 0.6
a

b<-ggscatter(dfOT_HC, x = "ctq_physab", y = "OXTR_CpG_sum",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Abuse", ylab = "OXTR Methylation [%]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.66, 0.63, 0.60), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.7)  # Extend y-axis to 0.6
b

c<-ggscatter(dfOT_FND, x = "ctq_physab", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Abuse", ylab = "Oxytocin [pg/ml]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6
c

d<-ggscatter(dfOT_HC, x = "ctq_physab", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Physical Abuse", ylab = "Oxytocin [pg/ml]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6

d
ggarrange(a,b,c,d)

#Sexual abuse
a <- ggscatter(dfOT_FND, x = "ctq_sexab", y = "OXTR_CpG_sum", fill = "rs53576", color = "rs53576",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = FALSE, cor.method = "spearman",
               xlab = "CTQ - Sexual Abuse", ylab = "OXTR Methylation [%]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.59, 0.56, 0.53), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.6)  # Extend y-axis to 0.6
a

b<-ggscatter(dfOT_HC, x = "ctq_sexab", y = "OXTR_CpG_sum",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Sexual Abuse", ylab = "OXTR Methylation [%]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(0.66, 0.63, 0.60), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0.3, 0.7)  # Extend y-axis to 0.6
b

c<-ggscatter(dfOT_FND, x = "ctq_sexab", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Sexual Abuse", ylab = "Oxytocin [pg/ml]", title = "FND Patients") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6
c

d<-ggscatter(dfOT_HC, x = "ctq_sexab", y = "Oxytocin",fill = "rs53576", color = "rs53576",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = FALSE, cor.method = "spearman",
             xlab = "CTQ - Sexual Abuse", ylab = "Oxytocin [pg/ml]", title = "Healthy Controls") + 
  stat_cor(aes(color = rs53576), method = "spearman", 
           label.x = Inf, label.y = c(25, 23, 21), # Position labels higher
           hjust = 1.1, vjust = 1.1) +  # Use the same vjust to keep them close to the top
  scale_fill_manual(values = c("#440154FF", "#21908CFF", "grey")) + 
  scale_color_manual(values = c("#440154FF", "#21908CFF", "grey")) +
  ylim(0, 25)  # Extend y-axis to 0.6

d
ggarrange(a,b,c,d)



## SUPPLEMENT
plota <- ggplot(dfOT_FND, aes_string(x = "gender", y = "OXTR_CpG_sum", fill = "gender")) +
  geom_violin(show.legend = FALSE, alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", show.legend = FALSE, outlier.shape = NA) +
  ylim(0.2, 0.8) +
  scale_fill_manual(values = c("#440154FF", "#440154FF")) +
  labs(y = "OXTR Methylation [%]", x = "Sex", title = "<b>A</b> OXTR Methylation FND") +
  theme_classic() +
  theme(plot.title = element_markdown()) +
  stat_compare_means(label = "p.signif", comparisons = list(c("male", "female")), label.y = 0.56)
plota

plotb <- ggplot(dfOT_HC, aes_string(x = "gender", y = "OXTR_CpG_sum", fill = "gender")) +
  geom_violin(show.legend = FALSE, alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", show.legend = FALSE, outlier.shape = NA) +
  ylim(0.2, 0.8) +
  scale_fill_manual(values = c("#21908CFF", "#21908CFF")) +
  labs(y = "OXTR Methylation [%]", x = "Sex", title = "<b>B</b> OXTR Methylation HC") +
  theme_classic() +
  theme(plot.title = element_markdown()) +
  stat_compare_means(label = "p.signif", comparisons = list(c("male", "female")), label.y = 0.53)
plotb

ggarrange(plota,plotb)
