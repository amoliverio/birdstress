## Microbial abundance
library(tidyverse)
library(ggplot2)
library(gridExtra)

# data
meta = read.csv("inputs/Cultivation_Abundance_CFU_Data.csv")

# drop the sample that is missing data
meta = meta %>%
  filter(TotalCFUTSA != 'x') %>%
  mutate_at("TotalCFUTSA", as.numeric) %>%
  mutate_at("TotalCFUBA", as.numeric)
  
  
# counts
table(meta$Treatment_name) # by treatment
table(meta$Treatment_name, meta$Sex)
sum(table(meta$Treatment_name)) # total n

# by week
weekly = meta %>%
  group_by(Treatment_name) %>%
  count(TreatmentWeek)

# plots by CFU
tsa<-ggplot(meta, aes(fct_relevel(Treatment_name, c("Reference",
                                                    "Stressed",
                                                    "Recovery")), TotalCFUTSA))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ylab("Log10 CFU- TSA") +
  xlab("")
pairwise.t.test(meta$TotalCFUTSA, meta$Treatment_name, p.adjust.method = 'hochberg')

pda<-ggplot(meta, aes(fct_relevel(Treatment_name, c("Reference",
                                                    "Stressed",
                                                    "Recovery")), TotalCFUPDA))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ylab("Log10 CFU- PDA")+
  xlab("")
pairwise.t.test(meta$TotalCFUPDA, meta$Treatment_name, p.adjust.method = 'hochberg')

mac<-ggplot(meta, aes(fct_relevel(Treatment_name, c("Reference",
                                                    "Stressed",
                                                    "Recovery")), TotalCFUMAC))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ylab("Log10 CFU- Mackonkey")+
  xlab("")
pairwise.t.test(meta$TotalCFUMAC, meta$Treatment_name, p.adjust.method = 'hochberg')

myco<-ggplot(meta, aes(fct_relevel(Treatment_name, c("Reference",
                                                     "Stressed",
                                                     "Recovery")), TotalCFUMYCO))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ylab("Log10 CFU- Myco")+
  xlab("")
pairwise.t.test(meta$TotalCFUMYCO, meta$Treatment_name, p.adjust.method = 'hochberg')

ba<-ggplot(meta, aes(fct_relevel(Treatment_name, c("Reference",
                                                   "Stressed",
                                                   "Recovery")), meta$TotalCFUBA))+
  geom_boxplot()+
  scale_y_log10()+
  theme_bw()+
  ylab("Log10 CFU- BA")+
  xlab("")

grid.arrange(tsa, mac, pda, myco, ba)

#Test significance
pairwise.t.test(meta$TotalCFUBA, meta$Treatment_name, p.adjust.method = 'hochberg')
