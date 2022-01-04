## Microbial abundance
library(ggplot2)
library(gridExtra)

# data
load(file="input_16S_rar_birdstress.rda")
meta = read.delim("bird_meta_AMO.txt")

# clean, add treat cat
meta = meta %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(Treatment_name = case_when(Treatment == 'Control' ~ 'Reference',
                                    Treatment == 'Experimental' ~ 'Stressed',
                                    Treatment =='No Tx' ~ 'Wild',
                                    Treatment == 'Recovery' ~ 'Recovery')) %>%
  mutate(Treatment_name = as.factor(Treatment_name))

# subset (drop 'wild')
meta = subset(meta, Treatment_name == 'Reference' |
                Treatment_name == 'Stressed' |
                Treatment_name == 'Recovery')

# counts
table(meta$Treatment_name) # by treatment
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
