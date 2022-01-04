## Format data for analyses, rarefy, community composition

library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
library(mctoolsr)
library(metagenomeSeq)
library(gridExtra)


# color guide for figures:
# Reference = Blue, Recovery = Purple, Stressed = Red, and Wild = Black

# treatment group order
# Reference, Stressed, Recovery, Wild

# input data AMO
tax_table_fp = 'inputs/16S_table_tax_filt_AMO.txt' #formatted for input
map_fp = 'inputs/bird_meta_AMO.txt' #formatted for input
input_16S = load_taxa_table(tax_table_fp, map_fp)

# check total N
metadata = input_16S$map_loaded
table(metadata$Treatment) # by treatment
sum(table(metadata$Treatment)) # total n

# by week
weekly = metadata %>%
  group_by(Treatment) %>%
  count(TreatmentWeek)

# drop unassigned otus
input_16S <- filter_taxa_from_input(input_16S, at_spec_level = 1, taxa_to_remove = "Unassigned")
input_16S <- filter_taxa_from_input(input_16S, at_spec_level = 2, taxa_to_remove = "NA")

# clean up and drop low seq samples
sort(colSums(input_16S$data_loaded))
input_16s_filt = filter_samples_by_counts(input_16S, min_seqs = 4608)

# rarefy and normalize
set.seed(1234)
input_16S_rar = single_rarefy(input_16s_filt, 4608)
input_16S_rar = convert_to_relative_abundances(input_16S_rar) 
input_16S_rar$map_loaded = input_16S_rar$map_loaded %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(Treatment_name = case_when(Treatment == 'Control' ~ 'Reference',
                                    Treatment == 'Experimental' ~ 'Stressed',
                                    Treatment =='No Tx' ~ 'Wild',
                                    Treatment == 'Recovery' ~ 'Recovery')) %>%
  mutate(Treatment_name = as.factor(Treatment_name))
  
# change back to initial structure
s16 = as.data.frame(t(input_16S_rar$data_loaded))
meta = input_16S_rar$map_loaded %>%
  rownames_to_column("SampleID")

# save as R object for other analyses
save(input_16S_rar, s16, meta, file="input_16S_rar_birdstress.rda")

###############################################################################
###############################################################################

load(file="inputs/input_16S_rar_birdstress.rda")
tax_num = input_16S_rar$taxonomy_loaded

#calculate PCoA based on Bray-curtis and calculate the percent variation explained
set.seed(1234)
s16.pcoa<-capscale(s16 ~ 1, distance='bray')
100*round(s16.pcoa$CA$eig[1]/sum(s16.pcoa$CA$eig), 3)
100*round(s16.pcoa$CA$eig[2]/sum(s16.pcoa$CA$eig), 3)

#extract the coordinates
s16.scores<-scores(s16.pcoa)
s16.coords<-as.data.frame(s16.scores$sites)
s16.coords$SampleID<-rownames(s16.coords)
s16.coords<-merge(s16.coords, meta, by=c('SampleID'))

#plot it
ggplot(s16.coords, aes(MDS1, MDS2, fill=Treatment_name))+  
  geom_point(size=3, pch=21, stroke=0.75, color="#262626")+
  theme_classic()+
  xlab("PCoA 1 (31.9%)")+
  ylab("PCoA 2 (17.1%)")+
  scale_fill_manual(values = c("purple", "blue", "red", "black")) +
  theme(text = element_text(size=14),
        axis.text = element_text(size=14), legend.text=element_text(size=14)) +
  labs(fill='Treatment') 

# Reference = Blue, Recovery = Purple, Stressed = Red, and Wild = Black

## Adonis

# overall treatment
dm.16S = calc_dm(input_16S_rar$data_loaded, method = "bray_sq_trans")
adonis(dm.16S ~ Treatment, data=input_16S_rar$map_loaded, permutations = 10000)

# pairwise permanovas
set.seed(1234)
perm.pair.treat = calc_pairwise_permanovas(dm.16S, input_16S_rar$map_loaded, "Treatment_name", 10000)

# Captive: stress vs. reference
input_16S_rar_exp_notx = filter_data(input_16S_rar, filter_cat = "Treatment_name", keep_vals = c("Stressed", "Reference"))
dm.16S.expref = calc_dm(input_16S_rar_exp_notx$data_loaded, method = "bray_sq_trans")
adonis(dm.16S.expref ~ Treatment, data=input_16S_rar_exp_notx$map_loaded, permutations = 10000)

# Captive: recovery vs. reference
input_16S_rar_rec_notx = filter_data(input_16S_rar, filter_cat = "Treatment_name", keep_vals = c("Recovery", "Reference"))
dm.16S.recref = calc_dm(input_16S_rar_rec_notx$data_loaded, method = "bray_sq_trans")
adonis(dm.16S.recref ~ Treatment, data=input_16S_rar_rec_notx$map_loaded, permutations = 10000)

# Captive: recovery vs. stress
input_16S_rar_exp_rec = filter_data(input_16S_rar, filter_cat = "Treatment_name", keep_vals = c("Stressed", "Recovery"))
dm.16S.exprec = calc_dm(input_16S_rar_exp_rec$data_loaded, method = "bray_sq_trans")
adonis(dm.16S.exprec ~ Treatment, data=input_16S_rar_exp_rec$map_loaded, permutations = 10000)

## Alpha Diversity

#Shannon
s16.shan<-diversity(s16, index='shannon')

#OTUs observed
s16.otus<-rowSums(s16>0)

#Pielou's Evenness
s16.even<-(diversity(s16))/s16.otus

#catenate data into a single frame
s16.div<-as.data.frame(cbind(s16.shan, s16.otus, s16.even))

#add sample IDs
s16.div$SampleID<-row.names(s16.div)

#fix coloumn names
names(s16.div)<-c('Shannon', 'OTUs_Obs', 'Pielous_Even', 'SampleID')

#add metadata
s16.div<-merge(s16.div, meta, by='SampleID')


shan = ggplot(s16.div, aes(fct_relevel(Treatment_name, c("Reference",
                                                        "Stressed",
                                                        "Recovery",
                                                        "Wild")), Shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3, pch=21, stroke=0.75, color="#262626", aes(fill=Treatment_name))+
  scale_fill_manual(values = c("purple", "blue", "red", "black")) +
  theme_classic() + 
  labs(x="Treatment", y="Shannon Diversity", fill='Treatment')
  
even = ggplot(s16.div, aes(fct_relevel(Treatment_name, c("Reference",
                                                         "Stressed",
                                                         "Recovery",
                                                         "Wild")), Pielous_Even)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size=3, pch=21, stroke=0.75, color="#262626", aes(fill=Treatment_name))+
  scale_fill_manual(values = c("purple", "blue", "red", "black")) +
  theme_classic() + 
  labs(x="Treatment", y="Pielou's Evenness", fill='Treatment')

grid.arrange(shan, even, ncol=2)

## Test significance of alpha diversity

#t-test for significance
pairwise.t.test(s16.div$Shannon, s16.div$Treatment_name, p.adjust.method = 'hochberg')
pairwise.t.test(s16.div$Pielous_Even, s16.div$Treatment_name, p.adjust.method = 'hochberg')
