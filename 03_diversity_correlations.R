## Correlations between diversity and meta data
### Alpha diversity

#read in data
load(file="inputs/input_16S_rar_birdstress.rda")

# subsample map to only include experimental birds, with measurements
input_captive = filter_data(input_16S_rar, "CaptiveorNot", keep_vals = 'Captive')
input_captive = filter_data(input_captive, "WeightChange", filter_vals = 'Not_Measured')

s16_exp = as.data.frame(t(input_captive$data_loaded)) 

meta_exp = input_captive$map_loaded %>%
  rownames_to_column("SampleID") %>%
  select(SampleID, 
         DaysTotalCaptivity,
         RelativeSpleenWt,
         WeightChange,
         AcuteCORT2,
         BaseCORT2)
  
# Shannon
s16.shan<-diversity(s16_exp, index='shannon')

#OTUs observed
s16.otus<-rowSums(s16_exp>0)

#Pielou's Evenness
s16.even<-(diversity(s16_exp))/s16.otus

#concatenate data into a single frame
s16.div<-as.data.frame(cbind(s16.shan, s16.otus, s16.even))

#add sample IDs
s16.div$SampleID<-row.names(s16.div)

#fix coloumn names
names(s16.div)<-c('Shannon', 'OTUs_Obs', 'Pielous_Even', 'SampleID')

#add metadata
s16.div<-merge(s16.div, meta_exp, by='SampleID') 
s16.div$WeightChange = as.numeric(s16.div$WeightChange)

#test correlations between diversity and parameters
cor.test(s16.div$DaysTotalCaptivity, s16.div$Shannon, method='spearman')
cor.test(s16.div$RelativeSpleenWt, s16.div$Shannon, method='spearman')
cor.test(s16.div$WeightChange, s16.div$Shannon, method='spearman')
cor.test(s16.div$AcuteCORT2, s16.div$Shannon, method='spearman')
cor.test(s16.div$BaseCORT2, s16.div$Shannon, method='spearman')

#test between sexes
t.test(s16.div$Shannon ~ s16.div$Sex)

### Beta diversity

#calculate a PCoA
set.seed(1234)
s16.pcoa<-capscale(s16_exp ~ 1, distance='bray')

#extract the coordinates
s16.scores<-scores(s16.pcoa)
s16.coords<-as.data.frame(s16.scores$sites)
s16.coords$SampleID<-rownames(s16.coords)
s16.coords<-merge(s16.coords, meta, by=c('SampleID'))

#test correlations between PCOA and parameters  # note: default is drop NA values
s16.div$MDS1 = s16.coords$MDS1
cor.test(s16.div$DaysTotalCaptivity, s16.div$MDS1, method='spearman')
cor.test(s16.div$RelativeSpleenWt, s16.div$MDS1, method='spearman')
cor.test(s16.div$WeightChange, s16.div$MDS1, method='spearman')
cor.test(s16.div$AcuteCORT2, s16.div$MDS1, method='spearman') # SIG
cor.test(s16.div$BaseCORT2, s16.div$MDS1, method='spearman')

# s16 to plot
s16_to_plot = s16.div %>%
  select(SampleID,
         MDS1,
         Shannon,
         DaysTotalCaptivity,
         RelativeSpleenWt,
         WeightChange,
         AcuteCORT2,
         BaseCORT2) %>%
  pivot_longer(!c("SampleID", "MDS1", "Shannon"))

plt.pcoa_cor = ggplot(s16_to_plot, aes(value, MDS1)) + 
  geom_point(pch=21, fill="gray", color="black", alpha=1) +
  facet_wrap(~name, scales = "free") + 
  theme_classic() + 
  labs(x="Factor", y= "PCoA Axis 1")

plt.shannon_cor = ggplot(s16_to_plot, aes(value, Shannon)) + 
  geom_point(pch=21, fill="gray", color="black", alpha=1) +
  facet_wrap(~name, scales = "free") + 
  theme_classic() + 
  labs(x="Factor", y= "Shannon diversity")


