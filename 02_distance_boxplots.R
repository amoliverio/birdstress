# Calculate distances within and among groups
 
library(ggplot2)
library(mctoolsr)
library(vegan)

  ## Distance boxplots
load(file="inputs/input_16S_rar_birdstress.rda")

## Calculate distances, melt, assign categories
treatment_names = meta %>%
  select(SampleID, Treatment_name)

dm = as.matrix(calc_dm(input_16S_rar$data_loaded, method = "bray_sq_trans"))
dm[lower.tri(dm)] <- NA

dm.long = melt(dm) %>%
  filter(!is.na(value)) %>%
  left_join(treatment_names, by= c("X1"="SampleID")) %>%
  rename(Treatment_X1 = Treatment_name) %>%
  left_join(treatment_names, by= c("X2"="SampleID")) %>%
  rename(Treatment_X2 = Treatment_name) %>%
  filter(Treatment_X1 != "Wild",
         Treatment_X2 != "Wild") %>%
  unite(Treatment_Comparison, c("Treatment_X1", "Treatment_X2")) %>%
  mutate(Treatment_Comparison_Clean = case_when(Treatment_Comparison == 'Recovery_Reference' ~ 'Reference_Recovery',
                                                Treatment_Comparison == 'Stressed_Reference' ~ 'Reference_Stressed',
                                                Treatment_Comparison =='Stressed_Recovery' ~ 'Recovery_Stressed',
                                                Treatment_Comparison =='Reference_Recovery' ~ 'Reference_Recovery',
                                                Treatment_Comparison =='Reference_Stressed' ~ 'Reference_Stressed',
                                                Treatment_Comparison =='Recovery_Stressed' ~ 'Recovery_Stressed',
                                                Treatment_Comparison =='Reference_Reference' ~ 'Reference_Reference',
                                                Treatment_Comparison =='Stressed_Stressed' ~ 'Stressed_Stressed',
                                                Treatment_Comparison == 'Recovery_Recovery' ~ 'Recovery_Recovery')) %>%
  mutate(bray_similarity = 1-value) %>%
  filter(bray_similarity != 1)

# plot
ggplot(dm.long, aes(fct_reorder(Treatment_Comparison_Clean,
                                -bray_similarity,
                                .fun = mean),
                                bray_similarity)) +
  geom_jitter(alpha=0.5) + 
  geom_boxplot(outlier.shape=NA, alpha=0.75) +
  theme_classic()+
  coord_flip()+
  xlab("")+
  ylab("Bray-Curtis Similarity")

ggplot(dm.long, aes(fct_reorder(Treatment_Comparison_Clean,
                                -bray_similarity,
                                .fun = mean),
                    bray_similarity)) +
  geom_boxplot(outlier.shape=NA, alpha=0.75) +
  theme_classic()+
  coord_flip()+
  xlab("")+
  ylab("Bray-Curtis Similarity")

## Test for significance

#Welch's t-test between categories
pairwise.t.test(dm.long$bray_similarity, dm.long$Treatment_Comparison_Clean,
                p.adjust.method = 'hochberg')
