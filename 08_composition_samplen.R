# Other calculations and tables as needed

# data
load(file="inputs/input_16S_rar_birdstress.rda")

# samples used overall 
table(input_16S_rar$map_loaded$Treatment_name)
table(input_16S_rar$map_loaded$Treatment_name, 
      input_16S_rar$map_loaded$Sex)

# samples used by week
weekly = input_16S_rar$map_loaded %>%
  group_by(Treatment_name) %>%
  count(TreatmentWeek)

# taxonomy
tax = input_16S_rar$taxonomy_loaded

## major taxonomic groups

# Captive only
input_16S_rar_captive = filter_data(input_16S_rar, filter_cat = "Treatment_name",
                                    filter_vals = 'Wild')

# phyla
bird_phyla_cap = summarize_taxonomy(input_16S_rar_captive, level = 2, relative = T)
bird_phyla_cap_means = as.data.frame(sort(rowMeans(bird_phyla_cap)))

# family
bird_family_cap = summarize_taxonomy(input_16S_rar_captive, level = 5, relative = T)
bird_family_cap_means = as.data.frame(sort(rowMeans(bird_family_cap)))

## split further by treatment and sum 
input_ref = filter_data(input_16S_rar, filter_cat = "Treatment_name",
                        keep_vals = 'Reference')
ref_means = summarize_taxonomy(input_ref, level = 5, relative = T)
ref_means = as.data.frame(sort(rowMeans(ref_means)))

input_stress = filter_data(input_16S_rar, filter_cat = "Treatment_name",
                           keep_vals = 'Stressed')
stress_means = summarize_taxonomy(input_stress, level = 5, relative = T)
stress_means = as.data.frame(sort(rowMeans(stress_means)))

input_rec = filter_data(input_16S_rar, filter_cat = "Treatment_name",
                        keep_vals = 'Recovery')
rec_means = summarize_taxonomy(input_rec, level = 5, relative = T)
rec_means = as.data.frame(sort(rowMeans(rec_means)))

# Wild only
input_16S_rar_wild = filter_data(input_16S_rar, filter_cat = "Treatment_name",
                                    keep_vals = 'Wild')
bird_family_wild = summarize_taxonomy(input_16S_rar_wild, level = 5, relative = T)
bird_family_wild_means = as.data.frame(sort(rowMeans(bird_family_wild)))

# check weight diffs
weight_diffs = input_16S_rar_captive$map_loaded %>%
  filter(Treatment_name != "Recovery") %>%
  filter(WeightChange != "Not_Measured") %>%
 # filter(BaselineCORT != "Not_Measured") %>%
  mutate_at("WeightChange", as.numeric)
#  mutate_at("BaselineCORT", as.numeric)

diff_acute = aov(weight_diffs$AcuteCORT2 ~ weight_diffs$Treatment_name)
summary(diff_acute)

diff_base = aov(weight_diffs$BaseCORT2 ~ weight_diffs$Treatment_name)
summary(diff_base)

ggplot(weight_diffs, aes(Treatment_name, AcuteCORT2)) + geom_boxplot()
