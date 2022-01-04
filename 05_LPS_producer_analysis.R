## LPS producer analysis
library(reshape)
library(ggplot2)
library(ggplot2)
library(plyr)

# cleaned data
load(file="input_16S_rar_birdstress.rda")

# summarize at family level
bird_fams = summarize_taxonomy(input_16S_rar, level = 5, relative = T) %>%
  rownames_to_column("X")

# drop whitespace
bird_fams$X = str_remove_all(bird_fams$X, " ")

##########################
#read in LPS data
fams_lps =  read.csv('outputs/bird_fams_AM.csv') %>%
  select(-X.1)

# melt data
fams_lps_m<-melt(fams_lps)

#tack on meta data
fams_lps2<-merge(fams_lps_m, meta, by.x='variable', by.y='SampleID')

#plot 
ggplot(fams_lps2,  aes(Probable_LPS.producer, value, fill=Treatment_name))+
geom_boxplot()+
xlab("LPS producer status")+
ylab("Relative abundance")+
theme_bw()

#make column that is aggregate between treatment/LPS status
fams_lps2$new<-paste(fams_lps2$Probable_LPS.producer, fams_lps2$Treatment_name, sep=',')

#Welch's t-test between categories
pairwise.t.test(fams_lps2$value, fams_lps2$new, p.adjust.method = 'hochberg')

# get means for each category
fams_lps_mean = fams_lps2 %>%
  group_by(Treatment_name, Probable_LPS.producer) %>%
  summarise(mean(value)) %>%
  rename(mean = `mean(value)`) %>%
  ungroup()

fams_lps_n = fams_lps2 %>%
  group_by(Treatment_name, Probable_LPS.producer) %>%
  tally()

fams_lps_se = fams_lps2 %>%
  group_by(Treatment_name, Probable_LPS.producer) %>%
  summarise(sd(value)) %>%
  rename(sd = `sd(value)`) %>%
  left_join(fams_lps_n) %>%
  mutate(se = sd/sqrt(n))

scale_factor = fams_lps_mean %>%
  group_by(Treatment_name) %>%
  summarise(sum(mean)) %>%
  rename(factor = `sum(mean)`)

lps = fams_lps_mean %>%
  left_join(fams_lps_se) %>%
  left_join(scale_factor) %>%
  mutate(scaled_mean = mean/factor*100,
         scaled_se = se/factor*100)

# make figure
ggplot(lps, aes(fct_relevel(Probable_LPS.producer,
                            c("Yes", "No", "Unknown")), scaled_mean, fill=Probable_LPS.producer))+
  geom_bar(stat="identity", width = 1, color = "black")+
  scale_fill_manual(values=c("black","gray","red"))+
  geom_errorbar(aes(ymin=scaled_mean, ymax=scaled_mean+scaled_se), width=.2)+
  facet_wrap(~fct_relevel(Treatment_name, c("Reference",
                                                "Stressed",
                                                "Recovery",
                                                "Wild"))) +
  xlab("LPS producer status")+
  ylab("Mean % relative abundance")+
  theme_bw()+
  theme(
    plot.title = element_text(size=12), text = element_text(size=12), 
    axis.title.x = element_text(size=12, face="bold"),
    axis.title.y = element_text( size=12, face="bold"),
    legend.position = "none")

# write out lps data for supp table
write.csv(lps, "outputs/lps_summary.csv")
