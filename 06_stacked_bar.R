## Stacked barplot of the rel. abun. of bacterial families
library(reshape)
library(tidyr)
library(stringi)
library(ggplot2)
library(RColorBrewer)

# data
load(file="input_16S_rar_birdstress.rda")

# summarize at family level
bird_fams = summarize_taxonomy(input_16S_rar, level = 5, relative = T)

# drop whitespace
rownames(bird_fams) = str_remove_all(rownames(bird_fams), " ")

# get row sums
bird_fams$sum<-rowSums(bird_fams)
bird_fams<-bird_fams[order(bird_fams$sum, decreasing=T) , ]

# top 20 families only
bird_fams<-bird_fams[1:15,]

# remove sum column
bird_fams<-as.data.frame(bird_fams[,-grep('sum', names(bird_fams))])

# get 'others' category (things not in top20)
others<-1-colSums(bird_fams)
bird_fams<-rbind(bird_fams, others)
rownames(bird_fams)[16]<-"zOthers;zOthers;zOthers;zOthers;zOthers;zOthers;zOthers"

# add taxonomy back
bird_fams$taxonomy<-row.names(bird_fams)

# melt data
fam_m<-melt(bird_fams)
names(fam_m)<-c("Taxonomy", "SampleID", "Rel_abun")

#change to percent
fam_m$Rel_abun<-fam_m$Rel_abun*100

#split taxonomy
fam_m_split<-separate(fam_m, Taxonomy, sep=";", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
fam_m_split$Family<-gsub("D_4__", "", fam_m_split$Family)

#bind metadata
fam_m_split<-merge(fam_m_split, meta, by='SampleID')

#reorder samples
fam_m_split<-fam_m_split[order(fam_m_split$SampleID, decreasing=T),]

#plot data as stacked barplot
table(meta$Treatment_name)

# palette
library(colorblindcheck)

pal = c("#eaae7f", "#d55e00", "#7f3800",
 "#e5bcd3", "#cc79a7", "#7a4864",
"#7fb8d8", "#0072b2", "#00446a",
"#f0e442", "#f7f1a0", "#908827",
"#66c4ab", "#009e73", "#005e45", 
"#f2f2f2")

colorblindcheck::palette_check(pal, plot=TRUE)

ggplot(fam_m_split, aes(SampleID, Rel_abun, fill=Family))+
  geom_bar(stat='identity')+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values=pal)+
  guides(fill=guide_legend(ncol=1))+
  xlab("")+
  ylab("% Relative Abundance")+
  theme_bw()+
  facet_wrap(~fct_relevel(Treatment_name, c("Reference",
                                            "Stressed",
                                            "Recovery",
                                            "Wild")), scales = 'free')+
  theme(text = element_text(size=14), axis.text.x = element_blank())

## Plot Phylum-level
ggplot(fam_m_split, aes(SampleID, Rel_abun, fill=Phylum))+
  geom_bar(stat='identity')+
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=pal)+
  guides(fill=guide_legend(ncol=1)) +
  xlab("") +
  ylab("% Relative Abundance") +
  theme_bw()+
  facet_wrap(~fct_relevel(Treatment_name, c("Reference",
                                            "Stressed",
                                            "Recovery",
                                            "Wild")), scales = 'free') +
  theme(text = element_text(size=14), axis.text.x = element_blank())

#identify rows with Proteobacteria (gram -) and Firmicutes (gram +)
firm<-which(fam_m_split$Phylum == "D_1__Firmicutes")
prot<- which(fam_m_split$Phylum == "D_1__Proteobacteria")

#create individual tables and then join them
firm_table<-fam_m_split[firm,]
prot_table<-fam_m_split[prot,]
firm_prot_table<-rbind(firm_table, prot_table)

#fix names
firm_prot_table$Phylum<-gsub('D_1__Proteobacteria', 'Proteobacteria', firm_prot_table$Phylum)
firm_prot_table$Phylum<-gsub('D_1__Firmicutes', 'Firmicutes', firm_prot_table$Phylum)

#plot
ggplot(firm_prot_table, aes(Treatment, Rel_abun, fill=Treatment))+
geom_boxplot()+
#geom_jitter()+
facet_wrap(~Phylum, scales='free')+
xlab("Treatment")+
ylab("Relative Abundance")+
  theme_bw()+
  theme(text = element_text(size=14), axis.text.x = element_blank())
  
#pairwise t-test for Proteobacteria
pairwise.t.test(prot_table$Rel_abun, prot_table$Treatment_name, p.adjust.method = 'hochberg')

#pairwise t-test for Firmicutes
pairwise.t.test(firm_table$Rel_abun, firm_table$Treatment_name, p.adjust.method = 'hochberg')
 