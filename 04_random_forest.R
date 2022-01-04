## Random forest of OTUs

library("randomForest")
library("plyr")
library("rfUtilities")
library("caret")
library("vegan")
library("reshape")
library("ggplot2")

# read in cleaned data
load(file="inputs/input_16S_rar_birdstress.rda")
s16 = as.data.frame(t(s16))

#how many nonzero counts?
otu_nonzero_counts<-apply(s16, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

#remove rare taxa
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=s16, cutoff_pro=0.1)
dim(otu_table_rare_removed)

#scale data
otu_table_scaled <- scale(otu_table_rare_removed, center = TRUE, scale = TRUE)  

#run RF model
otu_table_scaled_treatment <- data.frame(t(otu_table_scaled))
otu_table_scaled_treatment$SampleID<-row.names(otu_table_scaled_treatment)
meta_sub = meta %>%
  select(SampleID, Treatment_name)
names(meta_sub)<-c("SampleID", "Treatment_name")
otu_table_scaled_treatment<-merge(otu_table_scaled_treatment, meta_sub, by=c("SampleID"))
otu_table_scaled_treatment<-otu_table_scaled_treatment[,-1]

predictors= otu_table_scaled_treatment[,1:(ncol(otu_table_scaled_treatment)-1)]

set.seed(1234)
RF_treatment_classify<-randomForest(x=predictors,
                                    y=as.factor(otu_table_scaled_treatment$Treatment_name),
                                    ntree=501, importance=TRUE, proximities=TRUE)

#permutation test
RF_treatment_classify_sig<-rf.significance(x=RF_treatment_classify,
                                           xdata=predictors,
                                           nperm=1000,
                                           ntree=501)

#identifying important features
RF_state_classify_imp <- as.data.frame(RF_treatment_classify$importance)
RF_state_classify_imp$features <- rownames( RF_state_classify_imp )
RF_state_classify_imp_sorted <- arrange( RF_state_classify_imp  , desc(MeanDecreaseAccuracy)  )
barplot(RF_state_classify_imp_sorted$MeanDecreaseAccuracy,
        ylab="Mean Decrease in Accuracy (Variable Importance)",
        main="RF Classification Variable Importance Distribution")

#top 20 features
barplot(RF_state_classify_imp_sorted[1:20,"MeanDecreaseAccuracy"],
        las=2, names.arg=RF_state_classify_imp_sorted[1:20,"features"] ,
        ylab="Mean Decrease in Accuracy (Variable Importance)", main="Classification RF")  


## Plot mean abundance of top 20 most important OTUs
top20_feat<-as.data.frame(RF_state_classify_imp_sorted$features[1:20])
names(top20_feat)<-c("OTU")
str(top20_feat)

#read in OTU table, extract taxonomy
tax<-input_16S_rar$taxonomy_loaded
OTU<-row.names(s16)
s16_tax<-cbind(OTU, tax)

#extract top 20 from OTU table
s16.top20<-s16[rownames(s16) %in% top20_feat$OTU,]
dim(s16.top20)  
s16.top20$OTU<-row.names(s16.top20)

#get mean relative abundance of each OTU in the 4 cats
top20_m<-melt(s16.top20)
names(top20_m)<-c("OTU", "SampleID", 'Rel_abund')
top20_m<-merge(top20_m, meta, by='SampleID')
top20_m<-merge(top20_m, s16_tax, by='OTU')

#split taxonomy into groups
top20_m = top20_m %>%
  dplyr::rename(Kingdom = taxonomy1,
         Phylum = taxonomy2,
         Class = taxonomy3,
         Order = taxonomy4,
         Family = taxonomy5,
         Genus = taxonomy6,
         Species_bac = taxonomy7)

#fix the NAs and taxonomy
top20_m$Genus[is.na(top20_m$Genus)] <- "Unassigned"
top20_m$Family[is.na(top20_m$Family)] <- "Unassigned"
top20_m$Order[is.na(top20_m$Order)] <- "Unassigned"
top20_m$Class[is.na(top20_m$Class)] <- "Unassigned"
top20_m$Genus<-gsub('D_5__', '', top20_m$Genus)
top20_m$Family<-gsub('D_4__', '', top20_m$Family)
top20_m$Order<-gsub('D_3__', '', top20_m$Order)
top20_m$Class<-gsub('D_2__', '', top20_m$Class)
top20_m$Phylum<-gsub('D_1__', '', top20_m$Phylum)

top20_m_labels = top20_m %>%
  unite("genus_otu", Genus, OTU)

# graph labels
table(top20_m_labels$genus_otu)

otu_names <- c(
  `Acinetobacter_EF517956.1.1666` = "Acinetobacter 1",
  `Acinetobacter_JN082536.1.1536` = "Acinetobacter 2",
  `Acinetobacter_AB365066.1.1533` = "Acinetobacter 3",
  `Anaerostipes_DQ905930.1.1794` = "Anaerostipes",
  `Campylobacter_EU559331.1.1470` = "Campylobacter",
  `Catellicoccus_FJ192638.1.1515` = "Catellicoccus 1",
  `Catellicoccus_KF799139.1.1524` = "Catellicoccus 2",
  `Chryseobacterium_JPLY01000001.145690.147219` = "Chryseobacterium",
  `Clostridium sensu stricto 1_AF018036.1.1512` = "Clostridium",
  `Collinsella_DQ798456.1.1292` = "Collinsella",
  `Escherichia-Shigella_CCPS01000022.154.1916` = "Escherichia-Shigella",
  `Glutamicibacter_AF511517.1.1557` = "Glutamicibacter",
  `Lactobacillus_AF197125.1.1555` = "Lactobacillus 1",
  `Lactobacillus_KF178310.1.1559` = "Lactobacillus 2",
  `Massilia_CP012201.3677670.3679209` = "Massilia",
  `Pseudomonas_KJ161326.1.1708` = "Pseudomonas 1",
  `Pseudomonas_KJ535378.1.1545` = "Pseudomonas 2",
  `Serratia_KF625184.1.1787` = "Romboutsia",
  `Streptococcus_CDMW01000001.16532.18068` = "Streptococcus",
  `Unassigned_AMYT01000015.54.1603` = "Enterococcaceae")

# boxplots by OTU
ggplot(top20_m_labels, aes(fct_relevel(Treatment_name, c("Reference",
                                                         "Stressed",
                                                         "Recovery",
                                                         "Wild")), Rel_abund)) +
  geom_jitter(pch=21, aes(fill=Treatment_name), alpha=0.75) +
  geom_boxplot(alpha=0.25, outlier.shape = NA) + 
  facet_wrap(~genus_otu, scales="free",
             labeller = as_labeller(otu_names)) +
  scale_color_manual(values = c("purple", "blue", "red", "black")) +
  scale_fill_manual(values = c("purple", "blue", "red", "black")) +
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  labs(y= "Relative Abundance",
       fill="Treatment") +
  theme(axis.line = element_line(color='white'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

# select top 20 
ref = top20_m %>%
  select(OTU) %>%
  distinct()

# what % of comm do the top 20 comprise?
input_rf = filter_taxa_from_input(input_16S_rar,
                                  taxa_IDs_to_keep = ref$OTU)
mean(colSums(input_rf$data_loaded))

# what taxa?
tax_rf = input_rf$taxonomy_loaded
        
        
        
        