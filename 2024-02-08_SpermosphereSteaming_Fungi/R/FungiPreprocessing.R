################Data Preprocessing############

###### Libraries #####
library(phyloseq)
library(decontam)
library(vegan)
library(tidyverse)
library(metagenomeSeq)
library(ggplot2)
library(ggpubr)
library(Biostrings)
library(microbiome)

##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallets used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

########### Fungi #####

###### Read in data ####

#Loading the mapping file
samp_dat_fungi <- read.csv("2024-02-08_SpermosphereSteaming_Fungi/HPC_input_output/RawOutput/metadata/Spermosphere_Metadata_Fungi.csv", na.strings = "na")

rownames(samp_dat_fungi) <- samp_dat_fungi$Tube #row names must match OTU table headers
SAMP.fungi <- phyloseq::sample_data(samp_dat_fungi)

# OTU table 
otu <- read.csv("2024-02-08_SpermosphereSteaming_Fungi/HPC_input_output/RawOutput/otu/otu_table_ITS_Fungi.csv", na.strings = "na")
rownames(otu) <- otu$OTU
otu <- otu[,-1]
OTU.fungi <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)

any(is.na(otu)) # no NA in the OTU table

# Taxonomy
taxonomy.fungi <- read.csv("2024-02-08_SpermosphereSteaming_Fungi/HPC_input_output/RawOutput/taxonomy/NBCtaxonomy_comb.csv")
rownames(taxonomy.fungi) <- taxonomy.fungi$OTU

taxonomy.fungi$Mock_genus <- data.frame(str_split_fixed(taxonomy.fungi$Mock, "_", 4))$X3 
taxonomy.fungi$Mock_species <- data.frame(str_split_fixed(taxonomy.fungi$Mock, "_", 4))$X4 
taxonomy.fungi$Mock_LowestTaxonomicRank <- paste(data.frame(str_split_fixed(taxonomy.fungi$Mock, "_", 4))$X3, 
                                                 data.frame(str_split_fixed(taxonomy.fungi$Mock, "_", 4))$X4, sep = "_")

taxonomy.fungi$Kingdom <- ifelse(taxonomy.fungi$Mock == "", taxonomy.fungi$Kingdom, "Mocki")
taxonomy.fungi$Phylum <- ifelse(taxonomy.fungi$Mock == "", taxonomy.fungi$Phylum, "Mockimycota")
taxonomy.fungi$Class <- ifelse(taxonomy.fungi$Mock == "", taxonomy.fungi$Class, "Mockimycetes")
taxonomy.fungi$Order <- ifelse(taxonomy.fungi$Mock == "", taxonomy.fungi$Order, "Mockiales")
taxonomy.fungi$Family <- ifelse(taxonomy.fungi$Mock == "", taxonomy.fungi$Family, "Mockiaceae")
taxonomy.fungi$Genus <- ifelse(taxonomy.fungi$Mock_genus == "", taxonomy.fungi$Genus, "Mock")
taxonomy.fungi$Species <- ifelse(taxonomy.fungi$Mock_species == "", taxonomy.fungi$Species, taxonomy.fungi$Mock_species)
taxonomy.fungi$Label <- ifelse(taxonomy.fungi$Mock == "", taxonomy.fungi$Label, taxonomy.fungi$Mock)
taxonomy.fungi$Lowest_Taxonomic_Rank <- ifelse(taxonomy.fungi$Mock_LowestTaxonomicRank == "_", taxonomy.fungi$Lowest_Taxonomic_Rank, taxonomy.fungi$Mock_LowestTaxonomicRank)


taxonomy.fungi2 <- taxonomy.fungi %>%
  select(OTU, Kingdom:Label)

TAX.fungi <- phyloseq::tax_table(as.matrix(taxonomy.fungi2))

# Fasta
FASTA.fungi <- readDNAStringSet("2024-02-08_SpermosphereSteaming_Fungi/HPC_input_output/RawOutput/clustered/otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

###### Create Initial Phyloseq object ######
fungi.unedited <- phyloseq::phyloseq(OTU.fungi, TAX.fungi, FASTA.fungi, SAMP.fungi)

# Sanity check on the samples matching the metadata
samp_dat_fungi$Tube # 148 samples
colnames(otu) # 148 samples
sample_names(fungi.unedited) # 148 samples

setdiff(samp_dat_fungi$Tube, colnames(otu)) # what is in the metadata but not in the otu
setdiff(colnames(otu), samp_dat_fungi$Tube) # what is in the otu table but not in the metadata


###### Decontaminate  ######
fungi.unedited@sam_data$Sample_or_Control <- ifelse(fungi.unedited@sam_data$Type == "Negative Control", "Control Sample", "True Sample")
sample_data(fungi.unedited)$is.neg <- sample_data(fungi.unedited)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(fungi.unedited, method="prevalence", neg="is.neg", threshold = 0.1, normalize = TRUE)
badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

print(badTaxa)

ps.pa <- transform_sample_counts(fungi.unedited, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontaminate <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") + 
  ggtitle("Fungi") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette)

goodTaxa <- setdiff(taxa_names(fungi.unedited), badTaxa)
fungi_sub_no_bad <- prune_taxa(goodTaxa, fungi.unedited)

###### Mock Community #######
fungi_mock <- fungi_sub_no_bad %>% 
  subset_samples(Type == "Positive Control") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) 

mock2 <- microbiome::transform(fungi_mock, "compositional") # relative abundance transform

sequenced.mock.fungi <- mock2 %>%
  psmelt() %>% 
  mutate(Label_plot = ifelse(Kingdom == "Mocki", Lowest_Taxonomic_Rank, "Other")) %>%
  ggplot(aes(Sample, Abundance, fill = Label_plot)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb, "violet", "pink", "grey", "black", "blue", "green")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Fungi") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm')) 
sequenced.mock.fungi

mock.composition <- mock2 %>%
  psmelt() %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund) %>%
  ungroup() %>%
  group_by(Kingdom) %>%
  summarise(SumAbund = sum(MeanRelAbund)) 
mock.composition
 
# OTUs classified into the mock kingdom made up 99.9% of the reads in mock samples

fungi_not_mock <- fungi_sub_no_bad %>% 
  subset_samples(Type != "Positive Control") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) 

not_mock <- microbiome::transform(fungi_not_mock, "compositional") # relative abundance transform

not.mock.composition <- not_mock %>%
  psmelt() %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund) %>%
  ungroup() %>%
  group_by(Kingdom) %>%
  summarise(SumAbund = sum(MeanRelAbund)) %>%
  arrange(-SumAbund)
not.mock.composition
sum(not.mock.composition$SumAbund[c(2,4:9)])

# OTUs classified into the Fungal kingdom made up 98.1% of the reads in real samples

###### Taxonomy and Sample filtering #####
# remove OTUs that are not fungi, or unidentified at the kingdom level 
fungi_sperm <- fungi_sub_no_bad %>% 
  phyloseq::subset_taxa(Kingdom %in% c("Fungi")) %>%
  subset_samples(Type %in% c("Bulk.Soil", "Spermosphere", "Rhizosphere", "Soybean Epiphytes")) %>%
  prune_samples(sample_sums(.) > 1000, .) %>% # remove samples below 1,000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa not present in samples cut out

###### RDS of Non-normalized Fungi data #####
# Save an object to a file
saveRDS(fungi_sperm, file = "2024-02-08_SpermosphereSteaming_Fungi/RDS/Fungi_spermosphere_nonorm_2024-03-13.rds")
# Restore the object
fungi.unedited <- readRDS(file = "2024-02-08_SpermosphereSteaming_Fungi/RDS/Fungi_spermosphere_nonorm_2024-03-13.rds")

fungi.unedited@sam_data

###### READS PER SAMPLE #######
sample.sums <- data.frame(sample_sums(fungi.unedited))

read.dist <- ggplot(sample.sums, aes(x = sample_sums.fungi.unedited.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") +
  ylab("Number of samples") + 
  ggtitle("Fungi")

sum(taxa_sums(fungi.unedited)) # total reads = 8016866

mean(sample_sums(fungi.unedited)) # 57762.67
median(sample_sums(fungi.unedited)) # 51784 reads

######## Rarefaction analysis ######## 
sam.data <- data.frame(fungi.unedited@sam_data)
fOTU.table <- otu_table(fungi.unedited) %>%
  as.data.frame() %>%
  as.matrix()
rare.fun <- rarecurve(t(fOTU.table), step = 1000, sample = raremax, tidy = TRUE)

fungi.rare.curve.extract2 <- left_join(rare.fun, sam.data, by = c("Site" = "Tube"))

fungi.rare <- ggplot(fungi.rare.curve.extract2, aes(x = Sample, y = Species, group = Site, color = Trt)) + 
  #geom_point() +
  scale_color_manual(values = cbbPalette)+
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") + 
  ggtitle("Fungi")+
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(fungi.unedited)), linetype = "dashed") +
  ggtitle("") 

######### Metagenome CSS normalization #########
MGS <- phyloseq_to_metagenomeSeq(fungi.unedited)
p <- metagenomeSeq::cumNormStatFast(MGS)

MGS <- metagenomeSeq::cumNorm(MGS, p =p)

metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample

norm.fungi <- metagenomeSeq::MRcounts(MGS, norm = T)

norm.fungi.OTU <- phyloseq::otu_table(norm.fungi, taxa_are_rows = TRUE)

fungi.css.norm <- phyloseq::phyloseq(norm.fungi.OTU, TAX.fungi, FASTA.fungi, SAMP.fungi)

######## Save CSS object to a file ########
saveRDS(fungi.css.norm, file = "2024-02-08_SpermosphereSteaming_Fungi/RDS/Fungi_spermosphere_CSS_2024-03-13.rds")
# Restore the object
fungi.css.norm <- readRDS(file = "2024-02-08_SpermosphereSteaming_Fungi/RDS/Fungi_spermosphere_CSS_2024-03-13.rds")


###### Supplemental Figure 1 ######
supp.fig.1 <- ggpubr::ggarrange(sequenced.mock.fungi, 
                                fungi.rare, 
                                decontaminate, 
                                read.dist, nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))

fungi.css.norm.filt <- fungi.css.norm %>% 
  subset_samples(Type %in% c("Spermosphere", "Soybean Epiphytes", "Rhizosphere")) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa not present in samples cut out

fungi.css.norm.filt@sam_data$TimeXTrt <- interaction(fungi.css.norm.filt@sam_data$Time, fungi.css.norm.filt@sam_data$Trt)

bac.bray = phyloseq::distance(fungi.css.norm.filt, "bray") # create bray-curtis distance matrix

fungi.css.norm@sam_data
# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(bac.bray~Trt*Time, as(sample_data(fungi.css.norm.filt), "data.frame"))

global.ord.bray <- phyloseq::ordinate(fungi.css.norm.filt, "NMDS", "bray")
p1 = plot_ordination(fungi.css.norm.filt, global.ord.bray, type="samples", color="TimeXTrt", shape = "Type") + 
  stat_ellipse()
variation.axis1.bray <- round(100*global.ord.bray$values$Relative_eig[[1]], 2) # % variation on axis one
variation.axis2.bray <- round(100*global.ord.bray$values$Relative_eig[[2]], 2) # % variation on axis two

bac.ggplot.data.bray <- data.frame(global.ord.bray$vectors) # get the point values
bac.ggplot.data.bray$Tube <- rownames(bac.ggplot.data.bray) # add rownames
bac.ggplot.data.bray2 <- left_join(data.frame(fungi.css.norm.filt@sam_data), bac.ggplot.data.bray,  by = "Tube") # add rest of metadata

global.bray <- ggplot(bac.ggplot.data.bray2, aes(x = Axis.1, y = Axis.2, fill = Time, shape = Trt)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24), name = "Treatment")+
  scale_fill_manual(values=cbbPalette, name = "Environment")+
  xlab(paste("PcoA1 -", variation.axis1.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.bray, "%")) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  theme_bw() +
  ggtitle("Bray-Curtis")

# Rarefaction RDS Creation
ps.rarefied <- rarefy_even_depth(fungi.unedited, rngseed=12345, sample.size=0.9*min(sample_sums(fungi.unedited)), replace=T)
ps.rarefied

saveRDS(ps.rarefied, file = "2024-02-08_SpermosphereSteaming_Fungi/RDS/Fungi_spermosphere_Rarefied_2024-06-04.rds")












