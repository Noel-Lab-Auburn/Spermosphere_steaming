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


########### Bacteria #####

###### Read in data ####

# Metadata #
samp_dat_bac <- read.csv("metadata/Spermosphere_Metadata.csv", na.strings = "na")

rownames(samp_dat_bac) <- samp_dat_bac$Tube #row names must match OTU table headers
SAMP.bac <- phyloseq::sample_data(samp_dat_bac)

# OTU table #
otu_bac <- read.csv("otu_table/otu_table_16S_Bacteria.csv")
rownames(otu_bac) <- otu_bac$OTU
otu_bac <- otu_bac[,-1]
OTU.bac <- phyloseq::otu_table(otu_bac, taxa_are_rows = TRUE)

any(is.na(otu_bac)) # no NA in the OTU table

# Taxonomy #
taxonomy.bac <- read.csv("taxonomy/16s_taxonomy.csv")
rownames(taxonomy.bac) <- taxonomy.bac$OTU
TAX.bac <- phyloseq::tax_table(as.matrix(taxonomy.bac))

all.equal(rownames(samp_dat_bac), colnames(otu_bac))

# Fasta #
FASTA.bac <- readDNAStringSet("clustered/otus.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)

# Phylogentic tree #
tree <- phyloseq::read_tree("tree/otu_tree.tre")

###### Create Initial Phyloseq object #####
# Merge reads into Phyloseq object #
bac.unedited <- phyloseq::phyloseq(OTU.bac, TAX.bac, FASTA.bac, SAMP.bac, tree)

###### Skipped decontamination for now until I know what is what ####
###### Decontaminate #####
bac.unedited@sam_data$Sample_or_Control <- ifelse(bac.unedited@sam_data$Crop %in% c("NEC", "Water"), "Control Sample", "True Sample")
sample_data(bac.unedited)$is.neg <- sample_data(bac.unedited)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(bac.unedited, method="prevalence", neg="is.neg", threshold = 0.1, normalize = TRUE)
badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

print(badTaxa)

ps.pa <- transform_sample_counts(bac.unedited, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontaminate.bac <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") + 
  scale_color_manual(values = cbbPalette)+ 
  ggtitle("Prokaryote") +
  theme_classic()

goodTaxa <- setdiff(taxa_names(bac.unedited), badTaxa)
bac_sub_no_bad <- prune_taxa(goodTaxa, bac.unedited)


###### Taxonomy filtering #####
# remove OTUs that are mitochondria, chloroplast, or unidentified at the kingdom level 
bac_no_chloro <- bac.unedited %>% 
  phyloseq::subset_taxa(Order != "Chloroplast") %>%
  phyloseq::subset_taxa(Family != "Mitochondria") %>%
  phyloseq::subset_taxa(Kingdom != "unidentified")

###### Also not doing mock analsysi right now ####
###### Mock Community analysis ##### 
# positive controls
bac_mock <- bac_no_chloro %>% 
  subset_samples(Crop == "MOCK") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 2, TRUE) # filter OTUs to have more than 1 read in mock samples

mock2 <- microbiome::transform(bac_mock, "compositional") # relative abundance transform

sequenced.mock.bac <- mock2 %>%
  psmelt() %>% 
  ggplot(aes(Sample, Abundance, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Prokaryote") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm')) 
sequenced.mock.bac

# Adding in theoretical distribution - the last two are fungi and are not expected to be amplified with 16S
Label <- c("Pseudomonas aeruginosa", 
           "Escherichia coli",
           "Salmonella enterica", 
           "Lactobacillus fermentum", 
           "Enterococcus faecalis", 
           "Staphylococcus aureus", 
           "Listeria monocytogenes", 
           "Bacillus subtilis")

# theoretical species composition in the mock community
Abundance <- c(rep(0.125, 8))

th.mock <- data.frame(Label, Abundance)
th.mock$Sample <- "Theoretical"

th.mock$Label <- factor(th.mock$Label, levels = c("Lactobacillus fermentum", 
                                                  "Staphylococcus aureus", 
                                                  "Bacillus subtilis",
                                                  "Escherichia coli",
                                                  "Listeria monocytogenes",
                                                  "Enterococcus faecalis",
                                                  "Salmonella enterica",
                                                  "Pseudomonas aeruginosa"))


theory.mock <- ggplot(th.mock, aes(Sample, Abundance, fill = Label)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette[[1]], 
                              cbbPalette[[2]], 
                              cbbPalette[[3]], 
                              cbbPalette[[4]], 
                              cbbPalette[[5]],
                              cbbPalette[[6]],
                              cbbPalette[[8]],
                              "violet", "pink", "grey", "black", "blue")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative abundance (%)",
       title = "Theoretical composition") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic"),
        legend.title = element_blank())

# I think maybe the theoretical mock community can also be mentioned in the figure legend. 

mock.composition <- mock2 %>%
  psmelt() %>%
  group_by(Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  arrange(-MeanRelAbund)

# these 8 OTUs made up 99.9% of the mock composition. These OTUs also match the 8 supposed to be in the mock
sum(mock.composition[1:8,]$MeanRelAbund)

###### Data filtering #####
# remove samples with less than 1000 reads

bac_sperm <- bac_no_chloro %>% 
  subset_samples(Trt %in% c("Non.Steamed", "Steamed ", "Soybean Epiphytes")) %>%
  prune_samples(sample_sums(.) > 5000, .) %>% # remove samples below 10000 reads
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) # remove taxa with less than 1 reads

###### RDS of Non-normalized Prokaryote data ######
# Save an object to a file
saveRDS(bac_sperm, file = "Bacteria_spermosphere_nonnorm_2023-01-24.rds")
# Restore the object
bac_sperm <- readRDS(file = "Bacteria_spermosphere_nonnorm_2023-01-24.rds")

###### READS PER SAMPLE ######
sample.sums <- data.frame(sample_sums(bac_sperm))

read.dist.bac <- ggplot(sample.sums, aes(x = sample_sums.bac_sperm.)) +
  geom_histogram(color = "black", fill = cbbPalette[[4]]) + 
  theme_classic() +
  xlab("Read Depth") + 
  ggtitle("Prokaryote")

sum(sample_sums(bac_sperm)) # total reads = 4,389,027
median(sample_sums(bac_sperm)) # 35,522

###### Rarefaction analysis #####
sam.data <- data.frame(bac_sperm@sam_data)
bOTU.table <- otu_table(bac_sperm) %>%
  as.data.frame() %>%
  as.matrix()

raremax <- min(rowSums(t(bOTU.table)))
rare.fun <- rarecurve(t(bOTU.table), step = 1000, sample = raremax, tidy = T)

bac.rare.curve.extract2 <- left_join(sam.data, rare.fun, by = c("Tube" = "Site"))

bac.rare <- ggplot(bac.rare.curve.extract2, aes(x = Sample, y = Species, group = Tube, color = Trt)) + 
  #geom_point() +
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") +
  ggtitle("Prokaryote") +
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(bac_sperm)), linetype = "dashed") +
  scale_color_manual(values = cbbPalette)

###### Metagenome CSS normalization ######
MGS <- phyloseq_to_metagenomeSeq(bac_sperm) #converts to metagenomeseq format
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.bac <- metagenomeSeq::MRcounts(MGS, norm = T) 
norm.bac.OTU <- phyloseq::otu_table(norm.bac, taxa_are_rows = TRUE) #exports the new otu table
bac.css.norm <- phyloseq::phyloseq(norm.bac.OTU, FASTA.bac, SAMP.bac, TAX.bac, tree) #new otu table phyloseq object

saveRDS(bac.css.norm, file = "Bacteria_spermosphere_nonnorm_CSS_2023-01-24.rds")
# Restore the object
bac.css.norm <- readRDS(file = "Bacteria_spermosphere_nonnorm_CSS_2023-01-24.rds")

######## PROKARYOTE ########
###### Beta-diversity using bray-curtis ########

##### All Samples; Bray #####
bac.bray = phyloseq::distance(bac.css.norm, "bray") # create bray-curtis distance matrix

# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(bac.bray~Trt*Type, as(sample_data(bac.css.norm), "data.frame"))

global.ord.bray <- phyloseq::ordinate(bac.css.norm, "MDS", "bray")
variation.axis1.bray <- round(100*global.ord.bray$values$Relative_eig[[1]], 2) # % variation on axis one
variation.axis2.bray <- round(100*global.ord.bray$values$Relative_eig[[2]], 2) # % variation on axis two

bac.ggplot.data.bray <- data.frame(global.ord.bray$vectors) # get the point values
bac.ggplot.data.bray$Tube <- rownames(bac.ggplot.data.bray) # add rownames
bac.ggplot.data.bray2 <- left_join(data.frame(bac.css.norm@sam_data), bac.ggplot.data.bray,  by = "Tube") # add rest of metadata

global.bray <- ggplot(bac.ggplot.data.bray2, aes(x = Axis.1, y = Axis.2, fill = Trt, shape = Type)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24), name = "Treatment")+
  scale_fill_manual(values=cbbPalette, name = "Environment")+
  xlab(paste("PcoA1 -", variation.axis1.bray, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2.bray, "%")) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  theme_bw() +
  ggtitle("Bray-Curtis")

###### Beta-Diversity using Weighted Unifrac ####

##### All Samples ########
bac.unifrac = UniFrac(bac.css.norm, weighted = T) # create weighted unifrac distance matrix

# PERMANOVA - testing for differences in centroids
set.seed(12325)
adonis2(bac.unifrac~Trt*Type, as(sample_data(bac.css.norm), "data.frame"))

global.ord.uni <- phyloseq::ordinate(bac.css.norm, "MDS", "unifrac", weighted = TRUE)
variation.axis1 <- round(100*global.ord.uni$values$Relative_eig[[1]], 2)
variation.axis2 <- round(100*global.ord.uni$values$Relative_eig[[2]], 2)

bac.ggplot.data.uni <- data.frame(global.ord.uni$vectors)
bac.ggplot.data.uni$Tube <- rownames(bac.ggplot.data.uni)
bac.ggplot.data.uni2 <- left_join(data.frame(bac.css.norm@sam_data), bac.ggplot.data.uni,  by = "Tube")

global.uni <- ggplot(bac.ggplot.data.uni2, aes(x = Axis.1, y = Axis.2, fill = Trt, shape = Type)) +
  geom_point(size=4, alpha = 0.7)+
  scale_shape_manual(values=c(21,22, 23, 24))+
  scale_fill_manual(values=cbbPalette, name = "Environment")+
  xlab(paste("PcoA1 -", variation.axis1, "%")) + 
  ylab(paste("PcoA2 -", variation.axis2, "%")) +
  guides(fill=guide_legend(override.aes=list(shape= 21)),
         shape=guide_legend(override.aes=list(fill= "black")))+
  theme_bw() + 
  ggtitle("Weighted Unifrac")

### Figure - pcoa
figure1 <- ggpubr::ggarrange(global.uni, global.bray, common.legend = TRUE, labels = c("a", "b"), legend = "right")









