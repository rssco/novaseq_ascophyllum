## Table of content  
[1. Library](#library)  
  
[2. Barplots bacteria - SSU](#bacteria_barplots)  
- [Legend](#legend)   
- [Sites](#bacteria_sites)
- [Parts](#bacteria_parts)
- [Seasons](#bacteria_seasons)
- [Combined figures](#combined_bacteria)
  
[3. Barplot Eukaryotes - SSU](#Eukaryotes)   
- [Legend](#legend_euk)   
- [Sites](#sites_euk)  
- [Parts](#parts_euk)  
- [Seasons](#season_euk)  
- [Combined figures](#combined_euk)
  
[4. Heatmap fungi - ITS](#heatmap_fungi)  
- [Prepare tables for vsearch](#table_vsearch)  
- [Generate fasta file](#generate_fasta_file)  
- [Clustering OTU 99%](#vsearch)
- [Final table](#final_table)
- [Ps object with OTU 99%](#ps_object_ITS)
- [Generate Heatmap](#heatmap)


## 1. Libraries <a name="library"></a>
```r
library(phyloseq)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(magrittr)
library(microbiome)
library(RColorBrewer)
library(microDecon)
library(paletteer)
library(gridExtra)
library(vegan)
library(rstatix)
library(ggpubr)
library(ape)
library(ggrepel)
library(reshape2)
library(ggh4x)
library(colorspace)
library(multcompView)
library(ggridges)
```


## 2. Barplots bacteria <a name="bacteria_barplots"></a>

### Legend <a name="legend"></a>

```r
genus <- c('Other bacteria'="grey",
                                            'Other Phycisphaerae'="#FFCEF4",'Algisphaera'="#FFADC1",
                                            'Other Verrucomicrobiae'="#E2E6BD",'DEV007;NA'="#E5D75D",'Roseibacillus'="#E8C33C",'Luteolibacter'="#EAAF29",
                                            'Other Planctomycetes'="#E1C1F0",'Blastopirellula'="#645A9F",
                                            'Other Cyanobacteriia'="#b6d67e",'Phormidesmis ANT.LACV5.1'="#77914a",'Synechococcus IR11'="#238B45",'Pleurocapsa PCC-7319'="#00441B",
                                            'Other Bacteroidia'="#C35472",'Pricia'="#CB6CA2",'Maribacter'="#A93154",'Zobellia'="#A34261",'Croceitalea'="#853250",'Dokdonia'="#67223F",
                                            'Other Gammaproteobacteria'="#FBA453",'Alteromonadaceae;NA'="#F27E00",'Paraglaciecola'="#DE5600",'Cellvibrionaceae;NA'="#F16913",'Hahella'="#FB6A4A",'Granulosicoccus'="#B03F00",                                       
                                            'Other Alphaproteobacteria'="#7d9ac9",'Rhodobacteraceae;NA'="#5993f0",'Litorimonas'="#08519C",'Octadecabacter'="#6BAED6",'Erythrobacter'="#4292C6",'Jannaschia'="#2171B5" ,'Sedimentitalea'="#9ECAE1",'Alphaproteobacteria;NA;NA;NA'="#08306B")
```

```r
level_genus <- c(c('Other bacteria',
                                            'Other Phycisphaerae','Algisphaera',
                                            'Other Verrucomicrobiae','DEV007;NA','Roseibacillus','Luteolibacter',
                                            'Other Planctomycetes','Blastopirellula',
                                            'Other Cyanobacteriia','Phormidesmis ANT.LACV5.1','Synechococcus IR11','Pleurocapsa PCC-7319',
                                            'Other Bacteroidia','Pricia','Maribacter','Zobellia','Croceitalea','Dokdonia',
                                            'Other Gammaproteobacteria','Alteromonadaceae;NA','Paraglaciecola','Cellvibrionaceae;NA','Hahella','Granulosicoccus',                                       
                                            'Other Alphaproteobacteria','Rhodobacteraceae;NA','Litorimonas','Octadecabacter','Erythrobacter','Jannaschia' ,'Sedimentitalea','Alphaproteobacteria;NA;NA;NA'))
```


### Site <a name="bacteria_sites"></a>
On normalized data 

```r
otu_tax_ref <- read_tsv("01_Tables/05_otu_tax_ref_asco_norm.tsv")
otu_tax_ref %<>% tibble::column_to_rownames("ASV")

otu <- otu_tax_ref %>% dplyr::select(`16S_L1I1_MARS_BACT_S81_R1`:`16S_L3I4_NOV_BACT_S124_R1`) %>% as.matrix()
tax <- otu_tax_ref %>% dplyr::select(Kingdom:Genus) %>% as.matrix()

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)

ps_sites_asco_norm <- phyloseq(OTU, TAX, ps_sites_asco@refseq, ps_sites_asco@sam_data)

#otu_table()   OTU Table:         [ 6995 taxa and 24 samples ]
#sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
#tax_table()   Taxonomy Table:    [ 6995 taxa by 6 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 6995 reference sequences ]

# remove no abundant asvs
total=81035 
ps_sites_asco_norm_trans = filter_taxa(ps_sites_asco_norm, function(x) sum(x > total* 0.000050) > 0, TRUE)
ps_sites_asco_norm_trans

#otu_table()   OTU Table:         [ 3846 taxa and 24 samples ]
#sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
#tax_table()   Taxonomy Table:    [ 3846 taxa by 6 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 3846 reference sequences ]
```

```r
ps_agglomerate <- tax_glom(ps_sites_asco_norm_trans, taxrank = 'Genus', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Genus <- as.character(taxa$Genus)
write.table(taxa, "01_Tables/08_taxa_agglomerate_all_barplot_asco_norm.csv", quote=FALSE, sep=";") #Add ;NA by hand

taxa <- read.table("~/Documents/ownCloud/11_Rstudio/08_Novaseq/04_Bacteria_analysis/03_Sites/01_Tables/08_taxa_agglomerate_all_barplot_asco_norm.csv", sep=";", header=T)
taxa %<>% filter(Kingdom=="Bacteria")
taxa %<>% filter(!ID=="L1I1_MARCH_BACT") 

taxa$Month <- factor(taxa$Month, levels = c("November", "March")) 

taxa$Genus[!taxa$Genus %in% level_genus] <- "Other bacteria"

class <- c("Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia","Cyanobacteriia", "Planctomycetes", "Verrucomicrobiae","Phycisphaerae", "Ascomycota")

taxa$Class[!taxa$Class %in% class] <- "Others"

taxa$Genus[taxa$Class=="Alphaproteobacteria" & !taxa$Genus %in% c('Alphaproteobacteria;NA;NA;NA','Sedimentitalea', 'Jannaschia','Erythrobacter', 'Octadecabacter','Litorimonas', 'Rhodobacteraceae;NA' )] <- "Other Alphaproteobacteria" #9

taxa$Genus[taxa$Class=="Gammaproteobacteria" & !taxa$Genus %in% c('Granulosicoccus','Hahella', 'Cellvibrionaceae;NA','Paraglaciecola','Alteromonadaceae;NA')] <- "Other Gammaproteobacteria" #6

taxa$Genus[taxa$Class=="Bacteroidia" & !taxa$Genus %in% c('Dokdonia','Croceitalea', 'Zobellia', 'Maribacter','Pricia')] <- "Other Bacteroidia" #7

taxa$Genus[taxa$Class=="Cyanobacteriia" & !taxa$Genus %in% c('Pleurocapsa PCC-7319','Synechococcus IR11', 'Phormidesmis ANT.LACV5.1')] <- "Other Cyanobacteriia" #3

taxa$Genus[taxa$Class=="Planctomycetes" & !taxa$Genus %in% c('Blastopirellula')] <- "Other Planctomycetes" #2

taxa$Genus[taxa$Class=="Verrucomicrobiae" & !taxa$Genus %in% c('Luteolibacter','Roseibacillus','DEV007;NA')] <- "Other Verrucomicrobiae" #3

taxa$Genus[taxa$Class=="Phycisphaerae" & !taxa$Genus %in% c('Algisphaera')] <- "Other Phycisphaerae" #3

taxa$Genus <- factor(taxa$Genus, levels = level_genus)
                                            

site <- ggplot(data=taxa, aes(x=Algae, y=Abundance, fill=Genus))  +
  facet_nested(~Month+Site, scales = "free_x", space = "free_x")+
  geom_bar(stat="identity", position = "fill") +
  scale_color_manual(values = genus ) + scale_fill_manual(values = genus) +
   theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.50, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
 guides(fill=guide_legend(ncol = 1)) +
   ylab("Relative abundance")



taxa_summary <- taxa %>% group_by(Site, Month, ID) %>% summarise(Sum=sum(Abundance))

## STAT BOX 
aov_shannon <- stats::aov(Sum ~ Site+Month, taxa_summary)
summary(aov_shannon)
#Site         2 4.653e+09 2.326e+09   0.446  0.647
#Month        1 1.186e+09 1.186e+09   0.227  0.639
#Residuals   19 9.912e+10 5.217e+09      

palette_month <- c("March" = "#1b4d08", "November" = "#963414")
palette_site <- c("#F8766D", "#00BF7D", "#00B0F6")

box_sites <- ggplot(data=taxa_summary, aes(x=Site, y=Sum, color=Site))  +
  facet_nested(~Month, scales = "free_x", space = "free_x")+
geom_boxplot(aes(fill=Site),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) +
scale_color_manual(values = palette_site) +
  scale_fill_manual(values = palette_site) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11), 
    legend.position = "none") +
   ylab("Number of reads after normalization with A. nodosum reads")  + 
  ggtitle("A. Sites")

a <- box_sites + site + plot_layout(nrow=2)
```


### Parts <a name="bacteria_parts"></a>
On normalized data
```r
otu_tax_ref <- read_tsv("01_Table/11_otu_tax_ref_1asv_asco_norm.tsv")
otu_tax_ref %<>% tibble::column_to_rownames("ASV")

otu <- otu_tax_ref %>% dplyr::select(`16S_L1A1_JANV_BACT_S133_R1`:`16S_L1R5_AVRIL_BACT_S260_R1`) %>% as.matrix()
tax <- otu_tax_ref %>% dplyr::select(Kingdom:Genus) %>% as.matrix()

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)

ps_parts_asco_norm <- phyloseq(OTU, TAX, ps_parts_asco@refseq, ps_parts_asco@sam_data)

#otu_table()   OTU Table:         [ 5947 taxa and 31 samples ]
#sample_data() Sample Data:       [ 31 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 5947 taxa by 6 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 5947 reference sequences ]

# remove no abundant asvs
total=125943 
ps_parts_asco_norm_trans = filter_taxa(ps_parts_asco_norm, function(x) sum(x > total* 0.000050) > 0, TRUE)
ps_parts_asco_norm_trans

```

```r
ps_agglomerate <- tax_glom(ps_parts_asco_norm_trans, taxrank = 'Genus', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Genus <- as.character(taxa$Genus)

write.table(taxa, "01_Table/16_barplot_agglomerate_all.csv", quote=FALSE, sep=";") #add ;NA by hand

taxa <- read.table("~/Documents/ownCloud/11_Rstudio/08_Novaseq/04_Bacteria_analysis/04_Parts/01_Table/16_barplot_agglomerate_all.csv", sep=";", header=TRUE)
taxa %<>% filter(Kingdom=="Bacteria")
taxa %<>% filter(!Sample=="16S_L1B4_AVRIL_BACT_S159_R1")
taxa %<>% filter(!Sample=="16S_L1R3_AVRIL_BACT_S156_R1") #outlier past4

taxa$Algae_part <- factor(taxa$Algae_part, levels = c("Receptacle", "Apex", "Medium", "Base")) 
taxa$Month <- factor(taxa$Month, levels = c("January", "April")) 

taxa$Genus[!taxa$Genus %in% level_genus] <- "Other bacteria"

class <- c("Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia","Cyanobacteriia", "Planctomycetes", "Verrucomicrobiae","Phycisphaerae", "Ascomycota")

taxa$Class[!taxa$Class %in% class] <- "Others"

taxa$Genus[taxa$Class=="Alphaproteobacteria" & !taxa$Genus %in% c('Alphaproteobacteria;NA;NA;NA','Sedimentitalea', 'Jannaschia','Erythrobacter', 'Octadecabacter','Litorimonas', 'Rhodobacteraceae;NA' )] <- "Other Alphaproteobacteria" #9

taxa$Genus[taxa$Class=="Gammaproteobacteria" & !taxa$Genus %in% c('Granulosicoccus','Hahella', 'Cellvibrionaceae;NA','Paraglaciecola','Alteromonadaceae;NA')] <- "Other Gammaproteobacteria" #6

taxa$Genus[taxa$Class=="Bacteroidia" & !taxa$Genus %in% c('Dokdonia','Croceitalea', 'Zobellia', 'Maribacter','Pricia')] <- "Other Bacteroidia" #7

taxa$Genus[taxa$Class=="Cyanobacteriia" & !taxa$Genus %in% c('Pleurocapsa PCC-7319','Synechococcus IR11', 'Phormidesmis ANT.LACV5.1')] <- "Other Cyanobacteriia" #3

taxa$Genus[taxa$Class=="Planctomycetes" & !taxa$Genus %in% c('Blastopirellula')] <- "Other Planctomycetes" #2

taxa$Genus[taxa$Class=="Verrucomicrobiae" & !taxa$Genus %in% c('Luteolibacter','Roseibacillus','DEV007;NA')] <- "Other Verrucomicrobiae" #3

taxa$Genus[taxa$Class=="Phycisphaerae" & !taxa$Genus %in% c('Algisphaera')] <- "Other Phycisphaerae" #3

taxa$Genus <- factor(taxa$Genus, levels = level_genus)
     

parts <- ggplot(data=taxa, aes(x=Individual, y=Abundance, fill=Genus))  +
  facet_nested( ~ Month+Algae_part, space = "free", scales = "free") +
  geom_bar(stat="identity", position="fill") +
  scale_color_manual(values = genus ) + scale_fill_manual(values = genus) +
 theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.50, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
 guides(fill=guide_legend(ncol = 1)) +
   ylab("Relative abundance") 


palette_month <- c("#081D58","#73243C")
palette_parts <- c("#225EA8", "#D67500","#578B21", "#B696D6")


taxa_summary <- taxa %>% group_by(Month, Sample, Algae_part, ID) %>% summarise(Sum=sum(Abundance))

## STAT BOX 
aov_shannon <- stats::aov(Sum ~ Month+Algae_part, taxa_summary)
summary(aov_shannon)
#Month        1 3.073e+10 3.073e+10   3.887 0.0598 .
#Algae_part   3 1.586e+10 5.287e+09   0.669 0.5792  
#Residuals   25 1.976e+11 7.906e+09        

box_parts <- ggplot(data=taxa_summary, aes(x=Algae_part, y=Sum, color=Algae_part))  +
  facet_nested(~Month)+
geom_boxplot(aes(fill=Algae_part),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) +
  scale_color_manual(values = palette_parts) +
  scale_fill_manual(values = palette_parts) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
   ylab("Number of reads after normalization with A. nodosum reads") +
  ggtitle("B. Algal parts")

b <- box_parts + parts + plot_layout(nrow=2)
```

### Season <a name="bacteria_seasons"></a>
On normalized data
```r
otu_tax_ref <- read_tsv("01_Tables/05_otu_tax_ref_asco_norm.tsv")
otu_tax_ref %<>% tibble::column_to_rownames("ASV")

otu <- otu_tax_ref %>% dplyr::select(`16S_L1I1_AOUT_BACT_S101_R1`:`16S_L1I5_JANV_BACT_S76_R1`) %>% as.matrix()
tax <- otu_tax_ref %>% dplyr::select(Kingdom:Genus) %>% as.matrix()

OTU <- otu_table(otu, taxa_are_rows = TRUE)
TAX <- tax_table(tax)

ps_season_asco_norm <- phyloseq(OTU, TAX, ps_season_asco@refseq, ps_season_asco@sam_data)

#otu_table()   OTU Table:         [ 11194 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 8 sample variables ]
#tax_table()   Taxonomy Table:    [ 11194 taxa by 6 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 11194 reference sequences ]

# remove little asvs
total=115701 
ps_season_asco_norm_trans = filter_taxa(ps_season_asco_norm, function(x) sum(x > total* 0.000050) > 0, TRUE)
ps_season_asco_norm_trans

#otu_table()   OTU Table:         [ 5723 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 8 sample variables ]
#tax_table()   Taxonomy Table:    [ 5723 taxa by 6 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 5723 reference sequences ]
```

```r
ps_agglomerate <- tax_glom(ps_season_asco_norm_trans, taxrank = 'Genus', NArm=FALSE)
ps_agglomerate

taxa <- psmelt(ps_agglomerate) 
taxa$Genus <- as.character(taxa$Genus)
write.table(taxa, "01_Tables/06_barplot_agglomerate_all.csv", quote=FALSE, sep=";") #Add ;NA by hand

taxa <- read.table("~/Documents/ownCloud/11_Rstudio/08_Novaseq/04_Bacteria_analysis/06_Seasons/01_Tables/06_barplot_agglomerate_all.csv", sep=";", header=TRUE)
taxa %<>% filter(Kingdom=="Bacteria")
taxa %<>% filter(!ID=="L1I1_MARCH") 

taxa$Month <- factor(taxa$Month, levels = c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November21", "November22", "December")) 
taxa$Season <- factor(taxa$Season, levels = c("Winter", "Spring", "Summer", "Fall")) 
taxa$ID <- factor(taxa$ID, levels =c(
"L1I1_NOV21",
"L1I2_NOV21",
"L1I3_NOV21",
"L1I4_NOV21",
"L1I1_DEC",
"L1I3_DEC",
"L1I4_DEC",
"L1I5_DEC",
"L1I1_JAN",
"L1I2_JAN",
"L1I3_JAN",
"L1I5_JAN",
"L1I1_FEB",
"L1I2_FEB",
"L1I3_FEB",
"L1I4_FEB",
"L1I1_MARCH",
"L1I2_MARCH",
"L1I3_MARCH",
"L1I4_MARCH",
"L1I1_APRIL",
"L1I2_APRIL",
"L1I3_AVRIL",
"L1I4_AVRIL",
"L1I1_MAY",
"L1I2_MAY",
"L1I3_MAY",
"L1I4_MAY",
"L1I1_JUNE",
"L1I2_JUNE",
"L1I3_JUNE",
"L1I4_JUNE",
"L1I1_JUILY",
"L1I2_JUILY",
"L1I3_JUILY",
"L1I4_JUILY",
"L1I1_AUGUST",
"L1I2_AUGUST",
"L1I3_AUGUST",
"L1I4_AUGUST",
"L1I1_SEPT",
"L1I2_SEPT",
"L1I3_SEPT",
"L1I4_SEPT",
"L1I1_OCT",
"L1I2_OCT",
"L1I3_OCT",
"L1I4_OCT",
"L1I1_NOV22",
"L1I2_NOV22",
"L1I3_NOV22",
"L1I4_NOV22"))

taxa$Genus[!taxa$Genus %in% level_genus] <- "Other bacteria"

class <- c("Alphaproteobacteria", "Gammaproteobacteria", "Bacteroidia","Cyanobacteriia", "Planctomycetes", "Verrucomicrobiae","Phycisphaerae", "Ascomycota")

taxa$Class[!taxa$Class %in% class] <- "Others"

taxa$Genus[taxa$Class=="Alphaproteobacteria" & !taxa$Genus %in% c('Alphaproteobacteria;NA;NA;NA','Sedimentitalea', 'Jannaschia','Erythrobacter', 'Octadecabacter','Litorimonas', 'Rhodobacteraceae;NA' )] <- "Other Alphaproteobacteria" #9

taxa$Genus[taxa$Class=="Gammaproteobacteria" & !taxa$Genus %in% c('Granulosicoccus','Hahella', 'Cellvibrionaceae;NA','Paraglaciecola','Alteromonadaceae;NA')] <- "Other Gammaproteobacteria" #6

taxa$Genus[taxa$Class=="Bacteroidia" & !taxa$Genus %in% c('Dokdonia','Croceitalea', 'Zobellia', 'Maribacter','Pricia')] <- "Other Bacteroidia" #7

taxa$Genus[taxa$Class=="Cyanobacteriia" & !taxa$Genus %in% c('Pleurocapsa PCC-7319','Synechococcus IR11', 'Phormidesmis ANT.LACV5.1')] <- "Other Cyanobacteriia" #3

taxa$Genus[taxa$Class=="Planctomycetes" & !taxa$Genus %in% c('Blastopirellula')] <- "Other Planctomycetes" #2

taxa$Genus[taxa$Class=="Verrucomicrobiae" & !taxa$Genus %in% c('Luteolibacter','Roseibacillus','DEV007;NA')] <- "Other Verrucomicrobiae" #3

taxa$Genus[taxa$Class=="Phycisphaerae" & !taxa$Genus %in% c('Algisphaera')] <- "Other Phycisphaerae" #3

taxa$Genus <- factor(taxa$Genus, levels = level_genus)

season <- ggplot(data=taxa, aes(x=ID, y=Abundance, fill=Genus))  +
  geom_bar(stat="identity", position = "fill")  +
  scale_color_manual(values = genus ) + scale_fill_manual(values = genus) +
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
 guides(fill=guide_legend(ncol = 1)) +
   ylab("Relative abundance") 


palette_month <- c("January"="#3039b8", 
                   "February"="#9ECAE1",
                   "March" = "#1b4d08",
                   "April" ="#1a8a45",
                   "May"="#529139", 
                   "June"= "#f2b613", 
                   "July"= "#f56d05", 
                   "August"="#f20a2c",
                   "September"="#e0c32d",
                   "October"="#a17705",
                   "November22"="#963414",
                   "November21" = "#963414",
                   "December"="#733c05")


palette_season <- c("Winter"="#6BAED6",
                    "Spring"="#238B45",
                    "Summer"="#f2b613",
                    "Fall"="#963414")

taxa_summary <- taxa %>% group_by(Month, Sample) %>% summarise(Sum=sum(Abundance))

## STAT BOX 
aov_shannon <- stats::aov(Sum ~ Month, taxa_summary)
summary(aov_shannon)
#Month       12 1.438e+11 1.198e+10   3.817 0.000768 ***
#Residuals   38 1.193e+11 3.140e+09   

tukey <- TukeyHSD(aov_shannon)
tukey

cld <- multcompLetters4(aov_shannon, tukey)
cld


 TukeyHSD(aov_shannon) 
# diff         lwr           upr     p adj
#May-March             -204051.58 -355138.194 -52964.9727 0.0016742
#July-March            -173422.33 -324508.944 -22335.7227 0.0128970
#September-May          154962.25   15083.229 294841.2709 0.0188677
#May-December          -147164.75 -287043.771  -7285.7291 0.0315437
#May-December          -147164.75 -287043.771  -7285.7291 0.0315437

taxa_summary$Month <- factor(taxa_summary$Month, levels = c("November21", "December", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November22")) 

box_season <- ggplot(data=taxa_summary, aes(x=Month, y=Sum, color=Month))  +
geom_boxplot(aes(fill=Month),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) +
  scale_color_manual(values = palette_month) +
  scale_fill_manual(values = palette_month) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
   ylab("Number of reads after normalization with A. nodosum reads") +
  ggtitle("C. Seasons")

```

### Combined bacteria <a name="combined_bacteria"></a>
```{r}
library(patchwork)
#site + parts + plot_layout(guides = 'collect')

(box_sites / site) + theme(legend.position = "none") + plot_layout(heights = c(1,3)) | (box_parts / parts)+ theme(legend.position = "none") + plot_layout(heights = c(1,3)) | (box_season  / season) +  plot_layout(guides = 'collect', heights = c(1,3)) #+ pdf("02_Figures/19_combined_barplots_bacteria.pdf", height = 13, width = 25)
```

![Figure 2 | Barplot bacteria.](https://github.com/rssco/novaseq_ascophyllum/blob/main/01_Figures/Figure_2_Combined_barplots_bacteria.png)<!-- -->


## 2. Eukaryotes <a name="Eukaryota"></a>
### Legend <a name="legend_euk"></a>

```r
genus <- c("Other Eukaryota"="grey",
                                              "Eukaryota;NA;NA;NA;NA;NA"="maroon",
                                             "Other Arthropoda"="#FEEDDE","Hyale"="#E6B197","Pseudocalanus"="#FBA453",'Maxillopoda;NA'="#E6550D","Diarthrodes"="#A63603",
                                             "Other Phaeophyceae"="#9DBFEA","Stramenopiles;NA;NA;NA;NA"= "#1D91C0","Phaeophyceae_XX;NA"="#225EA8",
                                             "Other Ascomycota"="#BAE4B3", "Sordariomycetes;NA"="#74C476","Dothideomycetes;NA"="#238B45")

level_genus <- c("Other Eukaryota",
                                              "Eukaryota;NA;NA;NA;NA;NA",
                                             "Other Arthropoda","Hyale",'Pseudocalanus','Maxillopoda;NA',"Diarthrodes",
                                             "Other Phaeophyceae","Stramenopiles;NA;NA;NA;NA","Phaeophyceae_XX;NA",
                                             "Other Ascomycota", "Sordariomycetes;NA","Dothideomycetes;NA")
```

### Sites <a name="sites_euk"></a>
```r
taxa <- read.table("~/Documents/ownCloud/11_Rstudio/08_Novaseq/04_Bacteria_analysis/03_Sites/01_Tables/08_taxa_agglomerate_all_barplot_asco_norm.csv", sep=";", header=T)
taxa %<>% filter(Kingdom=="Eukaryota")
taxa %<>% filter(!Genus=="Ascophyllum")
taxa %<>% filter(!ID=="L1I1_MARCH_BACT") 

class <- c("Arthropoda","Eukaryota;NA;NA","Phaeophyceae","Ascomycota")

taxa$Class[!taxa$Class %in% class] <- "Others"

taxa$Genus[taxa$Class=="Ascomycota" & !taxa$Genus %in% c("Sordariomycetes;NA","Dothideomycetes;NA")] <- "Other Ascomycota" #9

taxa$Genus[taxa$Class=="Phaeophyceae" & !taxa$Genus %in% c("Saccharina","Phaeophyceae_XX;NA")] <- "Other Phaeophyceae" #6

taxa$Genus[taxa$Class=="Eukaryota;NA;NA" & !taxa$Genus %in% c("Eukaryota;NA;NA;NA;NA;NA")] <- "Other Eukaryota;NA;NA" #7

taxa$Genus[taxa$Class=="Arthropoda" & !taxa$Genus %in% c("Hyale")] <- "Other Arthropoda" #3

taxa$Genus[!taxa$Genus %in% level_genus] <- "Other Eukaryota"


taxa$Genus <- factor(taxa$Genus, levels = level_genus)
taxa$Month <- factor(taxa$Month, levels = c("November", "March"))

taxa$Sample <- factor(taxa$Sample, levels = c(
  "16S_L1I1_NOV21_BACT_S65_R1", 
  "16S_L1I2_NOV21_BACT_S66_R1",
  "16S_L1I3_NOV21_BACT_S67_R1",
  "16S_L1I4_NOV21_BACT_S68_R1",
  "16S_L2I1_NOV_BACT_S117_R1",
  "16S_L2I2_NOV_BACT_S118_R1", 
  "16S_L2I3_NOV_BACT_S119_R1",
  "16S_L2I4_NOV_BACT_S120_R1",
  "16S_L3I1_NOV_BACT_S121_R1",
  "16S_L3I2_NOV_BACT_S122_R1",
  "16S_L3I3_NOV_BACT_S123_R1",
  "16S_L3I4_NOV_BACT_S124_R1",
  "16S_L1I2_MARS_BACT_S82_R1",
  "16S_L1I3_MARS_BACT_S83_R1", 
  "16S_L1I4_MARS_BACT_S84_R1", 
  "16S_L2I1_MARS_BACT_S125_R1",
  "16S_L2I2_MARS_BACT_S126_R1",
  "16S_L2I3_MARS_BACT_S127_R1",
  "16S_L2I4_MARS_BACT_S128_R1",
  "16S_L3I1_MARS_BACT_S129_R1",
  "16S_L3I2_MARS_BACT_S130_R1", 
  "16S_L3I3_MARS_BACT_S131_R1",
  "16S_L3I4_MARS_BACT_S132_R1"))

sample_names <- c(
"16S_L1I1_NOV21_BACT_S65_R1"="L1I1_NOV",
"16S_L1I2_MARS_BACT_S82_R1"="L1I2_MARCH",
"16S_L1I2_NOV21_BACT_S66_R1"="L1I2_NOV",
"16S_L1I3_MARS_BACT_S83_R1"="L1I3_MARCH",
"16S_L1I3_NOV21_BACT_S67_R1"="L1I3_NOV",
"16S_L1I4_MARS_BACT_S84_R1"="L1I4_MARCH",
"16S_L1I4_NOV21_BACT_S68_R1"="L1I4_NOV",
"16S_L2I1_MARS_BACT_S125_R1"="L2I1_MARCH",
"16S_L2I1_NOV_BACT_S117_R1"="L2I1_NOV",
"16S_L2I2_MARS_BACT_S126_R1"="L2I2_MARCH",
"16S_L2I2_NOV_BACT_S118_R1"="L2I2_NOV",
"16S_L2I3_MARS_BACT_S127_R1"="L2I3_MARCH",
"16S_L2I3_NOV_BACT_S119_R1"="L2I3_NOV",
"16S_L2I4_MARS_BACT_S128_R1"="L2I4_MARCH",
"16S_L2I4_NOV_BACT_S120_R1"="L2I4_NOV",
"16S_L3I1_MARS_BACT_S129_R1"="L3I1_MARCH",
"16S_L3I1_NOV_BACT_S121_R1"="L3I1_NOV",
"16S_L3I2_MARS_BACT_S130_R1"="L3I2_MARCH",
"16S_L3I2_NOV_BACT_S122_R1"="L3I2_NOV",
"16S_L3I3_MARS_BACT_S131_R1"="L3I3_MARCH",
"16S_L3I3_NOV_BACT_S123_R1"="L3I3_NOV",
"16S_L3I4_MARS_BACT_S132_R1"="L3I4_MARCH",
"16S_L3I4_NOV_BACT_S124_R1"="L3I4_NOV")

sites <- ggplot(data=taxa, aes(x=Algae, y=Abundance, fill=Genus))  +
  facet_nested(~Month+Site, scales = "free_x", space = "free_x")+
  geom_bar(stat="identity", position = "fill") +
  scale_x_discrete(labels=sample_names)+ 
  scale_color_manual(values = genus ) + scale_fill_manual(values = genus) +
 theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.50, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
 guides(fill=guide_legend(ncol = 1)) +
   ylab("Relative abundance") 


palette_month <- c("March" = "#1b4d08", "November" = "#963414")
palette_site <- c("#F8766D", "#00BF7D", "#00B0F6")

taxa_summary <- taxa %>% group_by(Site, Month, ID) %>% summarise(Sum=sum(Abundance))

## STAT BOX 
aov_shannon <- stats::aov(Sum ~ Site+Month, taxa_summary)
summary(aov_shannon)
#Site         2  3557246 1778623   0.621  0.548
#Month        1   184963  184963   0.065  0.802
#Residuals   19 54461485 2866394     

box_sites <- ggplot(data=taxa_summary, aes(x=Site, y=Sum, color=Site))  +
  facet_nested(~Month, scales = "free_x", space = "free_x")+
geom_boxplot(aes(fill=Site),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) +
scale_color_manual(values = palette_site) +
  scale_fill_manual(values = palette_site) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11), 
    legend.position = "none") +
   ylab("Number of reads after normalization with A. nodosum reads")  + 
  ggtitle("A. Sites")
```

### Parts <a name="parts_euk"></a>

```r
taxa <- read.table("~/Documents/ownCloud/11_Rstudio/08_Novaseq/04_Bacteria_analysis/04_Parts/01_Table/16_barplot_agglomerate_all.csv", sep=";", header=T)
taxa %<>% filter(Kingdom=="Eukaryota")
taxa %<>% filter(!Genus=="Ascophyllum")
taxa %<>% filter(!Sample=="16S_L1B4_AVRIL_BACT_S159_R1")
taxa %<>% filter(!Sample=="16S_L1R3_AVRIL_BACT_S156_R1") #outlier past4

class <- c("Arthropoda", "Ascomycota","Eukaryota;NA;NA","Phaeophyceae", "Stramenopiles;NA", "Nematoda", "Platyhelminthes")

taxa$Class[!taxa$Class %in% class] <- "Others"


taxa$Genus[taxa$Class=="Arthropoda" & !taxa$Genus %in% c('Tisbe', 'Euterpina','Maxillopoda;NA', 'Pseudocalanus','Diarthrodes')] <- "Other Arthropoda" 

taxa$Genus[taxa$Class=="Phaeophyceae" & !taxa$Genus %in% c('Phaeophyceae_XX;NA')] <- "Other Phaeophyceae"


taxa$Genus[taxa$Class=="Ascomycota" & !taxa$Genus %in% c('Sordariomycetes;NA','Dothideomycetes;NA')] <- "Other Ascomycota" 

taxa$Genus[taxa$Class=="Stramenopiles;NA" & !taxa$Genus %in% c('Stramenopiles;NA;NA;NA;NA')] <- "Other Stramenopiles;NA" 

taxa$Genus[taxa$Class=="Nematoda" & !taxa$Genus %in% c('Pellioditis','Halomonhystera')] <- "Other Nematoda" 

taxa$Genus[taxa$Class=="Platyhelminthes" & !taxa$Genus %in% c('Trigonostomum')] <- "Other Platyhelminthes" 

taxa$Genus[!taxa$Genus %in% level_genus] <- "Other Eukaryota"


taxa$Algae_part <- factor(taxa$Algae_part, levels = c("Receptacle", "Apex", "Medium", "Base")) 
taxa$Month <- factor(taxa$Month, levels = c("January", "April")) 
taxa$ID <- factor(taxa$ID, levels =c(
  "L1R1_JANV_BACT","L1R2_JANV_BACT","L1R3_JANV_BACT","L1R4_JANV_BACT",
  "L1A1_JANV_BACT","L1A2_JANV_BACT","L1A3_JANV_BACT","L1A4_JANV_BACT", 
  "L1M1_JANV_BACT","L1M2_JANV_BACT","L1M3_JANV_BACT","L1M4_JANV_BACT",
  "L1B1_JANV_BACT","L1B2_JANV_BACT","L1B3_JANV_BACT","L1B4_JANV_BACT",
  "L1R2_AVRIL_BACT","L1R4_AVRIL_BACT","L1R5_AVRIL_BACT",
  "L1A2_AVRIL_BACT","L1A3_AVRIL_BACT","L1A4_AVRIL_BACT","L1A5_AVRIL_BACT", 
  "L1M2_AVRIL_BACT","L1M3_AVRIL_BACT","L1M4_AVRIL_BACT","L1M5_AVRIL_BACT",
  "L1B2_AVRIL_BACT","L1B3_AVRIL_BACT","L1B5_AVRIL_BACT"))

taxa$Genus <- factor(taxa$Genus, levels = level_genus)


parts <- ggplot(data=taxa, aes(x=ID, y=Abundance, fill=Genus))  +
  facet_nested(~Month, scales = "free_x", space = "free_x")+
  geom_bar(stat="identity", position = "fill") +
      scale_color_manual(values = genus ) + scale_fill_manual(values = genus) +
 theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.50, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
 guides(fill=guide_legend(ncol = 1)) +
   ylab("Relative abundance") 



palette_month <- c("#081D58","#73243C")
palette_parts <- c("#225EA8", "#D67500","#578B21", "#B696D6")

## STAT BOX 
taxa_summary <- taxa %>% group_by(Month, Sample, Algae_part, ID) %>% summarise(Sum=sum(Abundance))
taxa_summary$Algae_part <- factor(taxa_summary$Algae_part, levels = c("Receptacle", "Apex", "Medium", "Base")) 

aov_shannon <- stats::aov(Sum ~ Algae_part+Month, taxa_summary)
summary(aov_shannon)
#Algae_part   3 16239279 5413093   1.430  0.258
#Month        1   370467  370467   0.098  0.757
#Residuals   25 94620489 3784820             

box_parts <- ggplot(data=taxa_summary, aes(x=Algae_part, y=Sum, color=Algae_part))  +
  facet_nested(~Month)+
geom_boxplot(aes(fill=Algae_part),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) +
  scale_color_manual(values = palette_parts) +
  scale_fill_manual(values = palette_parts) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11), 
    legend.position = "none") +
   ylab("Number of reads after normalization with A. nodosum reads")  + 
  ggtitle("B. Algal parts")
```

### Season <a name="season_euk"></a>

```r
taxa <- read.table("~/Documents/ownCloud/11_Rstudio/08_Novaseq/04_Bacteria_analysis/06_Seasons/01_Tables/06_barplot_agglomerate_all.csv", sep=";", header=TRUE)
taxa %<>% filter(Kingdom=="Eukaryota")
taxa %<>% filter(!Genus=="Ascophyllum") 
taxa %<>% filter(!ID=="L1I1_MARCH") 

class <- c("Ascomycota","Phaeophyceae","Arthropoda")

taxa$Class[!taxa$Class %in% class] <- "Others"


taxa$Genus[taxa$Class=="Arthropoda" & !taxa$Genus %in% c('Tisbe','Hyale')] <- "Other Arthropoda" 

taxa$Genus[taxa$Class=="Phaeophyceae" & !taxa$Genus %in% c('Silvetia','Phaeophyceae_XX;NA')] <- "Other Phaeophyceae"

taxa$Genus[taxa$Class=="Ascomycota" & !taxa$Genus %in% c('Pezizomycotina;NA;NA','Sordariomycetes;NA','Dothideomycetes;NA')] <- "Other Ascomycota" 


Genus_others <- c(
                  "Other Arthropoda",'Tisbe','Hyale',
                  "Other Phaeophyceae",'Phaeophyceae_XX;NA',
                  "Other Ascomycota",'Sordariomycetes;NA','Dothideomycetes;NA') 

taxa$Genus[!taxa$Genus %in% level_genus] <- "Other Eukaryota"


taxa$Month <- factor(taxa$Month, levels = c("November21","December", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November22")) 
taxa$ID <- factor(taxa$ID, levels =c(
"L1I1_NOV21",
"L1I2_NOV21",
"L1I3_NOV21",
"L1I4_NOV21",
"L1I1_DEC",
"L1I3_DEC",
"L1I4_DEC",
"L1I5_DEC",
"L1I1_JAN",
"L1I2_JAN",
"L1I3_JAN",
"L1I5_JAN",
"L1I1_FEB",
"L1I2_FEB",
"L1I3_FEB",
"L1I4_FEB",
"L1I2_MARCH",
"L1I3_MARCH",
"L1I4_MARCH",
"L1I1_APRIL",
"L1I2_APRIL",
"L1I3_AVRIL",
"L1I4_AVRIL",
"L1I1_MAY",
"L1I2_MAY",
"L1I3_MAY",
"L1I4_MAY",
"L1I1_JUNE",
"L1I2_JUNE",
"L1I3_JUNE",
"L1I4_JUNE",
"L1I1_JUILY",
"L1I2_JUILY",
"L1I3_JUILY",
"L1I4_JUILY",
"L1I1_AUGUST",
"L1I2_AUGUST",
"L1I3_AUGUST",
"L1I4_AUGUST",
"L1I1_SEPT",
"L1I2_SEPT",
"L1I3_SEPT",
"L1I4_SEPT",
"L1I1_OCT",
"L1I2_OCT",
"L1I3_OCT",
"L1I4_OCT",
"L1I1_NOV22",
"L1I2_NOV22",
"L1I3_NOV22",
"L1I4_NOV22"))

taxa$Genus <- factor(taxa$Genus, levels =level_genus) 

seasons <- ggplot(data=taxa, aes(x=ID, y=Abundance, fill=Genus))  +
  geom_bar(stat="identity", position = "fill")  +
  scale_color_manual(values = genus ) + scale_fill_manual(values = genus) +
 theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.50, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11)) +
 guides(fill=guide_legend(ncol = 1)) +
   ylab("Relative abundance") 

palette_month <- c("January"="#3039b8", 
                   "February"="#9ECAE1",
                   "March" = "#1b4d08",
                   "April" ="#1a8a45",
                   "May"="#529139", 
                   "June"= "#f2b613", 
                   "July"= "#f56d05", 
                   "August"="#f20a2c",
                   "September"="#e0c32d",
                   "October"="#a17705",
                   "November22"="#963414",
                   "November21" = "#963414",
                   "December"="#733c05")

taxa_summary <- taxa %>% group_by(Month, Sample, ID) %>% summarise(Sum=sum(Abundance))
taxa_summary$Month <- factor(taxa_summary$Month, levels = c("November21","December", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November22")) 

aov_shannon <- stats::aov(Sum ~ Month, taxa_summary)
summary(aov_shannon)
#Month       12 133540910 11128409   1.028  0.444
#Residuals   38 411373173 10825610               

box_season <- ggplot(data=taxa_summary, aes(x=Month, y=Sum, color=Month))  +
geom_boxplot(aes(fill=Month),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) +
  scale_color_manual(values = palette_month) +
  scale_fill_manual(values = palette_month) +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
         strip.text=element_text(size=12),
        legend.title = element_text(size=12),
        legend.key.size = unit(0.65, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title=element_text(size=10),
    legend.text = element_text(size=11), 
    legend.position = "none") +
   ylab("Number of reads after normalization with A. nodosum reads")  + 
  ggtitle("C. Seasons")

```

### Combined Eukaryotes <a name="combined_euk"></a>
```r
library(patchwork)
#site + parts + plot_layout(guides = 'collect')

(box_sites / sites) + theme(legend.position = "none") + plot_layout(heights = c(1,3)) | (box_parts / parts)+ plot_layout(heights = c(1,3), guides = 'collect') | (box_season  / seasons) +  plot_layout(heights = c(1,3))+ theme(legend.position = "none")  + pdf("02_Figures/20_combined_barplots_eukaryotes2.pdf", height = 13, width = 25)
```

## 4. Heatmap fungi - ITS <a name="heatmap_fungi"></a>
### Prepare tables for vsearch <a name="table_vsearch"></a>
```r
#!/bin/Rscript
ps_fungi <- readRDS("../01_All/02_Phyloseq_objects/04_ps_ss_decon_wt_cc_mock.rds")

ref <- ps_fungi@refseq %>% as.data.frame()
colnames(ref)[1] <- "sequence"
ref <- cbind(rownames(ref), ref)
rownames(ref) <- NULL
colnames(ref)[1] <- "ASV"

tax <- ps_fungi@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_fungi@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

ref_tax <- merge(tax, ref, by="ASV")
otu_tax_ref <- merge(otu, ref_tax, by="ASV")
write.table(otu_tax_ref, "03_Tables/09_otu_tax_ref_for_vsearch.tsv", sep="\t", quote=FALSE)
```

### Generate fasta file with ASV_Species+sequence <a name="generate_fasta_file"></a>
```bash
#!/bin/bash
cat 09_otu_tax_ref_for_vsearch.tsv| sort -n -k 2 | cut -f1,130 | tr "\t" "\n" | sed "s/^ASV/>ASV/g" > fungi.fasta
```

### Clustering OTU 99% <a name="vsearch"></a>
```bash
#!/usr/bin/env bash
#SBATCH --job-name=vsearch
#SBATCH -p fast
#SBATCH --mem=1G
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL

module load vsearch/2.21.1


## NOVASEQ CLUSTERING 0.99%
input=../vsearch_ns/fungi.fasta
output=../fungi_consensus99.fasta
name=cluster

vsearch --cluster_smallmem ${input} --usersort --id 0.99 --sizeout --sizein --consout ${output} --clusters /scratch2/seabioz/tmp/${name}_
```

```bash 
#!/bin/bash
for i in {0..905}; do grep ">" cluster_${i} | sed "s/;size=1//g" | sed "s/>//g" | sed -e 's/,$//g' -e 's/ //g' > vector_${i}.txt; done 

for i in *.txt; do awk '{print $0, FILENAME}' $i |sed 's/\.txt//g' | awk 'NR==1 {line=$0} NR>=1 {$3=line" "$3} 1'  |sed 's/ /\t/g' | cut -f1,2,3 >  ${i%.txt}_v.txt; done

cat *v.txt > finale.txt
```


## Generate final table <a name="final_table"></a>
```r
#!/bin/Rscript
#merge 2 tables by hand because doesn't work well with R. 
#otu_tax_ref <- read.table("03_Tables/09_otu_tax_ref_for_vsearch.tsv", header=TRUE, sep="\t")
#vector <- read.table("03_Tables/10_asv_vector_representant_finale.txt", header=TRUE)

otu_tax_ref_vector <- read.table("03_Tables/11_otu_tax_ref_for_vsearch_vector_representant_vsearch99.tsv", header=TRUE, sep="\t")
tax_seq <- otu_tax_ref_vector %>% dplyr::select(Kingdom:Species, sequence, ASV)

otu_tax_ref_vector_agg <- rowsum(otu_tax_ref_vector[,c(4:123)],otu_tax_ref_vector$Representant,na.rm=T)
otu_tax_ref_vector_agg$ASV <- rownames(otu_tax_ref_vector_agg)


otu_tax_ref_agg_fungi <- merge(otu_tax_ref_vector_agg, tax_seq, by="ASV")
write.table(otu_tax_ref_agg_fungi,"03_Tables/12_otu_tax_ref_agg_fungi.tsv", sep="\t", quote=FALSE)

#remove Cecile samples + all ASV with sum = 0 
```

## ps object with OTU 99% <a name="ps_object_ITS"></a>
```r
#!/bin/Rscript
ps_old <- readRDS("../01_All/02_Phyloseq_objects/04_ps_ss_decon_wt_cc_mock.rds")
otu_tax_ref_99 <- read.table("03_Tables/14_otu_tax_ref_agg_fungi_wt_cecile_receptacles_vsearch99_blast_complete.csv", sep=";", header=TRUE, row.names = 1)

otu <- otu_tax_ref_99 %>% select(L1A1_JANV_FUNGI_S229_R1:L3I4_NOV_FUNGI_S220_R1) %>% as.matrix()
tax <- otu_tax_ref_99 %>% select(Kingdom:Species) %>% as.matrix()
sample <- read.table("03_Tables/15_samples_all.csv", sep=";", header=TRUE, row.names = 1)

ps_all_99 <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(sample), ps_old@refseq)
#otu_table()   OTU Table:         [ 381 taxa and 100 samples ]
#sample_data() Sample Data:       [ 100 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 381 taxa by 8 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 381 reference sequences ]
```

### Filtered and mediane normalization <a name="normalization"></a>
```r
#!/bin/Rscript
total = median(sample_sums(ps_all_99)) # Median= 294093
standf = function(x, t=total) round(t * (x / sum(x))) # Standardize abundances to the median sequencing depth
ps_all_filt = transform_sample_counts(ps_all_99, standf)

ps_all_trans_99 = filter_taxa(ps_all_filt, function(x) sum(x > total* 0.00010) > 0, TRUE)
ps_all_trans_99

#otu_table()   OTU Table:         [ 276 taxa and 100 samples ]
#sample_data() Sample Data:       [ 100 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 276 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 276 reference sequences ]
```

### SAVE
```r
#!/bin/Rscript
saveRDS(ps_all_99, "02_Phyloseq_objects/05_ps_all_vsearch99.rds")
saveRDS(ps_all_trans_99, "02_Phyloseq_objects/06_ps_all_trans_99.rds")

ref <- ps_all_99@refseq %>% as.data.frame()
colnames(ref)[1] <- "sequence"
ref <- cbind(rownames(ref), ref)
rownames(ref) <- NULL
colnames(ref)[1] <- "ASV"

tax <- ps_all_99@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_all_99@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

ref_tax <- merge(tax, ref, by="ASV")
otu_tax_ref <- merge(otu, ref_tax, by="ASV")
write.table(otu_tax_ref,"03_Tables/16_otu_tax_ref_vsearch99.tsv", sep="\t", quote=FALSE)

ref <- ps_all_trans_99@refseq %>% as.data.frame()
colnames(ref)[1] <- "sequence"
ref <- cbind(rownames(ref), ref)
rownames(ref) <- NULL
colnames(ref)[1] <- "ASV"

tax <- ps_all_trans_99@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_all_trans_99@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

ref_tax <- merge(tax, ref, by="ASV")
otu_tax_ref <- merge(otu, ref_tax, by="ASV")
write.table(otu_tax_ref,"03_Tables/17_otu_tax_ref_trans_vsearch99.tsv", sep="\t", quote=FALSE)
```

### Generate heatmap <a name="heatmap"></a>
```r
#!/bin/Rscript
ps_all_99 <- readRDS("02_Phyloseq_objects/06_ps_all_trans_99.rds")

tree <- read.tree("04_Figures/01_tree_asv_myco_mohei_lul_other.tree")

tree <- ggtree(tree) + 
        geom_tiplab(size=2, align = TRUE) +
        geom_treescale(color="black", fontsize = 0.5) 
tree


label_tree <- tree$data$label %>% as.data.frame()
#write.table(label_tree, "03_Tables/18_label_tree_asv.csv", quote=FALSE, sep=";")
asv <- read.table("03_Tables/18_label_tree_asv.csv", sep=";", header=TRUE)

#from heatmap part
otu_tax_ref <- read.table("03_Tables/19_label_myco_mohei_lul.tsv", sep="\t", header=TRUE, row.names = 1)

otu <- otu_tax_ref %>% dplyr::select(L1A1_JANV_FUNGI_S229_R1:L3I4_NOV_FUNGI_S220_R1) %>% as.matrix()
tax <- otu_tax_ref %>% dplyr::select(Kingdom:Species) %>% as.matrix()

ps_myco <- phyloseq(refseq(ps_all_99@refseq), otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(ps_all_99@sam_data))

tax <- ps_myco@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

sam <- ps_myco@sam_data %>% as.data.frame()
sam <- cbind(rownames(sam), sam)
rownames(sam) <- NULL
colnames(sam)[1] <- "variable"

otu <- ps_myco@otu_table %>% as.data.frame(stringsAsFactors = FALSE)
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"
otu %<>% melt(id.vars = "ASV")

merged <- merge(otu, tax, by="ASV")
merged <- merge(merged,sam, by="variable")
#write.table(merged, "03_Tables/20_merged_heatmap_tree_myco_mohei_lul_other.tsv", quote=FALSE, sep="\t")
merged <- read.table("03_Tables/20_merged_heatmap_tree_myco_mohei_lul_other.tsv", sep="\t", header=TRUE)

merged %<>% group_by(Month, Individual, ASV, Genus, Season, Dataset, Type) %>% summarise(Sum=sum(value))

label_tree_heatmap <- merge(asv, merged, by="ASV")
label_tree_heatmap %<>% dplyr::select("label_tree", "Sum", "Month", "Season", "Genus", "Dataset", "Individual", "Type")

label_tree_heatmap$Month <- factor(label_tree_heatmap$Month, levels = c("November21","December","January", "February", "March", "April", "May", "June", "July", "August", "September", "October",  "November22")) 
label_tree_heatmap$Season <- factor(label_tree_heatmap$Season, levels = c("Fall","Winter","Spring", "Summer")) 
label_tree_heatmap$Dataset <- factor(label_tree_heatmap$Dataset, levels = c("Site","Parts","Season")) 
label_tree_heatmap$Type <- factor(label_tree_heatmap$Type, levels = c("Site 1", "Site 2", "Site 3","Receptacle","Apex","Medium","Base", "Pool")) 


colors1 <- brewer.pal(11, "RdYlBu")

heatmap <- ggplot(label_tree_heatmap, aes(x=Individual,y=label_tree,fill=Sum)) + geom_tile() + 
  facet_nested(.~Dataset+Month+Type,space = "free", scales = "free")+
  scale_fill_gradientn(colours =rev(colors1),name="Abundance", trans='sqrt') +
  theme(axis.text.x = element_blank(), 
        strip.text=element_text(size=15),
        axis.text.y = element_text(size=8), 
        panel.background = element_blank(), 
        legend.title = element_text(size=14, face = "bold"),
        strip.text.y = element_blank(),
        panel.spacing.x  = unit(0.1, "lines"),
        panel.spacing.y  = unit(0.05, "lines")) +
  xlab(NULL) + ylab(NULL)


g <- heatmap %>% insert_left(tree, width = 0.2)
ggsave("04_Figures/02_tree_heatmap_myco_mohei_lul_other.pdf", g, height = 10, width = 25)
```

![Figure 3 | Barplot 18S.](https://github.com/rssco/novaseq_ascophyllum/blob/main/01_Figures/Figure_3_Combined_barplots_eukaryotes.png)<!-- -->