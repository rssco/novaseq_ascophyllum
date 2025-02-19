## Table of content 
[1. Library](#library)  

[2. Venn Diagrams](#venn)  
- [Ps object](#ps_object)  
- [Venn parts](#venn_parts)
- [Venn sites and seasons](#venn_sites_seasons)
- [Ps core](#ps_core)
- [Venn core](#venn_core)
- [Combined Venn](#combined_venn)

[3. Specifities](#specificities)   
- [Ps specificity parts](#ps_specificity_parts)   
- [Ps specificity season](#ps_specificity_seasons)  
- [Bubble plot](#bubble_plot)  


## 1. Library  <a name="library"></a>
```r
#!/bin/Rscript
library(tidyverse)
library(magrittr)
library(ggplot2)
library(aplot)
library(phyloseq)
library(microbiome)
library(gridExtra)
library(MicEco)
library(multcompView)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(ggh4x)
```


## 2. Venn Diagrams <a name="venn"></a>
### ps object <a name="ps_object"></a>
```r
#!/bin/Rscript
s_season_asco_norm <- readRDS("../06_Seasons/00_PHYLOSEQ_OBJECTS/03_ps_season_asco_norm.rds")
ps_season_asco_norm <- subset_samples(ps_season_asco_norm, ID !="L1I1_MARCH")
tax <- ps_season_asco_norm@tax_table %>% as.data.frame()
asco <- tax %>% as.data.frame() %>% filter(Genus=="Ascophyllum") #1
otu_asco <- rownames(asco[1])
all_ASV = taxa_names(ps_season_asco_norm)
all_ASV <- all_ASV[!(all_ASV %in% otu_asco)]
ps_season_asco_norm=prune_taxa(all_ASV, ps_season_asco_norm)
ps_season_asco_norm

ps_sites_asco_norm <- readRDS("../03_Sites/00_PHYLOSEQ_OBJECTS/05_ps_sites_asco_norm.rds")
ps_sites_asco_norm <- subset_samples(ps_sites_asco_norm, ID !="L1I1_MARCH_BACT")
tax <- ps_sites_asco_norm@tax_table %>% as.data.frame()
asco <- tax %>% as.data.frame() %>% filter(Genus=="Ascophyllum") #1
otu_asco <- rownames(asco[1])
all_ASV = taxa_names(ps_sites_asco_norm)
all_ASV <- all_ASV[!(all_ASV %in% otu_asco)]
ps_sites_asco_norm=prune_taxa(all_ASV, ps_sites_asco_norm)
ps_sites_asco_norm

ps_parts_asco_norm <- readRDS("../04_Parts/00_PHYLOSEQ_OBJECTS/06_ps_parts_asco_norm.rds")
ps_parts_asco_norm <- subset_samples(ps_parts_asco_norm, ID !="L1R3_AVRIL_BACT")
ps_parts_asco_norm <- subset_samples(ps_parts_asco_norm, ID !="L1B4_AVRIL_BACT")
tax <- ps_parts_asco_norm@tax_table %>% as.data.frame()
asco <- tax %>% as.data.frame() %>% filter(Genus=="Ascophyllum") #1
otu_asco <- rownames(asco[1])
all_ASV = taxa_names(ps_parts_asco_norm)
all_ASV <- all_ASV[!(all_ASV %in% otu_asco)]
ps_parts_asco_norm=prune_taxa(all_ASV, ps_parts_asco_norm)
ps_parts_asco_norm

ps_season_site <- merge_phyloseq(ps_season_asco_norm,ps_sites_asco_norm) 
sample <- read.table("01_Tables/11_samples_sites_seasons.csv", sep=";", header=TRUE, row.names = 1)
ps_season_site <- phyloseq(otu_table(ps_season_site@otu_table), tax_table(ps_season_site@tax_table), sample_data(sample))
```

### Venn Parts <a name="venn_parts"></a>
```r
#!/bin/Rscript
plot_parts_abund <- MicEco::ps_venn(ps_parts_asco_norm, "Algae_part", 
                                     fraction = 1,
                                     quantities = list(type="counts", font = 0.5, cex=1.5), 
                                     labels = FALSE, col = "grey50", 
                                     fill = c("Apex"="#D67500","Medium"="#578B21","Base"="#B696D6","Receptacle"="#225EA8"),alpha=0.9,
                                     legend=list(cex=1.5, col="black", side="bottom"), 
                                     main=list(cex=2, labels='Parts', col="black"))

#Lists
ps_parts_list <- ps_venn(ps_parts_asco_norm, "Algae_part", quantities = list(type=c("percent","counts"), font = 1.5), labels = list(cex = 1.5), col = "white", fill = c("#69000C", "#006027", "#D35000", "#023FA5"),plot = FALSE, fraction = 1)
```

### Venns Sites & Seasons <a name="venn_sites_seasons"></a>
```r
#!/bin/Rscript
ps_season_site@sam_data$Season <- factor(ps_season_site@sam_data$Season, levels = c("Fall", "Winter", "Spring", "Summer")) 

plot_season_site_abund <- MicEco::ps_venn(ps_season_site, "Season", 
                                     fraction = 1,
                                     quantities = list(type="counts", font = 0.5, cex=1.5), 
                                     labels = FALSE, col = "grey50", 
                                     fill = c("Fall"="#963414",
                                              "Winter"="#6BAED6","Spring"="#238B45","Summer"="#f2b613"), alpha=0.9, 
                                     legend=list(cex=1.5, col="black", side="bottom"), 
                                     main=list(cex=2, labels='Summer', col="black"))



#Lists
ps_season_site_list <- ps_venn(ps_season_site, "Season", quantities = list(type=c("percent","counts"), font = 1.5), labels = list(cex = 1.5), col = "white", fill = c("#69000C", "#006027", "#D35000", "#023FA5"),plot = FALSE, fraction = 1)
```

### ps_core <a name="ps_core"></a>
```r
#!/bin/Rscript
parts_core <- ps_parts_list$Apex__Base__Medium__Receptacle #28 ASV
season_site_core <- ps_season_site_list$Fall__Winter__Spring__Summer #46 ASV

ps_part_core <- subset(otu_table(ps_parts_asco_norm), rownames(otu_table(ps_parts_asco_norm)) %in% parts_core)
ps_season_site_core <- subset(otu_table(ps_season_site), rownames(otu_table(ps_season_site)) %in% season_site_core)

ps_part_core <- merge_phyloseq(ps_part_core, tax_table(ps_parts_asco_norm), sample_data(ps_parts_asco_norm)) #22
ps_season_site_core <- merge_phyloseq(ps_season_site_core, tax_table(ps_season_asco_norm),  sample_data(ps_season_site)) #39

ps_core <- merge_phyloseq(ps_part_core,ps_season_site_core) #43
sample <- read.table("01_Tables/12_samples_sites_seasons_parts.csv", sep=";", header=TRUE, row.names = 1)
ps_core <- phyloseq(otu_table(ps_core@otu_table), tax_table(ps_core@tax_table), sample_data(sample)) #52
```

### Venn core <a name="venn_core"></a>
```r
#!/bin/Rscript
ps_core@sam_data$Parts <- factor(ps_core@sam_data$Parts, levels = c("Receptacle", "Apex", "Medium", "Base", "Pool")) 

plot_core <- MicEco::ps_venn(ps_core, "Parts", 
                                     fraction = 1,
                                     quantities = list(type="counts", font = 0.5, cex=1.5), 
                                     labels = FALSE, col = "grey50", 
                                     fill = c(  "Receptacle"="#225EA8","Apex"="#D67500", "Medium"="#578B21","Base"="#B696D6","Pool"="grey30"), alpha=0.9, 
                                     legend=list(cex=1.5, col="black", side="bottom"), 
                                     main=list(cex=2, labels='Core', col="black"))

ps_core_list <- ps_venn(ps_core, "Parts", quantities = list(type=c("percent","counts"), font = 1.5), labels = list(cex = 1.5), col = "white", fill = c("#69000C", "#006027", "#D35000", "#023FA5"),plot = FALSE, fraction = 1)

core <- ps_core_list$Receptacle__Apex__Medium__Base__Pool #18 ASV

ps_core_2 <- subset(otu_table(ps_core), rownames(otu_table(ps_core)) %in% core)
ps_core_2 <- merge_phyloseq(ps_core_2, tax_table(ps_core), sample_data(ps_core)) #18

```

## Combined Venn <a name="combined_venn"></a>
```r
#!/bin/Rscript
grid.arrange(arrangeGrob(arrangeGrob(plot_parts_abund, plot_season_site_abund, nrow = 1), plot_core), ncol=2) 
```

## 3. Specificities <a name="specificities"></a>
### ps specificity parts <a name="ps_specificity_parts"></a>
```r
#!/bin/Rscript
parts_apex <- ps_parts_list$Apex
parts_base <- ps_parts_list$Base 
parts_receptacle <- ps_parts_list$Receptacle 
parts_medium <- ps_parts_list$Medium 

ps_parts_asco_norm_R <- subset_samples(ps_parts_asco_norm, Algae_part =="Receptacle")
ps_parts_asco_norm_A <- subset_samples(ps_parts_asco_norm, Algae_part =="Apex")
ps_parts_asco_norm_M <- subset_samples(ps_parts_asco_norm, Algae_part =="Medium")
ps_parts_asco_norm_B <- subset_samples(ps_parts_asco_norm, Algae_part =="Base")


ps_part_apex <- subset(otu_table(ps_parts_asco_norm_A), rownames(otu_table(ps_parts_asco_norm_A)) %in% parts_apex)
ps_part_base <- subset(otu_table(ps_parts_asco_norm_B), rownames(otu_table(ps_parts_asco_norm_B)) %in% parts_base)
ps_part_receptacle <- subset(otu_table(ps_parts_asco_norm_R), rownames(otu_table(ps_parts_asco_norm_R)) %in% parts_receptacle)
ps_part_medium <- subset(otu_table(ps_parts_asco_norm_M), rownames(otu_table(ps_parts_asco_norm_M)) %in% parts_medium)


ps_part_apex <- merge_phyloseq(ps_part_apex, tax_table(ps_parts_asco_norm), sample_data(ps_parts_asco_norm)) #22
ps_part_base <- merge_phyloseq(ps_part_base, tax_table(ps_parts_asco_norm), sample_data(ps_parts_asco_norm)) #22
ps_part_receptacle <- merge_phyloseq(ps_part_receptacle, tax_table(ps_parts_asco_norm), sample_data(ps_parts_asco_norm)) #22
ps_part_medium <- merge_phyloseq(ps_part_medium, tax_table(ps_parts_asco_norm), sample_data(ps_parts_asco_norm)) #22
```


## ps specificity seasons <a name="ps_specificity_seasons"></a>

```r
#!/bin/Rscript
winter <- ps_season_site_list$Winter
fall <-  ps_season_site_list$Fall
spring <- ps_season_site_list$Spring 
summer <- ps_season_site_list$Summer 

ps_season_asco_norm_W <- subset_samples(ps_season_asco_norm, Season =="Winter")
ps_season_asco_norm_F <- subset_samples(ps_season_asco_norm, Season =="Fall")
ps_season_asco_norm_Sp <- subset_samples(ps_season_asco_norm, Season =="Spring")
ps_season_asco_norm_Sm <- subset_samples(ps_season_asco_norm, Season =="Summer")

ps_winter <- subset(otu_table(ps_season_asco_norm_W), rownames(otu_table(ps_season_asco_norm_W)) %in% winter)
ps_fall <- subset(otu_table(ps_season_asco_norm_F), rownames(otu_table(ps_season_asco_norm_F)) %in% fall)
ps_spring <- subset(otu_table(ps_season_asco_norm_Sp), rownames(otu_table(ps_season_asco_norm_Sp)) %in% spring)
ps_summer <- subset(otu_table(ps_season_asco_norm_Sm), rownames(otu_table(ps_season_asco_norm_Sm)) %in% summer)


ps_winter <- merge_phyloseq(ps_winter, tax_table(ps_season_asco_norm), sample_data(ps_season_asco_norm))
ps_fall <- merge_phyloseq(ps_fall, tax_table(ps_season_asco_norm), sample_data(ps_season_asco_norm))
ps_spring <- merge_phyloseq(ps_spring, tax_table(ps_season_asco_norm), sample_data(ps_season_asco_norm)) 
ps_summer <- merge_phyloseq(ps_summer, tax_table(ps_season_asco_norm), sample_data(ps_season_asco_norm)) 
```

## Bubble plot specificity part & season <a name="bubble_plot"></a>

```r
#!/bin/Rscript
## Part core
parts_core <- ps_parts_list$Apex__Base__Medium__Receptacle #28 ASV
ps_part_core <- subset(otu_table(ps_parts_asco_norm), rownames(otu_table(ps_parts_asco_norm)) %in% parts_core)
ps_part_core <- merge_phyloseq(ps_part_core, tax_table(ps_parts_asco_norm), sample_data(ps_parts_asco_norm)) #22

tax <- ps_part_core@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

sam <- ps_part_core@sam_data %>% as.data.frame()
sam <- cbind(rownames(sam), sam)
rownames(sam) <- NULL
colnames(sam)[1] <- "variable"

otu <- ps_part_core@otu_table %>% as.data.frame(stringsAsFactors = FALSE)
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"
otu %<>% melt(id.vars = "ASV")

merged <- merge(otu, tax, by="ASV")
merged <- merge(merged,sam, by="variable")
#write.table(merged, "01_Tables/15_heatmap_core_parts.csv", sep=",", quote=FALSE)

## part specificity
ps_specificity <- merge_phyloseq(ps_part_apex,ps_part_medium, ps_part_receptacle, ps_part_base) #43

tax <- ps_specificity@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

sam <- ps_specificity@sam_data %>% as.data.frame()
sam <- cbind(rownames(sam), sam)
rownames(sam) <- NULL
colnames(sam)[1] <- "variable"

otu <- ps_specificity@otu_table %>% as.data.frame(stringsAsFactors = FALSE)
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"
otu %<>% melt(id.vars = "ASV")

merged <- merge(otu, tax, by="ASV")
merged <- merge(merged,sam, by="variable")
#write.table(merged, "01_Tables/14_heatmap_specificity_parts.csv", sep=",", quote=FALSE)


## season core
season_core <- ps_season_site_list$Fall__Winter__Spring__Summer #28 ASV
ps_season_core <- subset(otu_table(ps_season_asco_norm), rownames(otu_table(ps_season_asco_norm)) %in% season_core)
ps_season_core <- merge_phyloseq(ps_season_core, tax_table(ps_season_asco_norm), sample_data(ps_season_asco_norm)) #22

tax <- ps_season_core@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

sam <- ps_season_core@sam_data %>% as.data.frame()
sam <- cbind(rownames(sam), sam)
rownames(sam) <- NULL
colnames(sam)[1] <- "variable"

otu <- ps_season_core@otu_table %>% as.data.frame(stringsAsFactors = FALSE)
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"
otu %<>% melt(id.vars = "ASV")

merged <- merge(otu, tax, by="ASV")
merged <- merge(merged,sam, by="variable")
#write.table(merged, "01_Tables/15_heatmap_core_season.csv", sep=",", quote=FALSE)

## season specificity
ps_specificity <- merge_phyloseq(ps_winter,ps_fall, ps_summer, ps_spring) #43

tax <- ps_specificity@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

sam <- ps_specificity@sam_data %>% as.data.frame()
sam <- cbind(rownames(sam), sam)
rownames(sam) <- NULL
colnames(sam)[1] <- "variable"

otu <- ps_specificity@otu_table %>% as.data.frame(stringsAsFactors = FALSE)
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"
otu %<>% melt(id.vars = "ASV")

merged <- merge(otu, tax, by="ASV")
merged <- merge(merged,sam, by="variable")
#write.table(merged, "01_Tables/14_heatmap_specificity_season.csv", sep=",", quote=FALSE)

```

```r
#!/bin/Rscript
merged <- read.table("01_Tables/15_heatmap_specificity_parts_core_season.csv", sep=";", header=TRUE)

merged$Specificity <- factor(merged$Specificity, levels = c("Core_algae_part","Receptacle", "Apex","Medium","Base", "Core_season", "Fall", "Winter", "Spring", "Summer")) 
#merged$Type <- factor(merged$Type, levels = c("Core_algae_part","Core_season", "Specificity_algae_part", "Specificity_season")) 
merged$Class <- factor(merged$Class, levels = c("Alphaproteobacteria","Gammaproteobacteria","Verrucomicrobiae", "Bacteroidia","Planctomycetes", "Phycisphaerae",  "Ascomycota", "Bacteria,NA,NA","Bdellovibrionia","Cyanobacteriia",
  "Phaeophyceae",
  "Planctomycetacia",
  "OM190",
  "Oligoflexia",
  "Blastocatellia",
  "SAR324 clade(Marine group B),NA")) 

merged$Genus <- factor(merged$Genus, levels =c("Robiginitomaculum", "Litorimonas", "Hyphomonadaceae,NA", "Rhodobacteraceae,NA", "Yoonia-Loktanella", "Octadecabacter", "Roseovarius", "Erythrobacter", "Sphingorhabdus", "Altererythrobacter", "Hellea", "Alphaproteobacteria,NA,NA,NA", "Litoreibacter", "Sulfitobacter", "Brevundimonas", "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Silicimonas", "Arenicella", "Thalassotalea", "Granulosicoccus", "Gammaproteobacteria,NA,NA,NA", "Woeseia", "Acinetobacter", "pItb-vmat-80,NA,NA", "Pseudomonas", "Alteromonadaceae,NA", "Marinomonas", "Rubritalea", "DEV007,NA", "Luteolibacter", "Roseibacillus", "Lewinella", "Portibacter", "Saprospiraceae,NA", "Rubidimonas", "Pricia", "Pibocella", "Nonlabens", "Flavobacteriaceae,NA", "Marixanthomonas", "Chitinophagales,NA,NA", "Croceitalea", "Maribacter", "Dokdonia", "Winogradskyella", "Rubripirellula", "Bythopirellula", "Blastopirellula", "Planctomycetales,NA,NA", "Rhodopirellula", "Rubinisphaeraceae,NA", "Fuerstia", "Algisphaera", "Dothideomycetes,NA", "Bacteria,NA,NA,NA,NA,NA", "Bdellovibrio", "Peredibacter", "Pleurocapsa PCC-7319", "OM190,NA,NA,NA", "Oligoflexaceae,NA", "Blastocatella", "SAR324 clade(Marine group B),NA,NA,NA,NA"))
                         
palette_class=c(
  "Gammaproteobacteria"="#B03F00",
  "Alphaproteobacteria"="#08306B",
  "Planctomycetes"="#645A9F",
  "Bacteroidia"="#853250",
  "Bdellovibrionia"="#299B8A",
  "Cyanobacteriia"="#00441B",
  "Ascomycota"="#238B45",
  "Verrucomicrobiae"="#EAAF29",
  "Phaeophyceae"="#225EA8",
  "Phycisphaerae"="#FFADC1",
  "Planctomycetacia"="#645A9F",
  "OM190"="#9EBCDA",
  "Oligoflexia"="grey50",
  "Blastocatellia"="#769F3A",
  "Bacteria,NA,NA"="grey",
  "SAR324 clade(Marine group B),NA"="#E39879")

ggplot(merged, aes(x = Individual, y = Genus, fill=Class)) + 
  facet_nested(~Specificity, space = "free", scales = "free")+
  geom_point(aes(size = value, fill = Class), alpha = 1, shape = 21) +   
    scale_size_continuous(limits = c(0.000001, 190000), range = c(3)) +
   scale_color_manual(values = palette_class) + scale_fill_manual(values = palette_class) +
  theme(legend.key=element_blank(), 
        strip.text.y = element_text(size=10,angle = 0),
  axis.text.x = element_blank(), 
  axis.text.y = element_text(colour = "black", size = 9), 
  legend.text = element_text(size = 10, colour ="black"), 
  legend.title = element_text(size = 10), 
  panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
  legend.position = "right")#+ 
 # pdf("02_Figures/16_core_specificity_genus.pdf", height = 8, width = 15)
```


![Figure 4 | Core microbiome ](https://github.com/rssco/novaseq_ascophyllum/blob/main/01_Figures/Figure_4_Combined_specificity_venn_genus.png)<!-- -->