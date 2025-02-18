# Table of content  
[1. Bacteria](#miseq_bacteria)  
- [DADA2](#miseq_bacteria_dada2)   
- [Phyloseq](#miseq_phyloseq)  
- [Decontamination with microdecon](#miseq_decontamination)  
- [Remove mock](#miseq_remove_mock)  
- [Remove plastid reads](#miseq_remove_plastid_reads)  
  
[2. Fungi](#miseq_fungi)  
- [Cutadapt](#miseq_fungi_cutadapt)  
- [DADA2](#miseq_fungi_dada2)  
- [Phyloseq](#miseq_fungi_ps)
- [Decontamination with Micrdecon](#miseq_fungi_microdecon)
- [Remove mock](#miseq_fungi_remove_mock)


# 1. Library <a name="library"></a>
```r
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

# 2. Palette colors <a name="colors"></a>
```r
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
                   "November" = "#963414",
                   "December"="#733c05")
palette_season <- c("Winter"="#6BAED6","Spring"="#238B45","Summer"="#f2b613","Fall"="#963414")
palette_site <- c("#F8766D", "#00BF7D", "#00B0F6")
palette_parts <- c(  "Apex"="#D67500",
                     "Medium"="#578B21",
                     "Base"="#B696D6",
                     "Receptacle"="#225EA8")

palette_Parts_nmds <- c(
  "Apex"="#D67500", 
  "Medium"="#578B21",
  "Base"="#B696D6",
  "Receptacle"="#225EA8",
  "Pool"="grey30",
  "Pool_Site_1"="#F8766D",
  "Pool_Site_2"="#00BF7D",
  "Pool_Site_3"="#00B0F6"
)
```

# 3. NMDS <a name="nmds"></a>
## ps object
```r
ps_season <- readRDS("../06_Seasons/00_PHYLOSEQ_OBJECTS/02_ps_season_asco.rds")
ps_parts <- readRDS("../04_Parts/00_PHYLOSEQ_OBJECTS/05_ps_parts_asco.rds") 
ps_sites <- readRDS("../03_Sites/00_PHYLOSEQ_OBJECTS/04_ps_sites_asco.rds") 
```

## theme 
```r
theme_nmds <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_blank(), 
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(), 
                    panel.grid.major.x = element_blank(), 
                    panel.border = element_rect(color="black", fill=NA, size=1), 
                    aspect.ratio = 1)
```

## season
```r
ps_season <- prune_samples(sample_names(ps_season) != "16S_L1I1_MARS_BACT_S81_R1", ps_season) #L1I1MARS

# Distance
iDist <- distance(ps_season, method = "bray")
pn.nmds = ordinate(ps_season, 
                method = "NMDS", 
                distance = iDist)
pn.nmds
#Stress=0.1517905

# plot
plot.pn.nmds = plot_ordination(ps_season, pn.nmds, justDF = TRUE)
plot.pn.nmds$Season <- factor(plot.pn.nmds$Season, levels = c("Winter", "Spring", "Summer", "Fall")) 

season <- ggplot(plot.pn.nmds, aes(x = NMDS1, y = NMDS2, color = Season)) + 
        geom_point(size = 3, alpha = 0.75) + scale_color_manual(values = palette_season) + 
  annotate(geom="text", x=-0.7, y=-0.6, label="Stress=15%", size=3) +
        theme_nmds + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin()) + ggtitle("C. Season")


# permanova
p.df = as(sample_data(ps_season), "data.frame")
p.d = distance(ps_season, method = "bray")
p.adonis = adonis(p.d ~ Season+ Month , p.df)
p.adonis$aov.tab

#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Season     3    0.8141 0.271369  4.5595 0.18771  0.001 ***
#Month      9    1.2612 0.140137  2.3545 0.29081  0.001 ***
#Residuals 38    2.2617 0.059518         0.52148           
#Total     50    4.3370                  1.00000       
```
## parts
```r
ps_parts_wt_LIB4 <- prune_samples(sample_names(ps_parts) != "16S_L1R3_AVRIL_BACT_S156_R1", ps_parts)
ps_parts_wt_LIB4 <- prune_samples(sample_names(ps_parts_wt_LIB4) != "16S_L1B4_AVRIL_BACT_S159_R1", ps_parts_wt_LIB4)


# Distance
iDist <- distance(ps_parts_wt_LIB4, method = "bray")
pn.nmds = ordinate(ps_parts_wt_LIB4, 
                method = "NMDS", 
                distance = iDist)
pn.nmds
#Stress=0.1366201

ps_parts_wt_LIB4@sam_data$Algae_part <- factor(ps_parts_wt_LIB4@sam_data$Algae_part, levels = c("Receptacle", "Apex","Medium", "Base"))


# plot
plot.pn.nmds = plot_ordination(ps_parts_wt_LIB4, pn.nmds, justDF = TRUE)
parts <- ggplot(plot.pn.nmds, aes(x = NMDS1, y = NMDS2, color = Algae_part, shape=Month))+
        geom_point(size = 4) + scale_color_manual(values = palette_parts)+ 
    annotate(geom="text", x=-0.5, y=-0.8, label="Stress=14%", size=4) +
        theme_nmds + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin(),
                           legend.title = element_text(size=15), 
                           legend.text = element_text(size=12)) + ggtitle("B. Algae parts")
        
# permanova
p.df = as(sample_data(ps_parts_wt_LIB4), "data.frame")
p.d = distance(ps_parts_wt_LIB4, method = "bray")
p.adonis = adonis(p.d ~ Algae_part + Month + Age, p.df)
p.adonis$aov.tab

#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Algae_part  3   0.90589 0.30196  4.6197 0.30726  0.001 ***
#Month       1   0.34489 0.34489  5.2764 0.11698  0.001 ***
#Age         1   0.12880 0.12880  1.9705 0.04369  0.092 .  
#Residuals  24   1.56872 0.06536         0.53208           
#Total      29   2.94830                 1.00000    
```


## sites 
```r
ps_sites <- prune_samples(sample_names(ps_sites) != "16S_L1I1_MARS_BACT_S81_R1", ps_sites) #L1I1MARS

# Distance
iDist <- distance(ps_sites, method = "bray")
pn.nmds = ordinate(ps_sites, 
                method = "NMDS", 
                distance = iDist)
pn.nmds
#Stress=0.1242916

# plot
plot.pn.nmds = plot_ordination(ps_sites, pn.nmds, justDF = TRUE)
sites <- ggplot(plot.pn.nmds, aes(x = NMDS1, y = NMDS2, color = Site, shape=Month)) + 
        geom_point(size = 4) + scale_color_manual(values = palette_site) + 
  annotate(geom="text", x=-0.5, y=-0.4, label="Stress=12%", size=4) +
        theme_nmds + theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin(), 
                           legend.title = element_text(size=15), 
                           legend.text = element_text(size=12)) +ggtitle("A. Sites")

ggplot(plot.pn.nmds, aes(x = NMDS1, y = NMDS2, color = Site, shape=Month)) + 
        geom_point(size = 3, alpha = 0.75) + scale_color_manual(values = palette_site) + 
  annotate(geom="text", x=-0.5, y=-0.4, label="Stress=11%", size=3) + theme(legend.position = c(1,0.5))

# permanova
p.df = as(sample_data(ps_sites), "data.frame")
p.d = distance(ps_sites, method = "bray")
p.adonis = adonis(p.d ~ Site + Month, p.df)
p.adonis$aov.tab

#          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Site       2   0.25480 0.127400  1.6097 0.12891  0.072 .
#Month      1   0.21809 0.218093  2.7557 0.11034  0.015 *
#Residuals 19   1.50373 0.079144         0.76076         
#Total     22   1.97662                  1.00000     

#Dispersion
otu <- t(otu_table(ps_sites))
sample_info <- data.frame(ps_sites@sam_data)

bd <- betadisper(iDist, sample_info$Month)
anova(bd)

#Response: Distances
#          Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.023247 0.0232467  2.8128 0.1083
#Residuals 21 0.173559 0.0082647 

```

## NMDS figures together 

```r
figures_combined <- grid.arrange(sites, parts, season, nrow=1) 
ggsave("02_Figures/05_figures_combined_nmds2.pdf", figures_combined,width = 15, height = 10)
```

# 4. Alpha diversity
## ps objects
```r
ps_season_asco_norm <- readRDS("../06_Seasons/00_PHYLOSEQ_OBJECTS/03_ps_season_asco_norm.rds")
ps_season_asco_norm <- subset_samples(ps_season_asco_norm, ID !="L1I1_MARCH")
ps_sites_asco_norm <- readRDS("../03_Sites/00_PHYLOSEQ_OBJECTS/05_ps_sites_asco_norm.rds")
ps_sites_asco_norm <- subset_samples(ps_sites_asco_norm, ID !="L1I1_MARCH_BACT")
ps_parts_asco_norm <- readRDS("../04_Parts/00_PHYLOSEQ_OBJECTS/06_ps_parts_asco_norm.rds")
ps_parts_asco_norm <- subset_samples(ps_parts_asco_norm, ID !="L1R3_AVRIL_BACT")
ps_parts_asco_norm <- subset_samples(ps_parts_asco_norm, ID !="L1B4_AVRIL_BACT")
```

```r
source("https://raw.githubusercontent.com/microbiome/microbiome/master/R/boxplot_alpha.R")
```

## season

```r
tab <-microbiome::alpha(ps_season_asco_norm, index = "all")
ps1_meta <- meta(ps_season_asco_norm)

ps1_meta$Shannon <- tab$diversity_shannon 
ps1_meta$InverseSimpson <- tab$diversity_inverse_simpson
ps1_meta <- phyloseq::sample_data(ps1_meta)

physeq_alpha <- phyloseq::merge_phyloseq(ps_season_asco_norm, ps1_meta)
metadata <- data.frame(sample_data(physeq_alpha))

# create a list of pairwise comparisons
bmi <- levels(ps1_meta$Month) # get the variables
bmi <- unique(ps1_meta$Month)
print(bmi)

# make a pairwise list that we want to compare.
bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])

print(bmi.pairs)
```

```r
# MONTH

aov_shannon <- stats::aov(Shannon ~ Month , metadata)
summary(aov_shannon)
#Month       12  4.217  0.3514   3.254 0.0027 **
#Residuals   38  4.103  0.1080    

pairwise_test <- rstatix::tukey_hsd(metadata, Shannon ~ Month) 
pairwise_test %<>% filter(!p.adj.signif=="ns")
pairwise_test

#1 Month June   May                 0   -0.828   -1.65   -0.00778 0.046    *           
#2 Month March  May                 0   -1.02    -1.90   -0.130   0.0131   *           
#3 Month May    November21          0    1.00     0.181   1.82    0.00634  **          
#4 Month May    October             0    0.967    0.147   1.79    0.00951  **          
#5 Month May    September           0    1.15     0.333   1.97    0.000933 ***  
```

```r
tukey <- TukeyHSD(aov_shannon)
tukey

cld <- multcompLetters4(aov_shannon, tukey)
cld

metadata$Month <- factor(metadata$Month, levels = c("November21", "December", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November22")) 

TK <- group_by(metadata,Month) %>% summarise(mean=mean(Shannon), quant=quantile(Shannon, probs=0.75)) %>% arrange(desc(mean))

cld <- as.data.frame.list(cld$Month)
TK$cld <- cld$Letters
TK

season <- ggplot(metadata, aes(x = Month, y = Shannon, color=Month)) + 
  geom_boxplot(aes(fill=Month),alpha=0.6, lwd=0.2) + 
  geom_jitter(shape=16, position=position_jitter(0)) + 
    scale_color_manual(values = palette_month) +
  scale_fill_manual(values = palette_month) +
  geom_text(data = TK, aes(label=cld, x=Month, y=quant), vjust=-1, hjust=-0.6, size=3) + 
  theme_bw(base_line_size = 0) +
  theme(
axis.title.x = element_blank(), 
                    axis.title.y = element_blank(), 
                    axis.text.x=element_blank(),
                    axis.ticks=element_blank()) 
```

## parts

```r
tab <-microbiome::alpha(ps_parts_asco_norm, index = "all")
ps1_meta <- meta(ps_parts_asco_norm)

ps1_meta$Shannon <- tab$diversity_shannon 
ps1_meta$InverseSimpson <- tab$diversity_inverse_simpson
ps1_meta <- phyloseq::sample_data(ps1_meta)

physeq_alpha <- phyloseq::merge_phyloseq(ps_parts_asco_norm, ps1_meta)
metadata <- data.frame(sample_data(physeq_alpha))

# create a list of pairwise comparisons
bmi <- levels(ps1_meta$Algae_part) # get the variables
bmi <- unique(ps1_meta$Algae_part)
print(bmi)

# make a pairwise list that we want to compare.
bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])

print(bmi.pairs)
```

```r
factor(metadata$Month)

aov_shannon <- stats::aov(Shannon ~ Month + Algae_part , metadata)
summary(aov_shannon)
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#Month        1  0.008  0.0082   0.046  0.832
#Algae_part   3  0.825  0.2749   1.536  0.230
#Residuals   25  4.475  0.1790   
```

```r
metadata$Algae_part <- factor(metadata$Algae_part, levels = c("Receptacle","Apex", "Medium", "Base")) 
metadata$Month <- factor(metadata$Month, levels = c("January", "April")) 

parts <- ggplot(metadata, aes(x = Algae_part, y = Shannon, color=Algae_part)) + 
geom_boxplot(aes(fill=Algae_part),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) + scale_color_manual(values = palette_parts) +
  scale_fill_manual(values = palette_parts) +
theme_bw(base_line_size = 0) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 1, legend.position="none")
```

## sites

```r
tab <-microbiome::alpha(ps_sites_asco_norm, index = "all")
ps1_meta <- meta(ps_sites_asco_norm)

ps1_meta$Shannon <- tab$diversity_shannon 
ps1_meta$InverseSimpson <- tab$diversity_inverse_simpson
ps1_meta <- phyloseq::sample_data(ps1_meta)

physeq_alpha <- phyloseq::merge_phyloseq(ps_sites_asco_norm, ps1_meta)
metadata <- data.frame(sample_data(physeq_alpha))

# create a list of pairwise comparisons
bmi <- levels(ps1_meta$Site) # get the variables
bmi <- unique(ps1_meta$Site)
print(bmi)

# make a pairwise list that we want to compare.
bmi.pairs <- combn(seq_along(bmi), 2, simplify = FALSE, FUN = function(i)bmi[i])

print(bmi.pairs)
```

```r
factor(metadata$Site)

aov_shannon <- stats::aov(Shannon ~ Site + Month, metadata)
summary(aov_shannon)
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#Site         2  0.101  0.0507   0.270 0.7666  
#Month        1  0.810  0.8098   4.304 0.0519 .
#Residuals   19  3.575  0.1882  
```

## Theme
```r
theme_alpha <- theme(axis.line = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(), 
                    axis.title.x = element_blank(), 
                    axis.text.x=element_blank(),
                    axis.ticks=element_blank(), 
                    panel.grid.major.x = element_blank(), 
                    panel.border = element_rect(color="black", fill=NA, size=0.3), 
                    aspect.ratio = 1)
```


```r
metadata$Month <- factor(metadata$Month, levels = c("March","November")) 

site <- ggplot(metadata, aes(x = Month, y = Shannon, color=Month)) + 
geom_boxplot(aes(fill=Month),alpha=0.4, lwd=0.2) + geom_jitter(shape=16, position=position_jitter(0)) + scale_color_manual(values = palette_month) +
  scale_fill_manual(values = palette_month) +
theme_bw(base_line_size = 0) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 1, 
         legend.position="none")

```

![Figure 1 | Variation of SSU diversity and communities across the three datasets.](https://github.com/rssco/novaseq_ascophyllum/blob/main/01_Figures/Figure_1_Combined_nmds_alpha.pdf)<!-- -->

