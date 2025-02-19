## 1. Library
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

## 2. ASV order with pearson

```r
#!/bin/Rscript
nodes <- read.csv("../../07_Coocurrence_networks/03_Networks/07_nodes_modules_plus_easi.csv", sep=",", header=TRUE, dec=".")
ps_filt <- readRDS("../../07_Coocurrence_networks/00_Phyloseq_objects/04_ps_bact_fungi_45_filt_ss18Sasco.rds") # 45 reads in 3samples
tax <- nodes %>% select(Id, module,kingdom:species) %>% tibble::column_to_rownames("Id") %>% as.matrix

ps_net_plus <- phyloseq(otu_table(ps_filt@otu_table), tax_table(tax), sample_data(ps_filt@sam_data))
ps_net_plus <- subset_samples(ps_net_plus, Dataset_types == "Season")
ps_net_plus <- prune_taxa(taxa_sums(ps_net_plus) > 0, ps_net_plus)

tax <- ps_net_plus@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

sam <- ps_net_plus@sam_data %>% as.data.frame()
sam <- cbind(rownames(sam), sam)
rownames(sam) <- NULL
colnames(sam)[1] <- "variable"

otu <- ps_net_plus@otu_table %>% as.data.frame() 
#write.table(otu, "~/Desktop/otu_genus_core.csv", sep=";", quote=FALSE)
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"
otu %<>% melt(id.vars = "ASV")

#find ASV order
otu <- read.table("~/Desktop/otu_genus_core.tsv", sep="\t", header=TRUE, row.names = 1) %>% t()
cor_matrix <- cor(otu, method = "pearson", use = "pairwise.complete.obs")
dist_matrix <- as.dist(1 - cor_matrix) # Convertir la corrélation en distance
hc <- hclust(dist_matrix, method = "average")
asv_order <- rownames(cor_matrix)[hc$order]


merged <- merge(otu, tax, by = "ASV")
merged <- merge(merged, sam, by = "variable")

#write.table(merged, "01_Tables/16_core_heatmap_genus.tsv", sep="\t", quote = FALSE)
merged <- read.table("01_Tables/16_core_heatmap_genus.tsv", sep="\t", header=TRUE, dec=".")
```

## 3. Heatmap 
```r
#!/bin/Rscript
merged$Month <- factor(merged$Month, levels = c("November21","December","January", "February", "March", "April", "May", "June", "July", "August", "September", "October",  "November22")) 
merged$Dataset_types <- factor(merged$Dataset_types, levels = c("Sites","Algal_parts","Season")) 
merged$Type <- factor(merged$Type, levels = c("Site_2", "Site_3", "Receptacle","Apex","Medium","Base", "Individual")) 
merged$OTU <- factor(merged$OTU, levels = asv_order)

colors1 <- brewer.pal(11, "RdYlBu")

ggplot(merged, aes(x=variable,y=OTU,fill=value)) + geom_tile() + 
    facet_nested(genus~Month,space = "free", scales = "free")+
  scale_fill_gradientn(colours =rev(colors1),name="Abundance", trans='pseudo_log') +
  theme(axis.text.x = element_text(size=8, angle=90,hjust=0.95,vjust=0.5), 
        axis.text.y = element_text(size=8), 
        panel.background = element_blank(), 
        legend.title = element_text(size=9, face = "bold"),
         strip.text.y = element_text(size=9,angle = 0, face = "bold.italic"),
        panel.spacing.x  = unit(0.1, "lines"),
        panel.spacing.y  = unit(0.05, "lines")) +
  xlab(NULL) + ylab(NULL) + 
  pdf("02_Figures/21_heatmap_genus_core.pdf", height = 10, width = 10)

ggplot(merged, aes(x=variable,y=OTU,fill=value)) + geom_tile() + 
    facet_nested(~Month,space = "free", scales = "free")+
  scale_fill_gradientn(colours =rev(colors1),name="Abundance", trans='pseudo_log') +
  theme(axis.text.x = element_text(size=8, angle=90,hjust=0.95,vjust=0.5), 
        axis.text.y = element_text(size=7), 
        panel.background = element_blank(), 
        legend.title = element_text(size=9, face = "bold"),
         strip.text.y = element_text(size=9,angle = 0, face = "bold.italic"),
        panel.spacing.x  = unit(0.1, "lines"),
        panel.spacing.y  = unit(0.05, "lines")) +
  xlab(NULL) + ylab(NULL) + 
  pdf("02_Figures/21_heatmap_genus_core2.pdf", height = 8.5, width = 8)
```

![Figure 6 | Heatmap genus core microbiome ](https://github.com/rssco/novaseq_ascophyllum/blob/main/01_Figures/Figure_6_Heatmap_genus_core.png)<!-- -->