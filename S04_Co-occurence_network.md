## Table of content 
[1. Library](#library)  

[2. Fungi OTU 99%](#Fungi_OTU_99%)  
- [Prepare tables for vsearch](#Prepare_tables_for_vsearch)  
- [Generate fasta](#Generate_fasta)
- [Launch Vsearch](#Launch_Vsearch)
- [Final table](#Final_table)

[3. Bacteria OTU 99%](#Bacteria_OTU_99%)   
- [Prepare tables for vsearch](#Prepare_tables_for_vsearch_16S)   
- [Launch Vsearch](#Launch_Vsearch_16S)  
- [Final table](#Final_table_16S)  
  
[4. Phyloseq 16+18S and ITS](#ps_16+18S+ITS)  
  
[5. Filtering steps](#Filtering_steps)  
  
[6. Covariance estimation using SPIEC-EASI](#SPIEC_EASI)  

[7. Positive edges](#Positive_edges)      
- [Clustering algorithm compairison (modules)](#Clustering)  
- [Leiden or Spinglass or Louvain ?](#Leidein_Springlass_Louvain)    
- [Leiden algorithm for modularity](#Leiden)  
- [Observe network plot](#plot)  
  
[8. Negative edges](#Negative_edges)  
  
[9. Modules characteristics](#Modules_characteristics)
- [Add metadata to OTU and ps_gephi](#metadata)
- [Module palette](#palettes)
- [Plot seasons](#plot_seasons)
- [Plot parts](#plot_parts)
- [Module composition](#module_composition)
- [Combined figures](#combined_figures)





## 1. Library <a name="library"></a>
```r
#!/bin/Rscript
library(tidyverse)
library(magrittr)
library(phyloseq)
library(SpiecEasi)
        #installation not works 
        #mkdir .R (home directory)
        #echo "FLIBS=-L/opt/local/lib/gcc48/" > .R/Makevars
library(igraph) #network
library(gephi) #network
library(here) #network
library(ggdag) #network
library(tidygraph)
library(Matrix)
library(reshape2)
library(gridExtra)
library(ggh4x)
library(WGCNA)
library(Hmisc) #corr
library(ggtree) #tree
library(Biostrings) #tree
library(RColorBrewer)
library(aplot)
library(multcompView)
```

<https://loimai.github.io/BBobs/Lemonnier_Ushant_front_16S.html#Matrix_filtration>


## 2. Fungi OTU 99% <a name="Fungi_OTU_99%"></a>
### Prepare tables for vsearch <a name="Prepare_tables_for_vsearch"></a>

```r
#!/bin/Rscript
ps_fungi <- readRDS("../05_Fungi_analysis/01_All/02_Phyloseq_objects/04_ps_ss_decon_wt_cc_mock.rds")

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
write.table(otu_tax_ref, "../05_Fungi_analysis/01_All/03_Tables/09_otu_tax_ref_for_vsearch.tsv", sep="\t", quote=FALSE)
```

### Generate fasta file with ASV_Species+sequence <a name="Generate_fasta"></a>
```bash
#!/bin/bash
cat 09_otu_tax_ref_for_vsearch.tsv| sort -n -k 2 | cut -f1,130 | tr "\t" "\n" | sed "s/^ASV/>ASV/g" > fungi.fasta
```

### Launch Vsearch <a name="Launch_Vsearch"></a>
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
output=../vsearch_ns/fungi_consensus99.fasta
name=cluster

vsearch --cluster_smallmem ${input} --usersort --id 0.99 --sizeout --sizein --consout ${output} --clusters /scratch2/seabioz/tmp/${name}_ 
```

```bash
#!/bin/bash
for i in {0..905}; do grep ">" cluster_${i} | sed "s/;size=1//g" | sed "s/>//g" | sed -e 's/,$//g' -e 's/ //g' > vector_${i}.txt; done 

for i in *.txt; do awk '{print $0, FILENAME}' $i |sed 's/\.txt//g' | awk 'NR==1 {line=$0} NR>=1 {$3=line" "$3} 1'  |sed 's/ /\t/g' | cut -f1,2,3 >  ${i%.txt}_v.txt; done

cat *v.txt > finale.txt
```

### Generate final table <a name="Final_table"></a>
```r
#!/bin/Rscript
#merge 2 tables by hand because doesn't work well with R. 
#otu_tax_ref <- read.table("03_Tables/09_otu_tax_ref_for_vsearch.tsv", header=TRUE, sep="\t")
#vector <- read.table("03_Tables/10_asv_vector_representant_finale.txt", header=TRUE)

otu_tax_ref_vector <- read.table("../05_Fungi_analysis/01_All/03_Tables/11_otu_tax_ref_for_vsearch_vector_representant_vsearch99.tsv", header=TRUE, sep="\t")
tax_seq <- otu_tax_ref_vector %>% dplyr::select(Kingdom:Species, sequence, ASV)

otu_tax_ref_vector_agg <- rowsum(otu_tax_ref_vector[,c(4:123)],otu_tax_ref_vector$Representant,na.rm=T)
otu_tax_ref_vector_agg$ASV <- rownames(otu_tax_ref_vector_agg)


otu_tax_ref_agg_fungi <- merge(otu_tax_ref_vector_agg, tax_seq, by="ASV")
write.table(otu_tax_ref_agg_fungi,"../05_Fungi_analysis/01_All/03_Tables/12_otu_tax_ref_agg_fungi.tsv", sep="\t", quote=FALSE)

#cp 12_otu_tax_ref_agg_fungi.tsv in 07_Coocurrence_networks/01_Tables/01_otu_tax_ref_agg_fungi.csv
```

## 3. Bacteria otu 99% <a name="Bacteria_OTU_99"></a>
### Prepare tables for vsearch <a name="Prepare_tables_for_vsearch_16S"></a>

```r
#!/bin/Rscript
ps_season_asco_norm <- readRDS("../04_Bacteria_analysis/06_Seasons/00_PHYLOSEQ_OBJECTS/02_ps_season_asco.rds")
ps_season_asco_norm <- subset_samples(ps_season_asco_norm, ID !="L1I1_MARCH")

ps_sites_asco_norm <- readRDS("../04_Bacteria_analysis/03_Sites/00_PHYLOSEQ_OBJECTS/04_ps_sites_asco.rds")
ps_sites_asco_norm <- subset_samples(ps_sites_asco_norm, ID !="L1I1_MARCH_BACT")

ps_parts_asco_norm <- readRDS("../04_Bacteria_analysis/04_Parts/00_PHYLOSEQ_OBJECTS/05_ps_parts_asco.rds")
ps_parts_asco_norm <- subset_samples(ps_parts_asco_norm, ID !="L1R3_AVRIL_BACT")
ps_parts_asco_norm <- subset_samples(ps_parts_asco_norm, ID !="L1B4_AVRIL_BACT")

ps_season_site <- merge_phyloseq(ps_season_asco_norm,ps_sites_asco_norm) 
sample <- read.table("../04_Bacteria_analysis/05_All/01_Tables/11_samples_sites_seasons.csv", sep=";", header=TRUE, row.names = 1)
ps_season_site <- phyloseq(otu_table(ps_season_site@otu_table), tax_table(ps_season_site@tax_table), sample_data(sample))

ps <- merge_phyloseq(ps_season_asco_norm,ps_sites_asco_norm, ps_parts_asco_norm) 
#otu_table()   OTU Table:         [ 16894 taxa and 97 samples ]
#sample_data() Sample Data:       [ 97 samples by 10 sample variables ]
#tax_table()   Taxonomy Table:    [ 16894 taxa by 6 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 16894 reference sequences ]

ref <- ps@refseq %>% as.data.frame()
colnames(ref)[1] <- "sequence"
ref <- cbind(rownames(ref), ref)
rownames(ref) <- NULL
colnames(ref)[1] <- "ASV"

tax <- ps@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

ref_tax <- merge(tax, ref, by="ASV")
otu_tax_ref <- merge(otu, ref_tax, by="ASV")
write.table(otu_tax_ref, "01_Tables/02_otu_tax_ref_agg_bacteria.csv", sep=";", quote=FALSE)
```

### launch Vsearch <a name="Launch_Vsearch_16S"></a>
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
input=/shared/projects/seabioz/finalresult/test/vsearch_ns/06_ALL_BACTERIA/bact.fasta
output=/shared/projects/seabioz/finalresult/test/vsearch_ns/06_ALL_BACTERIA/bact_consensus99.fasta
name=cluster
vsearch --cluster_smallmem ${input} --usersort --id 0.99 --sizeout --sizein --consout ${output} --clusters /scratch2/seabioz/tmp/${name}_ 
```

```bash
#!/bin/bash
for i in {0..11846}; do grep ">" cluster_${i} | sed "s/;size=1//g" | sed "s/>//g" | sed -e 's/,$//g' -e 's/ //g' > vector_${i}.txt; done
```

```bash
#!/bin/bash
for i in *.txt; do awk '{print $0, FILENAME}' $i |sed 's/\.txt//g' | awk 'NR==1 {line=$0} NR>=1 {$3=line" "$3} 1'  |sed 's/ /\t/g' | cut -f1,2,3 >  ${i%.txt}_v.txt; done

cat *v.txt > finale.txt
```

### Generate final table <a name="Final_table_16S"></a>
```bash
#!/bin/bash
#merge 2 tables by hand because doesn't work well with R (02_otu_tax_ref_bacteria.csv+final.txt) 
otu_tax_ref_vector <- read.table("01_Tables/02_otu_tax_ref_bacteria.tsv", header=TRUE, sep="\t")
tax_seq <- otu_tax_ref_vector %>% dplyr::select(Kingdom:Genus, sequence, ASV)

otu_tax_ref_vector_agg <- rowsum(otu_tax_ref_vector[,c(4:100)],otu_tax_ref_vector$ASV_representant,na.rm=T)
otu_tax_ref_vector_agg$ASV <- rownames(otu_tax_ref_vector_agg)


otu_tax_ref_agg_bact <- merge(otu_tax_ref_vector_agg, tax_seq, by="ASV")
write.table(otu_tax_ref_agg_bact,"01_Tables/03_otu_tax_ref_agg_bacteria.tsv", sep="\t", quote=FALSE)

```

## 4. Phyloseq 16+18S+ITS <a name="ps_16+18S+ITS"></a>

```r
#!/bin/Rscript
otu_tax_ref <- read.table("01_Tables/04_otu_tax_ref_agg_bacteria_fungi_ssasco18S_norm.tsv", sep="\t", header=TRUE, row.names = 1)
sample <- read.table("01_Tables/05_sample.csv", sep=";",header=TRUE, row.names = 1)

otu <- otu_tax_ref %>% dplyr::select(L1I1_AOUT:L1R5_AVRIL) %>% as.matrix()
tax <- otu_tax_ref %>% dplyr::select(Kingdom:Species) %>% as.matrix()

ps <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(sample))

#otu_table()   OTU Table:         [ 12090 taxa and 96 samples ]
#sample_data() Sample Data:       [ 96 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 12090 taxa by 7 taxonomic ranks ]
```

## 5. Filtering steps <a name="Filtering_steps"></a>

```r
#!/bin/Rscript
#Filter OTUs that have 45 sequences in at least 3 samples
subset_3samples <- genefilter_sample(ps, filterfun_sample(function(x) x >= 45), A=3)
ps_filt <- prune_taxa(subset_3samples, ps)

#otu_table()   OTU Table:         [ 563 taxa and 96 samples ]
#sample_data() Sample Data:       [ 96 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 563 taxa by 7 taxonomic ranks ]
```

### SAVE

```r
#!/bin/Rscript
#saveRDS(ps, "00_Phyloseq_objects/01_ps_bact_fungi.rds")
#saveRDS(ps_filt, "00_Phyloseq_objects/02_ps_bact_fungi_45_filt.rds")

saveRDS(ps, "00_Phyloseq_objects/03_ps_bact_fungi_ss18Sasco.rds")
saveRDS(ps_filt, "00_Phyloseq_objects/04_ps_bact_fungi_45_filt_ss18Sasco.rds")

tax <- ps_filt@tax_table %>% as.data.frame()
tax <- cbind(rownames(tax), tax)
rownames(tax) <- NULL
colnames(tax)[1] <- "ASV"

otu <- ps_filt@otu_table %>% as.data.frame() 
otu <- cbind(rownames(otu), otu)
rownames(otu) <- NULL
colnames(otu)[1] <- "ASV"

otu_tax <- merge(otu, tax, by="ASV")
#write.table(otu_tax, "01_Tables/06_otu_tax_45_filt_ss18Sasco_norm.tsv", sep="\t", quote=FALSE)
#write.table(otu, "01_Tables/07_otu_45_filt_ss18Sasco_norm.tsv", sep="\t", quote=FALSE)
```

## 6. Covariance estimation using SPIEC-EASI <a name="SPIEC-EASI"></a>

Perfom SPIEC-EASI on cluster
```bash
#!/usr/bin/env bash
#SBATCH --job-name=spiec-easi
#SBATCH -p long
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --ntasks-per-node=1
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL

module load r

Rscript spiec-easi.R spiec.out 

```

Spiec-easi.R:
```r
#!/bin/Rscript
library(SpiecEasi, lib.loc = "/shared/software/miniconda/envs/r-4.2.3/lib/R/library")
library(phyloseq)
library(magrittr)

setwd("finalresult/03_NOVASEQ_METAB/05_SPIEC_EASI/")

ps_filt <- readRDS("04_ps_bact_fungi_45_filt_ss18Sasco.rds") # 45 reads in 3samples


ps_glasso <- spiec.easi(otu_table(ps_filt), method="glasso", lambda.min.ratio=0.01, nlambda=30, icov.select.params = list(rep.num=30, ncores=16))

save.image("spiec_easi_ps_ss18S_asco.RData")

```

```r
#!/bin/Rscript
load("00_Phyloseq_objects/spiec_easi_ps_ss18S_asco.RData")

ig2.mb <- adj2igraph(getRefit(ps_glasso),  vertex.attr=list(name=taxa_names(ps_filt)))
basic <- plot_network(ig2.mb, ps_filt, type='taxa', color="Phylum")
ggsave("02_Figures/00_basic.pdf", basic, height = 10, width = 15)


# Export significant covariance (Lemonnier et al)
secor <- cov2cor(forceSymmetric(getOptCov(ps_glasso), ifelse(sum(Matrix::tril(getOptCov(ps_glasso)))>sum(Matrix::triu(getOptCov(ps_glasso))), 'L', 'U')))
refit <- forceSymmetric(getRefit(ps_glasso), ifelse(sum(Matrix::tril(getRefit(ps_glasso)))>sum(Matrix::triu(getRefit(ps_glasso))), 'L', 'U'))
ps_glasso_net <- adj2igraph(as.matrix(secor*refit),  vertex.attr=list(name=taxa_names(ps_filt)))

E(ps_glasso_net)$weight %>% hist

graph.density(ps_glasso_net) #0.04037949
ps_glasso$est$sparsity[getOptInd(ps_glasso)]

g <- ps_glasso_net %>% as_tbl_graph()
E(g)$weight %>% hist(main="All weights", col="#009967", xlab = "edge weights")
summary(g) #531 nodes and 5682 egdes 
```

### SAVE

```r
#!/bin/Rscript
#write_graph(g, "03_Networks/01_spiec_easi_all.gml", format = "gml")

g <- read.graph(file = ("03_Networks/01_spiec_easi_all.gml"), format = c("gml"))
```

## 7. Positive edges <a name="Positive_edges"></a>

```r
#!/bin/Rscript
g_network_plus <- as_tbl_graph(g, directed=FALSE) %>% activate(edges) %>% filter(weight > 0)

E(g_network_plus)$weight %>% hist(main="Distribution of network weights After remove negative values", col="#009967", xlab = "edge weights")

edges_plus <- as_data_frame(g_network_plus, what="edges") 
vertice_plus <- as_data_frame(g_network_plus, what="vertices")

# keep large node
components <- igraph::clusters(g_network_plus, mode="weak")
print("Connected comp sizes:")
components$csize
#Keep largest component
biggest_comp_id <- which.max(components$csize)
vert_ids <- V(g_network_plus)[components$membership == biggest_comp_id]
g_network_plus <-igraph::induced_subgraph(g_network_plus, vert_ids)
```

### SAVE
```r
#!/bin/Rscript
write_graph(g_network_plus, "03_Networks/02_spiec_easi_plus.gml", format = "gml")
```

<https://github.com/AvantiShri/deterministic_louvain> louvain modif garde negative nodes


### Clustering algorithm comparaison <a name="Clustering"></a>

<https://github.com/hirotokaneko/plankton-from-satellite/blob/main/scripts/314_CommunityDetection.R>

```r
#!/bin/Rscript
g <- g_network_plus

print("Greedy algorithm")
imc <- cluster_fast_greedy(g)
mem_fast_greedy <- membership(imc)
print(modularity(g, membership(imc)))
#0.536742

print("Infomap")
m_max <- 0.0
for (i in 1:100) {
    set.seed(i)
    imc <- cluster_infomap(g)
    m_i <- modularity(g, membership(imc))
    if (m_max < m_i) {
        m_max <- m_i
        mem_infomap <- membership(imc)
    }
}
print(m_max)
#0.5561212


print("Label propagation")
m_max <- 0.0
for (i in 1:100) {
    set.seed(i)
    imc <- cluster_label_prop(g)
    m_i <- modularity(g, membership(imc))
    if (m_max < m_i) {
        m_max <- m_i
        mem_label_prop <- membership(imc)
    }
}
print(m_max)
#0.581


print("Eigenvector")
imc <- cluster_leading_eigen(g, options=list(maxiter=100000))
mem_leading_eigen <- membership(imc)
print(modularity(g, membership(imc)))
#0.5426426

print("Leiden algorithm")
m_max <- 0.0
for (i in 1:100) {
    set.seed(i)
    imc <- cluster_leiden(g, objective_function="modularity")
    m_i <- modularity(g, membership(imc))
    if (m_max < m_i) {
        m_max <- m_i
        mem_leiden <- membership(imc)
    }
}
print(m_max)
#0.6020787

print("Louvain algorithm")
imc <- cluster_louvain(g)
mem_louvain <- membership(imc)
print(modularity(g, membership(imc)))
#0.5956367

print("Spinglass algorithm")
m_max <- 0.0
for (i in 1:100) {
    set.seed(i)
    imc <- cluster_spinglass(g)
    m_i <- modularity(g, membership(imc))
    if (m_max < m_i) {
        m_max <- m_i
        mem_spinglass <- membership(imc)
    }
}
print(m_max)
#0.5950466

print("Random walk")
imc <- cluster_walktrap(g)
mem_walktrap <- membership(imc)
print(modularity(g, membership(imc)))
#0.5658616

DFmems <- as.data.frame(cbind(mem_fast_greedy,mem_infomap,mem_label_prop,mem_leading_eigen,
	mem_leiden,mem_louvain,mem_walktrap, mem_spinglass), stringsAsFactors=FALSE)

spin_leiden_louvain <- DFmems %>% dplyr::select(mem_leiden, mem_spinglass, mem_louvain) %>% as.data.frame()
#write.table(spin_leiden_louvain, "01_Tables/09_spinglass_leiden_louvain.csv", quote=FALSE, sep=";")
#write.table(DFmems, "01_Tables/08_all_algos_modularity.csv", quote=FALSE, sep=";")
```

### Leiden or Spinglass or Louvain ? <a name="Leidein_Springlass_Louvain"></a>

```r
#!/bin/Rscript
#NOT RUN
#Leiden, Louvain and Spinglass best modularity index
spin_leiden_lou <- read.table("01_Tables/16_spinglass_leiden_louvain_wt_roscoff_asv_number.csv", sep=";", header=TRUE)

spin <- spin_leiden_lou %>% filter(Algorithms=="mem_spinglass") %>% as.data.frame()
len <- spin_leiden_lou %>% filter(Algorithms=="mem_leiden") %>% as.data.frame()
lou <- spin_leiden_lou %>% filter(Algorithms=="mem_louvain") %>% as.data.frame()

spin <- ggplot(spin, aes(x=reorder(Module, -Number_ASVs), y=Number_ASVs)) + 
  facet_nested(~Algorithms, scale="free", space="free")+ geom_bar(stat="identity") +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank())+
   ylab("ASV count per modules")

len <- ggplot(len, aes(x=reorder(Module, -Number_ASVs), y=Number_ASVs)) + 
  facet_nested(~Algorithms, scale="free", space="free")+ geom_bar(stat="identity") +
  scale_y_continuous(limits=c(0,300))+
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank())+
   ylab("ASV count per modules")

lou <- ggplot(lou, aes(x=reorder(Module, -Number_ASVs), y=Number_ASVs)) + 
  facet_nested(~Algorithms, scale="free", space="free")+ geom_bar(stat="identity") +
  theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank())+
   ylab("ASV count per modules")

grid.arrange(spin, len, lou,nrow=1) 
```

### Leiden algorithm for modularity <a name="Leiden"></a>

```r
#!/bin/Rscript
print("Leiden algorithm")
m_max <- 0.0
for (i in 1:100) {
    set.seed(i)
    imc <- cluster_leiden(g_network_plus, objective_function="modularity")
    m_i <- modularity(g_network_plus, membership(imc))
    if (m_max < m_i) {
        m_max <- m_i
        mem_leiden <- membership(imc)
    }
}

print(m_max)
cat("Leiden identified", length(unique(imc$membership)), "modules.") #9modules
sizes(imc) 
#0.6020787
#1   2   3   4   5   6   7   8   
#78 106  34 122  53  79   4  44 

module.character <- as.character(imc$membership)
g_network_plus <- set_vertex_attr(g_network_plus, "module", index = V(g_network_plus), module.character)
```

```r
#!/bin/Rscript
edges_plus <- as_data_frame(g_network_plus, what="edges")
vertice_plus <- as_data_frame(g_network_plus, what="vertices")
summary(g_network_plus) #520 nodes & 4130 edges

otu_tax <- read.table("01_Tables/06_otu_tax_45_filt_ss18Sasco_norm.tsv", sep="\t", header=TRUE, row.names = 1)
otu <- otu_tax %>% dplyr::select(L1I1_AOUT:L1R5_AVRIL)
tax <- otu_tax %>% dplyr::select(Kingdom:Species)
names_asv <- vertice_plus$name 
tax <- tax[rownames(tax) %in% names_asv,]
otu <- otu[rownames(otu) %in% names_asv,]
tax %<>% arrange(factor(rownames(tax), levels = names_asv))
otu %<>% arrange(factor(rownames(otu), levels = names_asv))
otu %<>% mutate(Abundance=rowSums(across(where(is.numeric))))
#write.table(tax, "01_Tables/08_tax_g_network_plus.tsv", sep="\t", quote=FALSE)

V(g_network_plus)$Kingdom <- as.character(tax[,1])
V(g_network_plus)$Phylum <- as.character(tax[,2])
V(g_network_plus)$Class <- as.character(tax[,3])
V(g_network_plus)$Order <- as.character(tax[,4])
V(g_network_plus)$Family <- as.character(tax[,5])
V(g_network_plus)$Genus <- as.character(tax[,6])
V(g_network_plus)$Species <- as.character(tax[,7])
V(g_network_plus)$Abundance <- as.character(otu[,97])
```

### Observe network plot <a name="plot"></a>

```r
#!/bin/Rscript
g_network_iso <- which(degree(g_network_plus)==0)
g_network_iso_trim <- delete.vertices(g_network_plus, g_network_iso)
summary(g_network_iso_trim) # 520 nodes  4130 edges 

E(g_network_iso_trim)$weight %>% hist(main="Distribution of network weights After remove negative values", col="#009967", xlab = "edge weights")

plot(g_network_iso_trim, 
     layout = layout_(g_network_iso_trim, nicely()), 
     vertex.color=imc$membership, 
     vertex.label = NA,
     vertex.size= 5,
     edge.arrow.size= 5,
     edge.width=E(g_network_iso_trim)$weight) 
```

### SAVE

```r
#!/bin/Rscript
#write.graph(g_network_iso_trim, "03_Networks/02_spiec_easi_plus.gml", format = "gml")
#g_network_iso_trim <- read.graph(file = ("03_Networks/02_spiec_easi_plus.gml"), format = c("gml"))

edges_path <- here("~/Documents/ownCloud/11_Rstudio/08_Novaseq/07_Coocurrence_networks/03_Networks/03_edges_plus_spiec_easi.csv")
nodes_path <- here("~/Documents/ownCloud/11_Rstudio/08_Novaseq/07_Coocurrence_networks/03_Networks/04_nodes_plus_spiec_easi.csv")
gephi_write_both(g_network_iso_trim, edges_path, nodes_path, na = "", verbose = TRUE)
```

## 8. Negative edges <a name="Negative_edges"></a>

```r
#!/bin/Rscript
g_network_neg <- as_tbl_graph(g, directed=FALSE) %>% activate(edges) %>% filter(weight < 0) #1594 ASVs

E(g_network_neg)$weight %>% hist(main="Distribution of network weights After remove positive values", col="#009967", xlab = "edge weights")

edges_neg <- as_data_frame(g_network_neg, what="edges") #1594 edges
vertice_neg <- as_data_frame(g_network_neg, what="vertices") 
```

```r
#!/bin/Rscript
plot(g_network_neg, 
     layout = layout_(g_network_neg, nicely()), 
     vertex.label = NA,
     vertex.size= 5,
     edge.arrow.size= 5,
     edge.width=E(g_network_neg)$weight) 
```

### SAVE

```r
#!/bin/Rscript
otu_tax <- read.table("01_Tables/14_otu_tax_filt_3samples45reads.tsv", sep="\t", header=TRUE)
tax <- otu_tax %>% dplyr::select(Kingdom:Species)

names_asv <- vertice_neg$name 
tax <- tax[rownames(tax) %in% names_asv,]

V(g_network_neg)$Kingdom <- as.character(tax[,1])
V(g_network_neg)$Phylum <- as.character(tax[,2])
V(g_network_neg)$Class <- as.character(tax[,3])
V(g_network_neg)$Order <- as.character(tax[,4])
V(g_network_neg)$Family <- as.character(tax[,5])
V(g_network_neg)$Genus <- as.character(tax[,6])
V(g_network_neg)$Species <- as.character(tax[,7])

write.graph(g_network_neg, "03_Networks/06_spiec_easi_neg.gml", format = "gml")


## Table for GePhi
edges_path <- here("~/Documents/ownCloud/11_Rstudio/08_Novaseq/07_Coocurrence_networks/03_Networks/06_edges_neg_spiec_easi.csv")
nodes_path <- here("~/Documents/ownCloud/11_Rstudio/08_Novaseq/07_Coocurrence_networks/03_Networks/06_nodes_neg_spiec_easi.csv")
gephi_write_both(g_network_neg, edges_path, nodes_path, na = "", verbose = TRUE)
```


## 9. Modules characteristics <a name="Modules_characteristics"></a>
### Add metadata to OTU and ps_gephi <a name="metadata"></a>

```r
#!/bin/Rscript
nodes <- read.csv("01_Tables/08_nodes_modules_plus_easi.csv", sep=";", header=TRUE, dec=".")
otu <- read.csv("01_Tables/07_otu_45_filt_ss18Sasco_norm_for_module.csv", sep=";", header=TRUE, dec=",", row.names = 1)
otu %<>% as.matrix
tax <- nodes %>% select(Id, module,kingdom:species) %>% tibble::column_to_rownames("Id") %>% as.matrix

ps_net_plus <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tax), sample_data(ps_filt@sam_data))

taxa <- psmelt(ps_net_plus)
#write.table(taxa,"01_Tables/20_taxa_network_plus.tsv", sep="\t", quote=FALSE)

```

### Module palette <a name="palettes"></a>
```r
#!/bin/Rscript
color_modules <- c(
  "Module_1"="#3580C1",
  "Module_2"="#D966E7",
  "Module_3"="#35AB00",
  "Module_4"="#0D9165",
  "Module_5"="#9D7800",
  "Module_6"="#FD521D",
  "Module_7"="#D34971",
  "Module_8"="#FFFC2E")
```

### Season <a name="plot_seasons"></a>
```r
#!/bin/Rscript
taxa$Dataset_types[taxa$Month == "November21" & taxa$Type=="Site_1"] <- "Season"

mod1 <- taxa %>% filter(module=="Module_1")
#mod1 %<>% mutate(rel = Abundance / sum(Abundance))
mod1 %<>% filter(Dataset_types=="Season")
mod1 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod2 <- taxa %>% filter(module=="Module_2")
#mod2 %<>% mutate(rel = Abundance / sum(Abundance))
mod2 %<>% filter(Dataset_types=="Season")
mod2 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod3 <- taxa %>% filter(module=="Module_3")
#mod3 %<>% mutate(rel = Abundance / sum(Abundance))
mod3 %<>% filter(Dataset_types=="Season")
mod3 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod4 <- taxa %>% filter(module=="Module_4")
#mod4 %<>% mutate(rel = Abundance / sum(Abundance))
mod4 %<>% filter(Dataset_types=="Season")
mod4 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod5 <- taxa %>% filter(module=="Module_5")
#mod5 %<>% mutate(rel = Abundance / sum(Abundance))
mod5 %<>% filter(Dataset_types=="Season")
mod5 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod6 <- taxa %>% filter(module=="Module_6")
#mod6 %<>% mutate(rel = Abundance / sum(Abundance))
mod6 %<>% filter(Dataset_types=="Season")
mod6 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod7 <- taxa %>% filter(module=="Module_7")
#mod7 %<>% mutate(rel = Abundance / sum(Abundance))
mod7 %<>% filter(Dataset_types=="Season")
mod7 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

mod8 <- taxa %>% filter(module=="Module_8")
#mod8 %<>% mutate(rel = Abundance / sum(Abundance))
mod8 %<>% filter(Dataset_types=="Season")
mod8 %<>% group_by(Month, module) %>% summarise(Sum=sum(Abundance))

## mod 1 

mod1$Month <- factor(mod1$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod1$Month_numeric <- as.numeric(mod1$Month)

mod1 <- ggplot() +
  geom_area(data = subset(mod1, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod1, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod1, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod1, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod1, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod1, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod1$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 1")

##mod2
mod2$Month <- factor(mod2$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod2$Month_numeric <- as.numeric(mod2$Month)

mod2 <- ggplot() +
  geom_area(data = subset(mod2, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod2, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod2, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod2, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod2, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod2, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod2$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 2")

## Mod3
mod3$Month <- factor(mod3$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod3$Month_numeric <- as.numeric(mod3$Month)

mod3 <- ggplot() +
  geom_area(data = subset(mod3, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod3, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod3, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod3, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod3, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod3, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod3$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 3")

## Mod4 
mod4$Month <- factor(mod4$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod4$Month_numeric <- as.numeric(mod4$Month)

mod4 <- ggplot() +
  geom_area(data = subset(mod4, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod4, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod4, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod4, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod4, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod4, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod4$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 4")

## Mod 5
mod5$Month <- factor(mod5$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod5$Month_numeric <- as.numeric(mod5$Month)

mod5 <- ggplot() +
  geom_area(data = subset(mod5, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod5, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod5, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod5, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod5, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod5, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod5$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 5")


## Mod 6 
mod6$Month <- factor(mod6$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod6$Month_numeric <- as.numeric(mod6$Month)

mod6 <- ggplot() +
  geom_area(data = subset(mod6, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod6, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod6, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod6, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod6, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod6, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod6$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 6")


## Mod 7
mod7$Month <- factor(mod7$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod7$Month_numeric <- as.numeric(mod7$Month)

mod7 <- ggplot() +
  geom_area(data = subset(mod7, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod7, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod7, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod7, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod7, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod7, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod7$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 7")

 
## mod 8 

mod8$Month <- factor(mod8$Month, levels = c("November21", "December", "January", "February", "March", 
                                             "April", "May", "June", 
                                             "July", "August", "September", 
                                             "October", "November22"))
 
mod8$Month_numeric <- as.numeric(mod8$Month)

mod8 <- ggplot() +
  geom_area(data = subset(mod8, Month_numeric <= 3), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_area(data = subset(mod8, Month_numeric >= 3 & Month_numeric <=5), aes(x = Month_numeric, y = Sum), fill="#6BAED6", alpha=0.5) +
  geom_area(data = subset(mod8, Month_numeric >= 5 & Month_numeric <=8), aes(x = Month_numeric, y = Sum), fill="#238B45", alpha=0.5) +
  geom_area(data = subset(mod8, Month_numeric >= 8 & Month_numeric <=12), aes(x = Month_numeric, y = Sum), fill="#f2b613", alpha=0.5) +
  geom_area(data = subset(mod8, Month_numeric >= 12 & Month_numeric <=13), aes(x = Month_numeric, y = Sum), fill="#963414", alpha=0.5) +
  geom_line(data = mod8, aes(x = Month_numeric, y = Sum), color="black") +
  scale_x_continuous(breaks = 1:13, labels = levels(mod8$Month)) +  # Mettre des labels des mois
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    legend.title = element_text(size=9, face="bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines")
  ) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Module 8")

```

```r
#!/bin/Rscript
modules_ab <- grid.arrange(mod1, mod2, mod3, mod5, mod8, mod4, mod6, mod7, ncol=1)

ggsave("02_Figures/09_modules_abundances.pdf",modules_ab, width = 5, height = 10 )
```

### Parts <a name="plot_parts"></a>
```r
#!/bin/Rscript
p_mod1 <- taxa %>% filter(module=="Module_1")
#p_mod1 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod1 %<>% filter(Dataset_types=="Algal_parts")
p_mod1 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod2 <- taxa %>% filter(module=="Module_2")
#p_mod2 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod2 %<>% filter(Dataset_types=="Algal_parts")
p_mod2 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod3 <- taxa %>% filter(module=="Module_3")
#p_mod3 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod3 %<>% filter(Dataset_types=="Algal_parts")
p_mod3 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod4 <- taxa %>% filter(module=="Module_4")
#p_mod4 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod4 %<>% filter(Dataset_types=="Algal_parts")
p_mod4 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod5 <- taxa %>% filter(module=="Module_5")
#p_mod5 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod5 %<>% filter(Dataset_types=="Algal_parts")
p_mod5 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod6 <- taxa %>% filter(module=="Module_6")
#p_mod6 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod6 %<>% filter(Dataset_types=="Algal_parts")
p_mod6 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod7 <- taxa %>% filter(module=="Module_7")
#p_mod7 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod7 %<>% filter(Dataset_types=="Algal_parts")
p_mod7 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))

p_mod8 <- taxa %>% filter(module=="Module_8")
#p_mod8 %<>% mutate(rel = Abundance / sum(Abundance))
p_mod8 %<>% filter(Dataset_types=="Algal_parts")
p_mod8 %<>% group_by(Type, module) %>% summarise(Sum=sum(Abundance))


palette_parts <- c("Receptacle"="#225EA8", 
                   "Apex"="#D67500",
                   "Medium"="#578B21", 
                   "Base"="#B696D6")

## mod 1 

p_mod1$Type <- factor(p_mod1$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod1 <- ggplot(p_mod1, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    axis.text.x = element_blank(),
    panel.background = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 1")


## mod 2 

p_mod2$Type <- factor(p_mod2$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod2 <- ggplot(p_mod2, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 2")

## mod3 


p_mod3$Type <- factor(p_mod3$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod3 <- ggplot(p_mod3, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 3")

## mod4 


p_mod4$Type <- factor(p_mod4$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod4 <- ggplot(p_mod4, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 4")

## mod 5


p_mod5$Type <- factor(p_mod5$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod5 <- ggplot(p_mod5, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 5")

## mod 6 

p_mod6$Type <- factor(p_mod6$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod6 <- ggplot(p_mod6, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
     axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 6")

## mod 7 


p_mod7$Type <- factor(p_mod7$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod7 <- ggplot(p_mod7, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 7")

## mod 8 


p_mod8$Type <- factor(p_mod8$Type, levels = c("Receptacle", "Apex","Medium","Base")) 

p_mod8 <- ggplot(p_mod8, aes(x=Type, y=Sum, group=1, color=Type)) + geom_point(size=3)+ geom_line(color="black", alpha=0.9)+
  scale_color_manual(values=palette_parts)+
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.05, "lines"),
    legend.position = "none"
  ) +
  xlab(NULL)+
  ylab(NULL)+
  ggtitle("Module 8")
```

```r
#!/bin/Rscript
modules_part <- grid.arrange(p_mod1, p_mod2,p_mod3,p_mod5 , p_mod8, p_mod4, p_mod6, p_mod7, ncol=1)
ggsave("02_Figures/10_modules_abundances_parts.pdf",modules_ab, width = 5, height = 10 )
```


```r
#!/bin/Rscript
taxa3 <- rbind(mod1,mod2, mod3, mod4, mod5, mod6, mod7, mod8)

colors1 <- brewer.pal(11, "RdYlBu")

taxa3$Month <- factor(taxa3$Month, levels = c("November21","December","January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November22")) 

ggplot(taxa3, aes(x=Month,y=as.factor(module),fill=Sum)) + geom_tile() + 
  scale_fill_gradientn(colours =rev(colors1),
                       name="Abundance",
                        values=scales::rescale(c(0, 0.10,0.30, 0.4)), 
                       limits=c(0,0.45)) +
  theme(axis.text.x = element_text(size=8, angle=90,hjust=0.95,vjust=0.5), 
        axis.text.y = element_text(size=8), 
        panel.background = element_blank(), 
        legend.title = element_text(size=9, face = "bold"),
        strip.text.y = element_blank(),
        panel.spacing.x  = unit(0.1, "lines"),
        panel.spacing.y  = unit(0.05, "lines")) +
  xlab(NULL) + ylab(NULL)
```

### Modules compositions <a name="module_composition"></a>

```r
#!/bin/Rscript
otu_tax <- read.table("01_Tables/08_nodes_modules_plus_easi.csv", sep=";", row.names = 1, header = TRUE)
taxa2 <- otu_tax %>% group_by(module) %>% 
  mutate(pct = abundance/sum(abundance)*100) # pct for the pie - adds to 100

c_mod1 <- taxa2 %>% filter(module=="Module_1")
c_mod1 %<>% mutate(rel = abundance / sum(abundance))
c_mod1$class[c_mod1$rel < 0.005] <- "abund. < 0.5%"

c_mod2 <- taxa2 %>% filter(module=="Module_2")
c_mod2 %<>% mutate(rel = abundance / sum(abundance))
c_mod2$class[c_mod2$rel < 0.005] <- "abund. < 0.5%"

c_mod3 <- taxa2 %>% filter(module=="Module_3")
c_mod3 %<>% mutate(rel = abundance / sum(abundance))
c_mod3$class[c_mod3$rel < 0.005] <- "abund. < 0.5%"

c_mod4 <- taxa2 %>% filter(module=="Module_4")
c_mod4 %<>% mutate(rel = abundance / sum(abundance))
c_mod4$class[c_mod4$rel < 0.005] <- "abund. < 0.5%"

c_mod5 <- taxa2 %>% filter(module=="Module_5")
c_mod5 %<>% mutate(rel = abundance / sum(abundance))
c_mod5$class[c_mod5$rel < 0.005] <- "abund. < 0.5%"

c_mod6 <- taxa2 %>% filter(module=="Module_6")
c_mod6 %<>% mutate(rel = abundance / sum(abundance))
c_mod6$class[c_mod6$rel < 0.005] <- "abund. < 0.5%"

c_mod7 <- taxa2 %>% filter(module=="Module_7")
c_mod7 %<>% mutate(rel = abundance / sum(abundance))
c_mod7$class[c_mod7$rel < 0.005] <- "abund. < 0.5%"

c_mod8 <- taxa2 %>% filter(module=="Module_8")
c_mod8 %<>% mutate(rel = abundance / sum(abundance))
c_mod8$class[c_mod8$rel < 0.005] <- "abund. < 0.5%"

taxa3 <- rbind(c_mod1,c_mod2, c_mod3, c_mod4, c_mod5, c_mod6, c_mod7, c_mod8)

### abund < 0.05
palette_class <- c("Thermoanaerobaculia"= "#FE698E",
                   "NB1-j;NA"="#E4CDE2",
                   "Polyangia"="#FE7125",
                   "Verrucomicrobiae"="#EAAF29",
                   "Acidimicrobiia"="black", 
                   "Bacteroidia"="#853250",
                   "OM190"="#9EBCDA",
                   "Phycisphaerae"="#FFADC1",
                   "Planctomycetes"="#645A9F",
                   "Cyanobacteriia"="#00441B",
                   "Gammaproteobacteria"="#B03F00",
                   "Alphaproteobacteria"="#08306B", 
                   "Saccharomycetes"="#466147",
                   "Sordariomycetes"="#74C476",
                   "Dothideomycetes"="#238B45",
                   "Fungi;NA;NA"="#ADF1B0",
                   "Arthropoda"="#FEEDDE",
                   "Phaeophyceae"="#225EA8",
                   "abund. < 0.5%"="grey")



taxa3$class <- factor(taxa3$class, levels = c("Thermoanaerobaculia",
                                              "NB1-j;NA",
                                              "Polyangia",
                                              "Verrucomicrobiae",
                                             "Acidimicrobiia",
                                              "Bacteroidia",
                                              "OM190","Phycisphaerae","Planctomycetes",
                                              "Cyanobacteriia",
                                              "Gammaproteobacteria",
                                             "Alphaproteobacteria", 
                                              "Fungi;NA;NA","Saccharomycetes","Sordariomycetes","Dothideomycetes",
                                             "Arthropoda",
                                              "Phaeophyceae",
                                             "abund. < 0.5%"))

## mod 1 
c_mod1 <- ggplot(c_mod1, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

## mod 2
c_mod2 <- ggplot(c_mod2, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

#mod 3
c_mod3 <- ggplot(c_mod3, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))


# mod 4
c_mod4 <- ggplot(c_mod4, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

#mod 5
c_mod5 <- ggplot(c_mod5, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

# mod 6
c_mod6 <- ggplot(c_mod6, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

#mod 7
c_mod7 <- ggplot(c_mod7, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

#mod 8
c_mod8 <- ggplot(c_mod8, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(~module)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_text(size=30), 
        legend.position = "bottom")+
   guides(fill=guide_legend(ncol = 5))

# all 

taxa3$module <- factor(taxa3$module, levels = c( "Module_1", "Module_2", "Module_3", "Module_5", "Module_8", "Module_4", "Module_6", "Module_7"))

taxo <- ggplot(taxa3, aes(x="", y=rel, fill=class, color=class)) +
  facet_grid(module~.)+
     scale_color_manual(values = palette_class ) + scale_fill_manual(values = palette_class) +
  geom_bar(stat="identity") +
  coord_polar("y", start=0) +theme_void() +
  theme(legend.text = element_text(size=20), 
        strip.text=element_blank(), 
        legend.position = "none")+
   guides(fill=guide_legend(ncol = 5))# + 
  #pdf("02_Figures/05_module_compositions2.pdf", height = 10, width = 30)

```

### Combined figures <a name="combined_figures"></a>
```r
#!/bin/Rscript
taxo + modules_ab + modules_part
ggsave("02_Figures/11_module_taxo_season_parts.pdf", height = 10, width = 8)
```

![Figure 5 | Co-occurence network ](https://github.com/rssco/novaseq_ascophyllum/blob/main/01_Figures/Figure_5_Combined_network_degree_taxo_mod.png)<!-- -->