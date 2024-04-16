#BiocManager::install("microbiome")
#install.packages("hrbrthemes")
#install.packages("gcookbook")
#install.packages('rgl')
#BiocManager::install('DESeq2')

#import library ----
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microbiome)
library(vegan)
library(MASS)
library(agricolae)
library(DESeq2)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
#import data from RDS ----
phylo<-readRDS("./phylo_ASV.rds")
phylo.keep<-readRDS("./phylo_keep.RDS")
phylo.keep.norm<-readRDS('./phylo_keep_norm.RDS')
#set theme color and order ----
theme_color <- c("#0066cc","#ffbc0d","#FF69B4")
desired_order <- c("HC", "KS7", "KS9", "SS", "HBS", "TN3", "TN5", "TN11")
set.seed(314)
#import data ----
raw.data<-rio::import("./ASV_table_final.xlsx", which = 1)
row.names(raw.data)<-raw.data[,1]
sample.meta<-rio::import("./ASV_table_final.xlsx", which = 2)
row.names(sample.meta)<-sample.meta[,1]
sample.meta$Variety<-factor(sample.meta$Variety,levels=desired_order)
otu.table<-raw.data[,9:48]
tax.table<-raw.data[,2:8]
tree.file<-read_tree("./ASV_MAFFT.treefile")

#remove sample less than 2000 reads ----
sample.count<-colSums(otu.table)
sample.count<-sample.count>2000
otu.table.prun <- otu.table[,sample.count]
sample.meta.prun <- sample.meta[sample.count,]

#import as phyloseq file ----
OTU <- otu_table(as.matrix(otu.table.prun), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax.table))
samples <- sample_data(sample.meta.prun)
phylo <- phyloseq(OTU, TAX, samples, tree.file)

saveRDS(phylo,"./phylo_ASV.rds")
phylo
#otu_table()   OTU Table:         [ 560 taxa and 36 samples ]
#sample_data() Sample Data:       [ 36 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 560 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 560 tips and 554 internal nodes ]

#filtering OTUs with mean RA >= 0.001 and occupancy >= 5% ----

abund_occ <- function(otu.table) {
  otu_PA <- 1*((otu.table>0)==1)  # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA) # Occupancy
  otu_rel <- apply(decostand(otu.table, method="total", MARGIN=2),1, mean)     # mean relative abundance
  occ_abun <- rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),'otu') # combining occupancy and abundance data frame
  return(occ_abun)
}
phylo.keep <- phylo %>%
  otu_table() %>% # selecting the OTU table 
  as("matrix") %>% # turning the OTU table into a matrix for the function
  abund_occ() %>% # function defined above
  subset(otu_occ >= 0.05) %>% # filtering
  dplyr::select(otu) %>%
  as.matrix() %>%
  as.character()

phylo.keep <- prune_taxa(phylo.keep, phylo) # remove the taxa below 0.0001 and less than 5% occupancy
phylo.keep
saveRDS(phylo.keep, './phylo_keep.RDS')
#phyloseq-class experiment-level object
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 145 taxa and 36 samples ]
#sample_data() Sample Data:       [ 36 samples by 6 sample variables ]
#tax_table()   Taxonomy Table:    [ 145 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 145 tips and 144 internal nodes ]

#normalize by median of sequencing depth ----
total <- median(sample_sums(phylo.keep))
standf <- function(x, t=total) round(t * (x / sum(x)))
phylo.keep.norm <- transform_sample_counts(phylo.keep, standf)
saveRDS(phylo.keep.norm,'./phylo_keep_norm.RDS')

#Taxonomy properties ----
#phylum
tax.aggre<-tax_glom(phylo.keep.norm, "phylum")
otu.tax.aggre<-as.data.frame(otu_table(tax.aggre))
otu.tax.aggre<-rowSums(otu.tax.aggre)
tax.tax.aggre<-tax_table(tax.aggre)
pro.phylum.all<-cbind(tax.tax.aggre,otu.tax.aggre/sum(otu.tax.aggre)*100)
write.csv(pro.phylum.all,'./result/pro.phylum.all.csv')

#family
tax.aggre<-tax_glom(phylo.keep.norm, "family")
otu.tax.aggre<-as.data.frame(otu_table(tax.aggre))
otu.tax.aggre<-rowSums(otu.tax.aggre)
tax.tax.aggre<-tax_table(tax.aggre)
pro.family.all<-as.data.frame(cbind(tax.tax.aggre,otu.tax.aggre/sum(otu.tax.aggre)*100))
pro.family.all$V8<-as.numeric(pro.family.all$V8)
pro.family.all<-pro.family.all[order(pro.family.all$V8,decreasing=TRUE),]
write.csv(pro.family.all,'./result/pro.family.all.csv')

#genus
tax.aggre<-tax_glom(phylo.keep.norm, "genus")
otu.tax.aggre<-as.data.frame(otu_table(tax.aggre))
otu.tax.aggre<-rowSums(otu.tax.aggre)
tax.tax.aggre<-tax_table(tax.aggre)
pro.genus.all<-as.data.frame(cbind(tax.tax.aggre,otu.tax.aggre/sum(otu.tax.aggre)*100))
pro.genus.all$V8<-as.numeric(pro.genus.all$V8)
#pro.genus.all[,"V8">0.01]
pro.genus.all<-pro.genus.all[order(pro.genus.all$V8,decreasing=TRUE),]
write.csv(pro.genus.all,'./result/pro.genus.all.csv')

#species
tax.aggre<-tax_glom(phylo.keep.norm, "species")
otu.tax.aggre<-as.data.frame(otu_table(tax.aggre))
otu.tax.aggre<-rowSums(otu.tax.aggre)
tax.tax.aggre<-tax_table(tax.aggre)
pro.species.all<-as.data.frame(cbind(tax.tax.aggre,otu.tax.aggre/sum(otu.tax.aggre)*100))
pro.species.all$V8<-as.numeric(pro.species.all$V8)
#pro.species.all[,"V8">0.01]
pro.species.all<-pro.species.all[order(pro.species.all$V8,decreasing=TRUE),]
write.csv(pro.species.all,'./result/pro.species.all.csv')

#Alpha diversity ----
phylo.r<-rarefy_even_depth(phylo, sample.size = min(sample_sums(phylo)),rngseed=314)
adiv<-estimate_richness(phylo.r, split = TRUE) %>%
  dplyr::select(Observed,Shannon,InvSimpson)
adiv$Phenotype<-phylo.r@sam_data$Phenotype
adiv$Faith<- picante::pd(t(phylo.r@otu_table),phylo.r@phy_tree,include.root = F)[,1]
adiv$Pielou<-as.vector(evenness(phylo.r, index = "Pielou"))$pielou
adiv$Variety <-phylo.r@sam_data$Variety
adiv$Variety<-factor(adiv$Variety,levels=desired_order)
colnames(adiv)<-c("Richness","Shannon","InviSimpson","Phenotype","FaithPD", "Pielou", "Variety")
adiv.long<-melt(adiv,id = c("Phenotype","Variety"))

#ANOVA Resistant vs suceptible
adiv.S<-subset(adiv,Phenotype=='S')
adiv.R<-subset(adiv,Phenotype=='R')
with(adiv.S, shapiro.test(Richness))#0.6598
with(adiv.R, shapiro.test(Richness))#0.01508
wilcox.test(adiv.S$Richness,adiv.R$Richness)#0.1807

with(adiv.S, shapiro.test(Shannon))#0.5316
with(adiv.R, shapiro.test(Shannon))#0.1567
var.test(adiv.S$Shannon,adiv.R$Shannon)#0.7952
t.test(adiv.S$Shannon,adiv.R$Shannon,var.equal=TRUE) #0.082

with(adiv.S, shapiro.test(InviSimpson))#0.3978
with(adiv.R, shapiro.test(InviSimpson))#0.07353
var.test(adiv.S$InviSimpson,adiv.R$InviSimpson)#0.2514
t.test(adiv.S$InviSimpson,adiv.R$InviSimpson,var.equal=TRUE) #0.092

with(adiv.S, shapiro.test(FaithPD))#0.079
with(adiv.R, shapiro.test(FaithPD))#0.2988
var.test(adiv.S$FaithPD,adiv.R$FaithPD)#0.005294
t.test(adiv.S$Faith,adiv.R$Faith,var.equal=FALSE) #0.3055

with(adiv.S, shapiro.test(Pielou))#0.01204
with(adiv.R, shapiro.test(Pielou))#0.6261
wilcox.test(adiv.S$Pielou,adiv.R$Pielou)#0.168

##alpha diversity plot----
###RvsS
p <- ggplot(adiv.long, aes(x = Phenotype, y = value, fill=Phenotype)) + 
  geom_boxplot()+
  #geom_boxplot(aes(col=Phenotype)) +
  geom_jitter(alpha=0.2,size=0.5) +
  facet_wrap(~variable, scales = "free", nrow=1)+
  labs(x="Varieties",y="",fill="Phenotypes")+
  scale_fill_manual(values = theme_color)+
  #scale_color_manual(values = theme_color)+
  theme_bw()
p
ggsave("./Fig/alphadiv.png", height = 10, width = 18, unit='cm', dpi=600)
###Variety
p <- ggplot(adiv.long, aes(x = Variety, y = value, fill=Phenotype)) + 
  geom_boxplot() + 
  geom_jitter(alpha=0.2,size=0.5) +
  labs(x="Varieties",y="",fill="Phenotypes")+
  facet_wrap(~variable, scales = "free", nrow=1)+
  scale_fill_manual(values = theme_color)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("./Fig/alphadiv_var.png", height = 10, width = 22, unit='cm', dpi=600)

##Compare alpha-diveristy between cultivars
###observed
md <- aov(Richness~Variety,adiv)
shapiro.test(md$residuals) #0.305
car::leveneTest(adiv$Richness,adiv$Variety)#0.5145
summary(md) #0.000878
HSD.test(md,"Variety", console = TRUE)

#Richness groups
#HC   39.40000      a
#TN5  38.00000     ab
#SS   31.60000    abc
#KS9  21.66667    abc
#HBS  17.00000     bc
#TN11 16.80000      c
#TN3  15.80000      c
#KS7  13.00000      c

###Shannon
md <- aov(Shannon~Variety,adiv)
shapiro.test(md$residuals) #0.07
car::leveneTest(adiv$Shannon,adiv$Variety)#0.74
summary(md) #0.00412 **
HSD.test(md,"Variety", console = TRUE)

#Shannon groups
#KS9  2.429576      a
#SS   2.401701      a
#HC   2.338539      a
#TN5  2.277630      a
#TN11 2.118690     ab
#TN3  1.745505     ab
#KS7  1.427170     ab
#HBS  1.248663      b

###Inverse-simpson
md <- aov(InviSimpson~Variety,adiv)
shapiro.test(md$residuals) #0.044
md <- aov(log(InviSimpson)~Variety,adiv)
shapiro.test(md$residuals) #0.1069
car::leveneTest(log(adiv$InviSimpson),adiv$Variety)#0.5034
summary(md) #0.0273
HSD.test(md,"Variety", console = TRUE)
#log(InviSimpson) groups
#KS9         2.1319984      a
#SS          1.9598850      a
#TN11        1.9438676      a
#HC          1.8344557      a
#TN5         1.6231130      a
#TN3         1.4486082      a
#KS7         1.1361781      a
#HBS         0.9350226      a

###Faith
md <- aov(FaithPD~Variety,adiv)
shapiro.test(md$residuals) #0.8034
car::leveneTest(adiv$FaithPD,adiv$Variety)#0.5809
summary(md) #0.229

###Pielou
kruskal.test(Pielou~Variety,adiv)#0.016
md <- aov(Pielou~Variety,adiv)
shapiro.test(md$residuals) #0.0902
car::leveneTest(adiv$Pielou,adiv$Variety)#0.6413
summary(md) #0.00196
HSD.test(md,"Variety", console = TRUE)
#KS9  0.7928908      a
#TN11 0.7690404      a
#SS   0.7163800      a
#HC   0.6485489     ab
#TN3  0.6464401     ab
#TN5  0.6333343     ab
#KS7  0.5721562     ab
#HBS  0.4436598      b

#Beta diversity ----
#NMDS k=2
phylo.NMDS <- ordinate(phylo.keep.norm, "NMDS", "bray",k=2)
stress<-phylo.NMDS$stress #0.185
p<-plot_ordination(phylo.keep.norm, phylo.NMDS, type="samples", color="Phenotype", title=paste("NMDS (ASV)       ","stress = ", round(stress,3)),
                   shape='Variety') 
p + geom_point(size=3) + 
  scale_shape_manual(values = c(15,2,6,1,23,16,17,3)) +
  scale_color_manual(values = theme_color)+
  theme_bw()

ggsave("./Fig/NMDSplot_filter_norm.png", width=13, height=10, unit="cm",dpi=600)



##weighted unifrac ----
#conflict with vegan package
wunifrac_dist <- phyloseq::distance(phylo.keep.norm, method="unifrac", weighted=T)
phylo.PCoA.unifrac <- ordinate(phylo.keep.norm, method="PCoA", distance=wunifrac_dist)
p<-plot_ordination(phylo, phylo.PCoA.unifrac, color="Phenotype",shape="Variety",title="PCoA Unifrac")
p + geom_point(size=3) + 
  scale_shape_manual(values = c(1:4,7:10)) +
  scale_color_manual(values = theme_color)+
  theme_bw()
#looks like there are artifect

#eigenvalue
source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR1/evplot.R')
ev<- phylo.PCoA.unifrac$values$Eigenvalues
png("./ASV/PcoA_uni_broken.png", width=13, height=18, unit="cm",res=300)
evplot (ev)
dev.off()
PCoA.unifrac.scores<-phylo.PCoA.unifrac$vectors[,1:3]

##DCA ----
phylo.CA <- ordinate(phylo.keep.norm, "DCA")
#DCA1 length= 1.03 <3 linear model better
p<-plot_ordination(phylo.keep.norm, phylo.CA, type="samples", color="Phenotype", title="DCA (ASV)",
                   shape='Variety') 
p + geom_point(size=3) + 
  scale_shape_manual(values = c(1:4,7:10)) +
  scale_color_manual(values = theme_color)+
  theme_bw()
#not good

##jaccard ----
jc_dist <- phyloseq::distance(phylo.keep.norm, method="jaccard")
###PcOA
phylo.PCoA.jc <- ordinate(phylo.keep.norm, method="PCoA", distance=jc_dist)
p<-plot_ordination(phylo.keep.norm, phylo.PCoA.jc, color="Phenotype",shape="Variety",title="PCoA Unifrac")
p + geom_point(size=3) + 
  scale_shape_manual(values = c(1:4,7:10)) +
  scale_color_manual(values = theme_color)+
  theme_bw()
#looks like there are artifect

#NMDS
phylo.jc.NMDS <- ordinate(phylo.keep.norm, method="NMDS", distance=jc_dist ,k=2)
stress<-phylo.jc.NMDS$stress #0.158
p<-plot_ordination(phylo.keep.norm, phylo.jc.NMDS, type="samples", color="Phenotype", title=paste("NMDS (ASV)       ","stress = ", round(stress,3)),
                   shape='Variety') 
p + geom_point(size=3) + 
  scale_shape_manual(values = c(15,2,6,1,23,16,17,3)) +
  scale_color_manual(values = theme_color)+
  theme_bw()

#PERMANOVA ---- 

adonis2(t(otu_table(phylo.keep.norm))~phylo.keep.norm@sam_data$Phenotype+phylo.keep.norm@sam_data$Variety)
#                                   Df SumOfSqs      R2      F Pr(>F)    
#phylo.keep.norm@sam_data$Phenotype  1   1.0509 0.07123 4.3053  0.001 ***
#phylo.keep.norm@sam_data$Variety    6   6.8673 0.46548 4.6888  0.001 ***
#Residual                           28   6.8349 0.46328                  
#Total                              35  14.7531 1.00000       

adonis2(t(otu_table(phylo.keep.norm))~phylo.keep.norm@sam_data$Phenotype+phylo.keep.norm@sam_data$Site+phylo.keep.norm@sam_data$Variety)
#                                   Df SumOfSqs      R2      F Pr(>F)    
#phylo.keep.norm@sam_data$Phenotype  1   1.0509 0.07123 4.3053  0.002 ** 
#phylo.keep.norm@sam_data$Site       4   5.8597 0.39718 6.0012  0.001 ***
#phylo.keep.norm@sam_data$Variety    2   1.0076 0.06830 2.0639  0.003 ** 
#Residual                           28   6.8349 0.46328                  
#Total                              35  14.7531 1.00000      

#Relative abundance plot of each sample----
##color----
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
##transform to composition
phylo.filt<-microbiome::transform(phylo.keep.norm, "compositional")
phylo.filt <- aggregate_rare(phylo.filt, level = "species", detection = 1/100, prevalence = 10/100)

##species
colourCount = 13
#phylo.aggr<-phylo.filt%>%aggregate_taxa(level = "species")
p <- plot_composition(phylo.filt,
                      sample.sort = "Phenotype",
                      x.label = "Sample") +
  scale_fill_manual("Species",values = getPalette(colourCount)) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  xlab("Varieties")+ylab("Relative abundance")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic"))
p
png(filename = "./Fig/rabund_S.png", width=30 ,height=16, units='cm',res=600)
print(p)  
dev.off()

#Relative abundance plot group by resistance phenotype----

#Average by group-species
colourCount = 13
p <- plot_composition(phylo.filt,
                      average_by = "Variety", 
                      transform = "compositional",
                      group_by = "Phenotype") +
  scale_fill_manual("Species",values = getPalette(colourCount)) +
  scale_y_percent() +
  xlab("Varieties")+ylab("Relative abundance")+
  theme_bw(base_size = 20)+
  theme(#axis.text.x = element_text(angle=90, hjust=1),
        legend.text = element_text(face = "italic",size=10))
p
png(filename = "./Fig/rabund_aggre_S.png", width=28 ,height=16, units='cm',res=600)
print(p)
dev.off()


#Deseq2 by ASVs----
sample_data(phylo.keep)$Phenotype <- as.factor(sample_data(phylo.keep)$Phenotype)
dds <- phyloseq_to_deseq2(phylo.keep, ~ Phenotype)
###normalized
#https://github.com/joey711/phyloseq/issues/445
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(dds), 1, gm_mean)
###Deseq
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds)
dds <- DESeq(dds, fitType="local")
###result
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
#write.csv(res,file="./result/Deseq_RvsS_ASV.csv")
###select sig ASVs
alpha <- 0.05
res.sig <- res[(res$padj < alpha), ]
res.sig <- cbind(as(res.sig, "data.frame"), as(tax_table(phylo.keep)[rownames(res.sig), ], "matrix"))
res.sig
res.sig$log2FoldChange<-(res.sig$log2FoldChange)*(-1)
write.csv(res.sig,file="./result/Deseq_RvsS_ASV_sig95.csv")
res.sig.Rmore <- subset(res.sig, log2FoldChange < 0)

##Deseq2 plot----
###manhattan----
res <- cbind(as(res, "data.frame"), as(tax_table(phylo.keep)[rownames(res), ], "matrix"))
res$logpadj<--log10(res$padj)
range(res$log2FoldChange)
#lineartransformation
# Your vector of numbers
log2f.value <- res$log2FoldChange
# Define the old and new range
old_min <- 0
old_max <- 30
new_min <- 0
new_max <- 4

# Apply the linear transformation
transformed_l2f <- (log2f.value - old_min) / (old_max - old_min) * (new_max - new_min) + new_min

# Print the transformed numbers
res$tfl2f<-abs(transformed_l2f)

res <- res %>%
  mutate(classification = ifelse(padj > 0.05, 'No sig',
                                 ifelse(padj < 0.05 & log2FoldChange < 0, 'ENRICH',
                                        ifelse(padj < 0.05 & log2FoldChange > 0, 'DEPLETE', NA))))

getPalette =colorRampPalette(c("#fadc44","#FF5D1A","#ffbc0d","#00C891","#0097D1","#0066cc","#FC0049"))
cbbPalette<-getPalette(16)

ggplot(res, aes(x = genus, y = logpadj, color = genus, shape=classification))+
  geom_jitter(size = abs(transformed_l2f)+3)+
  #geom_jitter(size = 4.5)+
  scale_color_manual(values=cbbPalette) +
  scale_shape_manual(values = c(6,17,1)) +
  geom_hline(yintercept = 1.301, linetype = "dashed", color = "grey")+
  #scale_color_viridis_c(option = "C", direction = -1) +
  #scale_y_continuous(trans = "reverse") +
  ylim(0,25)+
  labs(x = "", y = expression(paste(-log[10],"(p-value)",sep="")),color="Genus",shape="Classification") +
  guides(color = guide_legend(ncol = 2))+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
        legend.text = element_text(face = 'italic',size=10))

ggsave("./Fig/DEG_Mahhaten.png",  width=15, height=8, dpi=600)
#ggsave("./Fig/DEG_Mahhaten_legend.png",  width=15, height=7, dpi=600)


###plot significant
ggdotchart(res.sig, x = "species", y = "log2FoldChange",
           color = "family",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#FC4E07",'#ce9eff'), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           #add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "order",                                # Order by groups
           dot.size = 3,                                 # Large dot size
           #label = round(dfm$mpg),                        # Add mpg values as dot labels
           #font.label = list(color = "white", size = 9, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_bw(base_size = 15, base_line_size = 1)                       # ggplot2 theme
) +geom_hline(aes(yintercept=0),color="dark grey") +
  theme(axis.text.y = element_text(face = 'italic'))

ggsave("./Fig/DEG_ASV.png",  width=10, height=5, dpi=300)


#DEseq by species ----
aggr.spe <- aggregate_taxa(phylo, level="species")
sample_data(aggr.spe)$Phenotype <- as.factor(sample_data(aggr.spe)$Phenotype)
dds <- phyloseq_to_deseq2(aggr.spe, ~ Phenotype)
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType="local")
#plotDispEsts(dds)
res <- results(dds)
res <- res[order(res$padj, na.last=NA), ]
alpha <- 0.05
res.sig <- res[(res$padj < alpha), ]
res.sig <- cbind(as(res.sig, "data.frame"), as(tax_table(aggr.spe)[rownames(res.sig), ], "matrix"))
res.sig$log2FoldChange<-(res.sig$log2FoldChange)*(-1)
res.sig
#B. altitudinis still there

###plot significant
ggdotchart(res.sig, x = "species", y = "log2FoldChange",
           color = "family",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#FC4E07",'#ce9eff'), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           #add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "order",                                # Order by groups
           dot.size = 3,                                 # Large dot size
           #label = round(dfm$mpg),                        # Add mpg values as dot labels
           #font.label = list(color = "white", size = 9, vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_bw(base_size = 15, base_line_size = 1)                       # ggplot2 theme
) +geom_hline(aes(yintercept=0),color="dark grey")+
  labs(x = "Species", y = "log2foldchange")+
  theme(axis.text.y = element_text(face='italic'),
        axis.title.x = element_text(vjust=0.3))

ggsave("./Fig/DEG_ASV_spe.png",  width=7, height=3, dpi=300)

#Shared ASV ----
KS7<-rowSums(otu_table(phylo)[,1:3])
KS9<-rowSums(otu_table(phylo)[,4:6])
HBS<-rowSums(otu_table(phylo)[,7:11])
HC<-rowSums(otu_table(phylo)[,12:16])
SS<-rowSums(otu_table(phylo)[,17:21])
TN11<-rowSums(otu_table(phylo)[,22:26])
TN3<-rowSums(otu_table(phylo)[,27:31])
TN5<-rowSums(otu_table(phylo)[,32:36])
aggrASV.df<-data.frame(HC,KS7,KS9,SS,HBS,TN3,TN5,TN11)
aggrASV.df.pa<-decostand(aggrASV.df,MARGIN = 1,method='pa')
aggrASV.df.pa$sums<-rowSums(aggrASV.df.pa)
aggrASV.df.pa<-cbind(as(aggrASV.df.pa, "data.frame"), as(tax_table(phylo)[rownames(aggrASV.df.pa), ], "matrix"))
write.csv(aggrASV.df.pa,'./result/aggregateASV.csv')
filter(aggrASV.df.pa,sums==2)

#Composition of ASV agrregated by cultivars
#Supplement table S3
df_ratio <- apply(aggrASV.df, 2, function(column) column/sum(column))
write.csv(df_ratio,'./result/cultivar_ratio.csv')

#Shared species ----
aggr.spe <- aggregate_taxa(phylo.keep, level="species")
aggr.spe.raw <- aggregate_taxa(phylo.keep, level="species")
KS7<-rowSums(otu_table(aggr.spe.raw)[,1:3])
KS9<-rowSums(otu_table(aggr.spe.raw)[,4:6])
HBS<-rowSums(otu_table(aggr.spe.raw)[,7:11])
HC<-rowSums(otu_table(aggr.spe.raw)[,12:16])
SS<-rowSums(otu_table(aggr.spe.raw)[,17:21])
TN11<-rowSums(otu_table(aggr.spe.raw)[,22:26])
TN3<-rowSums(otu_table(aggr.spe.raw)[,27:31])
TN5<-rowSums(otu_table(aggr.spe.raw)[,32:36])
aggr.spe_vari<-data.frame(KS7,KS9,HBS,HC,SS,TN11,TN3,TN5)
aggrspe.pa<-decostand(aggr.spe_vari,MARGIN = 1,method='pa')
aggrspe.pa$sums<-rowSums(aggrspe.pa)
###Species occur in 8 samples
filter(aggrspe.pa,sums==8) #none
###Species occur in 7 samples 
filter(aggrspe.pa,sums==7) #Priestia aryabhattai, Priestia megaterium
###Species occur in 6 samples
filter(aggrspe.pa,sums==6) #6 species
###Species occur in 5 samples
filter(aggrspe.pa,sums==5) #1 species
###Species occur in 4 samples
filter(aggrspe.pa,sums==4) #5 species
###Species occur less than 4 samples
filter(aggrspe.pa,sums<4) # 25 species

write.csv(aggrspe.pa,'./result/aggrspe_pa.csv')

#cluster analysis ----
spe<-t(otu_table(phylo.keep.norm))

##sorecent distance
spe.sor <- vegdist(spe, binary=TRUE)
###ward
spe.sor.cluster.ward <- cluster::agnes(spe.sor, method = 'ward')
png(filename = "./Fig/Cluster_spe_sor_ward.png", width=20 ,height=16, units='cm',res=300)
plot (spe.sor.cluster.ward, which.plot = 2)
dev.off()
###flexible
spe.sor.cluster.flex<- cluster::agnes(spe.sor, method ='flexible',par.method = 0.625)
png(filename = "./Fig/Cluster_spe_sor_flex.png", width=20 ,height=16, units='cm',res=300)
plot (spe.sor.cluster.flex, which.plot = 2)
dev.off()

#core microbiome ----
phylo.rel <- microbiome::transform(phylo.keep.norm, "compositional")
#Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(phylo.rel, detection = 1/100, sort = TRUE))
#Absolute population frequencies
head(prevalence(phylo.rel, detection = 1/100, sort = TRUE, count = TRUE))
core.taxa.standard <- core_members(phylo.rel, detection = 0, prevalence = 50/100)
# ASV3  ASV1 ASV18 ASV40  ASV9 ASV11  
pseq.core <- core(phylo.rel, detection = 0, prevalence = .5)
core.taxa <- taxa(pseq.core)
core.taxa
#"ASV1" "ASV3"


