##This code was adapted from B.J. Callahan "DADA2 + PacBio: Fecal Samples"
##https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html

#if (!requireNamespace("BiocManager", quietly = TRUE))
#BiocManager::install("dada2")
library(dada2)
library(ggpubr)
library(ggplot2)
library(Biostrings)
library(ShortRead)
library(reshape2)
library(RColorBrewer)
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(reshape2)
library(vegan)

#Setup ----
path <- "./Final/data"
list.files(path)
fns <- list.files(path, pattern="fastq.gz", full.names=TRUE)
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
theme_set(theme_bw())

#Remove Primers and Filter ----

##Remove primer and orient reads
nops <- file.path(path, "noprimers", basename(fns))
prim <- removePrimers(fns, nops, max.mismatch = 5, primer.fwd=F27, primer.rev=dada2:::rc(R1492), 
                      orient=TRUE)

##Inspect length distribution
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
mean(lens) #1453.742
hist(lens, 100)

prim[,2] / prim[,1]
prim
apply(prim, 2, sum) #reads.in 803220 , reads.out 780739, 97.08%

##Filter
filts <- file.path(path, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ=3, minLen=1000, maxLen=1800, maxN=0, rm.phix=FALSE, maxEE=2)
track
apply(track, 2, sum) #reads.in 780739, reads.out 729864, 93.48%

#Run DADA2 ----

##Dereplicate
drp <- derepFastq(filts, verbose=TRUE)

##Learn errors
err <- learnErrors(drp, errorEstimationFunction=dada2:::PacBioErrfun, BAND_SIZE=32, multithread=TRUE)

##Inspect errors
plotErrors(err)

##Denoise
dd2 <- dada(drp, err=err, BAND_SIZE=32, multithread=TRUE)
cbind(ccs=apply(prim, 2, sum)[1], primers=apply(prim, 2, sum)[2], filtered=apply(track, 2, sum)[2], denoised=sum(sapply(dd2, function(x) sum(x$denoised))))
#         ccs     primers  filtered denoised
#reads.in 803220  780739   729864   726404

#Make Sequence table
st <- makeSequenceTable(dd2); dim(st)
#40 835

#Check chimera
bim <- isBimeraDenovo(st, minFoldParentOverAbundance=3.5, multithread=TRUE)
table(bim)
sum(st[,bim])/sum(st)
st.nochim <- removeBimeraDenovo(st, method="consensus", multithread=TRUE, verbose=TRUE)
sum(st.nochim)/sum(st)
cbind(ccs=apply(prim, 2, sum)[1], primers=apply(prim, 2, sum)[2], filtered=apply(track, 2, sum)[2], denoised=sum(sapply(dd2, function(x) sum(x$denoised))),chimera=sum(st.nochim))
#         ccs     primers filtered denoised chimera
#reads.in 803220  780739   729864   726404  720895
#saveRDS(st.nochim,'./st.nochim.rds')
#st.nochim <-readRDS('./st.nochim.rds')
#Assign taxonomy with silva database v138.1 ----
tax <- assignTaxonomy(st.nochim, "./silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
#write.csv(tax,'./taxtable_silva.csv')

#make final ASV table ----
st_o<-t(st.nochim)
tax_final<-as.data.frame(cbind(tax, st_o),)
write.csv(tax_final,'./ASVtable_silva138.csv')


#Remove chloroplast, mitochondria and Eukaryota ASV ----

##create phyloseq object ----
###read sample metadata
sample.meta<-read.csv('./sample_meta.csv')
rownames(sample.meta)<-sample.meta$Sample
colnames(st_o)<-rownames(sample.meta)
###phyloseq object
OTU <- otu_table(as.matrix(st_o), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
samples <- sample_data(sample.meta)

phylo <- phyloseq(OTU, TAX, samples)

##remove uneanted ASV
pseq_chlor<-subset_taxa(phylo, Order=='Chloroplast')
pseq_chlor<-merge_taxa(pseq_chlor, 1:nrow(tax_table(pseq_chlor)))
rownames(pseq_chlor@tax_table)<-"Chloroplast"
rownames(pseq_chlor@otu_table)<-"Chloroplast"
pseq_mito<-subset_taxa(phylo, Family=='Mitochondria')
pseq_mito<-merge_taxa(pseq_mito, 1:nrow(tax_table(pseq_mito)))
rownames(pseq_mito@tax_table)<-"Mitochondria"
rownames(pseq_mito@otu_table)<-"Mitochondria"
pseq_euka<-subset_taxa(phylo, Kingdom=='Eukaryota')
pseq_euka<-merge_taxa(pseq_euka, 1:nrow(tax_table(pseq_euka)))
rownames(pseq_euka@tax_table)<-"Eukaryota"
rownames(pseq_euka@otu_table)<-"Eukaryota"

pseq_sep<-subset_taxa(phylo, Kingdom!='Eukaryota')
pseq_sep<-subset_taxa(pseq_sep, Order!='Chloroplast')
pseq_sep<-subset_taxa(pseq_sep,Family!='Mitochondria')
#saveRDS(pseq_sep,'./pseq_sep.rds')
ntaxa<-nrow(pseq_sep@otu_table)
pseq_sep_2 <- merge_taxa(pseq_sep, 1:ntaxa)
rownames(pseq_sep_2@tax_table)<-"Others"
rownames(pseq_sep_2@otu_table)<-"Others"

pseq.f<-merge_phyloseq(pseq_chlor,pseq_mito,pseq_euka,pseq_sep_2)
pseq.f<- transform(pseq.f, "compositional")
#saveRDS(pseq.f,'./pseq.f.rds')

sum(otu_table(pseq_sep))
##reads with out mit, chlor, eukaryote= 588419 (81.6%)
min(colSums(otu_table(pseq_sep))) #min read sample 366
max(colSums(otu_table(pseq_sep))) #max read sample 24736

##save trim data
asv_vector <- paste0("ASV", 1:ntaxa)
tax.table_trim<-data.frame(asv_vector,tax_table(pseq_sep)[,1:7],otu_table(pseq_sep))
write.csv(tax.table_trim,'./trim_ASV_table.csv')

#make fasta ----
format_fasta <- function(id, sequence) {
  return(paste0(">", id, "\n", sequence))
}
ASV.fasta <- mapply(format_fasta, tax.table_trim$asv_vector, rownames(tax.table_trim))
writeLines(ASV.fasta, "ASV.fasta")

#Plot the composition of each sample (Fig. 3A) ----

#pseq.f<-readRDS('./pseq.f.rds')
pseq.f<- microbiome::transform(pseq.f, "compositional")

p <- plot_composition(pseq.f)+ 
  scale_fill_manual(values = c("#00C891", "#ff6f69", "#ffcc5c" ,"#e6e6ea")) +
  scale_y_continuous(position="right")+
  labs(x="", y = "Abundance (%)",fill="Classification") +
  theme_minimal() +
  theme(legend.position = "left",
        axis.text.x = element_text(size=10,angle=90,hjust=0.5),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(vjust=1),
        legend.text = element_text())
p
ggsave('./Fig/contaimination_all.png',height = 5, width = 15, dpi=600, unit='cm')


#Plot reads count (Fig. 3A) ----
#pseq_sep<-readRDS('./pseq_sep.rds')
cbbPalette <- c("#000000", "#fadc44","#ffbc0d","#00C891","#0097D1","#0066cc","#FF5D1A","#FC0049")
reads.count<-apply(otu_table(pseq_sep), 2, sum)
reads.count<-log10(reads.count)
adiv<-estimate_richness(pseq_sep, split = TRUE)
desired_order <- c("HC", "KS7", "KS9", "SS", "HBS", "TN3", "TN5", "TN11")
adiv$Variety <-pseq_sep@sam_data$Variety
adiv$Variety <-factor(adiv$Variety,levels=desired_order)
adiv$reads.count<-reads.count
adiv$Sample<-rownames(adiv)
adiv$Sample<-factor(adiv$Sample,levels=c("HCn1","HCn2","HCn3","HCn4","HCn5",
                                         "KS7n1","KS7n2","KS7n3","KS7n4","KS7n5",
                                         "KS9n1","KS9n2","KS9n3","KS9n4","KS9n5",
                                         "SSn1","SSn2","SSn3","SSn4","SSn5",
                                         "HBSn1","HBSn2","HBSn3","HBSn4","HBSn5",
                                         "TN3n1","TN3n2","TN3n3","TN3n4","TN3n5",
                                         "TN5n1","TN5n2","TN5n3","TN5n4","TN5n5",
                                         "TN11n1","TN11n2","TN11n3","TN11n4","TN11n5"))
str(adiv)
ggbarplot(adiv,  x='Sample',  y='reads.count', fill= "Variety", color="white",legend = "right") +
  scale_fill_manual(values = cbbPalette) +
  labs(x="",y="log10(Reads counts)",fill="Varieties") +
  rotate_x_text(angle = 90)+
  scale_y_continuous(position="right")+
  guides(fill=guide_legend(ncol=2))+
  theme_minimal()+
  theme(legend.position = "left",
    axis.text.y=element_text(size=8),
        axis.text.x=element_blank(),
        legend.text=element_text(size=8),
        #axis.title.y=element_text(size=8),
        panel.grid=element_blank()
        )
  
ggsave("./Fig/readscount_bw.png", height = 5, width = 20, unit='cm', dpi=600)

#Rarecurve (Fig. S2) ----
sam.data <- data.frame(pseq_sep@sam_data)
asv.table <- otu_table(pseq_sep) %>%
  as.data.frame() %>%
  as.matrix()
rare.fun <- rarecurve(t(asv.table), step = 1000, sample = raremax, tidy = TRUE)
rare.curve.extract <- dplyr::left_join(rare.fun, sam.data, by = c("Site" = "Sample"))
rare.curve.extract$Variety <-factor(rare.curve.extract$Variety,levels=desired_order)
rare.plot <- ggplot(rare.curve.extract, aes(x = Sample, y = Species, group = Site, color = Variety)) + 
  scale_color_manual(values = cbbPalette)+
  geom_line() + 
  labs(x="Reads",  y="Number of ASVs") + 
  ggtitle("Bateria")+
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(pseq_sep)), linetype = "dashed") +
  ggtitle("") 
rare.plot

ggsave("./Fig/rarecurve.png", width=13, height=10, unit="cm",dpi=600)


