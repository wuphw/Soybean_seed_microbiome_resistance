#install.packages("taxize")
library(taxize)
library(dplyr)
library(ggpubr)
#https://vimeo.com/92883063
#o<-rio::import("./ASV/ASV_dada0802.xlsx",which=)
#read BLAST result ----
o<-read.delim('./all.16S.1e-5.besthit.txt')
spe<-as.factor(o$staxids)
spe.l<-levels(spe)
#find taxonumber in NCBI database ----
#ids<-get_ids("Pachyrhizus erosus",db="ncbi")
#generate taxonomy dataframe ----
taxo.list<-classification(spe.l,db='ncbi')
#arrange the result ----
sel_class<-c('superkingdom','phylum','class','order','family','genus','species')
spe.filter<-function(df) {
  df<-df[,1:2]
  factor(df$rank,ordered=TRUE,levels=sel_class)
  rank_d<-df$rank %in% sel_class
  df[rank_d,]
}
df.all<-lapply(taxo.list, spe.filter)
custom_merge <- function(df1, df2) {
  merge(df1, df2, by = "rank", all = TRUE)
}
merged_df <- Reduce(custom_merge, df.all)
merged_df<-merged_df[order(match(merged_df$rank,sel_class)),]
merged_df<-as.data.frame(t(merged_df))
colnames(merged_df)<-merged_df[1,]
merged_df<-merged_df[-1,]
rownames(merged_df)<-c(1:nrow(merged_df))
merged_df <- merged_df %>% distinct(species,.keep_all = TRUE)
#saveRDS(merged_df,"./merged_df.rds")
#match OTU table
#o<-read.delim('./all.16S.1e-5.besthit.txt')
#merged_df<-readRDS("./merged_df.rds")
otu_feature_table<-select(o, qseqid,pident,sscinames)
colnames(otu_feature_table) <- c('qseqid','pident','species')

extract_first_two_words <- function(strings) {
  sapply(strsplit(strings, " "), function(x) paste(head(x, 2), collapse=" "))
}

otu_feature_table$species <- extract_first_two_words(otu_feature_table$species)

otu_feature_table<-otu_feature_table %>%
  left_join(merged_df,by='species')


#makeOTUtable
asv.df<-read.csv("./trim_ASV_table.csv")
asv.df<-asv.df[,-c(1,3:9)]
ASV_table_final<-otu_feature_table %>%
  left_join(asv.df,join_by(qseqid==asv_vector))

write.csv(ASV_table_final,"./ASV_table_final.csv")

#BLAST identity
BLAST_raw<-rio::import("./ASV/ASV_dada0802.xlsx", which = 3)
BLAST_filtered<-rio::import("./ASV/ASV_dada0802.xlsx", which = 4)

gghistogram(BLAST_raw, x = "pident",
            #add = "mean", 
            rug = TRUE,
            fill="#0097D1") +
  theme_bw(base_size = 10) +
  labs(x="Blast identity (%)", y= "No.of ASV")+
  theme(axis.title.x = element_text(angle=180),
        axis.text.x = element_text(angle=90,vjust=0.5),
        axis.text.y = element_text(angle=90,hjust=0.5))
ggsave("./Fig/Blastidt.png",width=5.5,height=3.5,unit="cm",dpi=600)

gghistogram(BLAST_filtered, x = "pident",
            #add = "mean", 
            fill="#0097D1",
            color='white') +
  scale_y_continuous(position="right") +
  theme_bw(base_size = 10) +
  labs(x="BLAST identity (%)", y= "No. of ASVs")+
  theme(axis.title.x = element_text(angle=180),
        axis.text.x = element_text(angle=270,vjust=0.5),
        axis.text.y = element_text(angle=270))
ggsave("./Fig/Blastidt135.png",width=5.5,height=3.5,unit="cm",dpi=600)

#palette = c("#00AFBB", "#E7B800"))
