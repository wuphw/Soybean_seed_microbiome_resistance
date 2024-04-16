library(ggtreeExtra)
library(ggtree)
#https://yulab-smu.top/treedata-book/chapter4.html
library(phyloseq)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyverse)
library(reshape)
library(microbiome)
library(vegan)
#install.packages("paletteer")
library(paletteer)

phylo <- readRDS("./phylo_keep.rds")

#prepare meta/ absense presense
phylo.v <- merge_samples(phylo, "Variety")
ASV_table<-as.data.frame(otu_table(phylo.v))
ASV_pa<-decostand(ASV_table,method='pa')
ASV_pa<-t(ASV_pa)
ASV_pa<-melt(ASV_pa, variable.name="Veriaty", value.name="Existance")
colnames(ASV_pa)<-c('ASV','Varty','Existance')
ASV_pa$Existance <- paste(ASV_pa$Varty, ASV_pa$Existance)
ASV_pa$Existance[grepl("0", ASV_pa$Existance)] = "Absence"
ASV_pa$Varty <- factor(ASV_pa$Varty, levels=c("HC", "KS7", "KS9","SS","HBS", "TN3", "TN5", "TN11"))
ASV_pa$Existance <- factor(ASV_pa$Existance,
                    levels=c("HC 1", "KS7 1", "KS9 1","SS 1", "HBS 1", "TN3 1", "TN5 1", "TN11 1", 
                             "Absence"))
ASV_table2<-as.data.frame(t(ASV_table))

p<-ggtree(phylo, branch.length='none', 
          layout='circular',size=0.1) + 
  geom_tippoint(mapping=aes(color=family), 
                size=2,
                show.legend=TRUE)
p
#add ASVname
p1 <- p %<+% ASV_table2 +
  geom_tippoint(aes(color=family),
                alpha=0,
                show.legend=TRUE) +
  geom_tiplab(aes(color=family),
              align=TRUE,
              linetype=3,
              size=3,
              linesize=0.2,
              offset=1,
              show.legend=FALSE
  ) +
  scale_color_paletteer_d("ggthemes::Tableau_20")
#"ggthemes::Classic_20"
#"cartography::pastel.pal", 20
#"ggthemes::Tableau_20"
p1
ggsave("./Fig/Phylotree.png",width=20,height = 15, unit="cm",dpi=600)


#add pa result ----
p2 <- p1 +
  geom_fruit(
    data=ASV_pa,
    geom=geom_tile,
    mapping=aes(x=Varty, y=ASV, fill=Existance),
    #width=0.1,
    color="white",
    pwidth=0.1,
    offset=0.2
  ) +
  scale_fill_manual(
    name="Variety",
    values=c("#586085FF","#C9CCEAFF","#DFA398FF","#9C6755FF","#659794FF","#EA967CFF",
    "#F5C98EFF","#D65B5AFF","white"),
    na.translate=FALSE,
    guide=guide_legend(keywidth=0.5,
                       keyheight=0.5,
                       order=3
    )
  )

p2
ggsave("./Fig/Phylotree_var.png",width=25,height = 20, unit="cm",dpi=600)

