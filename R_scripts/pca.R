library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)
variant=read.csv("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/processed_data/Variants.csv")

variant_indexer = data.frame(
  var_id = colnames(variant)[!(colnames(variant) %in% c("Lineage", "SampleID"))],
  ref = unlist(lapply(strsplit(colnames(variant)[!(colnames(variant) %in% c("Lineage", "SampleID"))], "_"), `[`, 1)),
  pos = as.numeric(unlist(lapply(strsplit(colnames(variant)[!(colnames(variant) %in% c("Lineage", "SampleID"))], "_"), `[`, 2))),
  alt = unlist(lapply(strsplit(colnames(variant)[!(colnames(variant) %in% c("Lineage", "SampleID"))], "_"), `[`, 3))
)

variant_indexer = variant_indexer%>%
  filter( pos<4056375 | pos>4057733 )

inputs_tab = variant[, variant_indexer$var_id]

pca = prcomp(variant[, variant_indexer$var_id], scale=F, center=F)
pcs = data.frame(pca$x)
row.names(pcs) = variant$SampleID
pcs$Lineage = paste0("L", variant$Lineage)
write.csv(pcs, "/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/processed_data/VariantPCs.csv")

variance_explained = data.frame(
  variance=as.vector(pca$sdev^2/sum(pca$sdev^2))[1:15],
  PC=colnames(pcs)[1:15]
)

plt = ggplot(variance_explained)+
  geom_point(aes(x=PC, y=variance))+
  geom_text(aes(x=PC, y=variance+0.01, label=round(cumsum(variance), 2)))+
  geom_line(aes(x=PC, y=variance, group=1))+
  scale_x_discrete(limits=colnames(pcs)[1:15])+
  theme_bw()

ggsave("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/results/figures/scree.png", plt, width=8, height=6)


plt_pc1_pc2 = ggplot(pcs)+
  geom_point(aes(x=PC1, y=PC2, fill=Lineage), size=5, shape=21, color='white')+
  scale_fill_brewer(palette = "Set2")+
  theme_bw()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none")

plt_pc3_pc4 = ggplot(pcs)+
  geom_point(aes(x=PC3, y=PC4, fill=Lineage), size=5, shape=21, color='white')+
  scale_fill_brewer(palette = "Set2")+
  theme_bw()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = "none")

plt_pc5_pc6 = ggplot(pcs)+
  geom_point(aes(x=PC5, y=PC6, fill=Lineage), size=5, shape=21, color='white')+
  scale_fill_brewer(palette = "Set2")+
  theme_bw()+
  theme(axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.height = unit(2, "cm"),
        legend.position = c(1.6,0.45),
        legend.key.width = unit(1,"cm"))

plt = plot_grid(plt_pc1_pc2, plt_pc3_pc4, plt_pc5_pc6, ncol=2, align='hv')

ggsave("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/results/figures/lineage_pcs.png", plt, width=15, height=15)
