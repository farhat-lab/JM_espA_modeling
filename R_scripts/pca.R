library(ggplot2)
library(tidyverse)
library(gridExtra)
variant=read.csv("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/processed_data/Variants.csv")

variant_indexer = data.frame(
  var_id = colnames(variant)[2:ncol(variant)],
  ref = unlist(lapply(strsplit(colnames(variant)[2:ncol(variant)], "_"), `[`, 1)),
  pos = as.numeric(unlist(lapply(strsplit(colnames(variant)[2:ncol(variant)], "_"), `[`, 2))),
  alt = unlist(lapply(strsplit(colnames(variant)[2:ncol(variant)], "_"), `[`, 3))
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
  geom_line(aes(x=PC, y=variance, group=1))+
  scale_x_discrete(limits=colnames(pcs)[1:15])+
  theme_bw()

ggsave("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/results/figures/scree.png", plt, width=8, height=6)


plt_pc1_pc2 = ggplot(pcs)+
  geom_point(aes(x=PC1, y=PC2, color=Lineage))+
  scale_color_brewer(palette = "Set2")+
  theme_bw()

plt_pc3_pc4 = ggplot(pcs)+
  geom_point(aes(x=PC3, y=PC4, color=Lineage))+
  scale_color_brewer(palette = "Set2")+
  theme_bw()

plt_pc5_pc6 = ggplot(pcs)+
  geom_point(aes(x=PC5, y=PC6, color=Lineage))+
  scale_color_brewer(palette = "Set2")+
  theme_bw()

plt_pc7_pc8 = ggplot(pcs)+
  geom_point(aes(x=PC7, y=PC8, color=Lineage))+
  scale_color_brewer(palette = "Set2")+
  theme_bw()

plt = grid.arrange(plt_pc1_pc2, plt_pc3_pc4, plt_pc5_pc6, plt_pc7_pc8, ncol=2)

ggsave("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/results/figures/lineage_pcs.png", plt, width=10, height=10)
