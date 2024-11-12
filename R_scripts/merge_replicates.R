library(ggplot2)
library(tidyverse)
variant=read.csv("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/uncombined_data/Variants.csv")
rd=read.csv("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/uncombined_data/RD236a_MeanDepth.csv")
expression=read.csv("/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/uncombined_data/GeneExpressionLogFPKM.csv", check.names=F)


rd = rd %>%
  group_by(SampleID, Lineage)%>%
  summarise(RD236aMeanDepth = mean(RD236aMeanDepth))%>%
  mutate(HasDeletion = ifelse(RD236aMeanDepth>0.1, "No Deletion", "With Deletion"))

SampleID2Lineage = variant[,2:3]
SampleID2Lineage = SampleID2Lineage[duplicated(SampleID2Lineage)==F,]

tmp = aggregate(variant[,4:ncol(variant)], 
                by=list(SampleID = variant$SampleID), 
                FUN=max)
variant = left_join(tmp, SampleID2Lineage, by="SampleID")

tmp = aggregate(expression[,4:ncol(expression)], 
                by=list(SampleID = expression$SampleID), 
                FUN=mean)
expression = left_join(tmp, SampleID2Lineage, by="SampleID")

write.csv(rd, "/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/processed_data/RD236a_MeanDepth.csv", row.names = FALSE)
write.csv(variant, "/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/processed_data/Variants.csv", row.names = FALSE)
write.csv(expression, "/Users/miaojiazheng/Desktop/Harvard/Projects/Mycobacterium/eQTL/data/processed_data/GeneExpressionLogFPKM.csv", row.names = FALSE)

