#! /usr/bin/env Rscript
library("ggplot2")

genes<-list("Gag", "gp41")

for (gene in genes) {
  fileName = sprintf("/app/result_data/%sNaiveSeqsStat.csv", gene)
  data<-read.table(fileName, sep=",", header=TRUE, as.is=TRUE)
  png(filename=sprintf("/app/report/%s-naive-aachanges-dist.png", tolower(gene)),
      width=6, height=4, units = "in", res=300)
  diffAggs = aggregate(cbind(NumSequences = data$Accession) ~ NumAAChanges, data, FUN = length)
  print(ggplot(diffAggs, aes(NumAAChanges, NumSequences)) +
    geom_point() + ggtitle(sprintf(
      "AA Changes Distribution for %s Naive Sequences", gene
    )))
  dev.off()
  png(filename=sprintf("/app/report/%s-naive-stopcodons-dist.png", tolower(gene)),
      width=6, height=4.4, units = "in", res=300)
  stopAggs = aggregate(cbind(NumSequences = data$Accession) ~ NumStopCodons, data, FUN = length)
  print(ggplot(stopAggs, aes(NumStopCodons, NumSequences)) +
    geom_point() + ggtitle(sprintf(
        "Stop Codons Distribution for %s Naive Sequences", gene
      )))
  dev.off()
}
