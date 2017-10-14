#! /usr/bin/env Rscript
library("ggplot2")

ROOT = "/app"

genes<-list("gag", "gp41")

for (gene in genes) {
  fileName = sprintf("%s/data/naiveStudies/%sStatBySeq.csv", ROOT, gene)
  data<-read.table(fileName, sep=",", header=TRUE, as.is=TRUE)
  png(filename=sprintf("%s/report/%s-naive-aachanges-dist.png", ROOT, tolower(gene)),
      width=6, height=4, units = "in", res=300)
  diffAggs = aggregate(cbind(NumSequences = data$Accession) ~ NumAAChanges, data, FUN = length)
  print(ggplot(diffAggs, aes(NumAAChanges, NumSequences)) +
    geom_point() + ggtitle(sprintf(
      "AA Changes Distribution for %s Naive Sequences", gene
    )))
  dev.off()

  png(filename=sprintf("%s/report/%s-naive-stopcodons-dist.png", ROOT, tolower(gene)),
      width=6, height=4.4, units = "in", res=300)
  stopAggs = aggregate(cbind(NumSequences = data$Accession) ~ NumStopCodons, data, FUN = length)
  print(ggplot(stopAggs, aes(NumStopCodons, NumSequences)) +
    geom_point() + scale_y_log10(breaks=c(1,2,3,5,10,20,40,80,160,320,640,1280,2560,5120)) + ggtitle(sprintf(
      "Stop Codons Distribution for %s Naive Sequences", gene
    )))
  dev.off()

  adIndexFName = sprintf("%s/data/naiveStudies/apobec/%sNaiveADIndex.csv", ROOT, tolower(gene))
  data<-read.table(adIndexFName, sep=",", header=TRUE, as.is=TRUE)
  conserved <- data[1,]$NumConservedAPOBECSites
  numseqs <- nrow(data)
  png(filename=sprintf("%s/report/%s-adindex.png", ROOT, tolower(gene)),
      width=6, height=2, units = "in", res=300)
  print(ggplot(data, aes(NumAPOBECs)) + geom_histogram(binwidth=1) +
    scale_y_log10(breaks=c(1,4,16,64,256,1024,4096)) + ylab('# Sequences') + xlab('# APOBEC Signature Mutations'))
  dev.off()

}
