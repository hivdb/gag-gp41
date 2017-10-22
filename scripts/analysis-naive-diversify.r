#! /usr/bin/env Rscript
library("reshape2")
library("ggplot2")
ROOT = "/app"

genes<-list("gag", "gp41")

for (gene in genes) {
  fileName2 = sprintf("%s/data/naiveStudies/%sStatBySeq.csv", ROOT, gene)
  subtypes<-read.table(fileName2, sep=",", header=TRUE, as.is=TRUE)
  subtypes = aggregate(cbind(NumSequences = subtypes$Accession) ~ Subtype, subtypes, FUN=length)
  subtypes = subtypes[subtypes$NumSequences > 50,]
  subtypes = subtypes[order(-subtypes$NumSequences),]
  subtypes = subtypes$Subtype
  fileName = sprintf("%s/local/hyphyOutput/%sNaiveDistanceMatrix.txt", ROOT, gene)
  data<-read.table(fileName, sep="\t", header=TRUE, as.is=TRUE)
  data$Sequence1 = sub("^[^_]+_", "", data$Sequence1, perl=TRUE)
  data$Sequence2 = sub("^[^_]+_", "", data$Sequence2, perl=TRUE)
  colnames(data) <- c("Subtype", "Subtype2", "Distance")
  filtered = data[data$Subtype %in% subtypes & data$Subtype2 %in% subtypes,]
  aggs = aggregate(
    Distance ~ Subtype + Subtype2, filtered,
    FUN=function(x) sprintf("%.1f (%.1f-%.1f,%.1f-%.1f)", median(x) * 100, min(x) * 100, max(x) * 100, quantile(x, 1/4) * 100, quantile(x, 3/4) * 100)
  )
  aggs = dcast(aggs, Subtype ~ Subtype2)
  fileName = sprintf("%s/local/%sNaiveDiversifyMatrix.csv", ROOT, gene)
  write.csv(aggs, fileName, row.names=FALSE)

  xlimmax = round(max(data$Distance * 100) / 10) * 10
  pdf(sprintf("%s/local/%sNaiveDiversify.pdf", ROOT, tolower(gene)), width=20, height=10)
  filtered = filtered[filtered$Subtype >= filtered$Subtype2,]
  filtered$SubtypePair <- paste(filtered$Subtype, filtered$Subtype2, sep="-")
  print(ggplot(data, aes(Distance * 100)) + geom_freqpoly(binwidth=0.1) +
    geom_freqpoly(aes(Distance * 100, colour=SubtypePair), data=filtered, binwidth=0.1) +
    scale_y_log10(breaks=c(1,8,64,512,4096,32768,262144)) + ylab('# Sequences') + xlab('# Distance') + xlim(0, xlimmax) +
    ggtitle(sprintf("Distance distribution of %s Naive Sequences (All, n=%d)", gene, nrow(data))))
  for (subtype in subtypes) {
    for (subtype2 in subtypes) {
      if (subtype < subtype2) {
          next
      }
      partial = filtered[filtered$Subtype == subtype & filtered$Subtype2 == subtype2,]
      print(ggplot(partial, aes(Distance * 100)) + geom_freqpoly(binwidth=0.1) +
        scale_y_log10(breaks=c(1,8,64,512,4096,32768,262144)) + ylab('# Sequences') + xlab('# Distance') + xlim(0, xlimmax) +
        ggtitle(sprintf("Distance distribution of %s Naive Sequences (%s-%s, n=%d)", gene, subtype, subtype2, nrow(partial))))
    }
  }
  dev.off()
}
