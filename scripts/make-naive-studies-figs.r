#! /usr/bin/env Rscript
library("scales")
library("ggplot2")

ROOT = "/app"

genes<-list("gag", "gp41")

mylog_trans <- 
  function (base = exp(1), from = -0.3) 
  {
    trans <- function(x) log(x, base) - from
    inv <- function(x) base^(x + from)
    trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
  }


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
    geom_point() + scale_y_continuous(trans=mylog_trans(base=10), breaks=c(0,1,2,3,5,10,20,40,80,160,320,640,1280,2560,5120)) +
    ggtitle(sprintf(
      "Stop Codons Distribution for %s Naive Sequences", gene
    )))
  dev.off()

  adIndexFName = sprintf("%s/data/naiveStudies/apobec/%sNaiveADIndex.csv", ROOT, tolower(gene))
  data<-read.table(adIndexFName, sep=",", header=TRUE, as.is=TRUE)
  conserved <- data[1,]$NumConservedAPOBECSites
  numseqs <- nrow(data)
  png(filename=sprintf("%s/report/%s-adindex.png", ROOT, tolower(gene)),
      width=6, height=2, units = "in", res=300)
  print(ggplot(data, aes(NumAPOBECs, fill=NumAPOBECs >= 3)) + geom_histogram(binwidth=0.5, color="#595959", size=0.18) +
    scale_y_continuous(trans=mylog_trans(base=10), breaks=c(0,1,4,16,64,256,1024,4096)) +
    scale_fill_manual(values=c("white", "#595959")) + theme(legend.position="none") +
    ylab('# Sequences') + xlab('# APOBEC Signature Mutations'))
  dev.off()
# }

# for (gene in genes) {
  png(sprintf("%s/report/%s-naive-unusual-dist.png", ROOT, tolower(gene)), width=6, height=2, units="in", res=300)
  # pdf(sprintf("%s/report/%s-naive-unusual-dist.pdf", ROOT, tolower(gene)), width=6, height=4)
  fileName = sprintf("%s/data/naiveStudies/%sStatBySeq.csv", ROOT, gene)
  data<-read.table(fileName, sep=",", header=TRUE, as.is=TRUE)
  aggs = aggregate(cbind(All = data$Accession) ~ NumUnusuals, data, FUN=length)
  cutoff = 15
  if (identical(gene, 'gp41')) {
      cutoff = 10
  }
  print(ggplot(data, aes(NumUnusuals, fill=NumUnusuals >= cutoff)) + geom_histogram(binwidth=0.5, color="#595959", size=0.18) +
    xlim(-1, 40) + scale_fill_manual(values=c("white", "#595959")) + theme(legend.position="none") +
    scale_y_continuous(trans=mylog_trans(base=10), breaks=c(0,1,4,16,64,256,1024,4096)) + 
    ylab('# Sequences') + xlab('# Unusual Mutations'))
  #   ggtitle(sprintf("Unusual Mutations for %s Naive Sequences (All, n=%d)", gene, nrow(data))))

  # for (stype in unique(data$Subtype)) {
  #   partial <- data[data$Subtype == stype,]
  #   len = nrow(partial)
  #   n = 5
  #   if (len <= 5) {
  #       n = len
  #   }
  #   partialAggs = aggregate(cbind(NumSequence = partial$Accession) ~ NumUnusuals, partial, FUN=length)
  #   colnames(partialAggs)[2] <- stype

  #   aggs = merge(x=aggs, y=partialAggs, by="NumUnusuals", all=TRUE)
  #   print(ggplot(partial, aes(NumUnusuals)) + geom_histogram(binwidth=0.5) + xlim(-1, 40) +
  #     ylab('# Sequences') + xlab('# Unusual Mutations') + scale_y_continuous(breaks= pretty_breaks(n)) +
  #     ggtitle(sprintf("Unusual Mutations for %s Naive Sequences (%s, n=%d)", gene, stype, nrow(partial))))
  # }
  dev.off()
  # write.csv(aggs, sprintf("%s/report/%s-naive-unusual-dist.csv", ROOT, tolower(gene)), row.names=FALSE)
}
