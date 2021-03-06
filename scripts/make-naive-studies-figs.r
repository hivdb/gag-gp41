#! /usr/bin/env Rscript
library("scales")
library("ggplot2")
library("LocFDRPois")

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
  pdf(
    sprintf("%s/report/%s-naive-aachanges-dist.pdf", ROOT, tolower(gene)),
    width=6, height=4)
  diffAggs = aggregate(cbind(NumSequences = data$Accession) ~ NumAAChanges, data, FUN = length)
  print(ggplot(diffAggs, aes(NumAAChanges, NumSequences)) +
    geom_point() + ggtitle(sprintf(
      "AA Changes Distribution for %s Naive Sequences", gene
    )))
  dev.off()

  pdf(sprintf("%s/report/%s-naive-stopcodons-dist.pdf", ROOT, tolower(gene)),
      width=6, height=4.4)
  stopAggs = aggregate(cbind(NumSequences = data$Accession) ~ NumStopCodons, data, FUN = length)
  print(ggplot(stopAggs, aes(NumStopCodons, NumSequences)) +
    geom_point() + scale_y_continuous(trans=mylog_trans(base=10), breaks=c(0,1,2,3,5,10,20,40,80,160,320,640,1280,2560,5120)) +
    ggtitle(sprintf(
      "Stop Codons Distribution for %s Naive Sequences", gene
    )))
  dev.off()

  adIndexFName = sprintf("%s/internalFiles/naiveStudies/apobec/%sNaiveADIndex.csv", ROOT, tolower(gene))
  data<-read.table(adIndexFName, sep=",", header=TRUE, as.is=TRUE)
  conserved <- data[1,]$NumConservedAPOBECSites
  numseqs <- nrow(data)
  locfdrResult <- SummarizeLocfdr(data$NumAPOBECs)
  lambda0 <- locfdrResult$lambda0
  apobecCutoff = 0
  for (i in 1:length(locfdrResult$locfdr_res$fdr)) {
    print(sprintf("%s APOBEC: %d (fdr=%.3f)", gene, i - 1, locfdrResult$locfdr_res$fdr[[i]]))
    if (locfdrResult$locfdr_res$fdr[[i]] < 0.05) {  # <5%
        apobecCutoff = i - 1
        break
    }
  }
  pdf(sprintf("%s/report/%s-adindex-fdr.pdf", ROOT, tolower(gene)),
      width=6, height=2)
  print(locfdrResult$locfdr_fig)
  dev.off()
  pdf(sprintf("%s/report/%s-adindex.pdf", ROOT, tolower(gene)),
      width=6, height=2)
  print(ggplot(data, aes(NumAPOBECs, fill=NumAPOBECs >= apobecCutoff)) + geom_histogram(binwidth=0.5, color="#595959", size=0.18) +
    scale_y_continuous(trans=mylog_trans(base=10), breaks=c(0,1,4,16,64,256,1024,4096)) +
    scale_fill_manual(values=c("white", "#595959")) + theme(legend.position="none") +
    # ggtitle(sprintf("Distribution of %s APOBEC Signature Mutations (λ=%.3f)", gene, lambda0)) +
    ylab('# Sequences') + xlab('# APOBEC Signature Mutations'))
  dev.off()
# }

  fileName = sprintf("%s/internalFiles/naiveStudies/%sUnusuals.csv", ROOT, tolower(gene))
  data<-read.table(fileName, sep=",", header=TRUE, as.is=TRUE)
  # locfdrResult <- SummarizeLocfdr(data$NumUnusuals)
  # lambda0 <- locfdrResult$lambda0
  # unusualCutoff = 0
  # for (i in 1:length(locfdrResult$locfdr_res$fdr)) {
  #   print(sprintf("%s Unusuals: %d (fdr=%.3f)", gene, i - 1, locfdrResult$locfdr_res$fdr[[i]]))
  #   if (locfdrResult$locfdr_res$fdr[[i]] < 0.01) {  # <1%
  #       unusualCutoff = i - 1
  #       break
  #   }
  # }
  # print(sprintf("%s Unusuals lambda0 = %.3f", gene, lambda0))
  # pdf(sprintf("%s/report/%s-naive-unusual-dist-fdr.pdf", ROOT, tolower(gene)),
  #     width=6, height=2)
  # print(locfdrResult$locfdr_fig)
  # dev.off()

  pdf(sprintf("%s/report/%s-naive-unusual-dist.pdf", ROOT, tolower(gene)), width=6, height=2)
  unusualCutoff = 11
  if (identical(gene, 'gp41')) {
      unusualCutoff = 8
  }
  aggs = aggregate(cbind(All = data$Accession) ~ NumUnusuals, data, FUN=length)
  print(ggplot(data, aes(NumUnusuals, fill=NumUnusuals >= unusualCutoff)) + geom_histogram(binwidth=0.5, color="#595959", size=0.18) +
    scale_fill_manual(values=c("white", "#595959")) + theme(legend.position="none") +
    scale_y_continuous(trans=mylog_trans(base=10), breaks=c(0,1,4,16,64,256,1024,4096)) + 
    ylab('# Sequences') + xlab('# Unusual Mutations'))
  dev.off()
}
