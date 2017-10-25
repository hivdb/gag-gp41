#! /usr/bin/env Rscript
library("scales")
library("ggplot2")
library("grid")
library("gridExtra")
ROOT = "/app"

genes<-list("gag", "gp41")

for (gene in genes) {
  # fileName2 = sprintf("%s/data/naiveStudies/%sStatBySeq.csv", ROOT, gene)
  # subtypes<-read.table(fileName2, sep=",", header=TRUE, as.is=TRUE)
  # subtypes = aggregate(cbind(NumSequences = subtypes$Accession) ~ Subtype, subtypes, FUN=length)
  # subtypes = subtypes[subtypes$NumSequences > 100,]
  # subtypes = subtypes[order(-subtypes$NumSequences),]
  # subtypes = subtypes$Subtype
  subtypes = c('A1', 'B', 'C', '01_AE')
  fileName = sprintf("%s/local/naiveStudies/%sNaiveDistance.csv", ROOT, gene)
  data<-read.table(fileName, sep=",", header=TRUE, as.is=TRUE)
  data$Sequence1 = sub("^[^|]+\\|", "", data$Sequence1, perl=TRUE)
  data$Sequence2 = sub("^[^|]+\\|", "", data$Sequence2, perl=TRUE)
  colnames(data) <- c("Subtype", "Subtype2", "Distance")
  filtered = data[data$Subtype %in% subtypes & data$Subtype2 %in% subtypes,]
  # aggs = aggregate(
  #   Distance ~ Subtype + Subtype2, filtered,
  #   FUN=function(x) sprintf("%.1f (%.1f-%.1f,%.1f-%.1f)", median(x) * 100, min(x) * 100, max(x) * 100, quantile(x, 1/4) * 100, quantile(x, 3/4) * 100)
  # )
  # aggs = dcast(aggs, Subtype ~ Subtype2)
  # fileName = sprintf("%s/local/%sNaiveDiversifyMatrix.csv", ROOT, gene)
  # write.csv(aggs, fileName, row.names=FALSE)

  xmax = max(data$Distance) * 100
  xlimmax = ceiling(xmax / 5) * 5
  filtered = filtered[filtered$Subtype == filtered$Subtype2,]
  filtered$SubtypePair <- paste(filtered$Subtype, filtered$Subtype2, sep="-")
  plot = ggplot(data, aes(Distance * 100)) +
      geom_histogram(
        aes(y=(..count..)/sum(..count..)),
        binwidth=1, size=0.5) +
      geom_vline(xintercept=xmax, linetype="dashed") + scale_y_continuous(labels=percent, breaks=pretty_breaks(4)) +
      ylab('% Sequence Pairs') + xlab('% Distance') + xlim(0, xlimmax) +
      ggtitle("M") +
      theme(
        plot.title = element_text(hjust=0.02, margin=margin(10, 0, -22, unit = "pt")),
        axis.title.x=element_blank(),
        axis.title.y=element_text(color="white"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  plots = list(plot)

  for (subtype in subtypes) {
      partial = filtered[filtered$Subtype == subtype & filtered$Subtype2 == subtype,]
      xmax = max(partial$Distance) * 100
      plot = ggplot(partial, aes(Distance * 100)) +
        geom_histogram(
          aes(y=(..count..)/sum(..count..)),
          binwidth=1, size=0.5) +
        geom_vline(xintercept=xmax, linetype="dashed") + scale_y_continuous(labels=percent, breaks=pretty_breaks(4)) +
        ylab('% Sequence Pairs') + xlab('% Distance') + xlim(0, xlimmax) +
        ggtitle(subtype) + theme(plot.title = element_text(hjust=0.02, margin=margin(10, 0, -22, unit = "pt")))
      if (subtype != '01_AE') {
          plot = plot + 
            theme(
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      }
      if (subtype != 'B') {
          plot = plot + theme(axis.title.y=element_text(color="white"))
      }
      plots = append(plots, list(plot))
  }
  png(sprintf("%s/report/%s-naive-diversify.png", ROOT, tolower(gene)), width=6, height=8, units="in", res=300)
  ph = 457 / 300  # plot height
  bh = 115 / 300  # height of bottom x label
  grid.arrange(grobs=plots, nrow=length(plots), heights=unit(c(ph, ph, ph, ph, ph + bh), c('in', 'in', 'in', 'in', 'in')))
  dev.off()
}
