#! /usr/bin/env Rscript
library("scales")
library("ggplot2")
library("grid")
library("gridExtra")
ROOT = "/app"

genes<-list("gag", "gp41")


# Multiple plot function
#
# Copied from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  dev.set(file)
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


mylog_trans <- 
  function (base = exp(1), from = -0.3) 
  {
    trans <- function(x) log(x, base) - from
    inv <- function(x) base^(x + from)
    trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
  }

for (gene in genes) {
  # fileName2 = sprintf("%s/data/naiveStudies/%sStatBySeq.csv", ROOT, gene)
  # subtypes<-read.table(fileName2, sep=",", header=TRUE, as.is=TRUE)
  # subtypes = aggregate(cbind(NumSequences = subtypes$Accession) ~ Subtype, subtypes, FUN=length)
  # subtypes = subtypes[subtypes$NumSequences > 100,]
  # subtypes = subtypes[order(-subtypes$NumSequences),]
  # subtypes = subtypes$Subtype
  subtypes = c('A1', 'B', 'C', '01_AE')
  fileName = sprintf("%s/local/hyphyOutput/%sNaiveDistanceMatrix.txt", ROOT, gene)
  data<-read.table(fileName, sep="\t", header=TRUE, as.is=TRUE)
  data$Sequence1 = sub("^[^_]+_", "", data$Sequence1, perl=TRUE)
  data$Sequence2 = sub("^[^_]+_", "", data$Sequence2, perl=TRUE)
  colnames(data) <- c("Subtype", "Subtype2", "Distance")
  filtered = data[data$Subtype %in% subtypes & data$Subtype2 %in% subtypes,]
  # aggs = aggregate(
  #   Distance ~ Subtype + Subtype2, filtered,
  #   FUN=function(x) sprintf("%.1f (%.1f-%.1f,%.1f-%.1f)", median(x) * 100, min(x) * 100, max(x) * 100, quantile(x, 1/4) * 100, quantile(x, 3/4) * 100)
  # )
  # aggs = dcast(aggs, Subtype ~ Subtype2)
  # fileName = sprintf("%s/local/%sNaiveDiversifyMatrix.csv", ROOT, gene)
  # write.csv(aggs, fileName, row.names=FALSE)

  xlimmax = ceiling(quantile(data$Distance, 0.999) * 100 / 5) * 5
  filtered = filtered[filtered$Subtype == filtered$Subtype2,]
  filtered$SubtypePair <- paste(filtered$Subtype, filtered$Subtype2, sep="-")
  plot = ggplot(data, aes(Distance * 100)) +
      geom_histogram(
        aes(y=(..count..)/sum(..count..)),
        binwidth=1, size=0.5) + scale_y_continuous(labels=percent) +
      ylab('  ') + xlab('% Distance') + xlim(0, xlimmax) +
      ggtitle("M") +
      theme(
        plot.title = element_text(hjust=0.02, margin=margin(14, 0, -26, unit = "pt")),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  plots = list(plot)

  for (subtype in subtypes) {
      partial = filtered[filtered$Subtype == subtype & filtered$Subtype2 == subtype,]
      plot = ggplot(partial, aes(Distance * 100)) +
        geom_histogram(
          aes(y=(..count..)/sum(..count..)),
          binwidth=1, size=0.5) + scale_y_continuous(labels=percent) +
        ylab('  ') + xlab('% Distance') + xlim(0, xlimmax) +
        ggtitle(subtype) + theme(plot.title = element_text(hjust=0.02, margin=margin(14, 0, -26, unit = "pt")))
      if (subtype != '01_AE') {
          plot = plot + 
            theme(
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      }
      if (subtype == 'B') {
          plot = plot + ylab('% Sequences')
      }
      plots = append(plots, list(plot))
  }
  png(sprintf("%s/report/%s-naive-diversify.png", ROOT, tolower(gene)), width=6, height=8, units="in", res=300)
  ph = 456 / 300  # plot height
  bh = 120 / 300    # height of bottom x label
  grid.arrange(grobs=plots, nrow=length(plots), heights=unit(c(ph, ph, ph, ph, ph + bh), c('in', 'in', 'in', 'in', 'in')))
  dev.off()
}
