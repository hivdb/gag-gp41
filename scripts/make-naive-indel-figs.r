library("scales")
library("ggplot2")
library("gridExtra")

# data = aggregate(
#   cbind(Number=Accession) ~ Gene + Position + IndelType + InsLen,
#   data, FUN=length)
# data = data[data$IndelType == 'ins',]
# zzz = ave(data$Number, data$Position, function(x) rank(x, ties.method='first'))
# data$JitterX = runif(nrow(data), 0, 10)

mylog_trans <- 
  function (base = exp(1), from = -0.3) 
  {
    trans <- function(x) log(x, base) - from
    inv <- function(x) base^(x + from)
    trans_new("mylog", trans, inv, log_breaks(base = base), domain = c(base^from, Inf))
  }

plot_patterns <- function(pdfpath, genelen, data) {
  plot = ggplot(data, mapping=aes(y=Position, x=reorder(Pattern, -Count),
                                  ymin=Position, ymax=Position + InsLen)) +
    geom_linerange(size=2) + # position=position_dodge(.5)) +
    geom_text(size=2, aes(y=genelen, label=Count)) + # , hjust=0, vjust=0) +
    scale_y_continuous(expand = c(0, 5), limits = c(1, genelen), breaks = seq(0, genelen, 5)) +
    coord_flip() + theme_bw() +
    theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0.1, 0.2, 0.1, 0.2, "in"),
      axis.line = element_line(colour = "black"))
  pdf(pdfpath, width=24, height=8)
  print(plot)
  dev.off()
}



plot_histogram <- function(title, genelen, data) {
  plot = ggplot(data, aes(x=Position, y=Count)) +
    geom_col() +
    scale_x_continuous(expand = c(0, 0), limits = c(1, genelen), breaks = seq(0, genelen, 5)) +
    scale_y_continuous(expand = c(0, 0), trans = mylog_trans(base=10), breaks=c(1, 10, 100, 1000, 10000)) +
    theme_bw() + ggtitle(title) +
    theme(#axis.title.y=element_blank(),
          #axis.text.y=element_blank(),
          #axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.1, 0.2, 0.1, 0.2, "in"),
          axis.line = element_line(colour = "black"))
  plot
}

plot_gene <- function(gene, genelen) {
  data = read.csv(sprintf("/app/data/naiveStudies/%sNaiveIndels.csv", gene))
  data$InsLen = apply(matrix(data$InsertedCodons), 1, function(x) nchar(x) / 3)
  data$IndelTypeInt = apply(matrix(data$IndelType), 1, function(x) if (x == 'ins') 1 else -1)
  aggdata = aggregate(
    cbind(Count=Accession) ~ Position + IndelTypeInt,
    data, FUN=length)
  # aggdata$Count = aggdata$Count * aggdata$IndelTypeInt

  pdf(sprintf("/app/report/%s-naive-indels.pdf", gene), width=24, height=16)
  plots = list(
    plot_histogram(sprintf("Insertions of naive %s sequences", gene), genelen, aggdata[aggdata$IndelTypeInt == 1,]),
    plot_histogram(sprintf("Deletions of naive %s sequences", gene), genelen, aggdata[aggdata$IndelTypeInt == -1,])
  )
  grid.arrange(grobs=plots, nrow=2, ncol=1)
  dev.off()

  patterndata = data[data$IndelTypeInt == 1,]
  patterndata = aggregate(
    cbind(Pattern=sprintf('%s:%s', Position, InsLen)) ~ Accession,
    patterndata, FUN=function(x) paste(x, collapse=';')
  )
  data = merge(data[data$IndelTypeInt == 1,], patterndata, by='Accession')
  patternaggs = aggregate(
    cbind(Count=Accession) ~ Pattern + Position + InsLen,
    data, FUN=length)
  patternaggs = patternaggs[patternaggs$Count > 1,]
  plot_patterns(sprintf("/app/report/%s-ins-pattern.pdf", gene), genelen, patternaggs)
}

pdata = plot_gene('gag', 500)
plot_gene('gp41', 345)
