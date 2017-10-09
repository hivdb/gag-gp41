#! /usr/bin/env Rscript
genes<-list("gag", "gp41")

xticksGp41<-seq(50,350,20)
xticksGag<-seq(0,480,30)
xticks<-list(xticksGag, xticksGp41)
yticks<-c(seq(-4, 4, 1))
ylabels<-c(seq(-4,4,1))

gag_domains<-c("MA", "CA", "NC", "P6")
gp41_domains<-c("Extracellular", "TM",  "Cytoplasmic Domain")
gene_domains<-list(gag_domains, gp41_domains)

#### Information for rectangles:
rect_text_args<-list(col = 'black', cex = 0.8) 
MA_rect_args<-list(col = 'green', lty = 'solid')
CA_rect_args<-list(col = 'grey', lty = 'solid')
NC_rect_args<-list(col = 'dodgerblue', lty = 'solid')
p6_rect_args<-list(col = 'pink', lty = 'solid')
MA_rect_coords<-c(0, -3.6, 132, -3.2)
CA_rect_coords<-c(133, -3.6, 363, -3.2)
NC_rect_coords<-c(378, -3.6, 432, -3.2)
p6_rect_coords<-c(449, -3.6, 500, -3.2)
FU_rect_args<-list(col = 'green', lty = 'solid')
TM_rect_args<-list(col = 'grey', lty = 'solid')
CD_rect_args<-list(col = 'dodgerblue', lty = 'solid')
FU_rect_coords<-c(55, -3.6, 173, -3.2)
TM_rect_coords<-c(174,-3.6, 194, -3.2)
CD_rect_coords<-c(195, -3.6, 345, -3.2)


#### Exclude positions with fold < FOLD_FILTER
FOLD_FILTER<-0

############################################### Functions ############################################### 
readFile <- function(fileName) {
  data<-read.table(fileName, sep=",", header=TRUE, as.is=TRUE)
  return(data)
}

filterByLogFold <- function(df, cutOff) {
  df <- subset(df, df$LogFold > cutOff | df$LogFold < -cutOff)
  return(df)
} 

lookUpValue <- function(df, col1, col2, value) {
  col1Vector<-as.vector(df[,col1])
  col2Vector<-df[,col2]
  col1Index<-match(c(value),col1Vector)
  cat(paste("col1Index: ", col1Index), "\n")
  lookUpValue<-col2Vector[col1Index]
  cat(paste("lookUpValue: ", lookUpValue), "\n")
  return(lookUpValue)
}

#### Place text inside a rectangle
recttext <- function(xl, yb, xr, yt, text, rectArgs = NULL, textArgs = NULL) {
  # cat(paste('xl:', xl, ' yb:', yb, ' xr:', xr, ' yt:', yt, sep=""),"\n")
  center <- c(mean(c(xl, xr)), mean(c(yb, yt)))
  do.call('rect', c(list(xleft = xl, ybottom = yb, xright = xr, ytop = yt), rectArgs))
  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
}

sqrtLines <- function(x, y, type, cex, ...) {
  lines(x, y, type=type, cex=sqrt(cex), ...)
}

###########################################################################################################
####                                      FOLD DATA - PIs
###########################################################################################################
for (i in 1:length(genes)) {
  gene<-genes[i]
  xt<-unlist(xticks[i])

  #### Read in Fold Changes
  inputFileAAs = sprintf("/app/resultData/aaChangesByPosWPrev/%s.csv", gene)
  foldData = readFile(inputFileAAs)
  
  #### Get PI fold changes
  PI = foldData[foldData$Group=="PIs",]
  PI<-PI[complete.cases(PI),]
  PILogFolds<-PI$LogFold
  positionsPI<-unique(PI$Pos)
  numMutPI<-nrow(PI)
  PI<-filterByLogFold(PI,FOLD_FILTER)

  #### Get Control fold changes
  Controls = foldData[foldData$Group=="NNRTIs",]
  Controls<-Controls[complete.cases(Controls),]
  
  #### Remove mutations lower than lowest one in PI group
  if (gene == 'gp41') {
    Controls<-Controls[Controls$Pos>50,]
  }
  ControlsLogFolds<-Controls$LogFold
  
  ############################################################################################
  ####                               FOLD DATA - PIs
  ############################################################################################
  ## tiff(filename=sprintf("/app/report/%s-mutations.tiff", tolower(gene)), width=9, height=6.6, units = "in", res=300, compression="lzw")
  png(filename=sprintf("/app/report/%s-mutations.png", tolower(gene)), width=9, height=6.6, units = "in", res=300)
  par(fig=c(0, 1.0, 0.5, 1.0))
  par(mar=c(2, 4, 1, 1))

  plot(PI$Pos, PI$LogFold, type ='h', xaxt="n", yaxt="n", axes=FALSE, 
       ann=FALSE, ylim=c(-3.5,3.5))
  axis(2, at=yticks, label=ylabels, lwd.ticks=1, cex.axis=1.0)
  axis(1, at=xt, cex.axis=1.0)
  title(ylab = "Selection Index", line = 2, cex.lab = 1.5)
  sqrtLines(PI$Pos, PI$LogFold, type = 'p', cex = (0.7*PI$NumPts))
  abline(h=c(0,0))
  abline(h=c(1,1), lty=2)
  
  if (gene == 'gag') {
    text(-5, 3.4, pos=4, "A. PI-Treated Individuals", cex = 1.25)
    recttext(MA_rect_coords[1], MA_rect_coords[2], MA_rect_coords[3], MA_rect_coords[4],
             "MA", rectArgs = MA_rect_args, textArgs = rect_text_args)
    recttext(CA_rect_coords[1], CA_rect_coords[2], CA_rect_coords[3], CA_rect_coords[4],
             "CA", rectArgs = CA_rect_args, textArgs = rect_text_args)
    recttext(NC_rect_coords[1], NC_rect_coords[2], NC_rect_coords[3], NC_rect_coords[4],
             "NC", rectArgs = NC_rect_args, textArgs = rect_text_args)
    recttext(p6_rect_coords[1], p6_rect_coords[2], p6_rect_coords[3], p6_rect_coords[4],
             "P6", rectArgs = p6_rect_args, textArgs = rect_text_args)
  } else {
    text(45, 3.4, pos=4, "A. PI-Treated Individuals", cex = 1.25)
    recttext(FU_rect_coords[1], FU_rect_coords[2], FU_rect_coords[3], FU_rect_coords[4],
             "Extracellular", rectArgs = FU_rect_args, textArgs = rect_text_args)
    recttext(TM_rect_coords[1], TM_rect_coords[2], TM_rect_coords[3], TM_rect_coords[4],
             "TM", rectArgs = TM_rect_args, textArgs = rect_text_args)
    recttext(CD_rect_coords[1], CD_rect_coords[2], CD_rect_coords[3], CD_rect_coords[4],
             "Cytoplasmic Domain", rectArgs = CD_rect_args, textArgs = rect_text_args)
  }
  
 
  par(fig=c(0, 1.0, 0, 0.5), new = TRUE)
  par(mar=c(5, 4, 1, 1))

  plot(Controls$Pos, Controls$LogFold, type ='h', xaxt="n", yaxt="n", axes=FALSE, 
       ann=FALSE, ylim=c(-3.4,3.4))
  axis(2, at=yticks, label=ylabels, lwd.ticks=1, cex.axis=1.0)
  axis(1, at=xt, cex.axis=1.0)
  title(ylab = "Selection Index", line = 2, cex.lab = 1.5)
  title(xlab = "Amino Acid Position", line = 3, cex.lab = 1.5)
  sqrtLines(Controls$Pos, Controls$LogFold, type = 'p', cex = (0.7*Controls$NumPts))
  abline(h=c(0,0))
  abline(h=c(1,1), lty=2)
  
  
  if (gene == 'gag') {
    text(-10, 3.3, pos=4, "B. Control Individuals", cex = 1.25)
    recttext(MA_rect_coords[1], MA_rect_coords[2], MA_rect_coords[3], MA_rect_coords[4],
             "MA", rectArgs = MA_rect_args, textArgs = rect_text_args)
    recttext(CA_rect_coords[1], CA_rect_coords[2], CA_rect_coords[3], CA_rect_coords[4],
             "CA", rectArgs = CA_rect_args, textArgs = rect_text_args)
    recttext(NC_rect_coords[1], NC_rect_coords[2], NC_rect_coords[3], NC_rect_coords[4],
             "NC", rectArgs = NC_rect_args, textArgs = rect_text_args)
    recttext(p6_rect_coords[1], p6_rect_coords[2], p6_rect_coords[3], p6_rect_coords[4],
             "P6", rectArgs = p6_rect_args, textArgs = rect_text_args)
  } else {
    text(45, 3.3, pos=4, "B. Control Individuals", cex = 1.25)
    recttext(FU_rect_coords[1], FU_rect_coords[2], FU_rect_coords[3], FU_rect_coords[4],
             "Extracellular", rectArgs = FU_rect_args, textArgs = rect_text_args)
    recttext(TM_rect_coords[1], TM_rect_coords[2], TM_rect_coords[3], TM_rect_coords[4],
             "TM", rectArgs = TM_rect_args, textArgs = rect_text_args)
    recttext(CD_rect_coords[1], CD_rect_coords[2], CD_rect_coords[3], CD_rect_coords[4],
             "Cytoplasmic Domain", rectArgs = CD_rect_args, textArgs = rect_text_args)
  }

  dev.off()
}





