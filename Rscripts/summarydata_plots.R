library(ggplot2)
library(grid)
library(gridExtra)
library(binom)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

make_df <- function(threshold = 14400, lchem = c(0,50,100,200,300,400,500,600,700), be = c(0,1), kt = 15, data = "motile_summary.txt") {
  data <- read.delim(file = data, sep = " ", header = FALSE, col.names = c("MCS","LCHEM","KT","BE","CTL","CTL.INF","KCELL","INF","DEAD","BEH"))
  data <- subset(data, data$KT == kt)
  df <- data.frame()
  for ( lc in lchem ) {
    data_lc <- subset(data, data$LCHEM == lc)
    for ( entry in be ) {
      data_be <- subset(data_lc, data$BE == entry)
      data_be <- data_be[complete.cases(data_be),]
      
      data_be$DAMAGE <- as.numeric(as.character(data_be$INF)) + as.numeric(as.character(data_be$DEAD))
      
      #data_be$DAMAGE <- data_be$INF + data_be$DEAD
      conf <- binom.confint(sum((data_be$MCS < threshold)&(data_be$INF == 0)), dim(data_be)[1], methods = 'prop.test')
      df <- rbind(df, c(lc, entry, conf[,4], mean(data_be$DAMAGE), sd(data_be$DAMAGE), conf[,5], conf[,6]))
      #print(dim(data_be)[1])
      conf <- binom.confint(sum((data_be$MCS < threshold)&(data_be$INF == 0)), dim(data_be)[1], methods = 'prop.test')
    }
  }
  #df <- df[complete.cases(df),]
  names(df) <- c("LambdaChemotaxis", "BiasedEntry", "FractionCleared", "TissueDamageMean","TissueDamageSD", "FCmin", "FCmax")
  df$BiasedEntry <- factor(df$BiasedEntry)
  return(df)
}

make_plot_fractioncleared <- function(threshold = 14400, kt = 30, be = 0, legend = FALSE, titletext = "title") {
  col = "#AA0000"
  if (be == 1) {
    col = "#00AA00"
  }
  dfsens <- make_df(threshold = threshold, kt = kt, be = c(be), data = "data/sensitive/sensitive_summary.txt")
  dfsens$Behavior = "sensitive"
  dfarr <- make_df(threshold = threshold, kt = kt, be = c(be), data = "data/arresting/arresting_summary.txt")
  dfarr$Behavior = "arresting"
  dfmot <- make_df(threshold = threshold, kt = kt, be = c(be), data = "data/motile/motile_summary.txt")
  dfmot$Behavior = "motile"
  dftot <- rbind(dfsens, dfarr, dfmot)
  
  if (legend) {
  ggplot(data = dftot,aes(x = LambdaChemotaxis, y = FractionCleared, lty = Behavior)) + 
    geom_errorbar(data = dftot, mapping=aes(x=LambdaChemotaxis, ymin=FCmin, ymax=FCmax), width=25, size=.5, color="grey50") + 
    geom_point(size = 1, colour = col) + 
    geom_line(size = 1, colour = col) +
    theme_classic() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "#000000")) +
    theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, colour = "#000000"), plot.margin = margin(20, 20, 20, 20)) +
    theme(axis.title.y = element_text(angle = 90, vjust = 5, colour = "#000000"), plot.margin = margin(20, 20, 20, 20)) +
    theme(axis.text.x = element_text(colour = "#000000")) +
    theme(text = element_text(size=20)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs( x = expression(lambda[chemotaxis]), y = paste0("Fraction cleared at ",threshold/3600," hours")) +
    ylim(0,1) + 
    xlim(-25,725) +
    scale_x_continuous(breaks=c(0,50,100,200,300,400,500,600,700)) +
    ggtitle(titletext) +
    scale_linetype_manual(values=c("longdash", "dotted", "solid")) +
    theme(legend.justification=c(1,0), legend.position=c(0.95,0.05)) + 
    theme(legend.key.width = unit(2,"cm"), legend.text=element_text(size=20))
  } else {
  ggplot(data = dftot,aes(x = LambdaChemotaxis, y = FractionCleared, lty = Behavior)) + 
    geom_errorbar(data = dftot, mapping=aes(x=LambdaChemotaxis, ymin=FCmin, ymax=FCmax), width=25, size=.5, color="grey50") + 
    geom_point(size = 1, colour = col) + 
    geom_line(size = 1, colour = col) +
    theme_classic() + 
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5, colour = "#000000")) +
    theme(axis.title.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, colour = "#000000"), plot.margin = margin(20, 20, 20, 20)) +
    theme(axis.title.y = element_text(angle = 90, vjust = 5, colour = "#000000"), plot.margin = margin(20, 20, 20, 20)) +
    theme(axis.text.x = element_text(colour = "#000000")) +
    theme(text = element_text(size=20)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs( x = expression(lambda[chemotaxis]), y = paste0("Fraction cleared at ",threshold/3600," hours")) +
    ylim(0,1) + 
    xlim(-25,725) +
    scale_x_continuous(breaks=c(0,50,100,200,300,400,500,600,700)) +
    ggtitle(titletext) +
    scale_linetype_manual(values=c("longdash", "dotted", "solid")) +
    theme(legend.position = "none")
  }
}


#args = c(7200,45,0,"Results/kt_45_be_0.pdf")
threshold = strtoi(args[1])
killingtime = strtoi(args[2])
biasedentry = strtoi(args[3])
savename = args[4]
#print(args)
entrytype = ""
if(biasedentry == 0){
  entrytype = "random entry"
} else{
  entrytype = "biased entry"
}
result <- make_plot_fractioncleared(threshold = threshold, kt = killingtime, be = biasedentry, legend = FALSE, titletext = paste0("Killing time ",killingtime,", ",entrytype))
ggsave(savename)
