# plot results of information scores

library("ggplot2")
library("plyr")

density_plot = function(score_file, pred) {
  scores = read.csv(score_file, header = FALSE)
  
  df = data.frame(scores)
  
  ggplot(data=df, aes(V1)) + 
    geom_density( alpha = 0.2) +
    geom_vline(xintercept = pred, color = "red") +
    theme_classic()
}

  
density_plot("PI_ATV.csv.ts_1", 94)
