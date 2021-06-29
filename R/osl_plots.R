#' plots a profile of OSL sample densities
#'
#' @param df data frame, contains sample age-depth grids
#' @param Xmin integer, youngest age in years before present for the plot
#' @param Xmax integer, oldest age in years before present for the plot
#' @param densScale integer, value to scale the OSL densities on the plot
#' @param Ymin integer, shallowest depth for the plot
#' @param Ymax integer, deepest depth for the plot
#' @param Ybreaks vector, breaks for depths
#' @param Xbreaks vector, breaks for ages
#' @param Xlabels vector, labels for ages
#' @param vLine integer, value to draw a vertical red line showing the age of interest
#' @param facetNo integer, number of facets
#'
#' @return ggplot of OSL profile
#' @export
#'
#' @examples
plot.OSLProfile <- function(df,
                            Xmin = 0,
                            Xmax = 100000,
                            densScale = 8,
                            Ymin = 0,
                            Ymax = 250,
                            Ybreaks = seq(0,300,20),
                            Xbreaks = c(seq(0,plotRange[2],20000)),
                            Xlabels = c(0,paste(seq(20000,plotRange[2],20000)/1000,"k",sep="")),
                            vLine = 0,
                            facetNo = 7){

  plotRange <- c(Xmin,Xmax)
  df$dx <- df$ageGrid # df is a sample density data frame output from sample.density function() above
  df$dx[which(df$dx<=plotRange[1])]<-0
  df$dx[which(df$dx>=plotRange[2])]<-plotRange[2]
  df$Site <- factor(df$Site)


  p <- ggplot2::ggplot() + ggplot2::geom_polygon(data = df, ggplot2::aes(x = dx, y = -(densScale*dy-Depth), group = sample, fill = Site),alpha = 0.4) +
    ggplot2::scale_fill_brewer(palette = 'Set1') +
    ggplot2::scale_y_reverse(breaks = Ybreaks,limits = c(Ymax,Ymin)) +
    ggplot2::scale_x_reverse(breaks = Xbreaks,labels = Xlabels) +
    ggplot2::labs(y="Depth (cm)",x="Years BP") +
    ggplot2::facet_grid(facets = ~ df[,facetNo]) +
    ggplot2::theme(legend.position="none")

  if(vLine > 0){
    myRect <- data.frame(xmax = vLine, xmin = -Inf, ymin = -Inf, ymax = Inf)
    p <- p + ggplot2::geom_rect(data = myRect, ggplot2::aes(xmax = xmax, xmin = xmin, ymin = ymin, ymax = ymax),fill = 'gray',alpha = 0.4) +
      ggplot2::geom_vline(xintercept = vLine, col = 'gray',size = 2)
  }

  return(p)

}
