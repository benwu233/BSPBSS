#' @title levelplot for 2D images.
#' @description The function plots 2D images for a data matrix.
#'
#' @param S Data matrix with q rows (sample) and p colums (pixel).
#' @param lim 2-dimensional numeric vector, specifying the limits for the data.
#' @param xlim 2-dimensional numeric vector, specifying the lower and upper limits of \code{x}.
#' @param ylim 2-dimensional numeric vector, specifying the lower and upper limits of \code{y}.
#' @param coords Coordinates matrix with p rows (pixel) and 2 columns (dimension), specifying the coordinates of the data points.
#' @param layout 2-dimensional numeric vector, specifying the number of rows and number of columns for the layout of components.
#' @param color Colorbar.
#' @param file Name of the file to be saved.
#'
#' @return No return value.
#'
#' @import gplots
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' sim = sim_2Dimage(length = 30, sigma = 5e-4, n = 30, smooth = 6)
#' levelplot2D(sim$S,lim = c(-0.04,0.04), sim$coords)
#'
levelplot2D = function(S, coords, lim = c(min(S),max(S)), xlim=c(0,max(coords[,1])), ylim=c(0,max(coords[,2])),
                       color = bluered(100),layout = c(1,nrow(S)),file=NULL){

  q = nrow(S)
  p = ncol(S)

  coordx = NULL
  coordy = NULL
  value = NULL
  comp = NULL
  for(j in 1:q){
    coordx = c(coordx,coords[,1])
    coordy = c(coordy,coords[,2])
    value = c(value,S[j,])
    comp = c(comp,rep(j,p))
    #data0 = rbind(data0, data.frame(x = coords[,1], y = coords[,2], value = S[j,], comp = j) )
  }
  data0 = data.frame(coordx = coordx, coordy = coordy, value = value, comp = comp)

  plot =  ggplot(data0,aes(x=coordx, y=coordy,fill=value)) +
    geom_tile( ) + xlim(xlim) + ylim(ylim) +
    facet_wrap(~comp, nrow=layout[1],ncol=layout[2]) +
    labs(x = NULL, y = NULL) +
    scale_fill_gradientn(name = NULL,lim = c(lim[1],lim[2]),na.value = "black",colors=color) +
    theme_dark() +
    coord_equal()

  print(plot)

  if(!is.null(file)){
    ggsave(file,plot)
  }

}
