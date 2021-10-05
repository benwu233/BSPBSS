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
#' @importFrom gplots colorpanel
#' @import ggplot2
#'
#' @export
#'
levelplot2D = function(S, coords, lim = c(min(S),max(S)), xlim=c(0,max(coords[,1])), ylim=c(0,max(coords[,2])),
                       color = bluered(100),layout = c(1,nrow(S)),file=NULL){

  q = nrow(S)

  data.list = list()
  data0 = NULL
  for(j in 1:q){
    data0 = rbind(data0, data.frame(x = coords[,1], y = coords[,2], value = S[j,], comp = j) )
  }

  plot =  ggplot(data0,aes(x=x, y=y,fill=value)) +
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
