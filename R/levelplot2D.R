#' @title levelplot for 2D images.
#' @description The function plots 2D images for a data matrix.
#'
#' @param S Data matrix with q rows (sample) and p colums (pixel).
#' @param lim 2-dimensional vector specifying the limits for the data.
#' @param xgrid Cordinate matrix with p0 rows (voxel) and d columns (dimension).p0 >= p.
#' @param mask Logical vector of length p0, with p elements being equal to 1 and p0-p elements being equal to 0.
#' @param color Colorbar.
#' @param path Path that images are saved at, ended with a "/" or "\".
#' @param name File names.
#'
#' @importFrom lattice levelplot
#' @importFrom grDevices colorRampPalette dev.cur dev.off png
#' @import gplots
#'
#' @return
#' @export
#'
#' @examples
levelplot2D = function(S,lim = c(min(S),max(S)), xgrid,mask = rep(1,ncol(S)),color = bluered(100),path,name){

  p = nrow(xgrid)
  S0 = matrix(0,nrow = nrow(S),ncol = p)

  tag = 1
  for(i in 1:ncol(S0)){
    if(mask[i]==1){
      S0[,i] = S[,tag]
      tag = tag + 1
    }
  }

  S1 = matrix(-1,sqrt(p),sqrt(p))

  for(q in 1:nrow(S)){
    for(v in 1:p){
      if(mask[v]==1){
        S1[xgrid[v,2], xgrid[v,1]] = S0[q,v]
      }
    }

    while (dev.cur()>1) dev.off()
    my_palette = color
    seq1 = seq(lim[1],lim[2], by = (lim[2] - lim[1])/4)
    seq2 = seq(lim[1],lim[2], by = (lim[2] - lim[1])/length(color))

    if(q == nrow(S)){
      png(paste0(path,"/",name,"_",q,".png"),width = 750, height = 600)

      plot = levelplot(S1,
                       scales=list(x=list(cex=3.5),y=list(cex=3.5)),
                       xlab = NULL,
                       ylab = NULL,
                       colorkey = list(labels=list(cex=3.5, at = seq1 )),
                       at = seq2 ,col.regions = my_palette,cuts=length(my_palette)-1)

    }
    else{
      png(paste0(path,"/",name,"_",q,".png"),width = 650, height = 600)
      plot = levelplot(S1,
                       scales=list(x=list(cex=3.5),y=list(cex=3.5)),
                       xlab = NULL,
                       ylab = NULL,
                       colorkey = NULL,
                       at = seq2 ,col.regions = my_palette,cuts=length(my_palette)-1)
    }
    print(plot)
    while (dev.cur()>1) dev.off()
  }

}
