#' @title Write a NIfTI file.
#' @description This function saves a data matrix into a NIfTI file.
#'
#' @param X Data matrix with n rows (sample) and p colums (pixel).
#' @param nii a reference NIfTI-class object, representing a image with p0 voxels.
#' @param xgrid Cordinate matrix with p0 rows (voxel) and d columns (dimension). p0 >= p.
#' @param mask Logical vector of length p0, with p elements being equal to 1 and p0-p elements being equal to 0.
#' @param file The name to the file to be saved.
#'
#' @return
#' @export
#'
#' @importFrom oro.nifti writeNIfTI
#' @importFrom neurobase copyNIfTIHeader
#'
#' @examples
output_nii = function(X,nii,xgrid, filename=NULL, std=TRUE, thres = 0){

  if(std){
    X0 = matrix(0, ncol = ncol(X), nrow = nrow(X))
    for(i in 1:nrow(X)){
      if(sd(X[i,])>0){
        X0[i,] = ( X[i,] )/max(abs(X[i,]) )
      }
    }
  }
  else{
    X0 = X
  }

  if(thres>0){
    for(i in 1:nrow(X)){
      if(sd(X[i,])>0){
        ind = ( X0[i,] > quantile(X0[i,],thres))
        ind1 = ( X0[i,] < quantile(X0[i,],1-thres))
        X0[i,] = X0[i,] * ( (ind +ind1)>0 )
      }
    }
  }

  out_S = array(0, dim = c(dim(nii)[1], dim(nii)[2], dim(nii)[3], nrow(X)))

  tag = 1

  for(i in 1:nrow(xgrid)){
    out_S[xgrid[i,1],xgrid[i,2],xgrid[i,3],] = X0[,i]
  }

  copyNIfTIHeader(img = nii, arr = out_S, drop = FALSE)

  if(!is.null(filename)){
    writeNIfTI(out_S,filename = filename)
  }

  return(out_S)

}
