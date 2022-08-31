#' @title Transforms NIfTI to matrix
#'
#' @description This function transforms a NIfTI-class object into a matrix.
#'
#' @param nii 4D NIfTI-class object with dimensions x,y,z and t. Can be read from NIfTI file with \code{readNIfTI} function from the package \code{oro.nifti}.
#' @param mask Mask variable, also in NIfTI format.
#'
#' @importFrom oro.nifti readNIfTI
#'
#' @return List containing the data matrix with t rows and x*y*z colums (voxels), and the coordinates of the voxels.
#' @export
pre_nii = function(nii,mask){

  dim_nii = dim(mask)

  xgrid = as.matrix(expand.grid(1:dim_nii[1],1:dim_nii[2],1:dim_nii[3]))

  mask_vec = as.numeric(mask)
  data0 = NULL
  for(i in 1:dim(nii)[4]){
    data0 = cbind(data0,as.numeric(nii[,,,i])[mask_vec!=0])
  }

  out = list()
  out$data = t(data0)
  out$coords = xgrid[mask_vec!=0,]

  return(out)

}
