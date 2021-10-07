#' @title Read a NIfTI file
#'
#' @description This function transforms NIfTI file(s) into a data matrix.
#'
#' @param nii_file_list A list of NIfTI file(s).
#' @param mask_file Path to the mask file.
#'
#' @importFrom oro.nifti readNIfTI
#' @return
#' @export
#'
#' @examples
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
