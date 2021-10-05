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
pre_nii = function(nii_file_list,mask_file){

  xgrid = as.matrix(expand.grid(1:dim_nii[1],1:dim_nii[2],1:dim_nii[3]))

  mask = readNIfTI(mask_file)
  dim_nii = dim(mask)
  p = dim_nii[1]*dim_nii[2]*dim_nii[3]

  tmp_mask = rep(0,p)
  for(j in 1:nrow(xgrid)){
    tmp_mask[j] = mask[xgrid[j,][1],xgrid[j,][2],xgrid[j,][3]]
  }


  data0 = NULL
  for(i in length(nii_file_list)){

    nii = readNIfTI(nii_file_list[[i]])

    dim_nii0 = dim(nii)
    if(length(dim_nii0)==4){
      n0 = dim_nii0[4]
      tmp = matrix(0, nrow = p, ncol = n0)

      for(j in 1:nrow(xgrid)){
        if(tmp_mask[j]!=0){
          tmp[j,] = nii[xgrid[j,][1],xgrid[j,][2],xgrid[j,][3],]
        }
      }

    }
    else if(length(dim_nii0)==3){
      n0 = 1
      tmp = matrix(0, nrow = p, ncol = n0)

      for(j in 1:nrow(xgrid)){
        if(tmp_mask[j]!=0){
          tmp[j,] = nii[xgrid[j,][1],xgrid[j,][2],xgrid[j,][3]]
        }
      }

    }

    data0 = cbind(data0,tmp)

  }

  out = list()
  out$data = data0
  out$mask = tmp_mask
  out$grid = xgrid

  return(out)

}
