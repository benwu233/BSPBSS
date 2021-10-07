## code to prepare `DATASET` dataset goes here


mask0 = oro.nifti::readNIfTI('/home/ben/data/AAL_90_3mm.nii')
mask = (mask0 > 0)

path0 = "/home/ben/work/bspbss/real_kernel2/"
res_sum = readRDS(paste0(path0,"res_sum2_gk002_120_30k_b5k_nf500_q30_d03_newtest_noise1_.rds"))

p = ncol(res_sum$S)
n = 40
S0 = res_sum$S[c(7,19,26),]
A0 = res_sum$A[1:n,c(7,19,26)]
sigma0 = mean(res_sum$sigma)
X0 = A0%*%S0 + matrix( rnorm(n*p,0,1)*sqrt(sigma0), n,p)

mask1 = rep(0, nrow(coord))
coord0 = as.matrix(coord)
for(i in 1:nrow(coord)){
  mask1[i] = aal_map[coord0[i,][1],coord0[i,][2],coord0[i,][3]]
}

nii = readNIfTI('/home/ben/data/AAL_90_3mm.nii')

file0 = paste0(path0 ,paste0("example_data") )
output_nii(X0,mask0,coord,mask1>0,file0,std=FALSE)


usethis::use_data(DATASET, overwrite = TRUE)
