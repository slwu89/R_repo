##################################################
######Packages to install on fresh R install######
######Sean Wu (no set date)#######################
##################################################

packages <- c("Rcpp","RcppArmadillo","tmvtnorm","MASS","microbenchmark","devtools","ape","rjson","RcppProgress","scatterplot3d")
install.packages(packages)

library(devtools)
devtools::install_github("bwlewis/rthreejs")
devtools::install_github("sdwfrost/epiwidgets")
devtools::install_github("sbfnk/fitR")
devtools::install_github('thomasp85/ggforce')
devtools::install_github('thomasp85/ggraph')

source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
