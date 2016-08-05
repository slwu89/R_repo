#################################
######R Debugging Functions######
######Sean Wu 5/25/2016##########
#################################

###Useful functions for debugging R code

#retrieve source file location for loaded function
get_source <- function(FUN){
  return(attr(attr(FUN,"srcref"),"srcfile"))
}

#return all objects pointing to same memory location in R (OBJECT_OF_INTEREST)
#TURN INTO FUNCTION LATER, NEED TO REVIEW SCOPING RULES
require(pryr)
obj_address <- pryr:::address2(quote(OBJECT_OF_INTEREST),env=.GlobalEnv)
address_list <- NULL
obj_names <- ls(envir=.GlobalEnv)
for(i in 1:length(obj_names)){
  address_list[i] <- pryr:::address2(as.name(obj_names[i]),env=.GlobalEnv)
}
obj_index <- which(address_list == obj_address)
obj_names[obj_index]