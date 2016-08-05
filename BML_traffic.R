#############################
######BML Traffic Model######
######Sean Wu 5/24/2016######
#############################


####################
###Create lattice###
####################

cell_init <- function(i,j,p){
  if(runif(1) <= p){
    return(sample(c(1,2),1))
  } else {
    return(0)
  }
}

cell_initV <- Vectorize(cell_init,c("i","j"))

#lattice is of dimension rXc and with initial density p
bml_init <- function(r,c,p){
  outer(1:r,1:c,cell_initV,p=p)
}


###############################
###Create movement functions###
###############################

#blue (1) goes north, red (2) goes east 

#northward movement function
move_north <- function(mat){
  #index indicating if [i,j] was to attempt to move north if it would be blocked or not
  block_ind <- mat[c(nrow(mat),1:(nrow(mat)-1)),]!=0
  blue_ind <- mat*(mat==1)
  #hold non blue cars constant + hold blocked blue cars constant + non-blocked blue cars resorted to "move" north
  return(mat*(mat!=1) + blue_ind*block_ind + (blue_ind*!block_ind)[c(2:nrow(mat),1),])
}

#eastward movement function
move_east <- function(mat){
  #index indicating if [i,j] was to attempt to move east if it would be blocked or not
  block_ind <- mat[,c(2:ncol(mat),1)]!=0
  red_ind <- mat*(mat==2)
  #hold non red cars constant + hold blocked red cars constant + non-blocked red cars resorted to "move" east
  return(mat*(mat!=2) + red_ind*block_ind + (red_ind*!block_ind)[,c(ncol(mat),1:ncol(mat)-1)])
}

#single step function
bml_step <- function(mat){
  new_mat <- move_east(move_north(mat))
  jam <- identical(new_mat,mat)
  return(list(mat=new_mat,jam=jam))
}


##########################
###Simulation functions###
##########################


out <- list()
jam <- FALSE
i <- 1
init_mat <- bml_init(100,100,0.35)
while(!jam){
  if(i == 1){
    mat <- bml_step(init_mat)  
  } else {
    mat <- bml_step(mat$mat)
  }
  jam <- mat$jam
  out[[i]] <- mat$mat
  image(mat$mat,col=c("white","blue","red"))
  print(paste("Currently on iteration:",i))
  i <- i +1
}

for(i in 1:length(out)){
  write.table(x=out[[i]],file=paste0("time",i,".csv"),col.names=FALSE,sep=",")
}