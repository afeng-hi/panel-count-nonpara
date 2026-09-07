
################RUNNING SIMULATION IN PARALLEL ###############

library(parallel)

workerFunc <- function(core){
  system(paste("Rscript --vanilla Bots-CLHLS-MonthCov.R", core))
}


combinX=list()
no_cores <- 25


for(i in 1:no_cores)
{
  combinX[[i]]<- i
}

c1 <-makeCluster(no_cores)

parLapply(c1,combinX,workerFunc)

stopCluster(c1)

##############################################################

