model = "GRDPG" #underlying probability model. RDPG or GRDPG.

n = 100000 #network size

alist = seq(2.75,3.25,by = 0.25) #choices for a indicating subsample sizes m = ceiling(log(n)^{1+a}) 

dlist = c(5,10,20) #choices of true rank d

rholist = c(0.01,0.04) #choices of sparsity factor 

r = 1000 #number of Monte Carlo replications

source("Functions.R")
source("Generation.R")
source("Iterations.R")
source("Tables.R")

View(tab) #Recreates the table presented in the paper. 

