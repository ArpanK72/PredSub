#Set the directory to be this Hypothesis_Testing folder

n        <- 30000
d        <- 10
rho_n    <- 0.01
r        <- 200
n_boot   <- 100
alist    <- seq(2.875, 3.125, by = 0.125)
epsilons <- c(0, 0.02, 0.04, 0.08, 0.16)


source("Testing_Functions.R")
source("Testing_Generation.R")
source("Testing_Iterations.R")
source("Testing_Tables.R")

View(tab)
