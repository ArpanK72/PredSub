#Set the directory to be this Real_Data_Analysis folder

d     <- 10
B     <- 200
alist <- seq(2.75, 3.25, by = 0.25)
olist <- tail(alist, 2)


source("RDA_Functions.R")
source("DBLP_Preprocessing.R")
source("RDA_H0.R")
source("RDA_H1.R")
source("RDA_Tables.R")

View(tab_rda)
