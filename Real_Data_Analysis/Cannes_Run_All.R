# ---- Parameters ----
d     <- 10
B     <- 200
alist <- seq(2.75, 3.25, by = 0.25)
olist <- tail(alist, 2)

# ---- Pipeline ----
source("RDA_Functions.R")
source("Cannes_Preprocessing.R")
source("RDA_H0.R")
source("RDA_H1.R")
source("RDA_Tables.R")

View(tab_rda)