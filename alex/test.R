setwd("~/git/")
library(testthat)

# Installation
install.packages("~/git/ars_2.0.tar.gz", type="source", repos = NULL)
library(ars)
test_package('ars','main')


source("ars.tar.gz")
install(ars)
system("ls")
