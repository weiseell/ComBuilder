#Script to test the Combuilder and FindNCBI functions

#libraries
library(maptools)

#load in functions
source("ComBuilder.R")
source("Find_NCBI.R")
#read in shape file
area = readShapePoly("Input/greatlakes_subbasins/greatlakes_subbasins")

ComBuilder(area = area,community = "Decapoda", rank = "Order", loci = "16S",setmax = 10)

sl <- read.table("Input/example.community.txt", header = T, stringsAsFactors = F)
sl <- data.frame(sci.name = paste(sl$sample_number,sl$sci.name,sep = " "))
FindNCBI(species_list = sl, loci = "16S", setmax = 10)
