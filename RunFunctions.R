#Script to test the Combuilder and FindNCBI functions
#Example will use gbif to build a database of all decapoda species in the Great Lakes
#basin, and query NCBI for 16S species for that list (either through ComBuilder or FindNCBI)
#Input folder contains the shape files required for the Great Lakes Database
#and a community file if running NCBI

#library required to read in shapefile
library(maptools)

#load in functions
source("ComBuilder.R")
source("Find_NCBI.R")

#read in shape file for ComBuilder
area = readShapePoly("Input/greatlakes_subbasins/greatlakes_subbasins")

#run ComBuilder
ComBuilder(area = area,community = "Decapoda", rank = "Order", loci = "16S",setmax = 10)

#read in species list for FindNCBI
sl <- read.table("Input/example.community.txt", header = T, stringsAsFactors = F)
sl <- data.frame(sci.name = paste(sl$sample_number,sl$sci.name,sep = " "))

#run FindNCBI
FindNCBI(species_list = sl, loci = "16S", setmax = 10)
