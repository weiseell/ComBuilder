
if(wiki_check == "T"){
  #additional packages required for this section
  require(tidyverse)
  #additional function required for this section
  cleanFun <- function(htmlString) {
    return(gsub("<.*?>", "", htmlString))
  }
  
  #set j and output for the loop
  j <- NULL
  output <- NULL
  for(j in 1:nrow(df)){
    #print j to keep track of where the loop is for troubleshooting
    print(j)
    #reading in the data from the html webpage
    wp <- readLines(paste0('http://en.wikipedia.org/wiki/Special:Search/',df$species[j]))
    #getting lines of hmtl code that has both scientific and common names in them and 
    #cleaning up code so taxonomy of speceis can be extracted
    sn <- wp[grepl(pattern = "<td>",x = wp)]
    sn
    sn <- cleanFun(htmlString = sn)
    sn <- matrix(sn,ncol = 2,byrow = T)
    sn <- as.data.frame(sn)
    
    sn$V1 <- gsub(pattern = ":",replacement = "",x = sn$V1)
    sn$V2 <- gsub(pattern = "^[^;]*;",replacement = "",x = sn$V2)
    sn
    #making a data frame with all taxonomy levels
    taxon.ranks <- data.frame(V1 = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), V2=NA)
    
    #!# not sure what this does yet, need to go line by line
    k<-1
    for(k in 1:nrow(taxon.ranks)){
      if(length(as.character(sn[sn$V1 == taxon.ranks$V1[k],2]))==1){
        taxon.ranks$V2[k] <- as.character(sn[sn$V1 == taxon.ranks$V1[k],2])
      } #end of if statement
    } # end of k for loop
    
    #manipulating occurance data so it can be compared to wikipedia data species
    sn1 <- data.frame(common.name=occ$com.name[i])
    sn1 <- gather(data = sn1,key = V1,V2)
    sn2 <- data.frame(sci.name=occ$species[i])
    sn2 <- gather(data = sn2,key = V1,V2)
    sn <- rbind(taxon.ranks,sn1,sn2)

    #now making a df of the wide form
    sn <- spread(data = sn,key = V1,value = V2)
    sn <- sn %>% select(Kingdom,Phylum,Class,Order,Family,Genus,Species)
    output <- rbind(output,sn)  
  } #end of j loop
  
  #checking that all of the species that were selected by gbif are actually
  #a part of the selected taxonomy
  out1 <- subset(output,output$Phylum == taxon)
  out1 <- out1 %>% mutate(sci.name = paste(Genus,Species,sep = " "))
}