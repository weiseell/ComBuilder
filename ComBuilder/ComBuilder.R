ComBuilder <- function(area,community,rank,loci,setmax = 100, add_list = F, wiki_search = T,group = T, acc_only = F, summ_only = F){
  #load dependancies
  require(rgdal)
  require(rgbif)
  require(sf)
  require(dplyr)
  require(tidyr)
  require(rentrez)
  #make search term for gbif
  key = name_backbone(name = community)$usageKey
  
  #make a grid to put over area to cut the size into queryable by gbif
  startgrid = st_make_grid(spobj, n = c(10,10))   #Can adjust the dimensions of the grid -- this makes a 10 x 10 grid, probably unnecessary to have this many cells
  spgrid = as_Spatial(startgrid)
  
  #overlay grid over area
  plot(spobj)
  plot(spgrid, add = TRUE)
  
  #shape grid to area
  outgrid = startgrid[!is.na(over(spgrid,spobj))]
  plot(spobj)
  plot(outgrid, add = TRUE)
  
  #make grid gbif compatable
  mygrid = format(outgrid, width = Inf)
  
  #input grid to gbif in a for loop
  occ = NULL
  i <- NULL
  print("Collecting gbif occurance data")
  for(i in (1:length(outgrid))){
    #search for occurances in gbif
    res = occ_search(taxonKey = key, geometry = mygrid[i], hasCoordinate = TRUE, limit = 200000)
    #grab taxonomy info for the occurances and save into a bigger list
    data <- as.data.frame(res$data)
    #select taxonomy
    data <- data %>% select(kingdom,phylum,order,family,genus,species,name,collectionCode,basisOfRecord)
    #getting rid of data that isn't contemporary data (historical data not accurate to the current community)
    data <- subset(data,data$basisOfRecord != "PRESERVED_SPECIMEN")
    data <- subset(data,data$basisOfRecord != "UNKNOWN")
    data <- subset(data,data$basisOfRecord != "FOSSIL_SPECIMEN")
    occ <- rbind(occ,data)
  }
  print("Done collecting gbif occurance data")
  #condense occurance data to a unique species list
  df <- occ %>% distinct()
  out1 <- df %>% mutate(sci.name = paste(Genus,Species,sep = " "))
  
  #wikipedia taxonomy check
  if(wiki_check == "T"){
    print("Performing wikipedia taxonomy check")
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
  
  #create a species list to query NCBI with
  species_list <- out1 %>% select(sci.name)
  
  require(rentrez)
  
  #make an output folder and a sequence folder to put all of the sequences in
  dir.create("Output/")
  dir.create("Output/seqs")
  
  #add a count column to species_list
  species_list$count <- 0
  
  #loop to get accession numbers and NCBI sequences for each species
  i <- NULL
  for(i in 1:nrow(df)){
    print(i)
    
    #creating the search term for NCBI for each sci.name
    my.term <- paste0("(",species_list$sci.name[i],"[Organism]) AND",loci,"[Title]")
    my.term
    
    #doing the serach
    x <- entrez_search(db="nucleotide",
                       term=my.term,
                       retmax=setmax)
    #get the length of x to see the number of seqs per search
    length(x$ids)
    
    #counting the IDS
    species_list$locus.count[i] <- length(x$ids)
    x$ids
    
    if(acc_only == F){
      if(length(x$ids) > 0){
        #getting the ids and making a file name
        my.ids <- x$ids
        
        #getting all the fasta files
        my.fetch <- entrez_fetch(db = 'nucleotide',id=my.ids,rettype = "fasta")
        
        my.fetch = strsplit(my.fetch, split = "\n\n")[[1]]
        
        my.fetch <- gsub(pattern = "$\n",replacement = "",x = my.fetch)
        #putting the sequences for the species into the seq folder 
        write.table(my.fetch, file = paste0("Output/seqs/",df$sci.name[i],"_",loci,".fasta"), quote = F, row.names = F, col.names = F)
      }
    }
    if(acc_only == T){
      if(length(x$ids > 0)){
        out1 <- as.data.frame(x$ids)
        colnames(out1) <- "ids"
        zeros <- as.data.frame(rep(0,(100-length(x$ids))))
        colnames(zeros) <- "ids"
        out2 <- rbind(out1,zeros)
        colnames(out2) <- df$sci.name[i]
        out <- cbind(out,out2)
      }
    }
  } #end of per species for loop
  
  #compile .fasta files into one large files 
  if(group == T){
    my.file <- list.files("Output/seqs/", pattern = loci)
    #loop to condense all fasta files into one file
    i = 1
    for(i in 1:length(my.file)) {
      print(i)
      temp <- readLines(paste0("Output/seqs/", my.file[i]))
      cat(temp,file = paste0("Output/allseqs.fasta"),sep = "\n",append = TRUE)
    }
  }
  
  #!# need to combine gbif taxonomy file with NCBI counts
  write.table(x = df,file = "Output/community.summary.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
}