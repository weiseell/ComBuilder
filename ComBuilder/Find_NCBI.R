FindNCBI <- function(species_list,loci,setmax = 100, group = T, acc_only = F, summ_only = F){
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
      df <- readLines(paste0("Output/seqs/", my.file[i]))
      head(df)
      cat(df,file = paste0("Output/allseqs.fasta"),sep = "\n",append = TRUE)
    }
  }
  #write summary file
  write.table(x = df,file = "Output/community.summary.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
}