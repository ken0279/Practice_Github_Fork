## ----message = FALSE, warning = FALSE, echo = TRUE------------------------------------------------
# Loading packages
library(tidyverse)
library(rentrez)


## ----message = FALSE, warning = FALSE, echo = TRUE------------------------------------------------

species_files <- list.files(pattern = "\\.list\\.csv$")
species_files_list <- lapply(species_files, read_csv)
species <- do.call(rbind,species_files_list)

no_species <- nrow(species) #number of fishes in the list

#Set a small dataset for test procedure
#no_species <- 20 # for test
#species <- species[1:no_species,]


## ----message = FALSE, warning = FALSE, echo = TRUE------------------------------------------------
#https://bioconnector.github.io/workshops/r-ncbi.html#introduction
#https://pediatricsurgery.hatenadiary.jp/entry/2018/01/10/205737
#entrez_db_searchable(db ="genome")

#-------------------------------------------------------------------------
#Scraping information of genome-assembly of a target species using rentrez
genome_info_sp_NCBI <- function(scientific_name){
  scientific_name_vec <- unlist(strsplit(scientific_name,"\\s+"))
  genus <- scientific_name_vec[1]
  species <- scientific_name_vec[2]
  species_lab <- paste0(genus," ",species)
  query_species_lab <- paste0(genus," ",species,"[Organism]")
  spp_genome_info <- entrez_search(db="assembly", term =query_species_lab)
  no_assembly <- length(spp_genome_info$ids)
  
  if(no_assembly > 1){
  
    spp_genome_info_list <- entrez_summary(db="assembly", id=spp_genome_info$ids)

      for(i in 1:no_assembly){
        represent_sp_genome <- NA
        
        if(spp_genome_info_list[[i]]$refseq_category=="representative genome"){
          represent_sp_genome <- spp_genome_info_list[[i]] #Primary genome-assembly
          represent_sp_genome_accession <- spp_genome_info_list[[i]]$assemblyaccession
    
          #extracting genome size info: the object "represent_sp_genome$meta" includes genome size etc.
          represent_meta <- unlist(strsplit(represent_sp_genome$meta,"> <"))
          represent_genome_size_raw <- represent_meta[which(str_detect(represent_meta, pattern="total_length"))]
          represent_genome_size <- as.numeric(str_extract_all(represent_genome_size_raw, "[0-9.]+"))
          represent_contigN50 <- represent_sp_genome$contign50
          represent_scaffoldN50 <- represent_sp_genome$scaffoldn50
        }
        
        if(!is.list(represent_sp_genome)){
          arbitary_represent_id <- max(as.numeric(spp_genome_info$ids))
          represent_sp_genome <- entrez_summary(db="assembly", id=arbitary_represent_id)
          represent_sp_genome_accession <- represent_sp_genome$assemblyaccession
          
          #extracting genome size info
          represent_meta <- unlist(strsplit(represent_sp_genome$meta,"> <"))
          represent_genome_size_raw <- represent_meta[which(str_detect(represent_meta, pattern="total_length"))]
          represent_genome_size <- as.numeric(str_extract_all(represent_genome_size_raw, "[0-9.]+"))
          represent_contigN50 <- represent_sp_genome$contign50
          represent_scaffoldN50 <- represent_sp_genome$scaffoldn50
        }
      }
    
  }else if(no_assembly == 1){
    represent_sp_genome <- entrez_summary(db="assembly", id=spp_genome_info$ids)
    represent_sp_genome_accession <- represent_sp_genome$assemblyaccession
    
    #extracting genome size info
    represent_meta <- unlist(strsplit(represent_sp_genome$meta,"> <"))
    represent_genome_size_raw <- represent_meta[which(str_detect(represent_meta, pattern="total_length"))]
    represent_genome_size <- as.numeric(str_extract_all(represent_genome_size_raw, "[0-9.]+"))
    represent_contigN50 <- represent_sp_genome$contign50
    represent_scaffoldN50 <- represent_sp_genome$scaffoldn50
    
  }else{#(no_assembly == 0)
    represent_sp_genome <- NA
    represent_sp_genome_accession <- NA
    represent_genome_size <- NA
    represent_contigN50 <- NA
    represent_scaffoldN50 <- NA

  }
  
  represent_genome_size_Mbp <- round(represent_genome_size/10^6, digits=0)
  output <- list(species_lab,
                 no_assembly,
                 represent_sp_genome,
                 represent_sp_genome_accession, 
                 represent_genome_size_Mbp,
                 represent_contigN50,
                 represent_scaffoldN50
                 )
  return(output)
  Sys.sleep(1) #
}

#-------------------------------------------------------------------------


## ----message = FALSE, warning = FALSE, echo = TRUE------------------------------------------------
scientific_name_vec <- vector()
japanese_name_vec <- vector()
target_sp_genome_assembly_exist <- vector()
target_sp_genome_size_vec <- vector() # genome size of the target species
related_spp_genome_size_vec <- vector() # mean genome size of related species


No_deposited_sp_level_genome_vec <- vector()
represent_assembly_ID_vec <- vector()
represent_assembly_status_vec <- vector()
represent_contigN50_vec <- vector()
represent_scaffoldN50_vec <- vector()
represent_genome_size_Mbp_vec <- vector()



## ----message = FALSE, warning = FALSE, echo = TRUE------------------------------------------------

for(i in 1:no_species){
  # split genus/species names
  japanese_name <- species[i,]$Japanese_name 
  scientific_name <- species[i,]$Scientific_name 

  #------------------------------------------
  #scraping represent genome using rentrez
  rentrez_scraping_out <- genome_info_sp_NCBI(scientific_name)
  
  if(rentrez_scraping_out[[2]]>=1){
    japanese_name_vec[i] <- japanese_name
    scientific_name_vec[i] <- scientific_name
    represent_assembly_ID_vec[i] <- rentrez_scraping_out[[4]]
    represent_assembly_status_vec[i] <- rentrez_scraping_out[[3]]$assemblystatus
    represent_contigN50_vec[i] <- rentrez_scraping_out[[6]]
    represent_scaffoldN50_vec[i] <- rentrez_scraping_out[[7]]
    represent_genome_size_Mbp_vec[i] <- rentrez_scraping_out[[5]]

  }else{
    japanese_name_vec[i] <- japanese_name
    scientific_name_vec[i] <- scientific_name
    represent_assembly_ID_vec[i] <- NA
    represent_assembly_status_vec[i] <- NA
    represent_contigN50_vec[i] <- NA
    represent_scaffoldN50_vec[i] <- NA
    represent_genome_size_Mbp_vec[i] <- NA
  }
  No_deposited_sp_level_genome_vec[i] <- rentrez_scraping_out[[2]]
}



species_genome <- data.frame(Japanese_name=japanese_name_vec, 
                             Scientific_name = scientific_name_vec) %>% 
  as_tibble() %>% 
  mutate(No_assembly_of_the_species = No_deposited_sp_level_genome_vec,
          represent_assembly = represent_assembly_ID_vec,
          Genome_size_of_the_species_Mbp = represent_genome_size_Mbp_vec,
          represent_assembly_status = represent_assembly_status_vec,
          Contig_N50 = represent_contigN50_vec,
          Scaffold_N50 = represent_scaffoldN50_vec
          )


## ----message = FALSE, warning = FALSE, echo = TRUE------------------------------------------------
write_csv(species_genome, "output.csv")

