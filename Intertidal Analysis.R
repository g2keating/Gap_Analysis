library(taxize)
library(magrittr)
library(dplyr)
library(tidyverse)
library(bold)    # API interface to BOLD
library(seqinr)  # for FASTA output
library(rentrez)
library(readxl)
library(usethis)


######## load data
taxa_table <- read_excel("Biodiversity Species lookup table.xlsx")
#taxa_table2<-  read_csv("Bight taxa frequency of occurance2.csv")


######### Check for spelling, synonyms, and accepted names using Global Names Resolver and TOL resolver
lst<-unique(as.list(taxa_table$species_name))
taxa_table_gnr <-gnr_resolve(sci = as.character(lst), data_source_ids = c(11), with_canonical_ranks=T)
Success_GNR <- taxa_table_gnr %>% 
  filter(score > .5)
#write.csv(result.long, "intertidal.spellcheck.csv")
#no.resolve<- as.list(setdiff(taxa_table$species_name, taxa_table_gnr$user_supplied_name))
#write.csv(no.resolve, "intertidal.noresolve.csv")
Failed_GNR <- taxa_table_gnr %>% 
  filter( score <= .5)
Fixed_GNR <- tol_resolve(unique(Failed_GNR$submitted_name))
Success_TOL <- tol_resolve(unique(Success_GNR$matched_name2))
taxa_table_tol <- bind_rows(Success_TOL, Fixed_GNR)
lst_tol <- unique(taxa_table_tol$unique_name)



#########BOLD Downstream
##Bold Check
# check<-bold_search(sci = "Aberranta")
# check<-bold_search(sci = "Jassa")
# Aeolididae
# Aeolidiidae
# 
# check<-bold_search(sci = "Amphiporus californicus")
# check2<-bold_search(sci = "Ampharetidae")
# check<-bold_search(sci = "Ancistrosyllis")
#bold_downstream(id = check$taxid,downto = "Species")

Taxa_check <- sapply(lst_tol, bold_search)
Taxa_final <- bind_rows(Taxa_check)
Taxa_working <- Taxa_final %>% 
  filter(!is.na(taxon)) %>% 
  filter(tax_rank == "genus" | tax_rank == "species")
Missing_Taxa <- Taxa_final %>% 
  filter(is.na(taxon)) %>% 
  select(input)

# Genus_fix <- Missing_Taxa %>% 
#   filter(grepl("genus in", Missing_Taxa$input))
# Genus_fix<- word(Genus_fix$input, 1)
# Genus_check <- sapply(Genus_fix, bold_search)
# Genus_check <- bind_rows(Genus_check)
# write.csv(Taxa_final, "out1.csv")





############Bold seqspec
# check4<-bold_seqspec("Eteone")
# check2<-bold_seqspec("Psychoda")
# check2<-bold_seqspec("Jassa", cleanData = T)
# check3<-bold_seqspec("Jassa")
# bold_stats("Jassa")
# check<-bold_seqspec("Bryopsis")
# check<-bold_seqspec("Aglaophenia latirostris")
# check3<-bold_seqspec("Blidingia subsalsa")
# Taxa_working <- Taxa_working[1:60,]

bold_df_summary<-  lapply(unique(Taxa_working$taxon), function(x){
  print(x)
  mydf=bold_seqspec(x)
  bdf=mydf
  bldlength = nrow(mydf)
  if(class(mydf)=="data.frame"){ 
    mydf<- mydf %>% 
      mutate(Locus = case_when(
        
        str_detect(marker_codes %>% tolower(),"coi")~"COI", #Mitochondrial
        str_detect(marker_codes %>% tolower(),"co1")~"COI", #Mitochondrial
        str_detect(markercode %>% tolower(),"coi")~"COI",#Mitochondrial
        str_detect(markercode %>% tolower(),"co1")~"COI",#Mitochondrial
        str_detect(markercode %>% tolower(),"coxiii")~"COI", #Mitochondrial
        str_detect(marker_codes %>% tolower(),"coxiii")~"COI", #Mitochondrial
        T~"Other"))
    
    xdf<-tibble(FinalID=x,
                phylum = first(mydf$phylum_name),
                class = first(mydf$class_name),
                order = first(mydf$order_name),
                family = first(mydf$family_name),
                genus = first(mydf$genus_name),
                BOLD_Records= bldlength,
                BOLD_Records_USA=sum(bdf$country=="United States", na.rm=T),
                BOLD_Records_CA=sum(bdf$province_state=="California", na.rm=T),
                
                
                #COI, mitochondrial
                BOLD_Records_COI=sum(mydf$Locus=="COI", na.rm = T),
                BOLD_Records_COI_USA=sum(mydf$Locus=="COI" & mydf$country=="United States", na.rm=T),
                BOLD_Records_COI_CA=sum(mydf$Locus=="COI" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_COI_SW = sum( mydf$Locus=="COI" &mydf$province_state %in% c("California","Arizona", "New Mexico", 
                                                                                         "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
    )
  }
  else  {
    xdf= tibble(FinalID=x,  
                BOLD_Records=0, BOLD_Records_USA=0, BOLD_Records_CA=0, BOLD_Records_SW = 0,
                BOLD_Records_ITS=0, BOLD_Records_ITS_USA=0, BOLD_Records_ITS_CA=0, BOLD_Records_ITS_SW = 0,
                BOLD_Records_COI=0, BOLD_Records_COI_USA=0, BOLD_Records_COI_CA=0, BOLD_Records_COI_SW = 0,
                BOLD_Records_CYTB=0, BOLD_Records_CYTB_USA=0, BOLD_Records_CYTB_CA=0, BOLD_Records_CYTB_SW = 0,
                BOLD_Records_12S=0, BOLD_Records_12S_USA=0, BOLD_Records_12S_CA=0, BOLD_Records_12S_SW = 0,
                BOLD_Records_16S=0, BOLD_Records_16S_USA=0, BOLD_Records_16S_CA=0, BOLD_Records_16S_SW = 0,
                BOLD_Records_18S=0, BOLD_Records_18S_USA=0, BOLD_Records_18S_CA=0, BOLD_Records_18S_SW = 0,
                phylum = NA,
                class = NA,
                order = NA,
                family = NA,
                genus = NA)
  }
  xdf
})  %>% bind_rows()


FinalDF <- bold_df_summary %>% 
  select(FinalID:BOLD_Records_COI_SW)

#write_csv(FinalDF, "BOLD_Species_Check_MARINeUCSC.csv")

