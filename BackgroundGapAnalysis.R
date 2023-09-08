library(taxize)
library(magrittr)
library(dplyr)
library(tidyverse)
library(bold)    # API interface to BOLD
library(seqinr)  # for FASTA output
library(rentrez)
library(readxl)
library(usethis)
library(stringr)
library(Hmisc)
library(clipr)
library(knitr)
library(ggsci)


######## load data from lookup list
taxa_table <- read_excel("JanLookups.xlsx")

## Select column with taxa names  
lst<-unique(as.list(taxa_table$FinalID))

## Remove extra spaces around names 
lst<- str_squish(lst)


## Global Name Resolver for spelling mistakes 
taxa_table_gnr <-gnr_resolve(sci = as.character(lst), best_match_only = T, with_canonical_ranks=T)

## Check for taxa which got dropped by gnr 
Missing_GNR <- setdiff(taxa_table$FinalID, taxa_table_gnr$user_supplied_name)

##Check taxa dropped by GNR using Tree of Life resolver 
Missing_GNR <- tol_resolve(Missing_GNR, context_name = "All life")


##Filter successful GNR by scores. Below .5 not accurate enough  
Success_GNR <- taxa_table_gnr %>% 
  filter(score > .5)
Failed_GNR <- taxa_table_gnr %>% 
   filter( score <= .5) 


##If no GNR scores below .5, return empty df. Other wise check taxa with TOL resolver
if (dim(Failed_GNR)[1] == 0) {
  Fixed_TOL <- Missing_GNR[0,]
}else{
  Fixed_TOL <- tol_resolve(unique(Failed_GNR$submitted_name), context_name = "All life")
}

## TOL resolve on successful GNR  
taxa_table_tol <- tol_resolve(unique(Success_GNR$matched_name2), context_name = "All life")


##Check for taxa missed by TOL and with poor scores  
Missing_TOL <- taxa_table_tol %>%
  filter(is.na(unique_name) | score < .9)


##Add missing taxa from GNR to missing TOL 
Missing_TOL <- bind_rows(
  Missing_TOL,
  Missing_GNR) %>% 
  filter(is.na(unique_name) | score < .9 )

## save missing GNR that got fixed by TOL 
Missing_GNR <- Missing_GNR %>% 
  filter(!is.na(unique_name))


##Filter for successful TOL with scores over .9 
Success_TOL <- taxa_table_tol %>% 
  filter(score >= .9)

##Filter for successful TOL from failed GNR's 
Fixed_TOL<- Fixed_TOL %>% 
   filter(score >= .9)



##Combine all fixed and succesful results 
taxa_table_tol <- bind_rows(Success_TOL, Missing_GNR, Fixed_TOL)



##Remove TOL output addidition of "(class in kingdom...)" select just first word which is correct taxa 

taxa_table_tol2 <- taxa_table_tol %>%
  filter(str_detect(unique_name, " in ")) %>% 
  mutate(unique_name = word(unique_name, 1))
taxa_table_tol <- taxa_table_tol %>% 
  filter(!str_detect(unique_name, " in ")) %>%
  rbind(taxa_table_tol2)


##Final list for BOLD check
lst_tol <- unique(taxa_table_tol$unique_name)


##Compiling spelling changes and synonym changes 
SpellChangeGNR <- Success_GNR %>%
  subset(submitted_name != matched_name2) %>%
  select(submitted_name, matched_name2) %>%
  rename("OriginalName" = submitted_name,
         "SpellChange" = matched_name2)

SpellChangeTOL <- taxa_table_tol %>%
  mutate(search_string = capitalize(search_string)) %>%
  subset(search_string != unique_name) %>%
  filter(is_synonym == FALSE) %>%
  select(search_string, unique_name) %>%
  rename("OriginalName" = search_string,
         "SpellChange" = unique_name)


SpellChange <- rbind(SpellChangeGNR, SpellChangeTOL)

TOL_Synonyms <- taxa_table_tol %>%
  filter(is_synonym == TRUE) %>%
  select(search_string, unique_name) %>%
  rename("OriginalName" = search_string,
         "Synonym" = unique_name)

##Taxa missing from GNR and TOL databases
Absent_taxa <- Missing_TOL 
  


#########BOLD Downstream check
Taxa_check <- sapply(lst_tol, bold_search)
Taxa_final <- bind_rows(Taxa_check)
#write.csv(Taxa_final, "Taxa_final.csv")

## Filter for successful downstream check and filter by species ********************************
Taxa_working <- Taxa_final %>% 
  filter(!is.na(taxon)) %>% 
  filter(tax_rank == "species")


##Check for seqspec data 
bold_df_summary<-  lapply(unique(Taxa_working$taxon), function(x){
  print(x)
  mydf=bold_seqspec(x)
  bdf=mydf
  bldlength = nrow(mydf)
  if(class(mydf)=="data.frame" && colnames(mydf[1]) == "processid"){ 
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
                BOLD_Records_SW = sum(bdf$province_state %in% c("California","Arizona", "New Mexico", 
                                                                "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
                
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

##Select important columns for final dataframe 
FinalDF <- bold_df_summary %>% 
  select(FinalID:BOLD_Records_COI_SW)



##Isolate taxa missing from downstream search 
Missing_Bold_Taxa <- Taxa_final %>%
  filter(is.na(taxon)) %>%
  select(input) %>%
  rename("FinalID" = input)

##Isolate taxa with zero BOLD records form seqspec search 
Missing_seqspec_Taxa <- FinalDF %>%
  filter(BOLD_Records == 0) %>%
  select(FinalID)

##Combine missing taxa 
Missing_taxa <- rbind(Missing_Bold_Taxa, Missing_seqspec_Taxa)


## Add taxa rank to Final data frame
Taxa_rank <- Taxa_working %>%
  select(taxon, tax_rank)
colnames(Taxa_rank)[1] = "FinalID"
FinalDF = FinalDF %>%
  inner_join(Taxa_rank, by = "FinalID")
FinalDF <- FinalDF %>%
  relocate(tax_rank, .after = FinalID)


## Add taxa rank to missing taxa using GBIF
missing_lst <- as.character(as.list(Missing_taxa$FinalID))
tax_missing <- tax_rank(missing_lst, db ="gbif", rows = 1)
tax_missing <- as.data.frame(do.call(rbind, tax_missing))
tax_missing <- rownames_to_column(tax_missing, "FinalID")
tax_missing <- tax_missing %>%
  rename("tax_rank" = V1)

Missing_upstream <- rbind(classification(missing_lst, db = 'gbif', rows = 1))
Missing_final <- Missing_upstream %>%
  select( query,rank,name) %>%
  pivot_wider(names_from = rank,
              values_from = name) %>%
  select(query, phylum, class, order, family, genus) %>%
  rename("FinalID" = query) %>%
  inner_join(tax_missing, by = "FinalID") %>%
  relocate(tax_rank, .after = FinalID) %>% 
  filter(tax_rank == "species")


## Add missing taxa with ranks to final data frame
FinalDF <- FinalDF %>%
  filter(BOLD_Records != 0)
Final_Gap_analysis <- bind_rows(FinalDF, Missing_final)

##Replace NA's with 0 
Final_Gap_analysis[, 8:15][is.na(Final_Gap_analysis[, 8:15])] <- 0

# write.csv(Final_Gap_analysis, "Final_Gap_analysis.csv")

