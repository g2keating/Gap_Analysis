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

######## load data
taxa_table <- read_excel("Biodiversity Species lookup table.xlsx")
taxa_table<-  read_csv("Bight taxa frequency of occurance2.csv")
Arthro_taxa<- read_excel("Arthro_molecular analysis taxa_08182023.xlsx")
taxa_table <- read_excel("IntertidalSummaryFinal.xlsx")


######### Check for spelling, synonyms, and accepted names using Global Names Resolver and TOL resolver
# lst<-unique(as.list(taxa_table$species_name))
lst<-unique(as.list(taxa_table$Species_ed))
lst<- str_squish(lst)


taxa_table_gnr <-gnr_resolve(sci = as.character(lst), best_match_only = T, with_canonical_ranks=T)


Missing_GNR <- setdiff(taxa_table$Species_ed, taxa_table_gnr$user_supplied_name)

Success_GNR <- taxa_table_gnr %>% 
  filter(score > .5)
Failed_GNR <- taxa_table_gnr %>% 
  filter( score <= .5) #%>% 
  #filter(matched_name2 != "Othera")

#write.csv(result.long, "intertidal.spellcheck.csv")
#no.resolve<- as.list(setdiff(taxa_table$species_name, taxa_table_gnr$user_supplied_name))
#write.csv(no.resolve, "intertidal.noresolve.csv")

Fixed_TOL <- tol_resolve(unique(Failed_GNR$submitted_name), context_name = "All life")

Success_TOL <- tol_resolve(unique(Success_GNR$matched_name2), context_name = "All life")

Missing_TOL <- tol_resolve(Missing_GNR, context_name = "All life")

Missing_TOL <- Missing_TOL %>% 
  filter(!is.na(unique_name))

Success_TOL <- Success_TOL %>% 
  filter(score >= .9) %>% 
  filter(unique_name != "Omaloplia")

Fixed_TOL<- Fixed_TOL %>% 
  filter(score >= .9)

taxa_table_tol1 <- bind_rows(Success_TOL, Fixed_TOL, Missing_TOL)

taxa_table_tol2 <- taxa_table_tol1 %>%
  subset(str_count(unique_name, '\\w+') > 3) %>% 
  mutate(unique_name = word(unique_name, 1))

taxa_table_tol <- taxa_table_tol1 %>% 
  subset(str_count(unique_name, '\\w+') <= 3) %>% 
  rbind(taxa_table_tol2)

lst_tol <- unique(taxa_table_tol$unique_name)

SpellChange1 <- Success_GNR %>% 
  subset(submitted_name != matched_name2) %>% 
  select(submitted_name, matched_name2) %>% 
  rename("OriginalName" = submitted_name,
         "SpellChange" = matched_name2)

SpellChange2 <- taxa_table_tol %>% 
  mutate(search_string = capitalize(search_string)) %>% 
  subset(search_string != unique_name) %>%
  filter(is_synonym == FALSE) %>% 
  select(search_string, unique_name) %>% 
  rename("OriginalName" = search_string,
         "SpellChange" = unique_name)


SpellChange <- rbind(SpellChange1, SpellChange2)

TOL_Synonyms <- taxa_table_tol %>%
  filter(is_synonym == TRUE) %>% 
  select(search_string, unique_name) %>% 
  rename("OriginalName" = search_string,
         "Synonym" = unique_name)


# x <- taxa_table_tol[duplicated(taxa_table_tol$unique_name), ]


#########BOLD Downstream
##Bold Check
# check<-bold_search(sci = "Aberranta")
# check<-bold_search(sci = "Jassa")
# Aeolididae
# Aeolidiidae
# 
#check<-bold_search(sci = "Platyarthrus aiasensis")
# check2<-bold_search(sci = "Ampharetidae")
# check<-bold_search(sci = "Ancistrosyllis")
#bold_downstream(id = check$taxid,downto = "Species")

Taxa_check <- sapply(lst_tol, bold_search)
Taxa_final <- bind_rows(Taxa_check)
write.csv(Taxa_final, "Taxa_final.csv")

# Taxa_working <- Taxa_final %>% 
#   filter(!is.na(taxon)) #%>% 
#filter(tax_rank == "genus" | tax_rank == "species")
# Missing_Bold_Taxa <- Taxa_final %>% 
#   filter(is.na(taxon)) %>% 
#select(input)

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


# write_csv(Taxa_final, "arthro_taxa_final.csv")
# arthro_taxa_final <- read_csv("arthro_taxa_final.csv")

Taxa_working <- Taxa_final %>% 
  filter(!is.na(taxon))


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


FinalDF <- bold_df_summary %>% 
  select(FinalID:BOLD_Records_COI_SW)

write_csv(FinalDF, "Arthro_Rafi.csv")


Missing_Bold_Taxa <- Taxa_final %>% 
  filter(is.na(taxon)) %>% 
  select(input) %>% 
  rename("FinalID" = input)

Missing_seqspec_Taxa <- Arthro_Rafi %>% 
  filter(BOLD_Records == 0) %>% 
  select(FinalID)

Missing_taxa <- rbind(Missing_Bold_Taxa, Missing_seqspec_Taxa)



Taxa_rank <- Taxa_working %>% 
  select(taxon, tax_rank)
colnames(Taxa_rank)[1] = "FinalID"

Arthro_Rafi = Arthro_Rafi %>% 
  inner_join(Taxa_rank, by = "FinalID")
Arthro_Rafi <- Arthro_Rafi %>% 
  relocate(tax_rank, .after = FinalID)


missing_lst <- as.character(as.list(Missing_taxa$FinalID))


tax_missing <- tax_rank(missing_lst, db ="gbif")
tax_missing <- as.data.frame(do.call(rbind, tax_missing))
tax_missing <- rownames_to_column(tax_missing, "FinalID")
tax_missing <- tax_missing %>% 
  rename("tax_rank" = V1)


Missing_upstream <- rbind(classification(missing_lst, db = 'gbif'))

Missing_final <- Missing_upstream %>% 
  select( query,rank,name) %>% 
  pivot_wider(names_from = rank,
              values_from = name) %>% 
  select(query, phylum, class, order, family, genus) %>% 
  rename("FinalID" = query) %>% 
  inner_join(tax_missing, by = "FinalID") %>% 
  relocate(tax_rank, .after = FinalID)

Arthro_Rafi <- Arthro_Rafi %>% 
  filter(BOLD_Records != 0)

Arthro_Gap_analysis <- bind_rows(Arthro_Rafi, Missing_final) 
Arthro_Gap_analysis[, 8:15][is.na(Arthro_Gap_analysis[, 8:15])] <- 0

write_csv(Arthro_Gap_analysis, "Arthro_Gap_analysis.csv")
write_csv(SpellChange, "SpellChange.csv")
write_csv(TOL_Synonyms, "Synonyms.csv")

