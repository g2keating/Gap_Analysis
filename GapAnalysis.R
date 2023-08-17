# Intertidal gap analysis

# install ------
#install.packages("taxize")
library(taxize)
library(magrittr)
library(dplyr)
library(tidyverse)
library(bold)    # API interface to BOLD
library(seqinr)  # for FASTA output
library(rentrez)
library(readxl)

options(ENTREZ_KEY = "b2ccdcd6619a29d3bf31de74e7cde9a1c209")

# import ------

bight<-read.csv("Bight taxa frequency of occurance2.csv")
pcode<-read.csv("Ref - Offshore P-code Taxa.csv")
pcode<-unique(pcode)
taxa<-merge(bight, pcode, by="Species", all.x=T, all.y=T, sort=F)
write.csv(taxa, "taxa.csv")

# Spellcheck ------

lst<-unique(as.list(taxa$Species_ed.x, taxa$Species_ed.y))

#result.long <-gnr_resolve(sci = as.character(lst), data_source_ids = c(11), with_canonical_ranks=T)
#write.csv(result.long, "taxa.spellcheck.csv")

result.long<-read.csv("taxa.spellcheck.csv", row.names=1)

no.match<-setdiff(unique(lst), unique(result.long$matched_name2))
write.csv(no.match, "taxa.nomatch.csv")

taxa<-sort(unique(result.long$matched_name2))

# BOLD downstream ----------------------

check<-bold_search(sci = "Ancistrosyllis")
bold_downstream(id = check$taxid,downto = "Species")

out<-list()
for (i in taxa) { # i<-"Ancistrosyllis"
  foo<-bold_search(sci=i)
  if(length(foo)>1) {out[i]<-foo$tax_rank}
  print(out[i])
  }

out1<-as.data.frame(t(t(out)))
write.csv(as.matrix(out1), "out1.csv")

#foo<-which(taxa %in% c("Chironomidae", "Diptera", "Insecta", "Copepoda", "Anthozoa"))
#taxa<-taxa[-foo,]

foo<-which(out1$V1 != c("species", "genus"))
out2<-as.data.frame(out1[-foo,])
out2<-as.data.frame(t(out2))

# Run BOLD check to get taxonomic rank ------------

out<-list()

for (i in sort(unique(result.long$matched_name2))) { 
  foo<-bold_seqspec(i) # i<-"Lepidozona radians"
  out[i]<-nrow(foo)
  print(paste(i, nrow(foo)))
  
  }

write.csv(as.matrix(t(t(out))), "taxa.rank.csv")

#### 

checkme<-as.data.frame(t(t(out)))
foo<-which(checkme$V1 %in% c("genus", "species"))
checkme<-checkme[foo,]
checkme<-as.data.frame(t(t(checkme)))
out2<-list()

for (i in sort(unique(rownames(checkme)))) { 
  foo<-bold_seqspec(i) # i<-"Lepidozona radians"
  out2[i]<-nrow(foo)
  print(paste(i, nrow(foo)))
  
}

write.csv(as.matrix(t(t(out2))), "taxa.rank.csv")


bold_df_summary<-  lapply(unique(result.long$matched_name2), function(x){
  print(x)
  mydf=bold_seqspec(x) 
  if(class(mydf)=="data.frame"){ 
    mydf<-mydf %>% 
      mutate(Locus = case_when(

        str_detect(marker_codes %>% tolower(),"coi")~"COI", #Mitochondrial
        str_detect(marker_codes %>% tolower(),"co1")~"COI", #Mitochondrial
        str_detect(markercode %>% tolower(),"coi")~"COI",#Mitochondrial
        str_detect(markercode %>% tolower(),"co1")~"COI",#Mitochondrial
        str_detect(markercode %>% tolower(),"coxiii")~"COI", #Mitochondrial
        str_detect(marker_codes %>% tolower(),"coxiii")~"COI", #Mitochondrial
        T~"Other"))%>%
      filter(Locus!="Other")
    
    xdf<-tibble(FinalID=x,
                BOLD_Records=nrow(mydf),
                BOLD_Records_USA=sum(mydf$country=="United States", na.rm=T),
                BOLD_Records_CA=sum(mydf$province_state=="California", na.rm=T),

      
                #COI, mitochondrial
                BOLD_Records_COI=sum(mydf$Locus=="COI", na.rm = T),
                BOLD_Records_COI_USA=sum(mydf$Locus=="COI" & mydf$country=="United States", na.rm=T),
                BOLD_Records_COI_CA=sum(mydf$Locus=="COI" &mydf$province_state=="California", na.rm=T),
                BOLD_Records_COI_SW = sum( mydf$Locus=="COI" &mydf$province_state %in% c("California","Arizona", "New Mexico", "Texas", "Nevada", "Utah", "Colorado"), na.rm=T),
                
               
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
                BOLD_Records_18S=0, BOLD_Records_18S_USA=0, BOLD_Records_18S_CA=0, BOLD_Records_18S_SW = 0)
  }
  xdf
})  %>% bind_rows()



