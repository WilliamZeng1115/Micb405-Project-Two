#install.packages("tibble")
#install.packages("readr")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ggplot2")
#install.packages("stringr")
#install.packages("purrr")
#install.packages("forcats")
#install.packages("cowplot")
#install.packages("gvlma")
#source("https://bioconductor.org/biocLite.R")

#biocLite("BiocInstaller")
#biocLite("pathview")
#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("pathview", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("gage", version = "3.8")

#Load Packages:
library(plyr)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(forcats)

library(cowplot)
library(pathview)

library(readr)

library(RColorBrewer)
library(knitr)

##############################################
#Graph Palette
pal<-c("#3874ce", 
       "#2ed0bd", 
       "#54f2dc", 
       "#59d15e", 
       "#93ef34", 
       "#f5f132",
       "#eda943")

##############################################
#Cruise ID "SI0XX"
Cruises <- c(42,48,72,73,74,75)

##############################################
#MAGs/Bin #
apb_id = c("155","192","120","65","84","262","247","293")
#33,235,378 are no good, did not pass quality check

##############################################
#Corresponding MAGS to Alphaproteobacteria Classes
#Rhizobiales = c(155,192)
#Caulobacterales = c(120)
#Rhodobacterales = c(65,84,262,247)
#Sphinogomonadales = c(293)

#################################################################
#Data Clean Up for Transcriptome RPKM (Aligned to MAG contigs)
#################################################################

#Load all the RPKM ".csv" for transcriptome data aligned to contigs
setwd("~/Desktop/MICB_405/Project2/RPKM")
for(i in 1:length(apb_id)){
  for(j in 1:length(Cruises)){
    c= as.character(Cruises[j])
    a= as.character(apb_id[i])
    print(a)
    assign((paste("SI0",c,"_",a,sep="")), read_csv(paste(a,
                                                         "_SI0",
                                                         c,
                                                         "_MAG_ORFs_RPKM.csv",
                                                         sep=""), 
                                                   col_names = FALSE,
                                                   cols(X1 = col_character(),
                                                        X2 = col_double())))
  }
}

#Clean up RPKMXXX/RPKMXX where XX = MAG ID/bin ID, merge all the cruises of each MAG to the same dataframe
for(i in 1:length(apb_id)){
  apb=as.character(apb_id[i])
  assign(paste("RPKM",apb, sep=""),eval(parse(text =  paste("SI042_",apb,sep=""))))
  
  for(j in 2:length(Cruises)){
    c= as.character(Cruises[j])
    temp=paste("SI0",c,"_",apb, sep="")
    assign(paste("RPKM",apb, sep=""),
           merge(eval(parse(text = paste("RPKM",apb, sep=""))),
                 eval(parse(text = temp)), 
                 by.x="X1", 
                 by.y="X1"))
  }
}

#Combine RPKM values of different MAGs that assigned to the same Class of Alphaproteobacteria
Rhi <- rbind (RPKM155,RPKM192)
Cau <- RPKM120
Rho <- rbind(RPKM65,RPKM84,RPKM262,RPKM247)
Sph <- RPKM293
#Set column names for RPKM values
colnames(Rhi)<- c("Prokka_ID","SI042","SI048","SI072","SI073","SI074","SI075")
colnames(Cau)<- c("Prokka_ID","SI042","SI048","SI072","SI073","SI074","SI075")
colnames(Rho)<- c("Prokka_ID","SI042","SI048","SI072","SI073","SI074","SI075")
colnames(Sph)<- c("Prokka_ID","SI042","SI048","SI072","SI073","SI074","SI075")
setwd("~")

#Load KEGG ID that have corresponding Locus Tag ID, and assign colnames
KEGG_ID_LIST <- read_tsv("Desktop/MICB_405/Project2/KEGG/SaanichInlet_MAGx_ORFs_ko.cleaned.tsv", 
                         col_names = FALSE)
colnames(KEGG_ID_LIST) = c('MAG','KEGG')

#testing codes:
#rawRPKM <- as.tibble(cbind(as.numeric(SI042$X2),as.numeric(SI048$X2),SI072$X2,SI073$X2,SI074$X2,SI075$X2))
#colnames(rawRPKM) <- c('SI042','SI048','SI072','SI073','SI074','SI075')
#rawRPKM$MAG <- SI042$X1
#RPKM$KEGG <- c("NA")

#Filter for the Locus Tag ID of RPKM with corresponding KEGG ID
Rhi <- dplyr::filter(Rhi, Rhi$Prokka_ID %in% KEGG_ID_LIST$MAG )
Cau <- dplyr::filter(Cau, Cau$Prokka_ID %in% KEGG_ID_LIST$MAG )
Rho <- dplyr::filter(Rho, Rho$Prokka_ID %in% KEGG_ID_LIST$MAG )
Sph <- dplyr::filter(Sph, Sph$Prokka_ID %in% KEGG_ID_LIST$MAG )
#Write KEGG ID (corresponding to Locus Tag ID) to each row of RPKM values
for (x in 1:length(Rhi$Prokka_ID)) {
  Rhi$KEGG[x]<-KEGG_ID_LIST$KEGG[Rhi$Prokka_ID[x] == KEGG_ID_LIST$MAG]
}
for (x in 1:length(Cau$Prokka_ID)) {
  Cau$KEGG[x]<-KEGG_ID_LIST$KEGG[Cau$Prokka_ID[x] == KEGG_ID_LIST$MAG]
}
for (x in 1:length(Rho$Prokka_ID)) {
  Rho$KEGG[x]<-KEGG_ID_LIST$KEGG[Rho$Prokka_ID[x] == KEGG_ID_LIST$MAG]
}
for (x in 1:length(Sph$Prokka_ID)) {
  Sph$KEGG[x]<-KEGG_ID_LIST$KEGG[Sph$Prokka_ID[x] == KEGG_ID_LIST$MAG]
}
#Sum all the RPKM corresponding to the same KEGG ID in each Cruise of each Class 
Rhi<-Rhi %>%
  group_by(KEGG)%>%
  summarise(SI042=sum(SI042, na.rm=TRUE),
            SI048=sum(SI048, na.rm=TRUE),
            SI072=sum(SI072, na.rm=TRUE),
            SI073=sum(SI073, na.rm=TRUE),
            SI074=sum(SI074, na.rm=TRUE),
            SI075=sum(SI075, na.rm=TRUE))
Cau<-Cau %>%
  group_by(KEGG)%>%
  summarise(SI042=sum(SI042, na.rm=TRUE),
            SI048=sum(SI048, na.rm=TRUE),
            SI072=sum(SI072, na.rm=TRUE),
            SI073=sum(SI073, na.rm=TRUE),
            SI074=sum(SI074, na.rm=TRUE),
            SI075=sum(SI075, na.rm=TRUE))
Rho<-Rho %>%
  group_by(KEGG)%>%
  summarise(SI042=sum(SI042, na.rm=TRUE),
            SI048=sum(SI048, na.rm=TRUE),
            SI072=sum(SI072, na.rm=TRUE),
            SI073=sum(SI073, na.rm=TRUE),
            SI074=sum(SI074, na.rm=TRUE),
            SI075=sum(SI075, na.rm=TRUE))
Sph<-Sph %>%
  group_by(KEGG)%>%
  summarise(SI042=sum(SI042, na.rm=TRUE),
            SI048=sum(SI048, na.rm=TRUE),
            SI072=sum(SI072, na.rm=TRUE),
            SI073=sum(SI073, na.rm=TRUE),
            SI074=sum(SI074, na.rm=TRUE),
            SI075=sum(SI075, na.rm=TRUE))

#Filter out low RPKM (If RPKM is less than average of 1 RPKM per Cruise per Class, then it will be removed)
Rhi_6 <- dplyr::filter(Rhi, rowSums(Rhi[,2:7])>6)
Rhi_6<-as.matrix(Rhi_6[,2:7])
rownames(Rhi_6) <- (dplyr::filter(Rhi, rowSums(Rhi[,2:7])>6))$KEGG
colnames(Rhi_6) <- colnames(Rhi[,2:7])

Cau_6 <- dplyr::filter(Cau, rowSums(Cau[,2:7])>6)
Cau_6<-as.matrix(Cau_6[,2:7])
rownames(Cau_6) <- (dplyr::filter(Cau, rowSums(Cau[,2:7])>6))$KEGG
colnames(Cau_6) <- colnames(Cau[,2:7])

Rho_6 <- dplyr::filter(Rho, rowSums(Rho[,2:7])>6)
Rho_6<-as.matrix(Rho_6[,2:7])
rownames(Rho_6) <- (dplyr::filter(Rho, rowSums(Rho[,2:7])>6))$KEGG
colnames(Rho_6) <- colnames(Rho[,2:7])

Sph_6 <- dplyr::filter(Sph, rowSums(Sph[,2:7])>6)
Sph_6<-as.matrix(Sph_6[,2:7])
rownames(Sph_6) <- (dplyr::filter(Sph, rowSums(Sph[,2:7])>6))$KEGG
colnames(Sph_6) <- colnames(Sph[,2:7])

#Write a pathview for certain pathway.ids:
#00910 for Nitrogen
#01100 for metabolic pathway
#00920 for Sulphur
#dynamically (manually) changed the dataset and out.suffix to generate all pathways for all Classes
#multi-layer
#fancy colors
#plot to same layer

setwd("~/Desktop/MICB_405/Project2/PathView")
pv.out <- pathview(gene.data = Rhi_6[,1:6], 
                   pathway.id = "01100", 
                   species="ko",  
                   out.suffix = "Rhi.data",
                   multi.state=T,
                   same.layer =T,
                   limit = list(gene = c(0,50)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"))

str(pv.out) #show information of pathview
pv.out$plot.data.gene #show genes used

#################################################################
#Chemical Cycle Bubble Plots
#################################################################
#Read data downloaded from KEGG pathway for Sulfur/Nitrogen, clean out special characters, changed strings to upper case
KO_log <- read_table2("Desktop/MICB_405/Project2/Nitrogen:Sulfur/KO_log.csv", 
                      col_names = TRUE) %>%
  .[,-2]%>%
  .[,-(4:5)] %>%
  mutate(Nitrogen_1 = str_replace_all(Nitrogen_1, ",", "")) %>%
  mutate(Nitrogen_1 = str_replace_all(Nitrogen_1, ";", "")) %>%
  mutate(Nitrogen_1 = toupper(Nitrogen_1))

#Combine all RPKM by KEGG ID of all Class to a single dataset, rename column names
mega_cycle<-merge(Rho, Rhi, by.x="KEGG",by.y="KEGG", all=T)%>%
  merge(., Sph, by.x="KEGG",by.y="KEGG",all=T) %>%
  merge(., Cau, by.x="KEGG",by.y="KEGG",all=T) 
colnames(mega_cycle) <- c("KEGG",
                          "Rho SI042","Rho SI048","Rho SI072","Rho SI073","Rho SI074","Rho SI075",
                          "Rhi SI042","Rhi SI048","Rhi SI072","Rhi SI073","Rhi SI074","Rhi SI075",
                          "Sph SI042","Sph SI048","Sph SI072","Sph SI073","Sph SI074","Sph SI075",
                          "Cau SI042","Cau SI048","Cau SI072","Cau SI073","Cau SI074","Cau SI075")
#Separate Cruise to Class and Cruise ID in separate columns, filtered to keep only data with KEGG ID in Nitrogen Metabolism 00910
N_cycle_cruise <- gather(mega_cycle, key="Cruise",value="RPKM",-KEGG) %>%
  mutate(Class = substring(Cruise, 1,3)) %>%
  mutate(Cruise=substring(Cruise,5)) %>%
  dplyr::filter( KEGG %in% (KO_log$`00910`[1:64]))
#Generate new Gene column, prefilled with NA
N_cycle_cruise$Gene = "NA"
#Add Gene ID to dataset corresponding to KEGG ID
for (x in 1:length(N_cycle_cruise$KEGG)) {
  N_cycle_cruise$Gene[x]=KO_log$Nitrogen_1[N_cycle_cruise$KEGG[x]==(KO_log$`00910`[1:64])]
}
#Plot Genes with RPKM bubble plot of Nitrogen Metabolism Pathway 
N_cycle_cruise_dotplot <- N_cycle_cruise %>%
  ggplot (aes(x=Gene, y= Cruise, color = Class, size = RPKM)) +
  geom_point () +
  scale_x_discrete(name ="Expressed Genes")+
  scale_y_discrete(name ="Cruises")+
  facet_wrap(~Class, scales="free_y", nrow = 4) +
  theme(legend.text = element_text(size=8), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position ="none") 

#Separate Cruise to Class and Cruise ID in separate columns, filtered to keep only data with KEGG ID in Sulfur Metabolism 00920
S_cycle_cruise <- gather(mega_cycle, key="Cruise",value="RPKM",-KEGG) %>%
  mutate(Class = substring(Cruise, 1,3)) %>%
  mutate(Cruise=substring(Cruise,5)) %>%
  dplyr::filter( KEGG %in% (KO_log$`00910`[1:64]))
#Generate new Gene column, prefilled with NA
S_cycle_cruise$Gene = "NA"
#Add Gene ID to dataset corresponding to KEGG ID
for (x in 1:length(S_cycle_cruise$KEGG)) {
  S_cycle_cruise$Gene[x]=KO_log$Nitrogen_1[S_cycle_cruise$KEGG[x]==(KO_log$`00910`[1:64])]
}
#Plot Genes with RPKM bubble plot of Sulfur Metabolism Pathway 
S_cycle_cruise_dotplot <- S_cycle_cruise %>%
  ggplot (aes(x=Gene, y= Cruise, color = Class, size = RPKM)) +
  geom_point () +
  scale_x_discrete(name ="Expressed Genes")+
  scale_y_discrete(name ="Cruises")+
  facet_wrap(~Class, scales="free_y", nrow = 4) +
  theme(legend.text = element_text(size=8), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "right") 

#Combine bubble plots of Nitrogen and Sulfur Metabolism Pathway
CruisevsGenes<-plot_grid(N_cycle_cruise_dotplot, S_cycle_cruise_dotplot, labels=c("A", "B"), 
                         align="h", axis="tb", rel_widths=c(2/5, 3/5))

CruisevsGenes

#testing:
#mega_cycle$Rho <- rowSums(mega_cycle[2:7])
#mega_cycle$Rhi <- rowSums(mega_cycle[8:13])
#mega_cycle$Sph <- rowSums(mega_cycle[14:19])
#mega_cycle$Cau <- rowSums(mega_cycle[20:25])
# 
# mega_cycle <- mega_cycle %>%
#   .[,-(2:25)] %>%
#   dplyr::select(KEGG,Rho,Rhi,Sph,Cau)
##############################################################
#Sulphur Cycle Bubble Plot
# 
# N_cycle <- mega_cycle %>%
#   dplyr::filter( KEGG %in% (KO_log$`00910`[1:64]))
# N_cycle$Gene = "NA"
# 
# for (x in 1:length(N_cycle$KEGG)) {
#   N_cycle$Gene[x]=KO_log$Nitrogen_1[N_cycle$KEGG[x]==(KO_log$`00910`[1:64])]
# }
# 
# N_cycle %>%
#   dplyr::mutate_if(~ any(is.na(.)),~ if_else(is.na(.),0,.)) %>%
#   filter(rowSums(.[2:5])>0) %>%
#   gather(key="Class",value="RPKM",-KEGG,-Gene) %>%
#   
#   ggplot (aes(x=Gene, y=Class, color = Class, size = RPKM)) +
#   geom_point () +
#   scale_x_discrete(name ="Expressed Genes")+
#   scale_y_discrete(name ="Class")+
#   theme(legend.text = element_text(size=8), 
#         axis.text.x = element_text(angle = 90, hjust = 1)) 
# ##############################################################
# #Sulphur Cycle Bubble Plot
# S_cycle <- mega_cycle %>%
#   dplyr::filter( KEGG %in% (KO_log$`00910`[65:166]))
# S_cycle$Gene = "NA"
# 
#  for (x in 1:length(S_cycle$KEGG)) {
#    S_cycle$Gene[x]=KO_log$Nitrogen_1[S_cycle$KEGG[x]==(KO_log$`00910`[65:166])]
#  }
#  
# S_cycle %>%
#    dplyr::mutate_if(~ any(is.na(.)),~ if_else(is.na(.),0,.)) %>%
#    filter(rowSums(.[2:5])>0) %>%
#    gather(key="Class",value="RPKM",-KEGG,-Gene) %>%
#    
#    ggplot (aes(x=Gene, y=Class, color = Class, size = RPKM)) +
#    geom_point () +
#    scale_x_discrete(name ="Expressed Genes")+
#    scale_y_discrete(name ="Class")+
#    theme(legend.text = element_text(size=8), 
#          axis.text.x = element_text(angle = 90, hjust = 1)) 

##############################################
#ggplot of Chemical Data, clean up
##############################################
rawChemDat <- read_csv("Desktop/MICB_405/Project2/Saanich_TimeSeries_Chemical_DATA.csv.txt", 
                       col_names = TRUE)
#Clean up NA, NAN, ND in Chemical Data, filtered for 10 m depth
ChemDat <- 
  rawChemDat %>%
  # filter(Cruise %in% Cruises) %>%
  filter(Depth == 10) %>%
  dplyr::select(Cruise, Date, CTD_O2, PO4, NO3, SI, Mean_NH4,Mean_N2O,Mean_NO2,Mean_CH4) %>%
  filter((!is.na(PO4)) & ("NAN" != PO4) & (PO4!="ND"))%>%
  filter((!is.na(NO3)) & ("NAN" != NO3) & (NO3!="ND"))%>%
  filter((!is.na(CTD_O2)) & ("NAN" != CTD_O2) & (CTD_O2!="ND"))%>%
  filter((!is.na(SI)) & ("NAN" != SI) & (SI!="ND"))%>%
  filter((!is.na(Mean_NH4)) & ("NAN" != Mean_NH4) & (Mean_NH4!="ND"))%>%
  filter((!is.na(Mean_N2O)) & ("NAN" != Mean_N2O) & (Mean_N2O!="ND"))%>%
  filter((!is.na(Mean_NO2)) & ("NAN" != Mean_NO2) & (Mean_NO2!="ND"))%>%
  filter((!is.na(Mean_CH4)) & ("NAN" != Mean_CH4) & (Mean_CH4!="ND"))%>%
  dplyr::rename(O2=CTD_O2)
#N2, O2, H2S, CO2 does not contain enough data to be plotted

#Store all chemical data as numeric data
ChemDat <-
  ChemDat %>%
  dplyr::mutate(O2=as.numeric(O2)) %>%
  dplyr::mutate(PO4=as.numeric(PO4)) %>%
  dplyr::mutate(NO3=as.numeric(NO3)) %>%
  dplyr::mutate(SI=as.numeric(SI)) %>%
  dplyr::mutate(Mean_NH4=as.numeric(Mean_NH4)) %>%
  dplyr::mutate(Mean_N2O=as.numeric(Mean_N2O)) %>%
  dplyr::mutate(Mean_NO2=as.numeric(Mean_NO2)) %>%
  dplyr::mutate(Mean_CH4=as.numeric(Mean_CH4))
#dplyr::select(Cruise, Date, Depth,O2, WS_NO3, WS_H2S) 
#O2,PO4,SI,NO3....Longitude, Latitude,Date, Cruise, Depth
#Depth = 10
#Longitude = -123.505
#Latitude = 48.59167

#Filter the Cruises for the Cruises that have transcriptome data
#Clean up the Data and rename the columns, also gather the data
ChemDat_cruise <- ChemDat %>%
  filter(Cruise %in% Cruises) %>%
  dplyr::select(Date, Cruise, O2, PO4, NO3, SI, Mean_NH4,Mean_N2O,Mean_NO2,Mean_CH4) %>% 
  dplyr::rename(  "NO3 (μM)"=NO3,
                  "PO4 (μM)"=PO4,
                  "O2 (μM)"=O2,
                  "SiO2 (μM)"=SI,
                  "NH4 (μM)"=Mean_NH4,
                  "N2O (nM)"=Mean_N2O,
                  "NO2 (μM)"=Mean_NO2,
                  "CH4 (nM)"=Mean_CH4)   %>%
  gather(key="Chemical", value="Concentration", -Date, -Cruise)

######################
#Chemicals vs Time/Date graph
######################
#Plot chemicals over time graph
TimeGRAPH<-ChemDat %>%
  dplyr::select(Date, O2, PO4, NO3, SI, Mean_NH4,Mean_N2O,Mean_NO2, Mean_CH4) %>% 
  dplyr::rename(  "NO3 (μM)"=NO3,
                  "PO4 (μM)"=PO4,
                  "O2 (μM)"=O2,
                  "SiO2 (μM)"=SI,
                  "NH4 (μM)"=Mean_NH4,
                  "N2O (nM)"=Mean_N2O,
                  "NO2 (μM)"=Mean_NO2,
                  "CH4 (nM)"=Mean_CH4)  %>%
  gather(key="Chemical", value="Concentration", -Date) %>% 
  
  ggplot(aes(x=Date, y=Concentration, color=Chemical)) +
  geom_point() +
  facet_wrap(~Chemical, scales="free_y", nrow = 4) +
  theme_bw() +
  theme(legend.position="none")+
  geom_text(data = ChemDat_cruise, aes(x = Date, 
                                       y = Concentration, 
                                       label = paste("SI0", as.character(Cruise),
                                                     sep="")), 
            colour="black",
            vjust=1.25,
            size=3) +
  geom_point(data = ChemDat_cruise, colour="black") 

TimeGRAPH

#Select Chemical Data filtered and not filtered for Cruise
ChemDat_cruise2 <- ChemDat %>%
  filter(Cruise %in% Cruises) %>%
  dplyr::select(Cruise, O2, PO4, NO3, SI, Mean_NH4,Mean_N2O,Mean_NO2, Mean_CH4) 
ChemDat2 <- ChemDat %>%
  dplyr::select(Cruise,O2, PO4, NO3, SI,Mean_NH4,Mean_N2O,Mean_NO2, Mean_CH4) 

#Function for plotting chemical X against chemical Y, and generate a single graph
#Not used as a final figure
plotting<-function(X_AXIS,Y_AXIS,COLOR_ID){
  ChemDat2%>% 
    ggplot(aes_string(x=X_AXIS, y=Y_AXIS)) +
    geom_point(color=pal[COLOR_ID]) +
    theme_bw() +
    theme(legend.position="none")+
    geom_text(data = ChemDat_cruise2, aes_string(x = X_AXIS, y = Y_AXIS, label="Cruise"),
              colour="black",
              vjust=-0.25,
              size=3) +
    #  aes(label = paste("SI0", as.character(ChemDat_cruise2$Cruise),sep="")) +
    geom_point(data = ChemDat_cruise2, colour="black") 
}

p1<-plotting("O2","PO4",1)
p2<-plotting("O2","NO3",2)
p3<-plotting("O2","SI",3)
p4<-plotting("O2","Mean_N2O",4)
p5<-plotting("O2","Mean_CH4",5)
p6<-plotting("O2","Mean_NH4",6)
p7<-plotting("O2","Mean_NO2",7)

p<-plot_grid(p1, p2, p3,p4,p5,p6,p7, labels=c("A", "B", "C","D","E","F","G"), align="h", axis="tb", ncol = 4)

#Function for plotting chemical X against all chemicals, and generate a single graph to save as .pdf
#Not used in final figures
plotting2<-function(X_AXIS){
  ChemDat_cruise2 <- ChemDat2 %>%
    filter(Cruise %in% Cruises) %>%
    dplyr::select(Cruise,X_AXIS,O2, PO4, NO3, SI, Mean_NH4,Mean_N2O,Mean_NO2, Mean_CH4) %>% 
    gather_(key="Chemical", value="Concentration", gather_cols=colnames(.[,3:8])) 
  
  ChemDat2%>% 
    dplyr::select(X_AXIS,O2, PO4, NO3, SI, Mean_NH4,Mean_N2O,Mean_NO2, Mean_CH4) %>% 
    gather_(key="Chemical", value="Concentration", gather_cols= colnames(.[,2:7])) %>% 
    
    ggplot(aes_string(x=X_AXIS, y="Concentration", color="Chemical")) +
    geom_point() +
    facet_wrap(~Chemical, scales="free_y", nrow = 4) +
    theme_bw() +
    theme(legend.position="none")+
    geom_smooth(method="lm") +
    geom_text(data = ChemDat_cruise2, 
              aes_string(x = X_AXIS, y = "Concentration", label ="Cruise"), 
              colour="black",
              vjust=-0.25,
              size=3) +
    geom_point(data = ChemDat_cruise2, colour="black")
  
}

Chemicals_ID <- c("O2","PO4","Mean_N2O","Mean_NO2","SI","NO3","Mean_NH4","Mean_CH4")
for(x in 1:length(Chemicals_ID)){
  X=Chemicals_ID[x]
  pdf(paste(X,"pdf",sep=".")) 
  print(plotting2(as.character(X)))
  dev.off()
}

############################################### 
#Heatmap of R2 correlations between filtered chemicals
###############################################
#Function for r square, m, b of given data set
lm_eqn = function(df){
  m = lm(y ~ x, df);
  # eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
  #                  list(a = format(coef(m)[1], digits = 2), 
  #                       b = format(coef(m)[2], digits = 2), 
  #                       r2 = format(summary(m)$r.squared, digits = 3)))
  # as.character(as.expression(eq));    
  
  eq <-  list(a = format(coef(m)[1], digits = 2),
              b = format(coef(m)[2], digits = 2),
              r2 = format(summary(m)$r.squared, digits = 3))
  as.character(eq);   
} 

#Chemicals Used/Have enough non-NA values to be used 
Chemicals_ID <- c("O2", "NO3", "PO4","Mean_N2O","Mean_NO2","SI","Mean_NH4","Mean_CH4")

#Generate linear equation, R square and store in a matrix, cleaned up the matrix
a <-c()
for(i in 1:length(Chemicals_ID)){ 
  for(j in 1:length(Chemicals_ID)){
    X = Chemicals_ID[i]
    Y = Chemicals_ID[j]
    if(X!=Y){
      a<-ChemDat2 %>%
        rename_(x=X,y=Y) %>%
        lm_eqn() %>%
        c(X, Y, .) %>%
        c(a,.)
    }
    else
      a<-c(a,X,Y,1,1,1)
  }
}
#Re-arrange Rsquare, y = mx + b values to data frame with corresponding column names
lineregression<-matrix(a, ncol=5, byrow=TRUE)
colnames(lineregression) = c("X", "Y","m","b","r^2")
lineregression <- as.data.frame(lineregression)
lineregression$`r^2`<- as.double(matrix(a, ncol=5, byrow=TRUE)[,5])

#Plot RSquare between chemicals as heatmap
Chemical_RSquare<-ggplot(lineregression, aes(X, Y )) +
  geom_tile(aes(fill = `r^2`), color = "white") +
  scale_fill_gradient(low = "white", high = pal[7]) +
  ylab("Chemicals") +
  xlab("Chemicals") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "r^2")

#Function to generate equation on graph based on previously calculated correlations and linear model
lm_eqn_graph = function(M,B,R){
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = M, 
                        b = B, 
                        r2 = R))
  as.character(as.expression(eq));    
} 

#Plot PO4 vs NO3 linear graph with linear equation model
ChemDat_cruise2 <- ChemDat2 %>%
  filter(Cruise %in% Cruises) %>%
  dplyr::select(Cruise, PO4, NO3) 
NO3vsPO4<-ChemDat2%>% 
  dplyr::select(PO4, NO3) %>% 
  ggplot(aes_string(x="PO4", y="NO3", color="NO3")) +
  geom_point() +
  theme(legend.position="none")+
  geom_smooth(method="lm") +
  ylab("NO3 (μM)") +
  xlab("PO4 (μM)") +
  geom_text(data = ChemDat_cruise2, 
            aes(x = PO4, y = NO3, label = paste("SI0",Cruise, sep="")), 
            colour="black",
            vjust=-0.25,
            size=3) +
  geom_point(data = ChemDat_cruise2, colour="black")
#Add linear equation model to the graph
NO3vsPO4 <- NO3vsPO4 + geom_text(x = 1, 
                              y = 30,
                              label=lm_eqn_graph(as.character(lineregression$m[18]),
                                                 as.character(lineregression$b[18]),
                                                 lineregression$`r^2`[18]), 
                              parse = TRUE)

#Combine chemical correlations and NO3 vs PO4 data graphs
Chemical_Figure <-plot_grid(Chemical_RSquare,NO3vsPO4, 
                            labels=c("A", "B"), 
                            align="h", 
                            axis="tb", 
                            ncol = 2)
Chemical_Figure

##############################################
#CheckM, Plot Class, RPKM, Completeness, Contamination
##############################################
checkm <- read_tsv("~/Desktop/MICB_405/Project2/MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv",
                   col_names = TRUE)
#Load backteria and archea classification by gtdb package
arc_class <- read.table("~/Desktop/MICB_405/Project2/gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("~/Desktop/MICB_405/Project2/gtdbtk.bac120.classification_pplacer.tsv", sep="\t")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
#Load total RPKM from MAG reads
binned_rpkm <- read_csv("~/Desktop/MICB_405/Project2/SaanichInlet_10m_binned.rpkm.csv",
                        col_names = TRUE)
#Clean-up RPKM, separate Bin_ID from Sequence, and summarize total RPKM for each Bin
binned_rpkm <- binned_rpkm %>%
  mutate(Bin_Id=sapply(strsplit(Sequence, split='_', fixed=TRUE), `[`, 3)) %>%
  mutate(Bin_Id=paste("SaanichInlet_10m.",Bin_Id,sep="")) %>%
  group_by(Bin_Id) %>%
  summarise(Total_RPKM=sum(RPKM, na.rm=TRUE))

#Combine gtdb data + checkM output + rpkm by Bin Id
mergedat <- merge(checkm, gtdb_dat, by.x="Bin Id", by.y="mag") %>%
  merge(binned_rpkm, by.x="Bin Id", by.y="Bin_Id")

#Plot Abundance(RPKM), Completeness, Contamination, Class 
circle_plot<- mergedat %>%
  mutate(Class=substring(Class,4)) %>%
  dplyr::rename(Abundance=Total_RPKM)%>%
  ggplot (aes(x=Completeness, y=Contamination, color = Class, size = Abundance)) +
  geom_point () +
  scale_x_continuous(name ="Completeness (%)")+
  scale_y_continuous(name ="Contamination (%)")+
  theme(legend.text = element_text(size=8)) +
  geom_text(data=mergedat_apb,aes(label=BinID),
            colour="black",
            size=3) 

circle_plot

#Plot Abundance(RPKM), Completeness, Contamination, Class for Alphaproteobacteria
circle_plot_apb<- mergedat %>%
  mutate(Class=substring(Class,4)) %>%
  dplyr::rename(Abundance=Total_RPKM)%>%
  mutate(BinID=substring(`Bin Id`,18)) %>%
  filter(BinID%in%apb_id) %>%
  ggplot (aes(x=Completeness, y=Contamination, color = Class, size = Abundance)) +
  geom_point () +
  theme(legend.position = "bottom")+
  geom_text(data=mergedat_apb,aes(label=BinID),
            colour="black",
            size=3) 

circle_plot_apb

#Visualize quality cutoff, Quality=Completeness-(5*Contamination) > 50
mergedat %>%
  dplyr::select(Phylum, Class, Order, Family,Genus, `Bin Id`, Total_RPKM,Completeness,Contamination)%>%
  filter(Class =="c__Alphaproteobacteria") %>%
  dplyr::mutate(Quality=Completeness-(5*Contamination)) %>%
  dplyr::arrange(desc(Quality)) %>%
  View()

##############################################
#Prokka Output
##############################################
#Clean up Prokka output
prokka_output <- read_tsv("Desktop/MICB_405/Project2/combined_prokka_num.tsv",
                          col_names=FALSE) 
for(x in 1:length(prokka_output$X1)){
  prokka_output$Value[x] <- unlist(strsplit(prokka_output$X1[x], ": "))[2]
  prokka_output$X1[x] <- unlist(strsplit(prokka_output$X1[x], ": "))[1]
}
write_csv(prokka_output, path="~/Desktop/MICB_405/Project2/prokka_output.csv" )
#Use excel to select and transpose prokka_output files to the selected MAGs 
prokka_output <- read_csv("Desktop/MICB_405/Project2/prokka_output.csv",
                          col_names=TRUE) 

#Fix the order of Prokka MAGs/Bins when plotting
prokka_output$X1 <- 
  factor(prokka_output$X1, levels = prokka_output$X1) 

#Plot a prokka output graph for each of: CDS,tRNA, tmRNA, rRNA,Contigs, Bases
prokka_CDS<-prokka_output %>%
  ggplot (aes(x=X1, y=CDS, fill = Class)) +
  geom_col () +
  scale_x_discrete(name ="MAGs")+
  scale_y_continuous(name ="CDS")+
  theme(legend.position = "none",
        legend.text = element_text(size=8))
prokka_tRNA<-prokka_output %>%
  ggplot (aes(x=X1, y=tRNA, fill = Class)) +
  geom_col () +
  scale_x_discrete(name ="MAGs")+
  scale_y_continuous(name ="tRNA")+
  theme(legend.position = "bottom",
        legend.text = element_text(size=8))
prokka_rRNA<-prokka_output %>%
  ggplot (aes(x=X1, y=rRNA, fill = Class)) +
  geom_col () +
  scale_x_discrete(name ="MAGs")+
  scale_y_continuous(name ="rRNA")+
  theme(legend.position = "none",
        legend.text = element_text(size=8))
prokka_tmRNA<-prokka_output %>%
  ggplot (aes(x=X1, y=tmRNA, fill = Class)) +
  geom_col () +
  scale_x_discrete(name ="MAGs")+
  scale_y_continuous(name ="tmRNA")+
  theme(legend.position = "none",
        legend.text = element_text(size=8))
prokka_contigs<-prokka_output %>%
  ggplot (aes(x=X1, y=contigs, fill = Class)) +
  geom_col () +
  scale_x_discrete(name ="MAGs")+
  scale_y_continuous(name ="Contigs")+
  theme(legend.position = "none",
        legend.text = element_text(size=8))
prokka_base<-prokka_output %>%
  ggplot (aes(x=X1, y=bases, fill = Class)) +
  geom_col () +
  scale_x_discrete(name ="MAGs")+
  scale_y_continuous(name ="Bases")+
  theme(legend.position = "none",
        legend.text = element_text(size=8))
#Generate a single plot for Prokka Output
prokka_plot<-plot_grid(prokka_base,
                       prokka_contigs,
                       prokka_CDS,
                       prokka_tRNA + theme(legend.position="none"),
                       prokka_tmRNA,
                       prokka_rRNA, 
                       labels=c("A", "B","C","D","E","F"), 
                       axis="tb",nrow=2) %>%
  plot_grid(., get_legend(prokka_tRNA),nrow=2, 
            rel_heights = c(9/10,1/10), 
            axis="lr")
prokka_plot