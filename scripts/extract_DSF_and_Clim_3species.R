# script cree le 16/12/2021 par julien ruffault (julien.ruff@gmail.com)
# selection de placettes DSF +  ajout de donnes climatiques annuelles et moyennes (SAFRAN)
# especes selectionnees : pinus halepensis, Quercus Pubescens, Quercus Ilex 
#  notes  1 : une placette est selectionne si au moins 10 individus ont été échantillonés au moins un annee 
#         2 : Les donnes DSF ont ete pretraitres cf (raw_tree_data/DSF_France/notes_extraction_DSF_02_2021.png)

# initialization ----------------------------------------------------------
rm(list=ls())
gc()
# libraries -------------------------------------------------------------
library(dplyr)
library(plyr)
library(rgdal)
library(sp)
library(raster)
library(dplyr)
library(lubridate)
library(SPEI)
# directories -------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change if necessary
getwd()



# load datasets -------------------------------------------------------------
# load spatial data for representation 
DEPS= readOGR(dsn=file.path(paste0('../aux_data/DEPS')))

# load spatial data / countries 
CTR = readOGR(dsn = file.path('../aux_data/ne_50m_admin_0_countries') ,layer= 'ne_50m_admin_0_countries')
OCEAN = readOGR(dsn=file.path('../aux_data/ne_50m_ocean'))

# SAFRAN directory (works only in INRAE Avignon) 
SAFRAN_DIR = file.path("/Volumes/SHARED_DATA/CLIMAT/SAFRAN") # directory to NAS PEF  /change if necessary
COORD_SAF <- read.table(file.path(SAFRAN_DIR, "COORD_SAFRAN.csv"), sep = ";", header = T)

# load raw data DSF 
DSF_L1 = read.csv("../raw_tree_data/DSF_France/ALL_Reseau_systematique_niveau1_02_2021.csv",dec=',',sep=';')
dim(DSF_L1)
DSF_L2  = read.csv("../raw_tree_data/DSF_France/ALL_Reseau_systematique_niveau2_02_2021.csv",dec=',',sep=';')
dim(DSF_L2)
DSF_ALL= rbind(DSF_L1,DSF_L2)
DSF_ALL <- DSF_ALL[-which(is.na(DSF_ALL$Deficit_foliaire)==T),]
#plot(DSF_ALL$X_Lambert,DSF_ALL$Y_Lambert,col='red') # 


# load climate utils 
source("utils/climate_utils.R")

# SELEC SITES  -----------------------------------------------------------
#placettes PH
io = aggregate(DSF_ALL$Essence,by =list(Annee=DSF_ALL$Annee,Code_placette=DSF_ALL$Code_placette),function(x) sum(x==125))
selec_placette_PH = unique(io[io$x>=10,"Code_placette"])
print( paste0('numbber of sites for pinus hal. : ',length(selec_placette_PH)))

# placettes QI
io = aggregate(DSF_ALL$Essence,by =list(Annee=DSF_ALL$Annee,Code_placette=DSF_ALL$Code_placette),function(x) sum(x==46))
selec_placette_QI = unique(io[io$x>=10,"Code_placette"])
print( paste0('numbber of sites for Quercus ilex : ',length(selec_placette_QI)))

#placettes QP
io = aggregate(DSF_ALL$Essence,by =list(Annee=DSF_ALL$Annee,Code_placette=DSF_ALL$Code_placette),function(x) sum(x==49))
selec_placette_QP = unique(io[io$x>=10,"Code_placette"])
print( paste0('numbber of sites for Quercus pubescens : ',length(selec_placette_QP)))

# DATASET with only QI, QP, and QI 
DSF = DSF_ALL[DSF_ALL$Code_placette %in% unique(c(selec_placette_PH,selec_placette_QI,selec_placette_QP)),]



# keeping only useful column
DSF = DSF[,c("Code_placette","Annee","X_Lambert","Y_Lambert",
             "Essence","Deficit_foliaire","Numero_darbre","Mortalite_branches",
             "Circonference","Classe_dage_mesuree","Statut_annee_suivante" ,
             "RNO_annee_suivante","Agent_dommages")]
head(DSF)









# On associe un numéro de maille Safran à placette DSF
io = DSF %>% group_by(Code_placette) %>% filter(row_number()==1)
for(i in 1:nrow(io)){
# print(io[i,])
print(i)
  Tmp <- pointDistance(COORD_SAF [,c('XL2e','YL2e')],io[i,c('X_Lambert','Y_Lambert')],lonlat=F)
  # On règle les conflits éventuels liés à deux mailles équidistantes en prenant le max des valeurs
  io$Code_safran[i] <- as.numeric(max(COORD_SAF[Tmp==min(Tmp),][1]))
}



# Pour chaque placette on recupere les donnes meteo : (PET-PPT) + SPEI3 et SPEI6 de chaque annee
DSF_CLIM = DSF
DSF_CLIM$WBmoy =NA
DSF_CLIM$SPEI3_7 = NA
DSF_CLIM$SPEI6_7 = NA

DSF_CLIM$SPEI3_8 = NA
DSF_CLIM$SPEI6_8 = NA

DSF_CLIM$SPEI3_9 = NA
DSF_CLIM$SPEI6_9 = NA

DSF_CLIM$Tmax3_9 = NA 
DSF_CLIM$VPDmax3_9 = NA

DSF_CLIM$Tmean3_9 = NA 
DSF_CLIM$VPDmean3_9 = NA 

for (i in 1:nrow(io)){
  
  tmp_SAF = read.table(file= paste0(SAFRAN_DIR,'/1960_2020/quotidiennes_1959_2020_maille_',io$Code_safran[i],'.csv.gz'),sep=';',h=T)
  tmp_SAF$date=as.Date(tmp_SAF$date,format='%d/%m/%Y')
  tmp_SAF$year <- year(tmp_SAF$date)
  tmp_SAF$month <-  month(tmp_SAF$date) 
  tmp_SAF=tmp_SAF[tmp_SAF$year>=1989,]
  
  tmp_SAF$PPT_sum     = tmp_SAF[,"preliq_q"]+tmp_SAF[,"prenei_q"]## precipitation (mm)
  tmp_SAF$RH          = RHfromSHandT(Tc=tmp_SAF[, "t_q" ],HS = tmp_SAF[, "q_q"] /1000)
  tmp_SAF$VPD         = compute.VPDfromRHandT(relative_humidity= tmp_SAF$RH,temperature=tmp_SAF$t_q) 
  
  tmp_SAF_Month <-  aggregate(tmp_SAF[,c('t_q')],by=list(month = tmp_SAF$month,year=tmp_SAF$year),mean)
  colnames(tmp_SAF_Month) <- c('month','year','tmean')
  tmp_SAF_Month$ppt = aggregate(tmp_SAF[,c('PPT_sum')],by=list(month = tmp_SAF$month,year=tmp_SAF$year),sum)[,3]
  tmp_SAF_Month$date= as.Date(paste('1',tmp_SAF_Month$month,tmp_SAF_Month$year,sep='/'),format='%d/%m/%Y')

  
  tmp_SAF_Month$pet= as.vector(thornthwaite(Tave=tmp_SAF_Month$tmean,lat=43.24))
  tmp_SAF_Month$wb = tmp_SAF_Month$ppt- tmp_SAF_Month$pet
  
  tmp_SAF_Month$spei3 = as.vector(spei(data =   tmp_SAF_Month$wb,scale=3)$fitted)
  tmp_SAF_Month$spei6 = as.vector(spei(data =   tmp_SAF_Month$wb,scale=6)$fitted)

  wb_SAF_year = aggregate(tmp_SAF_Month$wb,by=list(year=tmp_SAF_Month$year),sum)
  pet_SAF_year = aggregate(tmp_SAF_Month$pet,by=list(year=tmp_SAF_Month$year),sum)
  ppt_SAF_year = aggregate(tmp_SAF_Month$ppt,by=list(year=tmp_SAF_Month$year),sum)
  wb_ete_SAF_year = aggregate(tmp_SAF_Month$wb[tmp_SAF_Month$month %in% c(6,7,8)],by=list(year=tmp_SAF_Month$year[tmp_SAF_Month$month %in% c(6,7,8)]),sum)
 
  
  io_789 = tmp_SAF[tmp_SAF$month %in% c(7,8,9),]
  
  Tmean3_SAF_YEAR = aggregate(io_789$t_q,by=list(year=io_789$year),mean)
  Tmax3_SAF_YEAR= aggregate(io_789$t_q,by=list(year=io_789$year),max)
  
  VPDmean3_SAF_YEAR = aggregate(io_789$VPD,by=list(year=io_789$year),mean)
  VPDmax3_SAF_YEAR  = aggregate(io_789$VPD,by=list(year=io_789$year),max)
  
  
 # atrribut mean values 
 DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i],'WBmoy'] = mean(wb_SAF_year$x)
 DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i],'PPTmoy'] = mean(ppt_SAF_year$x)
 DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i],'PETmoy'] = mean(pet_SAF_year$x)
 DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i],'WBmoy_ete'] = mean(wb_ete_SAF_year$x)
 
 for (Y in 1989:2020)
 {
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'SPEI3_8'] = tmp_SAF_Month[tmp_SAF_Month$month==8 & tmp_SAF_Month$year==Y,'spei3']
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'SPEI6_8'] = tmp_SAF_Month[tmp_SAF_Month$month==8 & tmp_SAF_Month$year==Y,'spei6']
 
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'SPEI3_7'] = tmp_SAF_Month[tmp_SAF_Month$month==7 & tmp_SAF_Month$year==Y,'spei3']
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'SPEI6_7'] = tmp_SAF_Month[tmp_SAF_Month$month==7 & tmp_SAF_Month$year==Y,'spei6']
   
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'SPEI3_9'] = tmp_SAF_Month[tmp_SAF_Month$month==9 & tmp_SAF_Month$year==Y,'spei3']
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'SPEI6_9'] = tmp_SAF_Month[tmp_SAF_Month$month==9 & tmp_SAF_Month$year==Y,'spei6']
   
   
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'VPDmean3_9'] =   VPDmean3_SAF_YEAR[VPDmean3_SAF_YEAR$year==Y,'x']
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'VPDmax3_9'] = VPDmax3_SAF_YEAR[VPDmax3_SAF_YEAR$year==Y,'x']
   
   
   
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'Tmean3_9'] =  Tmean3_SAF_YEAR[Tmean3_SAF_YEAR$year==Y,'x']
   DSF_CLIM[DSF_CLIM$Code_placette==io$Code_placette[i] & DSF_CLIM$Annee==Y,'Tmax3_9']  =  Tmax3_SAF_YEAR[Tmax3_SAF_YEAR$year==Y,'x']
   
   
   
   
}
 }
 

# write_final_data_set for future use
write.csv(file= '../compiled_datasets/DSF_3species_CLIM.csv', DSF_CLIM)



# A=DSF%>% 
#   group_by(Code_placette,Annee) %>% 
#   dplyr::summarise(defol=sum(Deficit_foliaire==100)/n()) %>%
#   group_by(Annee) %>% dplyr::summarise(defol = mean(defol))
# plot(A)

#get SAFRAN neareast point for each pixel
#plot(COORD_SAF$XL2e,COORD_SAF$YL2e)
#points(DSF$X_Lambert,DSF$Y_Lambert,col="red")





A=DSF_CLIM %>% 
  group_by(Code_placette,Annee) %>% 
  dplyr::summarise(defol=mean(Deficit_foliaire)) %>%
  group_by(Annee) %>% dplyr::summarise(defol = mean(defol))
plot(A)


sYear=1990
eYear=2020
sp = 46
x=DSF_CLIM
pThres=10
# x    :  data.frame DSF_CLIM 
# sp   : species number  (49 : QP, 46: QI, 125 : PH) 
# sYear : starting year  ( take second year by default) 
# eYear  : ending Year  (lastYEAR)
# pThres : number of minimum individuals of the selected species to keep the site (default = T)
calculate.deficit.fol <- function(x,sp,sYear,eYear,pThres)
{

if (missing(sYear)==T){sYear = sort(unique(x$Annee))[2]}
if (missing(eYear)==T){eYear = sort(unique(x$Annee))[length(unique(x$Annee))]}
if (missing(pThres)==T){pThres=10} 
  
io = aggregate(DSF_CLIM$Essence,by =list(Annee=DSF_CLIM$Annee,Code_placette=DSF_CLIM$Code_placette),function(x) sum(x==sp))
selec_placette = unique(io[io$x>=pThres,"Code_placette"])
  


FFF = data.frame(matrix(nrow=length(selec_placette)*length(sYear:eYear),ncol=ncol(x)+3))
colnames(FFF) = c(colnames(x),'diff_def','mean_def','n_common_trees')
count=0
for (PP in selec_placette){
# add a loop for placettes 
tmp = DSF_CLIM[DSF_CLIM$Code_placette==PP,]

for  (YY in sYear:eYear)
{
  count=count+1
  print(YY)
  io = tmp[tmp$Annee==YY & tmp$Essence==sp,]
  #print(unique(io$Numero_darbre))
  
  ioM1 = tmp[tmp$Annee==(YY-1) & tmp$Essence==sp,]
  #print(unique(ioM1$Numero_darbre))
  
  common_trees <- Reduce(intersect, list(io$Numero_darbre, ioM1$Numero_darbre))
  

  aaa = mean(io[io$Numero_darbre %in% common_trees,'Deficit_foliaire']) - mean(ioM1[ioM1$Numero_darbre %in% common_trees,'Deficit_foliaire'])
  
  FFF[count,'Code_placette'] = PP
  FFF[count,'Annee'] = YY
  FFF[count,'X_Lambert'] = io[1,'X_Lambert']
  FFF[count,'Y_Lambert'] = io[1,'Y_Lambert']
  FFF[count,'WBmoy'] = io[1,'WBmoy'] 
  FFF[count,'SPEI3_7'] = io[1,'SPEI3_7'] 
  FFF[count,'SPEI6_7'] = io[1,'SPEI6_7'] 
  FFF[count,'SPEI3_8'] = io[1,'SPEI3_8'] 
  FFF[count,'SPEI6_8'] = io[1,'SPEI6_8'] 
  FFF[count,'SPEI3_9'] = io[1,'SPEI3_9'] 
  FFF[count,'SPEI6_9'] = io[1,'SPEI6_9'] 
  FFF[count,'PPTmoy'] = io[1,'PPTmoy'] 
  FFF[count,'PETmoy'] = io[1,'PETmoy']    
  FFF[count,'WBmoy_ete'] = io[1,'WBmoy_ete']   
  FFF[count,'Essence'] = sp
  
  FFF[count,'n_common_trees'] = length(common_trees)
  FFF[count,'diff_def'] = aaa
  FFF[count,'mean_def']= mean(io[io$Numero_darbre %in% common_trees,'Deficit_foliaire'])

}
    
  
}
  return(FFF)
}




QP = calculate.deficit.fol(DSF_CLIM,sp = 49,sYear=2000)
QI = calculate.deficit.fol(DSF_CLIM,sp = 46,sYear=2000)
PH = calculate.deficit.fol(DSF_CLIM,sp = 125,sYear=2000)
ALL = rbind(QP,QI,PH)
ALL$Essence=as.factor(ALL$Essence)


library(mgcv)

plot(QP$WBmoy,QP$diff_def)


A = gam(mean_def ~ s(SPEI6_9) + s(WBmoy) +  Essence ,data=QP)
summary(A)
gam.check(A)
plot(A)




AAA= aggregate(QP,by=list(QP$Annee),function(x) mean(x,na.rm=T))

plot(AAA$mean_def,type='l')
plot(AAA$diff_def,type='l')
lines(AAA$SPEI6_8,col='blue')

plot(AAA$SPEI6_8,AAA$diff_def)
cor(AAA$SPEI3_9,AAA$diff_def)



plot(aggregate(TEST$diff_def,by=list(TEST$Annee),function(x) mean(x,na.rm=T)),type='l')


abline(v=2003)



# Selec sites with a least 10 individual of the selected species were sampled during the period (note that this is not a strict selection but that was made for keeping the maximum sites as lower trees might be have sample during other years; We might need to change this criteria afterwards)
io = aggregate(DSF_CLIM$Essence,by =list(Annee=DSF_CLIM$Annee,Code_placette=DSF_CLIM$Code_placette),function(x) sum(x==49))
selec_placette_QP = unique(io[io$x>=10,"Code_placette"])

io = aggregate(DSF_CLIM$Essence,by =list(Annee=DSF_CLIM$Annee,Code_placette=DSF_CLIM$Code_placette),function(x) sum(x==46))
selec_placette_QI = unique(io[io$x>=10,"Code_placette"])

io = aggregate(DSF_CLIM$Essence,by =list(Annee=DSF_CLIM$Annee,Code_placette=DSF_CLIM$Code_placette),function(x) sum(x==125))
selec_placette_PH = unique(io[io$x>=10,"Code_placette"])

selec_placette_ALL = c(selec_placette_PH,selec_placette_QI,selec_placette_QP)




# calcul du deficit foliaire moyne par site annee pour l'epsece seletionne sur les sites avec au moins 10 indiviuds pour une annee
azer = DSF_CLIM[DSF_CLIM$Code_placette %in% selec_placette_QP,]
azer = azer[azer$Essence==49,] # keep only selected species 
azer2 = aggregate(azer$Deficit_foliaire,by=list(Code_placette=azer$Code_placette, Annee = azer$Annee),mean)
azer3 <-azer2[order(azer2$Code_placette,azer2$Annee),]
azer3$Species='QP'


io  = DSF_CLIM %>% group_by(Code_placette,Annee) %>% filter(row_number()==1)


azer4=  join(azer3,io,by=c("Code_placette","Annee"))

library(mgcv)

A = gam(x ~ s(SPEI3)+s(WBmoy) + s(SPEI3,WBmoy), data=azer4,method = "REML")
summary(A)

gam.check(A)
plot(A)
#A =DSF_CLIM %>% group_by(Code_placette) %>% filt






# for plot mean conditons for DSF points (mean water balance to check restuls)
A = DSF_CLIM %>% group_by(Code_placette) %>% filter(row_number()==1)
ggplot(data=A, aes(x=X_Lambert,y=Y_Lambert,colour=WBmoy_ete))+ 
  geom_point(alpha=.5, size = 6) + scale_color_gradientn(colours=rainbow(5)) 
graphics.off()
p1 = ggplot(data=A, xName='X_Lambert',yName='Y_Lambert')
plot(A$X_Lambert,A$Y_Lambert,colors=A$WBmoy)


 
 
 
 
 # a finir pour affecter les info climat  à chaque site*espece*anne  A note que le climat moyen lui, sera le meme pour chaque annee

 

 
 
# general PLOTS  / attention les plots si dessous sont tre gneerique . elle mélagnent l'ensmble ede site et ne differentcient pas les especes au sein des sites ! 
  
plot(aggregate(DSF$Deficit_foliaire,by=list(DSF$Annee),mean))
abline(v=2003) 
plot(aggregate(DSF$Deficit_foliaire,by=list(DSF$Annee),mean)) 

plot(aggregate(DSF$Deficit_foliaire,by=list(DSF$Annee),function(x) 100*sum(x>=90)/length(x))) 
abline(v=2003) 


 
 

#### PLOTS ####### (for furtheer use / to be moved in another scritps afterwards)
# make big maps showing the position of each station for the diufferet species 

library(PNWColors)
pal <- pnw_palette("Sailboat",3)


QI = DSF %>% filter(Code_placette %in% selec_placette_QI)  %>% group_by(Code_placette) %>% filter(row_number()==1)
coordinates(QI) <- c('X_Lambert','Y_Lambert')
crs(QI) =CRS("+init=EPSG:27572")
QI.WGS84= spTransform(QI,crs(OCEAN))

PH = DSF %>% filter(Code_placette %in% selec_placette_PH)  %>% group_by(Code_placette) %>% filter(row_number()==1)
coordinates(PH) <- c('X_Lambert','Y_Lambert')
crs(PH) =CRS("+init=EPSG:27572")
PH.WGS84= spTransform(PH,crs(OCEAN))

QP= DSF %>% filter(Code_placette %in% selec_placette_QP)  %>% group_by(Code_placette) %>% filter(row_number()==1)
coordinates(QP) <- c('X_Lambert','Y_Lambert')
crs(QP) =CRS("+init=EPSG:27572")
QP.WGS84= spTransform(QP,crs(OCEAN))


plot(OCEAN,xlim=c(-5,12),ylim=c(40,51),col=rgb(240,248,252,maxColorValue=256),lwd=0.5)
plot(CTR,xlim=c(-5,12),ylim=c(40,51),add=T,col="grey60",lwd=0.5)
plot(QI.WGS84,add=T,pch=21,bg=pal[1],cex=0.7)
plot(PH.WGS84,add=T,pch=21,bg=pal[2],cex=0.7)
plot(QP.WGS84,add=T,pch=21,bg=pal[3],cex=0.7)
legend('topright',pch=21,legend=c('Quercus pubescens','Quercus ilex', 'Pinus halepensis'),pt.bg=c(pal[3],pal[1],pal[2]))









