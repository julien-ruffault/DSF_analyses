RHfromSHandT<-function(Tc, HS, allowSaturated = FALSE){
  #-------------------------------------------------------------
  #D?claration des constantes pour le calcul de l'HR et de la T en ?C
  Mas=28.966 # masse molaire air sec (g/mol)
  Mh2o=18 # masse molaire H2O(g/mol)
  Rgz=8.314472 # %J/mol/K cste gaz parfait
  p_air=101325 #  %en Pa
  #-------------------------------------------------------------
  

  #Calcul de l'Humidit? relative bas? sur Hs et Tk
  Tk=Tc+273.15
  Dair=((p_air)/(Rgz*(Tk)))*Mas # Calcul de la masse volumique de l'air sec en g.m-3
  masshum=HS*Dair # Masse d'eau dans l'air g.m-3
  nhum=masshum/Mh2o # en mol.m-3
  ea=nhum*(Rgz*(Tk)) # Pression de vapeur r?elle (en Pa)
  es=6.108*exp(17.27*Tc/(237.2+Tc))*100 # Pression de vapeur saturante ? la Temp?rature T (en Pa)
  HR=ea/es*100 #Calcul de l'HR
  if(!allowSaturated) HR[HR>100]=100 # On ne passe pas audessus de 100
  return(HR)
}



compute.VPDfromRHandT <- function(relative_humidity, temperature, air_pressure = 101325) {
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # AUTHORS     :  Julien Ruffault (julien.ruff@gmail.com)
  #                Nicolas Martin-StPaul (nicolas.martin@inrae.fr)
  
  # Constants
  Mas <- 28.966    # molar weight dry air (g/mol)
  Mh2o <- 18       # molar weight H20 H2O(g/mol)
  Rgz <- 8.314472  # pPerfect gaz constant %J/mol/K
  
  Tk <- temperature + 273.15 # conversion of temperature in K
  Dair <- ((air_pressure) / (Rgz * (Tk))) * Mas
  es <- 6.108 * exp(17.27 * temperature / (237.2 + temperature)) * 100
  ea <- relative_humidity * es / 100
  vpd <- (es - ea) / 1000
  vpd[vpd < 0] <- 0
  return(vpd)
}


