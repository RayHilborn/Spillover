#This code runs the 9 cases for the base model

source("Spillover functions multistep v2.r")
source("BaseSetup.r")
Policy="effort"

DispList=c(3,10,20)   #dispersal rates
UMSYList=c(.5,1,2)

for (Dbase in DispList){
  for (UUMSY in UMSYList){

DispersalSD=Dbase
  M=BaseMPA()

   save(file=paste0("Model ",Dbase," ",UUMSY," ",Policy,".rdata"),M)
  
  }
}
  
 
  
  