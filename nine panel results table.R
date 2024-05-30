

DispList=c(3,10,20) #dispersal rates
UMSYList=c(.5,1,2)
Res=array(dim=c(9,7))
colnames(Res)=c("D","UUMSY","Cbefore","Cafter","Bbefore","Bafter","CPUEChange")
DispersalSD=3; UUMSY=1
row=0
Policy="effort"
for (DispersalSD in DispList){
  for (UUMSY in UMSYList){
row=row+1    
load(file=paste0("Model ",DispersalSD," ",UUMSY," ",Policy,".rdata"))
Res[row,"D"]=DispersalSD
Res[row,"UUMSY"]=UUMSY
Res[row,"Cbefore"]=M$Chist[M$MPAStartYear-1]
Res[row,"Cafter"]=M$Chist[M$Ntime]
Res[row,"Bbefore"]=sum(M$B[,M$MPAStartYear-1])
Res[row,"Bafter"]=sum(M$B[,M$Ntime])
Res[row,"CPUEChange"]=M$CPUE[100]/M$CPUE[49]#   M doesn't have year of MPA establishment

}} #end of loops over models
D=as.data.frame(Res)
D$CatchChange=D$Cafter/D$Cbefore
D$BChange=D$Bafter/D$Bbefore
D$CMSYbefore=D$Cbefore/(M$MSY*M$ncell)
D$CMSYafter=D$Cafter/(M$MSY*M$ncell)
write.csv(file=paste0("9panel ",Policy,".csv"),D)

