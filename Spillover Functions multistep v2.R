


MPAmodel=function(Data,mpaList,D,LD){
  if (Debug==TRUE) print("Start MPAModel")
  B=array(dim=c(Data$ncell,(Data$Ntime+1)));  HR=B ; Recruitment=B;BoatStore=B
  N=B #this will store the numbers
  Chist=array(dim=Data$Ntime)
  B[,1:10]=Data$Bzero ; HR[,1:10]=0 #initialize arrays
  N[,1:10]=Data$Rzero*(Data$Surv/(1-Data$Surv)) #initialize numbers
  Recruitment[,1:10]=Data$Rzero
  CPUE=array(dim=Data$Ntime); TotHR=CPUE; ActualEffort=CPUE
  StorempaList=mpaList;  mpaList=c(0,0) # no mpas untl MPAStartYear
  for (t in 10:Data$Ntime) { #this is a loop over time

    #now harvest and move multiple steps
    Ctot=as.numeric(array(dim=ncell,0)) #will store the total catch over all moves
    Bstart=B[,t]
    BeginBio=Bstart
    BeginN=N[,t] #store numbers
    for (step in 1:Nsteps){
   # if  (t==45) browser()  
    
        if (t>Data$MPAStartYear)mpaList=StorempaList
    if (sum(mpaList)<1) {MaxPop=max(Bstart) #get MaxPop to scale boats B1 is intermediary calculation
    }  else {MaxPop=max(Bstart[-mpaList]) } #maximum pop outside MPAs
        
        
    B1=exp(-Data$Concent*(1-Bstart/MaxPop))*EffScale #intermediate calculation of boats
    #now get pop in and pop out and rescale boats
    if (sum(mpaList>1)) {
      B1[mpaList]=0; PopIn=sum(Bstart[mpaList]); PopOut=sum(Bstart[-mpaList])
    } else { PopIn=0; PopOut=sum(Bstart)}
    Boats=Data$Nboat*B1/sum(B1); Tpop=PopIn+PopOut
    BoatStore[,t]=Boats
    ActualEffort[t]=sum(Boats*Data$q)
    #  harvest rate
    hr=1-exp(-Boats*Data$q) #this is the default harvest rate for each area
    # if it is a harvest rate policy then calculate outside hr
    if (Data$Policy=="hr") {
     # if(CallBrowser==T) {browser()}
     
      DesiredCatch=Tpop*Data$TargetHR  #what the catch should be
      DesiredHR=DesiredCatch/max(PopOut,.1)  #make sure this doesn't go infinite
      if (DesiredHR>.8) DesiredHR=.8
      Effort=(ncell)*-log(1-DesiredHR)/Data$q
      Firsthr=1-exp(-Boats*Data$q) #using base q
      FirstCatch=sum(Firsthr*Bstart) #see what the catch would be
      qScale=(DesiredCatch/FirstCatch)
      Secondhr=1-exp(-Boats*Data$q*qScale)
      SecondCatch=sum(Secondhr*Bstart) #see what the catch would be
      qScale=qScale*(DesiredCatch/SecondCatch)
      Thirdhr=1-exp(-Boats*Data$q*qScale)
      ThirdCatch=sum(Thirdhr*Bstart) #see what the catch would be
      hr=Thirdhr
      if (Debug==T) print(paste(t,"Desired",round(DesiredCatch,0), "2nd",round(SecondCatch,0),
                  "3rd",round(ThirdCatch,0)))
      ActualEffort[t]=sum(Boats*Data$q*qScale)
    }  #end of hr policy
    if (Data$Policy=="effort") {hr=1-exp(-Boats*Data$q)} #don't really need this but it makes logic clear
  
    hr[which(hr>.8)]=.8
  

      Cat=hr/Nsteps*Bstart #catch for each cell in biomass
      Bstart=Bstart-Cat  #take catch away from biomass
      CatN= hr/Nsteps*BeginN #catch in numbers

      BeginN=BeginN-CatN

      Ctot=Ctot+(Cat) #Ctot is the total catch for each area
      #now move the adults multiply times dispersal matrix
      newpop=D %*% Bstart;  Bstart=newpop
      newN=D %*% BeginN;  BeginN=newN
      if (Debug==TRUE) plot(Bstart,type="l")
    }#end of loop over steps
    
#############################################################    
    
    
    
     #now readjust the harvest rate for each cell
    if (is.nan(sum(B[,t]))==T) browser()
    ActualHR=Ctot/B[,t] #harvest rate is the catch/the starting pop
    ActualHR[which(ActualHR>0.8)]=0.8
     B[,t]=Bstart
     N[,t]=BeginN
    HR[,t]=ActualHR 
     
     #need to dispers larvae and or recruits
     Eggs=B[,t]
     if (Data$Movement=="before") {  #movement happens before density dependence
       EggsMoved=LD %*% Eggs
       Recruitment[,t]=EggsMoved/(Data$Alpha+Data$Beta*EggsMoved)
     }
     if (Data$Movement=="after") { #movement happens after density dependent
       EggsMoved=Eggs
       RecruitmentBefore=Eggs/(Data$Alpha+Data$Beta*Eggs)
       Recruitment[,t]= LD %*% RecruitmentBefore
     }
      
    #TSurv=Data$Surv*(1-ActualHR)
     TSurv=Data$Surv #we have already removed harvest!
    N[,t+1]=TSurv*N[,t]+Recruitment[,t-Data$Lag]
    NewBiomass=TSurv*(Brodyalpha*N[,t]+Data$Rho*B[,t])+Data$Wrec*Recruitment[,t-Data$Lag]
    NewBiomass[which(NewBiomass<0)]=0
    B[,t+1]=NewBiomass
    catch =Ctot                      #the catch is pop * hr
 
    TotHR[t]=sum(catch)/sum(BeginBio)
    TotBio=apply(B,2,FUN=sum)
    CPUE[t]=sum(catch)/ActualEffort[t]
    
    if (Debug==TRUE) print(paste(t,"Catch",sum(catch),"Effort",ActualEffort[t]))
    
    Chist[t]=sum(catch)
  }  #end of loop over time
  
  return(list(pop=B[,Ntime],TotBio=TotBio,catch=catch,HR=HR,boats=Boats,
              BoatStore=BoatStore,CPUE=CPUE,TotHR=TotHR,
              B=B,N=N,Chist=Chist,MSY=Data$MSY,ActualEffort=ActualEffort,
              Bzero=Data$Bzero,BMSY=Data$BMSY,ncell=Data$ncell,mpaList=mpaList,
              MPAStartYear=Data$MPAStartYear,Ntime=Data$Ntime,NMPA=Data$NMPA,
              Recruitment=Recruitment,EggsMoved=EggsMoved,HR=HR,q=Data$q,D=D,
              Cellwidth=Cellwidth,EffScale=EffScale,mpaList=mpaList))
} #end of MPA function

############################################################

BaseMPA=function(){
  #
  #  This function sets initial values,  calls FindMSYBzero to get those parameters then calls evaluate
  #
  
  Info=FindMSYBzero(Surv,Wpre,Wrec,Rho,Rzero,Steep,Lag)
  Bzero=Info$Bzero; BMSY=Info$BMSY; MSY=Info$MSY
  UMSY=MSY/BMSY
  Alpha=(Bzero/Rzero)*(1-(Steep-0.2)/(0.8*Steep))
  Beta=(Steep-0.2)/(0.8*Steep*Rzero)
  
  cellid=c(1:ncell)
  EffScale<<-rep(1,ncell)
  
  TargetHR=min(UUMSY*UMSY,.8)
  q=-log(1-TargetHR)/(Nboat/(ncell*mean(EffScale)))  #q value to achieve UUMSY 
  
  Data<<-list(Surv=Surv,Wpre=Wpre,Wrec=Wrec,Rho=Rho,Rzero=Rzero,Steep=Steep,Lag=Lag,Bzero=Bzero,
              BMSY=BMSY, MSY=MSY,UMSY=UMSY,Alpha=Alpha,Beta=Beta,cellid=cellid,Ntime=Ntime,ncell=ncell,Nboat=Nboat,
              UUMSY=UUMSY,TargetHR=TargetHR,q=q,Policy=Policy,MPAStartYear=MPAStartYear,
              Concent=Concent,mpaWidth=mpaWidth,NMPA=NMPA,DispersalSD=DispersalSD,LarvalDispersalSD=LarvalDispersalSD,
              Movement=Movement,PlotResult=PlotResult)
  M=evaluate(Data)
  return(M)
} #end BaseMPA

######################################################



evaluate=function(Data)  {  #this function sets up the mpaList and the dispersal matrices then calls the MPA function
  # set up MPA grid
  if (Data$NMPA>0 ){
    nIn=Data$mpaWidth*Data$NMPA
    nOut=Data$ncell-nIn
    gaps=round(nOut/(Data$NMPA+1),0)
    mpaseq1=c(rep(0,gaps),rep(1,Data$mpaWidth))
    mpaseq=c(rep(mpaseq1,Data$NMPA),rep(0,Data$ncell-(length(mpaseq1)*Data$NMPA)))
    mpaList=which(mpaseq==1)
  } else  { mpaList=c(0,0)}
  #build dispersal matrix
  #adult dispersal
  x=seq(-Data$ncell/2,(Data$ncell/2)-1)
  PDisp=dnorm(x,mean=0,sd=Data$DispersalSD/Nsteps)
  PDisp=PDisp/sum(PDisp)
 
  D=array(dim=c(Data$ncell,Data$ncell))
  for (ix in 1:Data$ncell){
    iy=ix+x
    LT1=which(iy<1)
    iy[LT1]=Data$ncell+iy[LT1]
    GTncell=which(iy>Data$ncell)
    iy[GTncell]=iy[GTncell]-Data$ncell
    D[ix,iy]=PDisp
  }
  #Larval dispersal
  x=seq(-Data$ncell/2,(Data$ncell/2)-1)
  PDisp=dnorm(x,mean=0,sd=Data$LarvalDispersalSD)
  PDisp=PDisp/sum(PDisp)
  LD=array(dim=c(Data$ncell,Data$ncell))
  for (ix in 1:Data$ncell){
    iy=ix+x
    LT1=which(iy<1)
    iy[LT1]=Data$ncell+iy[LT1]
    GTncell=which(iy>Data$ncell)
    iy[GTncell]=iy[GTncell]-Data$ncell
    LD[ix,iy]=PDisp
  }
  M=MPAmodel(Data,mpaList,D,LD)
  return(M)
}  #end of evaluate function

######################################################


FindMSYBzero=function(Surv,Wpre,Wrec,Rho,Rzero,Steep,Lag){
  HR=rep(0,50)
  B=array(dim=50); Recruitment=array(dim=50); N=B
  B[1:10]=60; Recruitment[1:10]=Rzero
  N[1:10]=Rzero*Surv/(1-Surv)
  #simulate forward with Rzero recruitment to get Bzero
  for (y in 10:40){
  
    Recruitment[y]=Rzero
     N[y+1]=Surv*N[y]+Recruitment[y-Lag]
     B[y+1]=Surv*(Brodyalpha*N[  y]+Rho*B[y]) + Wrec*Recruitment[y-Lag]
    
 
  } #end of loop to get Bzero
  Bzero=B[y+1]
  #now simulate forward
  Alpha = (Bzero/Rzero)*(1-(Steep-0.2)/(0.8*Steep))
  Beta= (Steep-0.2)/(0.8*Steep*Rzero)
  B[1:10]=Bzero; Recruitment[1:10]=Rzero;
  HRList=seq(.01,.9,.01); NHR=length(HRList)
  Catch=array(dim=NHR); Biomass=Catch
  for (i in 1:NHR ){
    HR=HRList[i]
    for (y in 10:100){
      #changed version to make equivalent to multistep removal of catch before growth equations
      B[y]=B[y]*(1-HR)
      N[y]=N[y]*(1-HR)
      TSurv=Surv
      #TSurv=Surv*(1-HR)
      Recruitment[y]=B[y]/(Alpha+Beta*B[y])
      N[y+1]=TSurv*N[y]+Recruitment[y-Lag]
      B[y+1]=TSurv*(Brodyalpha*N[y]+Rho*B[y])+Wrec*Recruitment[y-Lag]
      if (B[y+1]<0) B[y+1]=1e-6
    } #end of loop over y
    Catch[i]=B[y+1]*HR
    Biomass[i]=B[y+1]
    # plot(HRList,Catch)
    # plot(HRList,Biomass)
  } #end of loop over hr
  j=which(Catch==max(Catch))
  BMSY=Biomass[j]; MSY=Catch[j]
  return(list(Bzero=Bzero, BMSY=BMSY, MSY=MSY))
} # end of function Bzero
#########################################################
shadeMPA=function(mpaList,MaxVal){  #shade the MPA plots to show what areas are closed
  d=which(mpaList[2:length(mpaList)]-mpaList[1:(length(mpaList)-1)]!=1)
  if (length(d)<1) {shadePoly(mpaList,MaxVal)}
  Start=mpaList[1]
  for (i in 1:(length(d)+1))  {
    End=mpaList[d[i]]
    if(i==(length(d)+1)) End=max(mpaList)
    List=seq(Start,End)
    # print(paste(i,Start,End))
    #  print(List)
    shadePoly(List,MaxVal)
    Start=mpaList[d[i]+1]
  }
} #end of function
shadePoly=function(List,MaxVal){
  minx=min(List)
  maxx=max(List)
  xseq=c(minx,maxx,maxx,minx,minx)
  yseq=c(0,0,MaxVal,MaxVal,0)
  polygon(x=xseq,y=yseq,density=10,col="red")
}

shadePoly2=function(List,MaxVal){
  minx=min(List)*Cellwidth
  maxx=max(List)*Cellwidth
  xseq=c(minx,maxx,maxx,minx,minx)
  yseq=c(0,0,MaxVal,MaxVal,0)
  polygon(x=xseq,y=yseq,density=10,col="red")
}
##################################################################

#This is he basic plot of equilibrium distribution of B and E and time trend
TwoPanel=function(){
  
  par(mfcol=c(2,1),mar=c(4,4,1,2))
  Btot=apply(M$B,2,sum)
  par(mfcol=c(2,1),mar=c(4,4,1,2))
  plot(x=seq(1,M$ncell),M$pop/M$Bzero,ylab="Abundance",xlab="Area",type="l",lwd=3,ylim=c(0,1))
  points(x=seq(1,M$ncell),M$boats/max(M$boats),col="purple",pch=19,cex=.8)
  if (M$NMPA>0) shadePoly(M$mpaList,1)
  #text(x=500,y=0.9,"MPA area in red shading")
  
  plot(M$Chist/(M$MSY*M$ncell),ylab="Catch Biomass and CPUE",xlab="Year",type="l",lwd=3,ylim=c(0,2))
  lines(Btot/(M$Bzero*M$ncell),col="red",lwd=3)
  lines(2*M$CPUE/max(M$CPUE,na.rm=TRUE),col="purple",lwd=3)
  lines(x=c(M$MPAStartYear,M$MPAStartYear),y=c(0,2),lwd=2,col="orange")
  
  
  text(x=M$MPAStartYear,y=1.5,paste("MPA starts in year",M$MPAStartYear))
  legend(x=70,y=2,cex=.8,legend=c("Catch","Total Biomass","CPUE"),col=c("black","red","purple"),lwd=c(3,3,3))
  
  par(mfcol=c(1,1),mar=c(4,4,1,2))
}

################################################

MeanDistance=function(){
  #base on midpoint
  Mid=round(Data$ncell/2,0)
  distance=abs(seq(1,Data$ncell)-Mid)
  MeanDist=sum(distance*M$D[Mid,])*Cellwidth
  return(MeanDist)
}
###########################################

#This is he basic plot of equilibrium distribution of B and E and time trend
TwoPanelBigeye=function(){
  
  par(mfcol=c(2,1),mar=c(4,4,1,2))
  Btot=apply(M$B,2,sum)
  par(mfcol=c(2,1),mar=c(4,4,1,2))
  plot(x=seq(1,M$ncell)*Cellwidth,M$pop/M$Bzero,ylab="Abundance",xlab="Distance (km)",type="l",lwd=3,ylim=c(0,1))
  points(x=seq(1,M$ncell)*Cellwidth,M$boats/max(M$boats),col="purple",pch=19,cex=.7)
  shadePoly2(M$mpaList,1)
  #text(x=50,y=0.9,"MPA area in red shading",cex=.8)
  
  
  YearstoPlot=seq(40,100)
  plot(YearstoPlot,M$Chist[YearstoPlot]/(M$MSY*M$ncell),ylab="Catch Biomass and CPUE",xlab="Year",
       type="l",lwd=3,ylim=c(0,2),xlim=c(min(YearstoPlot),max(YearstoPlot)))
  lines(YearstoPlot,Btot[YearstoPlot]/(M$Bzero*M$ncell),col="red",lwd=3)
  lines(YearstoPlot,2*M$CPUE[YearstoPlot]/max(M$CPUE,na.rm=TRUE),col="purple",lwd=3)
  lines(x=c(M$MPAStartYear,M$MPAStartYear),y=c(0,2),lwd=2,col="orange")
  
  
  text(x=M$MPAStartYear,y=1.5,paste("MPA starts in year",M$MPAStartYear))
  legend(x=70,y=2,cex=.8,legend=c("Catch","Total Biomass","CPUE"),col=c("black","red","purple"),lwd=c(3,3,3))
  par(mfcol=c(1,1),mar=c(4,4,1,2))
  plot(M$HR[,100])
}

##########################################

GetEffortScalar=function(){
  E <- read_excel("BigeyeEffortvsDistance.xlsx", sheet = "Effort") 
  EffScale=array(dim=ncell)
  FirstMPA=ncell/2-mpaWidth/2
  LastMPA=ncell/2+mpaWidth/2
  InMPA=seq(FirstMPA,LastMPA)
  for (i in 1:ncell){
    if (i %in% InMPA){ EffScale[i]= E$RelativeEffort[1]
    } else {
      D=min(abs(i-InMPA)*Cellwidth)
      D2=abs(E$Distance-D)
      D3=min(D2)
      j=which(D2==D3)  #this will be the closest distance
      EffScale[i]=E$RelativeEffort[j]
      
    } #end of else
  } #end of for
  return(EffScale)
}

#end of function

CalcMoved=function(M,D){
  B=rep(0,ncell)
  B[ncell/2]=1000
  plot(B,type="l")
  x=seq(1,ncell)
  for (step in 1:Nsteps){
    
    newpop=M$D %*% B
    B=newpop
    points(x,B,type="l")
  }
  distance=abs(x-ncell/2)
  avg=sum(distance*B)/sum(B)
  #print(paste("average cells moved ",round(avg,2)))
  return(avg)
}


######################################################