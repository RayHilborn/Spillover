
jpeg(file="Figure 3.jpg",width=7,height=7,units="in",res=300)
par(mfrow=c(3,3),mar=(c(0,0,0,0)),oma=c(4,3.5,0.1,0.1))
DispList=c(3,10,20) #dispersal rates
UMSYList=c(.5,1,2)
Adisp=3; UUMSY=1
Policy="effort"
#Policy="hr"


for (Adisp in DispList){
  for (UUMSY in UMSYList){
    
load(file=paste0("Model ",Adisp," ",UUMSY," ", Policy, ".rdata"))

plot(x=seq(1,M$ncell),M$pop/M$Bzero,ylab="Abundance",xlab="Area",type="l",lwd=3,
     ylim=c(0,1.1),xaxt="n",yaxt="n")
#shadeMPA(M$mpaList,1)
points(x=seq(1,M$ncell),M$boats/max(M$boats),col="grey",pch=19,cex=.5)

lines(x=rep(M$mpaList[1],2)+.1,y=c(0,1),lwd=1,col="black")
lines(x=rep(M$mpaList[length(M$mpaList)],2)-.1,y=c(0,1),lwd=1,col="black")

#x=expression(u/u[MSY] )
x <- substitute(paste(u/u[MSY], "=", UUMSY), list(UUMSY=UUMSY)) 
text(140,1,x,adj=0,cex=1)
x <- substitute(paste(sigma, "=", Adisp), list(Adisp=Adisp)) 
 #x=(expression(sigma ~      " =" ))
 text(1,1,x,adj=0,cex=1)


}} #end of loops over models

mtext("Location", side = 1,  outer = TRUE,line=2,cex=1.5)
mtext("Abundance and boats", side = 2,  outer = TRUE,line=2,cex=1.5)

dev.off()
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1,oma=c(2,2,2,2))


