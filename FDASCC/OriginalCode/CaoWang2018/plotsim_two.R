read.csv("covAR1.csv", header=T)->ar1

read.csv("covtoep.csv", header=T)->toep

seq(0,0.27,0.03)->x


plot(cbind(x,ar1[1:10,3]),xlim=c(0,0.27),ylim=c(0,1), cex.axis=1.5,cex.lab=1.5,xlab=expression(delta),type="n")

lines(cbind(x,ar1[1:10,3]),col=3,lty=2,lwd=4)

lines(cbind(x,ar1[1:10,2]),col=2,lty=1,lwd=4)


plot(cbind(x,toep[1:10,3]),xlim=c(0,0.27),ylim=c(0,1), cex.axis=1.5,cex.lab=1.5,xlab=expression(delta),type="n")

lines(cbind(x,toep[1:10,3]),col=3,lty=2,lwd=4)

lines(cbind(x,toep[1:10,2]),col=2,lty=1,lwd=4)
#title("TOEP structure")