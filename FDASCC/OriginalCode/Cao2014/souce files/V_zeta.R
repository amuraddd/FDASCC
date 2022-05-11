V_zeta<- function(phiavgv,lamdaK,rho.est,J,structure){
      kappa=length(lamdaK)
      Ntest=dim(phiavgv)[2]
	Vjj= matrix(0,Ntest,1)
	for(ka1 in 1:kappa){
	  Vjj=Vjj + diag((as.matrix(c(phiavgv[ka1,]))%*%t(as.matrix(phiavgv[ka1,])))*lamdaK[ka1])
	  #Vjj=Vjj+matrix(diag(t(phiavgv[ka1,])%*%phiavgv[ka1,]*lamdaK[ka1]), Ntest,1)
       
	} 

	if(structure=="AR1"){a=(J+sum(2*seq((J-1),1)*rho.est^seq(1,(J-1)))) /J^2}
	if(structure=="EX") {a=(J+sum(2*seq((J-1),1)*rho.est))/J^2}
	if(structure=="BAND"){a=(J+rho.est*2*(J-1))/J^2}
	if(structure=="TOEP"){a=(J+sum(2*seq((J-1),1)*rho.est))/J^2}
	if(structure=="IND"){a=1/J}
	V.zeta=a*Vjj

     	return(Vjj)
}