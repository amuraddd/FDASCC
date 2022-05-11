var_zeta<- function(phiavg,lamdaK,rho.est,J,structure){
      kappa=length(lamdaK)
      Ntest=dim(phiavg)[1]
	Gjj=matrix(0,Ntest,1)
	for(ka1 in 1:kappa){
		Gjj=Gjj+diag((phiavg[,ka1]%*%t(phiavg[,ka1]))*lamdaK[ka1])
	} 
	if(structure=="AR1"){a=(J+sum(2*seq((J-1),1)*rho.est^seq(1,(J-1)))) /J^2}
	if(structure=="EX") {a=(J+sum(2*seq((J-1),1)*rho.est))/J^2}
	if(structure=="BAND"){a=(J+rho.est*2*(J-1))/J^2}
	if(structure=="TOEP"){a=(J+sum(2*seq((J-1),1)*rho.est))/J^2}
	if(structure=="IND"){a=1/J}
	var.zeta=a*Gjj
	return(var.zeta)
}