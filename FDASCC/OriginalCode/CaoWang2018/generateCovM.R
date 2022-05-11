  generateCovM<-function(J,rho.est,structure,phiavg,lamdaK,var.zeta){
	
      kappa=length(lamdaK)
      Ntest=dim(phiavg)[1]
	if(structure=="AR1"){
		covM=cor_matrix(J,J,rho.est,"AR1")
      } 
	if(structure=="EX"){
		covM=cor_matrix(J,J,rho.est,"EX")
	}
	if(structure=="BAND"){ 
		covM=cor_matrix(J,J,rho.est,"BAND")
	}
      if(structure=="TOEP"){
		covM=toeplitz(c(1,rho.est))
     }
	if(structure=="IND"){covM=diag(J)}



	list(covM=covM)
}