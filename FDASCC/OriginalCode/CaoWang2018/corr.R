########################################################
## define a covariance matrix ***AR(1)*** for epi, kxi  
########################################################
cor_matrix<-function(row, col, rho, structure){       
	cor_m=matrix(0,row, col)
	if(structure=="AR1")  
	cor_m= toeplitz(c(1,as.vector(rho)^seq(1,(row-1))))                                                    ### No.1   AR(1)
	if(structure=="EX")  
            cor_m= toeplitz(c(1,rep(rho,(row-1))))       ### No.2   Exchangeable
	if(structure=="BAND")         
		cor_m= toeplitz(c(1,rho,rep(0,(row-2))))     ###  No.3   banded
	if(structure=="TOEP")   
		cor_m= toeplitz(c(1,rho))                    ### No.4 Toeplitz
	if(structure=="IND")   
		cor_m= diag(row)                             ### No.5   Independent
  
    return(cor_m)
}



  