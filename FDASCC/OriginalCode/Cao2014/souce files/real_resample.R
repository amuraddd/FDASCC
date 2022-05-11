  
 rm(list=ls(all=TRUE))
source("corr.R")
source("xknots2.R")
source("Ghat.R")
source("generateQ.R")
source("cov_rate.R")
source("var_zeta.R")
source("rho_est.R")

  Kappa=kappa=2
  boot=1000	### generate supermum
  bootnum=1000 
  alpha=c(0.05,0.01)
  Order1=Order2=3
  struct=c( "IND","EX", "AR1","TOEP")  #c("IND", "EX", "AR1","TOEP")
  conboot=array(0,c(bootnum,2,length(struct)))
  choose_str=matrix(0,length(struct),2)
  index=matrix(0,1,1)
  covrate=array(0,c(length(struct),1,length(alpha)))



 library(MASS)
 library(splines)

   f <- file.path("E:/repeated FDA/bltper_1x1", 
  c("AUS.bltper_1x1.txt", "AUT.bltper_1x1.txt","BLR.bltper_1x1.txt",
    "BEL.bltper_1x1.txt", "BGR.bltper_1x1.txt", "BGR.bltper_1x1.txt",
    "CAN.bltper_1x1.txt", "CHL.bltper_1x1.txt", "CZE.bltper_1x1.txt",
    "DNK.bltper_1x1.txt", "EST.bltper_1x1.txt", "FIN.bltper_1x1.txt", 
    "FRATNP.bltper_1x1.txt", "DEUTE.bltper_1x1.txt", "DEUTW.bltper_1x1.txt",
    "HUN.bltper_1x1.txt","ISL.bltper_1x1.txt",  "IRL.bltper_1x1.txt",
    "ISR.bltper_1x1.txt", "ITA.bltper_1x1.txt", "JPN.bltper_1x1.txt",
    "LVA.bltper_1x1.txt", "LTU.bltper_1x1.txt", "LUX.bltper_1x1.txt",
    "NLD.bltper_1x1.txt",  "NZL_MA.bltper_1x1.txt", "NZL_NM.bltper_1x1.txt",
    "NOR.bltper_1x1.txt","POL.bltper_1x1.txt",  "PRT.bltper_1x1.txt",
    "RUS.bltper_1x1.txt", "SVK.bltper_1x1.txt",   "ESP.bltper_1x1.txt",
    "SWE.bltper_1x1.txt",  "TWN.bltper_1x1.txt",
    "GBRTENW.bltper_1x1.txt", "GBRCENW.bltper_1x1.txt",  "GBR_SCO.bltper_1x1.txt",
    "GBR_NIR.bltper_1x1.txt", "UKR.bltper_1x1.txt",  "USA.bltper_1x1.txt"))

  lapply(f, read.table)->d


  n=41
  N=Ntest=111
  J=14  ### repeated measurements number
  YY=array(0,c(J,N,n))

  d[[1]]->aus
  as.numeric(as.character(d[[1]][[3]])[which( as.numeric(as.character(d[[1]][[1]][-1]))<=2005 & as.numeric( as.character(d[[1]][[1]][-1]))>=1992 )])->new.aus
  matrix(new.aus, byrow=T,nrow=length(new.aus)/111,ncol=111)->start->YY[,,1]


     
 
  for (i in 2:n)

  {
  
   as.numeric(as.character(d[[i]][[3]][-1])[which( as.numeric(as.character(d[[i]][[1]][-1]))<=2005 & as.numeric( as.character(d[[i]][[1]][-1]))>=1992 )])->new.aus
  
  print(c(length(new.aus), i))

   matrix(new.aus, byrow=T,nrow=J,ncol=111)->new->YY[,,i]

   rbind(new,start)->start

  }

  plot(cbind(seq(1,111),colMeans(start)))

  begin=60
  end=100
  
  YY[,c(begin:end),]->YY
  N=Ntest=end-begin+1

   
  YY->Yplot
  
  par(mgp = c(2.5, 1, 0))
  plot(Yplot[1,,1],type="n",xlim=c(begin,end),ylim=c(-0.01,max(Yplot[1,,])), xlab="Age", ylab=paste("Mortality Rates in ",1992+2),cex.axis=1.1,cex.lab=1.1,cex.main=1.2)
  
   for (i in 1:n)
  lines(cbind(seq(begin,end,by=1),Yplot[1,,i]),col=4,lty=1,lwd=1)
   

  J=14
   for(j in 1:J)
   { pdf(paste("rate_",j,".pdf", sep=""))
   
   par(mgp = c(2.5, 1, 0))
   plot(Yplot[1,,1],type="n",xlim=c(begin,end),ylim=c(0,max(Yplot[,,1])), xlab="Age", ylab="",cex.axis=1.1,cex.lab=1.1,cex.main=1.2)
   
   title(paste("Mortality Rates from 1992 to", 1991+j), col.main=j)
   for (jj in 1:j)
    {lines(cbind(seq(begin,end,by=1),Yplot[jj,,1]),col=jj,lty=1,lwd=2)}

    dev.off()

   }
  

  { for(j in 1:J)
     {lines(cbind(seq(begin,end,by=1),Yplot[j,,1]),col=j,lty=1,lwd=1)
   }
  }


##############################################
####			start simulation		####
##############################################

		
	#############################################
	#####			knots matrix		###
	#############################################
	 
	knots=knots2(N,Ntest,order1,order2,n,J)
	XB1=knots$XB1
	beta=knots$beta
	Beta=knots$Beta
	xknots=knots$xknots
	xtestknots=knots$xtestknots
      XBtest1=knots$XBtest1

	#############################################
	###			generate data		###
	#############################################
         	      
              
               	Y=matrix(apply(YY,2,mean),N,1)

	 		#############################################
	 		###		start estimation mhat		###
	 		#############################################
			beta.est=beta%*%Y
			mhat.data=XB1%*%beta.est
			mhat.test=XBtest1%*%beta.est         

			#############################################
			###	start estimation Ghat (jj1,jj2)	###
			#############################################
			Gest=Ghat(YY,mhat.data,Beta,xknots,xtestknots,Ntest)
			phiavgY.data=Gest$phiavgY
			lamdaK.data=Gest$lamdaK
			phiavg.data=Gest$phiavg
                  
			
				Mtest<-mhat.test
				mhat.fix<- mhat.data
              		for(str in 1:length(struct)){
					ST=struct[str]
					rho.est=rho_est(YY,ST) 
					var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
					Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha)
					Q=Q.est$Q
					#covM=Q.est$covM

					###################################
					# bootstrapping generate new data #
					###################################
	                        YYnew=array(0, c(J,N,n))	
					for(new in 1:bootnum){                    
						
						boot.index=sample(1:n,n,replace=T)
						YYnew=YY[,,boot.index]
						
                        		Ynew=matrix(apply(YYnew,2,mean),N,1)
						beta.est=beta%*%Ynew
						mhat=XB1%*%beta.est
						mhat.test=XBtest1%*%beta.est                                
						Gest=Ghat(YYnew,mhat,Beta,xknots,xtestknots,Ntest)
						lamdaK=Gest$lamdaK
						phiavg=Gest$phiavg

						rho.est=rho_est(YYnew,ST)
						var.zeta=var_zeta(phiavg,lamdaK,rho.est,J,ST)
						cov.rate=cov_rate(mhat.test,Mtest,var.zeta,Q,ST,alpha,n)
						conboot[new,,str]=cov.rate
                    		}
                            
					choose_str[str,]=abs(colMeans(conboot[ , , str])-(1-alpha)) 
				} 
                          aa=which(rowMeans(choose_str)==min(rowMeans(choose_str))) 
                        index=min(aa)
                         
                        index=2
	                 
				ST=struct[index]
				rho.est=rho_est(YY,ST)
                        var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
				Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha)
				Q=Q.est$Q
				covM=Q.est$covM

				 mhat.test+Q[1]*(var.zeta)^0.5/(n)^0.5->uptest1
                         mhat.test-Q[1]*(var.zeta)^0.5/(n)^0.5->lowtest1

                         mhat.test+Q[2]*(var.zeta)^0.5/(n)^0.5->uptest2
                         mhat.test-Q[2]*(var.zeta)^0.5/(n)^0.5->lowtest2



 
  plot(mhat.test,type="n",xlim=c(begin,end), xlab="Age", ylab=expression(hat(mu)),cex.axis=1,cex.lab=1.1,cex.main=1.2)
  lines(cbind(seq(begin,end,by=1),lowtest2),col=2,lty=2,lwd=3)
  lines(cbind(seq(begin,end,by=1),uptest2),col=2,lty=2,lwd=3)
  lines(cbind(seq(begin,end,by=1),mhat.test),col=1,lty=1,lwd=1)

  
	                  ST=struct[1]
				rho.est=rho_est(YY,ST)
                        var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
				Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha)
				Q=Q.est$Q
				covM=Q.est$covM

				 mhat.test+Q[1]*(var.zeta)^0.5/(n)^0.5->uptest1
                         mhat.test-Q[1]*(var.zeta)^0.5/(n)^0.5->lowtest1

                         mhat.test+Q[2]*(var.zeta)^0.5/(n)^0.5->uptest2
                         mhat.test-Q[2]*(var.zeta)^0.5/(n)^0.5->lowtest2

  lines(cbind(seq(begin,end,by=1),lowtest2),col=3,lty=4,lwd=2)
  lines(cbind(seq(begin,end,by=1),uptest2),col=3,lty=4,lwd=2)






        