Dmatrix1<-function(n,Order1,Ns1)
 { 
  #Ns1=floor(1*log(n)*n^(1/(2*(Order1+1))))	## knots number (Ns1>>Ns2)
  #Ns2=floor(2*n^(1/(2*(Order2+1)))*log(log(n)))


  p=Order1+1
  D=diag(-1,Ns1+p)
  D[row(D) == col(D) + 1]<- 1
  D[,-(Ns1+p)]->D1 
  t(D1)->A1
  
  s=1
  for (i in 1:(p-s-1))
  {
   A1[i,]*(p-s)/(p-(p-s-i+1))->A1[i,]
   A1[(Ns1+p-i-s+1),]*(p-s)/(p-(p-s-i+1))->A1[(Ns1+p-i-s+1),]
  }
  t(A1)->D1


  return(D1=D1)
 }