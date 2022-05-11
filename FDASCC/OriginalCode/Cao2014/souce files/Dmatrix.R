Dmatrix<-function(n,Order1,Order2)
 { 
  Ns1=floor(1*log(n)*n^(1/(2*(Order1+1))))	## knots number (Ns1>>Ns2)
  Ns2=floor(2*n^(1/(2*(Order2+1)))*log(log(n)))


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

  p=Order2+1
  D=diag(-1,Ns2+p)
  D[row(D) == col(D) + 1]<- 1
  D[,-(Ns2+p)]->D2 
  t(D2)->A2
  

  for (i in 1:(p-s-1))
  {
   A2[i,]*(p-s)/(p-(p-s-i+1))->A2[i,]
   A2[(Ns2+p-i-s+1),]*(p-s)/(p-(p-s-i+1))->A2[(Ns2+p-i-s+1),]
  }
  t(A2)->D2
  list(D1=D1,D2=D2)
 }