Dmatrix2<-function(n,Order2)
 { 
  	             ## knots number (Ns1>>Ns2)
  Ns2=floor(2*n^(1/(2*(Order2+1)))*log(log(n)))
  s=1
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
 return(D2=D2)
 }