library("distribs")
library("smoother")
library("stats")
library("graphic") 
library("xplore")
library("multi")

seed=394678
alpha=0.05
n=200
coeff=0.5
d=2
nsim=500

proc(mlfhat,Ns,hs,x)=ml(data,fp,n,nmax,nt,nf)
x=data[,1]
y=data[,2]
d=2
a=min(x)
b=max(x)
length=b-a
Ns=floor(0.5*nt^(1/3)*log(n))
hs=length/(Ns+1)
knotsl=a+hs*(0:Ns)'
knotsu=a+hs*(1:Ns+1)'
xknots=(x.>=knotsl[,1:Ns]).*(x.<knotsu[,1:Ns])~(x.>=knotsl[,Ns+1]).*(x.<=knotsu[,Ns+1])
fpknots=(fp.>=knotsl[,1:Ns]).*(fp.<knotsu[,1:Ns])~(fp.>=knotsl[,Ns+1]).*(fp.<=knotsu[,Ns+1])

bx=kron(matrix(1,d),xknots)
T=kron(data[,3:d+2],matrix(1,Ns+1))
D=bx.*T
lamda=ginv(D'*D)*D'*y
lamdal=reshape(lamda,#(Ns+1,d))
;mlfhat=xknots*lamdal
mlfhat=fpknots*lamdal
;mlhat
endp


proc(Glfhat,ff,refhat)=Gl(data,fp,n,nt,ni)
x=data[,1]
y=data[,2]
d=2
a=min(x)
b=max(x)
length=b-a
;N=floor(nt^(1/3)*log(n))
N=floor(0.5*n^(1/3))
h=length/(N+1)
knotsl=a+h*(0:N)'
knotsu=a+h*(1:N+1)'
xknots=(x.>=knotsl[,1:N]).*(x.<knotsu[,1:N])~(x.>=knotsl[,N+1]).*(x.<=knotsu[,N+1])
fpknots=(fp.>=knotsl[,1:N]).*(fp.<knotsu[,1:N])~(fp.>=knotsl[,N+1]).*(fp.<=knotsu[,N+1])

;f=(1/(nt*h))*xknots*xknots'*matrix(nt,1)
ff=(1/(nt*h))*fpknots*xknots'*matrix(nt,1)

bx=kron(matrix(1,d),xknots)
T=kron(data[,3:d+2],matrix(1,N+1))
D=bx.*T
lamda=ginv(D'*D)*D'*y
resi=(y-D*lamda).^2
bbx=kron(matrix(1,d+1),xknots)
TT=kron(data[,3:d+2].^2~matrix(nt,1),matrix(1,N+1))
DD=bbx.*TT
rho=ginv(DD'*DD)*DD'*resi
rhol=reshape(rho,#(N+1,d+1))
Glrho=rhol[,1:d]
;Glhat=xknots*Glrho
Glfhat=fpknots*Glrho
;rehat=xknots*rhol[,d+1]
refhat=fpknots*rhol[,d+1]

endp


proc(sigmalf,sigmal2f)=sigf(tdata,Glfhat,refhat,ff,n,nt,nf,ni,hs)
d=2
Gt=Glfhat*(tdata'.^2)
sigmayhat=Gt+refhat
fsig=(nt*hs*ff)^(-1).*sigmayhat
fsig=reshape(fsig',#(n,1,nf))
tfsig=tdata.*fsig
tfsig=reshape(tfsig,#(n,d*nf))
ttfsig=n^(-1)*tdata'*tfsig
ttfsig=reshape(ttfsig,#(d,d,nf))

Qn=n^(-1)*tdata'*tdata
invQn=inv(Qn)
invQn=kron(matrix(1,nf),invQn)
invQn=reshape(invQn,#(d,d,nf))
sigma=invQn*ttfsig*invQn
sigmal=xdiag(sigma) 
sigmal=reshape(sigmal,#(d,nf))
sigmalf=sigmal'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      


fGt=(ff*hs).*Gt./sigmayhat
fGt=reshape(fGt',#(n,1,nf))
tfsig2=(fsig.*(1+(sum(ni.^2)/nt-1).*fGt)).*tdata
tfsig2=reshape(tfsig2,#(n,d*nf))
ttfsig2=n^(-1)*tdata'*tfsig2
ttfsig2=reshape(ttfsig2,#(d,d,nf))
sigma2=invQn*ttfsig2*invQn
sigma12=xdiag(sigma2)                                                       
sigma12=reshape(sigma12,#(d,nf))
sigmal2f=sigma12'

endp

proc(per1,per2)=simu(coeff,n,alpha,nsim,seed,d)

per1=matrix(nsim,1)
per2=matrix(nsim,1)
randomize(seed)

isim=1

while(isim<=nsim)
u=uniform(n,1)
ni=floor(13*u)+2
nt=sum(ni)
nmax=max(ni)

x=uniform(nt,1)
e=normal(nt,1)
t1n=normal(n,1)
t2n=randbin(1,0.5,n)
;t1n=(t1n-0.5)*2
;t2n=normal(n,1)
xi11=normal(n,1)
xi21=normal(n,1)
xi12=normal(n,1)
xi22=normal(n,1)
xi32=normal(n,1)
tdata=t1n~t2n
r=matrix(1,nmax)
xi11=xi11.*r
xi21=xi21.*r
xi12=xi12.*r
xi22=xi22.*r
xi32=xi32.*r
t1=t1n.*r
t2=(t2n+1).*r
;t2=t2n.*r
id=#(1:nmax)
id=matrix(n,1).*id'
id=id.<=ni
xi11=xi11.*id
xi21=xi21.*id
xi12=xi12.*id
xi22=xi22.*id
xi32=xi32.*id
t1=t1.*id
t2=t2.*id

xi11=reshape(xi11',#(nmax*n,1))
xi21=reshape(xi21',#(nmax*n,1))
xi12=reshape(xi12',#(nmax*n,1))
xi22=reshape(xi22',#(nmax*n,1))
xi32=reshape(xi32',#(nmax*n,1))
t1=reshape(t1',#(nmax*n,1))
t2=reshape(t2',#(nmax*n,1))

xi11=paf(xi11,abs(xi11)>0)
xi21=paf(xi21,abs(xi21)>0)
xi12=paf(xi12,abs(xi12)>0)
xi22=paf(xi22,abs(xi22)>0)
xi32=paf(xi32,abs(xi32)>0)
t1=paf(t1,abs(t1)>0)
;t1=t1-1
t2=paf(t2,abs(t2)>0)
t2=t2-1

k1=1/sqrt(5)
m1=sin(2*pi*(x-0.5))
ph11=k1*(-2*cos(pi*(x-0.5)))
ph21=k1*sin(pi*(x-0.5))

k2=1
m2=5*(x-0.6)^2
ph12=1
ph22=k2*sqrt(2)*sin(2*pi*x)
ph32=k2*sqrt(2)*cos(2*pi*x)

y=(m1+xi11.*ph11+xi21.*ph21).*t1+(m2+xi12.*ph12+xi22.*ph22+xi32.*ph32).*t2+e*coeff
;y=(m1+xi11.*ph11+xi21.*ph21).*t1+(m2+xi12.*ph12).*t2+e*coeff
data=x~y~t1~t2
m=m1~m2

a=min(x)
b=max(x)
length=b-a
int=#(0:100)
fp=a+int*length/100
nf=rows(fp)
m1fp=sin(2*pi*(fp-0.5))
m2fp=5*(fp-0.6)^2


{mlfhat,Ns,hs,x}=ml(data,fp,n,nmax,nt,nf)
{Glfhat,ff,refhat}=Gl(data,fp,n,nt,ni)
{sigmalf,sigmal2f}=sigf(tdata,Glfhat,refhat,ff,n,nt,nf,ni,hs)


/*confidence band*/
dn=0.5*(log(log(Ns+1))+log(4*pi))
dn=log(-0.5*log(1-alpha))+dn
dn=1-(2*log(Ns+1))^(-1)*dn
;lmf=dn*sqrt(2*log(Ns+1)).*(sigmalf.^(0.5))
lmf=dn*sqrt(2*log(Ns+1)).*(sigmal2f.^(0.5))

ma1=max(abs(mlfhat[,1]-m1fp)./lmf[,1])
per1[isim]=(ma1<=1)
ma2=max(abs(mlfhat[,2]-m2fp)./lmf[,2])
per2[isim]=(ma2<=1)
isim=isim+1
endo
endp

{per1,per2}=simu(coeff,n,alpha,nsim,seed,d)
mean(per1)
mean(per2)


