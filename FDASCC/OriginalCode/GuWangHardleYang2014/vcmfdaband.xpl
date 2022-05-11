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
nsim=1

;randomize(seed)

proc(ni)=ob(n)
u=uniform(n,1)
;r=1/10
;ni=2.*(u.<=r)+3.*(u.>r).*(u.<=2*r)+4.*(u.>2*r).*(u.<=3*r)+5.*(u.>3*r).*(u.<=4*r)+6.*(u.>4*r).*(u.<=5*r)+7.*(u.>5*r).*(u.<=6*r)+8.*(u.>6*r).*(u.<=7*r)+9.*(u.>7*r).*(u.<=8*r)+10.*(u.>9*r)
;ni=ni+3
ni=floor(13*u)+2
endp

proc(data,tdata,fp,nt,nmax,m,G,ni,sigmay)=dt(n,coeff)
ni=ob(n)
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
G1=ph11^2+ph21^2

k2=1
;m2=cos(2*pi*(x-0.5))
m2=5*(x-0.6)^2
ph12=1
ph22=k2*sqrt(2)*sin(2*pi*x)
ph32=k2*sqrt(2)*cos(2*pi*x)
G2=matrix(nt,1).*ph12^2+ph22^2+ph32^2
G=G1~G2
sigmay=G*(tdata'.^2)+(coeff^2)*matrix(nt,1)

;y=(m1+xi11.*ph11+xi21.*ph21).*t1+(m2+xi12.*ph12).*t2+e*coeff
y=(m1+xi11.*ph11+xi21.*ph21).*t1+(m2+xi12.*ph12+xi22.*ph22+xi32.*ph32).*t2+e*coeff
data=x~y~t1~t2
m=m1~m2


a=min(x)
b=max(x)
length=b-a
int=#(0:100)
fp=a+int*length/100
endp

proc(mlhat,resi,Ns,hs,x)=ml(data,n,nmax,nt)
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
;fpknots=(fp.>=knotsl[,1:Ns]).*(fp.<knotsu[,1:Ns])~(fp.>=knotsl[,Ns+1]).*(fp.<=knotsu[,Ns+1])

bx=kron(matrix(1,d),xknots)
T=kron(data[,3:d+2],matrix(1,Ns+1))
D=bx.*T
lamda=ginv(D'*D)*D'*y
lamdal=reshape(lamda,#(Ns+1,d))
mlhat=xknots*lamdal
;mlfphat=fpknots*lamdal
resi=(y-D*lamda).^2
;mlhat


endp


proc(Glhat,f,ff,sigmal,sigmal2)=Gl(data,tdata,fp,resi,n,nt,hs,ni,sigmay)
x=data[,1]
y=data[,2]
d=2
a=min(x)
b=max(x)
length=b-a
N=floor(0.5*n^(1/3))
h=length/(N+1)
knotsl=a+h*(0:N)'
knotsu=a+h*(1:N+1)'
xknots=(x.>=knotsl[,1:N]).*(x.<knotsu[,1:N])~(x.>=knotsl[,N+1]).*(x.<=knotsu[,N+1])
fpknots=(fp.>=knotsl[,1:N]).*(fp.<knotsu[,1:N])~(fp.>=knotsl[,N+1]).*(fp.<=knotsu[,N+1])

f=(1/(nt*h))*xknots*xknots'*matrix(nt,1)
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
Glhat=xknots*Glrho
rehat=xknots*rhol[,d+1]
;Glhat~rehat

Gt=Glhat*(tdata'.^2)
sigmayhat=Gt+rehat
ss=sigmay~sigmayhat
ss=reshape(ss,#(nt*n,2))
ss=corr(ss)
ss
fsig=(nt*hs*f)^(-1).*sigmayhat
fsig=reshape(fsig',#(n,1,nt))
tfsig=tdata.*fsig
tfsig=reshape(tfsig,#(n,d*nt))
ttfsig=n^(-1)*tdata'*tfsig
ttfsig=reshape(ttfsig,#(d,d,nt))

Qn=n^(-1)*tdata'*tdata
invQn=inv(Qn)
invQn=kron(matrix(1,nt),invQn)
invQn=reshape(invQn,#(d,d,nt))
sigma=invQn*ttfsig*invQn
sigmal=xdiag(sigma) 
sigmal=reshape(sigmal,#(d,nt))
sigmal=sigmal'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
;sigmal

fGt=(f*hs).*Gt./sigmayhat
fGt=reshape(fGt',#(n,1,nt))
tfsig2=(fsig.*(1+(sum(ni.^2)/nt-1).*fGt)).*tdata
tfsig2=reshape(tfsig2,#(n,d*nt))
ttfsig2=n^(-1)*tdata'*tfsig2
ttfsig2=reshape(ttfsig2,#(d,d,nt))
sigma2=invQn*ttfsig2*invQn
sigma12=xdiag(sigma2)                                                       
sigma12=reshape(sigma12,#(d,nt))
sigmal2=sigma12'
;sigmal2


endp

{data,tdata,fp,nt,nmax,m,G,ni,sigmay}=dt(n,coeff)
{mlhat,resi,Ns,hs,x}=ml(data,n,nmax,nt)
{Glhat,f,ff,sigmal,sigmal2}=Gl(data,tdata,fp,resi,n,nt,hs,ni,sigmay)


/*confidence band*/
dn=0.5*(log(log(Ns+1))+log(4*pi))
dn=log(-0.5*log(1-alpha))+dn
dn=1-(2*log(Ns+1))^(-1)*dn
;lm=dn*sqrt(2*log(Ns+1)).*(sigmal.^(0.5))
lm=dn*sqrt(2*log(Ns+1)).*(sigmal2.^(0.5))
l=mlhat-lm
u=mlhat+lm

lmp=qfn(1-alpha/2).*(sigmal2.^(0.5))
lp=mlhat-lmp
up=mlhat+lmp


dd=createdisplay(1,2)
;dd1=createdisplay(1,1)
lb1=setmask(sort(x~l[,1]), "line", "red", "size", 3)
ub1=setmask(sort(x~u[,1]), "line", "red" , "size", 3)
lpw1=setmask(sort(x~lp[,1]), "line","black", "dotted", "size", 2)
upw1=setmask(sort(x~up[,1]), "line","black", "dotted", "size", 2)
m1truep=setmask(sort(x~m[,1]), "line", "blue", "size", 6)
m1estp=setmask(sort(x~mlhat[,1]),"line", "magenta", "size", 3)
lb2=setmask(sort(x~l[,2]), "line", "red", "size", 3)
ub2=setmask(sort(x~u[,2]), "line", "red" , "size", 3)
lpw2=setmask(sort(x~lp[,2]), "line","black", "dotted", "size", 2)
upw2=setmask(sort(x~up[,2]), "line","black", "dotted", "size", 2)
m2truep=setmask(sort(x~m[,2]), "line", "blue", "size", 6)
m2estp=setmask(sort(x~mlhat[,2]),"line", "magenta", "size", 3)
data=setmask(data,"black","small", "circle")

G1true=setmask(sort(x~G[,1]), "line", "purple")
G1est=setmask(sort(x~Glhat[,1]),"line", "red", "dashed")
fest=setmask(sort(fp~ff),"line", "red")

show(dd,1,1,lb1,ub1,lpw1,upw1,m1truep,m1estp)
;show(dd1,1,1,lb1,ub1,lpw1,upw1,m1truep,m1estp)
show(dd,1,2,lb2,ub2,lpw2,upw2,m2truep,m2estp)
;show(dd1,1,1,G1true,G1est)
;show(dd1,1,2,G2true,G2est)
;show(dd,1,1,m1truep,m1estp)
;show(dd,1,2,m2truep,m2estp)
;show(dd1,1,1,fest)
setgopt(dd,1,1,"ylim",-4|4,"ymajor",1)
setgopt(dd,1,2,"ylim",-5|5,"ymajor",1)


