library("distribs")
library("smoother")
library("stats")
library("graphic") 
library("xplore")
library("multi")
library("plot")

seed=39467800
alpha=0.05


d=4
nsim=1

;randomize(seed)

;data=read("E:\vcmfdaband\cd4.txt")
data=read("E:\cd4.txt")
;data=read("E:\vcmfdaband\cd4_1.txt")
nt=rows(data)
x=data[,4]
t1=data[,5]
t2=data[,7]
t3=data[,6]
y=data[,8]./100
r=rank(x)./nt

xmin=min(x)
xmax=max(x)
rmin=min(r)
rmax=max(r)
lengthr=rmax-rmin
int=#(0:400)
rgrid=rmin+int.*lengthr/400
fp=(xmax-xmin)*rgrid
nf=rows(fp)


proc(dd)=ob(data)
id=data[,3]
nt=rows(data)
dd=matrix(nt,1)
i=1
while(i<=(nt-1))
dd[i]=(id[i+1]!=id[i])
i=i+1
endo
endp

dd=ob(data)
t1n=dd.*(t1+1)
t1n=paf(t1n,abs(t1n)>0)
t1n=t1n-1

t2n=dd.*t2
t2n=paf(t2n,abs(t2n)>0)
t2m=t2n-mean(t2n)

t3n=dd.*t3
t3n=paf(t3n,abs(t3n)>0)
t3m=t3n-mean(t3n)

;tdata=t1n~t2n~t3n
n=rows(t1n)
tdata=matrix(n,1)~t1n~t2m~t3m

tt2=t2-mean(t2n)
tt3=t3-mean(t3n)
;ttdata=t1~t2~t3
ttdata=matrix(nt,1)~t1~tt2~tt3

proc(ni)=nn(data,n)
dd=ob(data)
nt=rows(data)
num=#(1:nt)
dd1=dd.*num
dd1=paf(dd1,abs(dd1)>0)
ni=matrix(n,1)
i=1
while(i<=(n-1))
ni[1]=dd1[1]
ni[i+1]=dd1[i+1]-dd1[i]
i=i+1
endo
endp

ni=nn(data,n)

proc(mlhat,mlfphat,Ns,hs,x)=ml(data,ttdata,n,nt,d,rgrid)
x=data[,4]
r=rank(x)./nt
y=data[,8]./100
a=min(r)
b=max(r)
length=b-a
Ns=floor(0.3*nt^(1/3)*log(n))
Ns
hs=length/(Ns+1)
knotsl=a+hs*(0:Ns)'
knotsu=a+hs*(1:Ns+1)'
;xknots=(x.>=knotsl[,1:Ns]).*(x.<knotsu[,1:Ns])~(x.>=knotsl[,Ns+1]).*(x.<=knotsu[,Ns+1])
xknots=(r.>=knotsl[,1:Ns]).*(r.<knotsu[,1:Ns])~(r.>=knotsl[,Ns+1]).*(r.<=knotsu[,Ns+1])
rknots=(rgrid.>=knotsl[,1:Ns]).*(rgrid.<knotsu[,1:Ns])~(rgrid.>=knotsl[,Ns+1]).*(rgrid.<=knotsu[,Ns+1])

bx=kron(matrix(1,d),xknots)
T=kron(ttdata,matrix(1,Ns+1))
D=bx.*T
lamda=(ginv(D'*D))*D'*y
lamdal=reshape(lamda,#(Ns+1,d))
mlhat=xknots*lamdal
mlfphat=rknots*lamdal

endp

proc(Glhat,f,sigmalf,sigmal2f)=Gl(data,ttdata,tdata,n,nt,ni,hs,d,rgrid,nf)
x=data[,4]
r=rank(x)./nt
y=data[,8]./100
a=min(r)
b=max(r)
length=b-a
N=floor(0.4*n^(1/3))
N
h=length/(N+1)
knotsl=a+h*(0:N)'
knotsu=a+h*(1:N+1)'
;xknots=(x.>=knotsl[,1:N]).*(x.<knotsu[,1:N])~(x.>=knotsl[,N+1]).*(x.<=knotsu[,N+1])
xknots=(r.>=knotsl[,1:N]).*(r.<knotsu[,1:N])~(r.>=knotsl[,N+1]).*(r.<=knotsu[,N+1])
rknots=(rgrid.>=knotsl[,1:N]).*(rgrid.<knotsu[,1:N])~(rgrid.>=knotsl[,N+1]).*(rgrid.<=knotsu[,N+1])

;f=(1/(nt*h))*xknots*xknots'*matrix(nt,1)
ff=(1/(nt*h))*rknots*xknots'*matrix(nt,1)


bx=kron(matrix(1,d),xknots)
T=kron(ttdata,matrix(1,N+1))
D=bx.*T
lamda=(ginv(D'*D))*D'*y
resi=(y-D*lamda).^2
bbx=kron(matrix(1,d+1),xknots)
TT=kron(ttdata.^2~matrix(nt,1),matrix(1,N+1))
DD=bbx.*TT
rho=ginv(DD'*DD)*DD'*resi
rhol=reshape(rho,#(N+1,d+1))
Glrho=rhol[,1:d]
;Glhat=xknots*Glrho
;rehat=xknots*rhol[,d+1]
Glfhat=rknots*Glrho
refhat=rknots*rhol[,d+1]


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

{mlhat,mlfphat,Ns,hs,x}=ml(data,ttdata,n,nt,d,rgrid)
{Glhat,f,sigmalf,sigmal2f}=Gl(data,ttdata,tdata,n,nt,ni,hs,d,rgrid,nf)
Ns


/*confidence band*/
ans=(2*log(Ns+1))^(1/2)
bns=ans-log(2*pi*(ans^2))/(2*ans)
qns=bns-log(-0.5*log(1-alpha))/ans
;lm=qns.*(sigmalf.^(0.5))
lm=qns.*(sigmal2f.^(0.5))
l=mlfphat-lm
u=mlfphat+lm

lmp=qfn(1-alpha/2).*(sigmal2f.^(0.5))
lp=mlfphat-lmp
up=mlfphat+lmp

/*baseline cd4 effect test*/
bfp=matrix(nf,1)~fp
coehat=inv(bfp'*bfp)*bfp'*mlfphat[,1]
coehat
m0=bfp*coehat
supm0=max(abs(mlfphat[,1]-m0)./(sigmal2f[,1].^(0.5)))
alpha0=1-exp(-2*exp(-ans*(supm0-bns)))
alpha0
qns0=bns-log(-0.5*log(1-alpha0))/ans
lm0=qns0.*(sigmal2f[,1].^(0.5))
l0=mlfphat[,1]-lm0
u0=mlfphat[,1]+lm0
lmp0=qfn(1-alpha0/2).*(sigmal2f[,1].^(0.5))
lp0=mlfphat[,1]-lmp0
up0=mlfphat[,1]+lmp0


/*smoking effect test*/
m1=0*matrix(nf,1)
supm1=max(abs(mlfphat[,2]-m1)./(sigmal2f[,2].^(0.5)))
alpha1=1-exp(-2*exp(-ans*(supm1-bns)))
alpha1
qns1=bns-log(-0.5*log(1-alpha1))/ans
lm1=qns1.*(sigmal2f[,2].^(0.5))
l1=mlfphat[,2]-lm1
u1=mlfphat[,2]+lm1
lmp1=qfn(1-alpha1/2).*(sigmal2f[,2].^(0.5))
lp1=mlfphat[,2]-lmp1
up1=mlfphat[,2]+lmp1



/*precd4 effect test*/
xfp=matrix(nf,1)~fp
abhat=inv(xfp'*xfp)*xfp'*mlfphat[,3]
abhat
m2=xfp*abhat
m2mean=mean(mlfphat[,3])*matrix(nf,1)
supm2=max(abs(mlfphat[,3]-m2mean)./(sigmal2f[,3].^(0.5)))
alpha2=1-exp(-2*exp(-ans*(supm2-bns)))
alpha2
qns2=bns-log(-0.5*log(1-alpha2))/ans
lm2=qns2.*(sigmal2f[,3].^(0.5))
l2=mlfphat[,3]-lm2
u2=mlfphat[,3]+lm2
lmp2=qfn(1-alpha2/2).*(sigmal2f[,3].^(0.5))
lp2=mlfphat[,3]-lmp2
up2=mlfphat[,3]+lmp2


/*age effect test*/
m3=0*matrix(nf,1)
supm3=max(abs(mlfphat[,4]-m3)./(sigmal2f[,4].^(0.5)))
alpha3=1-exp(-2*exp(-ans*(supm3-bns)))
qns3=bns-log(-0.5*log(1-alpha3))/ans
lm3=qns3.*(sigmal2f[,4].^(0.5))
l3=mlfphat[,4]-lm3
u3=mlfphat[,4]+lm3
alpha3
lmp3=qfn(1-alpha3/2).*(sigmal2f[,4].^(0.5))
lp3=mlfphat[,4]-lmp3
up3=mlfphat[,4]+lmp3




ss=0*matrix(nf,1)

dd1=createdisplay(1,1)
dd2=createdisplay(1,1)
dd3=createdisplay(1,1)
dd4=createdisplay(1,1)
dd5=createdisplay(1,1)
dd6=createdisplay(1,1)
dd7=createdisplay(1,1)
dd8=createdisplay(1,1)
ss0=setmask(sort(fp~ss),"line", "green",  "size",3)

lb0=setmask(sort(fp~l[,1]), "line", "red", "size", 3)
ub0=setmask(sort(fp~u[,1]), "line", "red" , "size", 3)
lpw0=setmask(sort(fp~lp[,1]), "line","black", "dotted", "size", 2)
upw0=setmask(sort(fp~up[,1]), "line","black", "dotted", "size", 2)
m0estp=setmask(sort(fp~mlfphat[,1]),"line", "magenta", "size", 3)
lb00=setmask(sort(fp~l0), "line", "red", "size", 3)
ub00=setmask(sort(fp~u0), "line", "red" , "size", 3)
m0line=setmask(sort(fp~m0),"line", "blue", "size", 3)
lpw00=setmask(sort(fp~lp0), "line","black", "dotted", "size", 2)
upw00=setmask(sort(fp~up0), "line","black", "dotted", "size", 2)


lb1=setmask(sort(fp~l[,2]), "line", "red", "size", 3)
ub1=setmask(sort(fp~u[,2]), "line", "red" , "size", 3)
lpw1=setmask(sort(fp~lp[,2]), "line","black", "dotted", "size", 2)
upw1=setmask(sort(fp~up[,2]), "line","black", "dotted", "size", 2)
m1estp=setmask(sort(fp~mlfphat[,2]),"line", "magenta", "size", 3)
lb11=setmask(sort(fp~l1), "line", "red", "size", 3)
ub11=setmask(sort(fp~u1), "line", "red" , "size", 3)
m1line=setmask(sort(fp~m1),"line", "blue", "size", 3)
lpw11=setmask(sort(fp~lp1), "line","black", "dotted", "size", 2)
upw11=setmask(sort(fp~up1), "line","black", "dotted", "size", 2)




lb2=setmask(sort(fp~l[,3]), "line", "red", "size", 3)
ub2=setmask(sort(fp~u[,3]), "line", "red" , "size", 3)
lpw2=setmask(sort(fp~lp[,3]), "line","black", "dotted", "size", 2)
upw2=setmask(sort(fp~up[,3]), "line","black", "dotted", "size", 2)
m2estp=setmask(sort(fp~mlfphat[,3]),"line", "magenta", "size", 3)
m2line=setmask(sort(fp~m2mean),"line", "blue", "size", 3)
lb22=setmask(sort(fp~l2), "line", "red", "size", 3)
ub22=setmask(sort(fp~u2), "line", "red" , "size", 3)
lpw22=setmask(sort(fp~lp2), "line","black", "dotted", "size", 2)
upw22=setmask(sort(fp~up2), "line","black", "dotted", "size", 2)



lb3=setmask(sort(fp~l[,4]), "line", "red", "size", 3)
ub3=setmask(sort(fp~u[,4]), "line", "red" , "size", 3)
lpw3=setmask(sort(fp~lp[,4]), "line","black", "dotted", "size", 2)
upw3=setmask(sort(fp~up[,4]), "line","black", "dotted", "size", 2)
m3estp=setmask(sort(fp~mlfphat[,4]),"line", "magenta", "size", 3)
m3line=setmask(sort(fp~m3),"line", "blue", "size", 3)
lb33=setmask(sort(fp~l3), "line", "red", "size", 3)
ub33=setmask(sort(fp~u3), "line", "red" , "size", 3)
lpw33=setmask(sort(fp~lp3), "line","black", "dotted", "size", 2)
upw33=setmask(sort(fp~up3), "line","black", "dotted", "size", 2)




show(dd1,1,1,lb0,ub0,lpw0,upw0,m0estp,ss0)
show(dd2,1,1,lb00,ub00,lpw00,upw00,m0estp,m0line)
setgopt(dd1,1,1,"title","(a) Baseline CD4","ylim",-0.1|0.5,"xlabel","Year","ylabel","Baseline CD4","ymajor",0.1)
setgopt(dd2,1,1,"title","(b) Baseline CD4","ylim",-0.1|0.5 ,"yvalue",0|1,"xlabel","Year","ylabel","Baseline CD4","ymajor",0.1)


show(dd3,1,1,lb1,ub1,lpw1,upw1,m1estp,ss0)
show(dd4,1,1,lb11,ub11,lpw11,upw11,m1estp,m1line)
setgopt(dd3,1,1,"title","(a) Smoking Effect","ylim",-0.4|0.4, "yvalue",0|1, "xlabel","Year","ylabel","Coeff. of Smoke","ymajor",0.1)
setgopt(dd4,1,1,"title","(b) Smoking Effect","ylim",-0.4|0.4, "yvalue",0|1, "xlabel","Year","ylabel","Coeff. of Smoke","ymajor",0.1)

show(dd5,1,1,lb2,ub2,lpw2,upw2,m2estp,ss0)
show(dd6,1,1,lb22,ub22,lpw22,upw22,m2estp,m2line)
setgopt(dd5,1,1,"title","(a) PreCD4 Effect","ylim",-0.03|0.03, "yvalue",0|1,"xlabel","Year","ylabel","Coeff. of PreCD4","ymajor",0.01)
setgopt(dd6,1,1,"title","(b) PreCD4 Effect","ylim",-0.03|0.03, "yvalue",0|1,"xlabel","Year","ylabel","Coeff. of PreCD4","ymajor",0.01)


show(dd7,1,1,lb3,ub3,lpw3,upw3,m3estp,ss0)
show(dd8,1,1,lb33,ub33,lpw33,upw33,m3estp,m3line)
setgopt(dd7,1,1,"title","(a) Age Effect","ylim",-0.03|0.03, "yvalue",0|1,"xlabel","Year","ylabel","Coeff. of Age","ymajor",0.01)
setgopt(dd8,1,1,"title","(b) Age Effect","ylim",-0.03|0.03, "yvalue",0|1,"xlabel","Year","ylabel","Coeff. of Age","ymajor",0.01)


