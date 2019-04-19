function G=JacobiMatrix(ivp,ivs,iden,ang,wvletMatrix,Ks)

m=length(ivp);
G11=zeros(m,m);G12=zeros(m,m);G13=zeros(m,m);
G21=zeros(m,m);G22=zeros(m,m);G23=zeros(m,m);
G31=zeros(m,m);G32=zeros(m,m);G33=zeros(m,m);
          

for i=1:m-1
    vp1=ivp(i);vp2=ivp(i+1);
    vs1=ivs(i);vs2=ivs(i+1);
    p1=iden(i);p2=iden(i+1);
    thet=ang(1);
    thet2=asin(vp2/vp1*sin(thet));
    thet_=0.5*(thet+thet2);
    [dRdvp1,dRdvp2,dRdvs1,dRdvs2,dRdp1,dRdp2]=elementdiff(vp1,vp2,vs1,vs2,p1,p2,thet_);
    G11(i,i:i+1)=[dRdvp1,dRdvp2];
    G12(i,i:i+1)=[dRdvs1,dRdvs2];
    G13(i,i:i+1)=[dRdp1,dRdp2];
end
[dRdvp1,dRdvp2,dRdvs1,dRdvs2,dRdp1,dRdp2]=elementdiff(ivp(end),ivp(end),ivs(end),ivs(end),iden(end),iden(end),ang(1));
G11(m,m)=dRdvp1;G12(m,m)=dRdvs1;G13(m,m)=dRdp1;

for i=1:m-1
    vp1=ivp(i);vp2=ivp(i+1);
    vs1=ivs(i);vs2=ivs(i+1);
    p1=iden(i);p2=iden(i+1);
    thet=ang(2);
    thet2=asin(vp2/vp1*sin(thet));
    thet_=0.5*(thet+thet2);
    [dRdvp1,dRdvp2,dRdvs1,dRdvs2,dRdp1,dRdp2]=elementdiff(vp1,vp2,vs1,vs2,p1,p2,thet_);
    G21(i,i:i+1)=[dRdvp1,dRdvp2];
    G22(i,i:i+1)=[dRdvs1,dRdvs2];
    G23(i,i:i+1)=[dRdp1,dRdp2];
end
[dRdvp1,dRdvp2,dRdvs1,dRdvs2,dRdp1,dRdp2]=elementdiff(ivp(end),ivp(end),ivs(end),ivs(end),iden(end),iden(end),ang(2));
G21(m,m)=dRdvp1;G22(m,m)=dRdvs1;G23(m,m)=dRdp1;



for i=1:m-1
    vp1=ivp(i);vp2=ivp(i+1);
    vs1=ivs(i);vs2=ivs(i+1);
    p1=iden(i);p2=iden(i+1);
    thet=ang(3);
    thet2=asin(vp2/vp1*sin(thet));
    thet_=0.5*(thet+thet2);
    [dRdvp1,dRdvp2,dRdvs1,dRdvs2,dRdp1,dRdp2]=elementdiff(vp1,vp2,vs1,vs2,p1,p2,thet_);
    G31(i,i:i+1)=[dRdvp1,dRdvp2];
    G32(i,i:i+1)=[dRdvs1,dRdvs2];
    G33(i,i:i+1)=[dRdp1,dRdp2];
end
[dRdvp1,dRdvp2,dRdvs1,dRdvs2,dRdp1,dRdp2]=elementdiff(ivp(end),ivp(end),ivs(end),ivs(end),iden(end),iden(end),ang(3));
G31(m,m)=dRdvp1;G32(m,m)=dRdvs1;G33(m,m)=dRdp1;

G11=synrec(wvletMatrix,Ks,G11);
G12=synrec(wvletMatrix,Ks,G12);
G13=synrec(wvletMatrix,Ks,G13);

G21=synrec(wvletMatrix,Ks,G21);
G22=synrec(wvletMatrix,Ks,G22);
G23=synrec(wvletMatrix,Ks,G23);

G31=synrec(wvletMatrix,Ks,G31);
G32=synrec(wvletMatrix,Ks,G32);
G33=synrec(wvletMatrix,Ks,G33);

G=[G11,G12,G13;
   G21,G22,G23;
   G31,G32,G33;];




