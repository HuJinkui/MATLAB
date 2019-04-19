clc;clear;

% syms vp1 vp2 vs1 vs2 p1 p2 thet_ rthet
% rthet=0.5*( 1-4*( (vs1+vs2)*0.5 )^2*sin(thet_)^2/(0.5*(vp1+vp2))^2)*(p2-p1)/(0.5*(p1+p2))...
%        +0.5*(1+tan(thet_)^2)*(vp2-vp1)/(0.5*(vp1+vp2))...
%        -4*((vs1+vs2)*0.5)^2*sin(thet_)^2/(0.5*(vp1+vp2))^2*(vs2-vs1)/(0.5*(vs1+vs2));
% drdp2=diff(rthet,p2) 
% 
% vp1=20;vp2=20;
% vs1=30;vs2=30;
% p1=10;p2=10;
% thet_=1;

dRdvp1=((tan(thet_)^2/2 + 1/2)*(vp1 - vp2))/(2*(vp1/2 + vp2/2)^2)... 
        - (tan(thet_)^2/2 + 1/2)/(vp1/2 + vp2/2) ...
        - (4*sin(thet_)^2*(vs1/2 + vs2/2)*(vs1 - vs2))/(vp1/2 + vp2/2)^3 ...
        - (2*sin(thet_)^2*(vs1/2 + vs2/2)^2*(p1 - p2))/((vp1/2 + vp2/2)^3*(p1/2 + p2/2));
dRdvp2=(tan(thet_)^2/2 + 1/2)/(vp1/2 + vp2/2)...
       + ((tan(thet_)^2/2 + 1/2)*(vp1 - vp2))/(2*(vp1/2 + vp2/2)^2) ...
       - (4*sin(thet_)^2*(vs1/2 + vs2/2)*(vs1 - vs2))/(vp1/2 + vp2/2)^3 ...
       - (2*sin(thet_)^2*(vs1/2 + vs2/2)^2*(p1 - p2))/((vp1/2 + vp2/2)^3*(p1/2 + p2/2));
dRdvs1=(4*sin(thet_)^2*(vs1/2 + vs2/2))/(vp1/2 + vp2/2)^2 ...
       + (2*sin(thet_)^2*(vs1 - vs2))/(vp1/2 + vp2/2)^2 ...
       + (2*sin(thet_)^2*(vs1/2 + vs2/2)*(p1 - p2))/((vp1/2 + vp2/2)^2*(p1/2 + p2/2)) ;
dRdvs2=(2*sin(thet_)^2*(vs1 - vs2))/(vp1/2 + vp2/2)^2 ...
       - (4*sin(thet_)^2*(vs1/2 + vs2/2))/(vp1/2 + vp2/2)^2 ...
       + (2*sin(thet_)^2*(vs1/2 + vs2/2)*(p1 - p2))/((vp1/2 + vp2/2)^2*(p1/2 + p2/2));
dRdp1=((2*sin(thet_)^2*(vs1/2 + vs2/2)^2)/(vp1/2 + vp2/2)^2 - 1/2)/(p1/2 + p2/2) ...
       - ((p1 - p2)*((2*sin(thet_)^2*(vs1/2 + vs2/2)^2)/(vp1/2 + vp2/2)^2 - 1/2))/(2*(p1/2 + p2/2)^2);
dRdp2=- ((2*sin(thet_)^2*(vs1/2 + vs2/2)^2)/(vp1/2 + vp2/2)^2 - 1/2)/(p1/2 + p2/2) ...
      - ((p1 - p2)*((2*sin(thet_)^2*(vs1/2 + vs2/2)^2)/(vp1/2 + vp2/2)^2 - 1/2))/(2*(p1/2 + p2/2)^2);









% vp1=rand(1,1);vp2=rand(1,1);
% vs1=rand(1,1);vs2=rand(1,1);
% p1=rand(1,1);p2=rand(1,1);
% thet_=rand(1,1)+0.1;
% 
% vs_=(vs2+vs1)*0.5;dvs=vs2-vs1;
% vp_=(vp2+vp1)*0.5;dvp=vp2-vp1;
% p_=(p2+p1)*0.5;dp=p2-p1;
% 
% num=-2*vs_^2*sin(thet_)^2*(-1*vp_^-3)*dp/p_...
%     +0.5*(1+tan(thet_)^2)*(-vp_-0.5*dvp)/vp_^2 ...
%     -4*vs_*sin(thet_)^2*dvs*(-vp_^-3)
% drdvp1=subs(drdvp1)