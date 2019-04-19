clc;clear;clf
% load model.txt
% load initial.txt
load w01_1.txt
model=w01_1;


% *************************************************************************
%   理论模型、初始模型
% *************************************************************************

vp=model(:,2);vs=model(:,3);den=model(:,4);
% ivp=initial(:,2);ivs=initial(:,3);iden=initial(:,4);
ivp=smooth(vp,101);ivs=smooth(vs,101);iden=smooth(den,101);

% [vp,vs,den,ivp,ivs,iden]=well2model(well);  


modellen=length(vp);                      % 时间长度
dep=1:modellen;



% *************************************************************************
%   反射系数、子波矩阵、合成理论记录
% *************************************************************************
% ang=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70].*pi/180;  % 角度（弧度制）
ang=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50].*pi/180;  % 角度（弧度制）
% ang=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30].*pi/180;  % 角度（弧度制）
% ang=[10,20,30].*pi/180;  % 角度（弧度制）
anglen=length(ang);

r=rflecoe(vp,vs,den,ang,'shuey');     % 两种反射系数的计算  'zoeppritz'

Ks=60;f0=40;dt=0.001;   % 子波长度、主频、采样率
wvletMatrix=waveletMatrix(f0,Ks,dt,length(r));  % 子波矩阵

syndat=synrec(wvletMatrix,Ks,r);  % 合成理论记录


% *************************************************************************
%   迭代计算
% *************************************************************************


for i=1:anglen
    ir=rflecoe(ivp,ivs,iden,ang,'aki_richard');
    isyndat=synrec(wvletMatrix,Ks,ir);
    ds=syndat-isyndat;
    ds=ds(:);
    G=JacobiMatrix(ivp,ivs,iden,ang,wvletMatrix,Ks);
    dm=(G'*G+eps/0.001*eye(length(G)))\G'*ds;
    ivp=ivp+dm(1:modellen);
    ivs=ivs+dm(modellen+1:2*modellen);
    iden=iden+dm(2*modellen+1:end);
    i
end

% *************************************************************************
%   成图比较
% *************************************************************************

offset=200;
figure(1)
subplot(1,3,1),plot(vp,dep);
hold on;plot(ivp,dep,'r');hold off;
set(gca,'YDir','reverse','XAxisLocation','top'); 
axis([min(vp)-offset,max(vp)+offset,0,modellen+10])
title('Vp')
subplot(1,3,2),plot(vs,dep);
hold on;plot(ivs,dep,'r');hold off;
set(gca,'YDir','reverse','XAxisLocation','top'); 
axis([min(vs)-offset,max(vs)+offset,0,modellen+10])
title('Vs')
subplot(1,3,3),plot(den,dep);
hold on;plot(iden,dep,'r');hold off;
set(gca,'YDir','reverse','XAxisLocation','top');   
axis([min(den)-offset,max(den)+offset,0,modellen+10])
title('Density')
% figure(2)  
% wigb(syndat)
% figure(3)  
% wigb(isyndat)

