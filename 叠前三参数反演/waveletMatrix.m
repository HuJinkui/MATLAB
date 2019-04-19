function [wvletMatrix]=waveletMatrix(f0,Ks,dt,rlen)

% rlen---反射系数长度

source(1:2*Ks+1,1)=(1-2.*((pi*f0.*(-Ks:Ks).*dt).^2)).*exp(-((pi*f0.*(-Ks:Ks).*dt).^2));   %  ricker子波
wvletMatrix=zeros(rlen+2*Ks,rlen);
for i=1:rlen
    wvletMatrix(i:i+2*Ks,i)=source;
end