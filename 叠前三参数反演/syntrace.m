function [S,R] = syntrace(AI,tlen,f0,dt,Ks)

R=zeros(tlen,1);           % 反射系数
for i=1:tlen-1
    R(i)=(AI(i+1)-AI(i))/(AI(i+1)+AI(i));
end

source(1:2*Ks+1,1)=(1-2.*((pi*f0.*(-Ks:Ks).*dt).^2)).*exp(-((pi*f0.*(-Ks:Ks).*dt).^2));   %  子波


S=conv(source,R);   % 合成地震道
S=S(Ks+1:end-Ks);

if nargout>1 varargout{2}=R;end





