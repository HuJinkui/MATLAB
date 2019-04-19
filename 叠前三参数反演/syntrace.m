function [S,R] = syntrace(AI,tlen,f0,dt,Ks)

R=zeros(tlen,1);           % ����ϵ��
for i=1:tlen-1
    R(i)=(AI(i+1)-AI(i))/(AI(i+1)+AI(i));
end

source(1:2*Ks+1,1)=(1-2.*((pi*f0.*(-Ks:Ks).*dt).^2)).*exp(-((pi*f0.*(-Ks:Ks).*dt).^2));   %  �Ӳ�


S=conv(source,R);   % �ϳɵ����
S=S(Ks+1:end-Ks);

if nargout>1 varargout{2}=R;end





