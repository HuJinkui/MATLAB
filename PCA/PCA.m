clear;
a=load('result45.dat');
X=a(:,3:11);                  %3-11��Ϊ�������ݵ�������
x=zscore(X);                       %��׼��
[coef,score,eig,t]=pca(x);   %����princomp�������
t ;                                %ÿһ���������������µ�ԭ��ľ���
s=0;
i=1;
while s/sum(eig)<0.974        %����ۼƹ����ʴ���xx%��������
    s=s+eig(i);
    i=i+1;
end                              
NEW=x*coef(:,1:i-1)  ;            %����µ�����
[R,P] = corrcoef(NEW); % ���� R = corrcoef(x,y)  ����Է���
pca(:,1:2)=a(:,1:2);        %1��2��Ϊ������������cdp��
pca(:,3:9)=NEW;             %3-9��Ϊ����������ɷ�
% xlswrite('pca89.xlsx',pca);   %����µ����ɷ����ݵ�Excel
figure
pareto(eig/sum(eig));          %���������ֱ��ͼ
figure(2)
plot(eig,'r+');
hold on
plot(eig,'b-');