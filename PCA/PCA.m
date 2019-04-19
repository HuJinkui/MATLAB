clear;
a=load('result45.dat');
X=a(:,3:11);                  %3-11列为输入数据的特征列
x=zscore(X);                       %标准化
[coef,score,eig,t]=pca(x);   %利用princomp处理矩阵
t ;                                %每一组数据在新坐标下到原点的距离
s=0;
i=1;
while s/sum(eig)<0.974        %获得累计贡献率大于xx%几组数据
    s=s+eig(i);
    i=i+1;
end                              
NEW=x*coef(:,1:i-1)  ;            %输出新的数据
[R,P] = corrcoef(NEW); % 或者 R = corrcoef(x,y)  相关性分析
pca(:,1:2)=a(:,1:2);        %1、2列为坐标或地震数据cdp号
pca(:,3:9)=NEW;             %3-9列为输出的新主成分
% xlswrite('pca89.xlsx',pca);   %输出新的主成分数据到Excel
figure
pareto(eig/sum(eig));          %输出贡献率直方图
figure(2)
plot(eig,'r+');
hold on
plot(eig,'b-');