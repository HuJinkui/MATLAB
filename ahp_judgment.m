% A=[1,1.2,1.9;0.83,1,1.5;0.53,0.67,1];
% A=[1,1.2,1.9;0.83,1,1.5;0.53,0.67,1];
% A=[1,1.6;0.63,1];
a11=1.0;a12=0.67;a13=1.25;
a21=1.0/a12;a22=1.0;a23=1.67;
a31=1.0/a13;a32=1.0/a23;a33=1.0;
A=[a11,a12,a13;a21,a22,a23;a31,a32,a33];
m=3;
RI=0.52;
% [w,CR]=mycom(A,m,RI);
[x,lumda]=eig(A); 
r=abs(sum(lumda)); 
n=find(r==max(r)); 
max_lumda_A=lumda(n,n); 
max_x_A=x(:,n);
w=A/sum(A);
CR=(max_lumda_A-m)/(m-1)/RI;