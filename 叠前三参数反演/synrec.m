function syndat=synrec(wvletMatrix,Ks,r)

[m,n]=size(r);
syndat=zeros(m,n);

for i=1:n
    temp=wvletMatrix*r(:,i);
    syndat(:,i)=temp(Ks+1:end-Ks);
end