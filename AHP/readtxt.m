clc;
clear;
txt_rms=importdata('rms45.txt');
tmp=txt_rms;
inline_min=min(tmp(:,1));
inline_max=max(tmp(:,1));
xline_min=min(tmp(:,2));
xline_max=max(tmp(:,2));

dem_x=inline_max-inline_min+1;
dem_y=xline_max-xline_min+1;

B=zeros(dem_x,dem_y);
c1=zeros(dem_x,1);
c2=zeros(dem_y,1);

for i=inline_min:inline_max
    j=i-inline_min+1;
    c1(i,1)=tmp(j,1)-inline_min+1;
   
end

for i=xline_min:xline_max
    j=i-xline_min+1;
    c2(i,1)=tmp(j,1)-xline_min+1;
   
end

% txt_rms89=importdata('rms89.txt');
% 
% tmp45=txt_rms45;
% tmp89=txt_rms89;
% inline_min45=min(tmp45(:,1));
% inline_max45=max(tmp45(:,1));
% inline_min89=min(tmp89(:,1));
% inline_max89=max(tmp89(:,1));
% 
% xline_min45=min(tmp45(:,2));
% xline_max45=max(tmp45(:,2));
% xline_min89=min(tmp89(:,2));
% xline_max89=max(tmp89(:,2));
% 
% 
% 
% B=zeros(inline_max45,xline_max45);
% c1=zeros(inline_max45);
% c2=zeros(xline_max45);
% 
% for i=inline_min45:inline_max45
%     j=i-inline_min45+1;
%     c1(i)=tmp45(j,1);
%     
% end
