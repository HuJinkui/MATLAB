function s=seiswavelet(fm,t0,dt,halflen)
% This function for seismic source (zero phase rick wavelet)
% IN   fm:        major frequency(Hz)
%      dt:        time sample interval
%      halflen:   half length of seismic source
% OUT  s:         a vector of seismic source      

% s=(1-2.*(pi*fm.*((-halflen:halflen).*dt-t0)).^2).*exp(-((pi*fm.*((-halflen:halflen).*dt-t0).^2));         
s=(1-2.*((pi*fm.*((-halflen:halflen).*dt-t0)).^2)).*exp(-((pi*fm.*((-halflen:halflen).*dt-t0)).^2));    
s=s';