
function [yy] = funFID(ampDet,freqDet,phaseDet, t, N)

offset = centerIdx(N)-1;

%modulation = exp(-t./N./tau) ; 
%harmonic = exp(-1i.*   (-2*pi*(freqDet-centerIdx(N))./N.*(t-offset) + phaseDet)   );
%harmonic = exp(-1i*   (-2*pi*(freqDet-centerIdx(N))./N.*(t-offset) + phaseDet)   );
harmonic =  exp(-1i*   (2*pi*freqDet/N*(t-offset) + phaseDet)   );
%harmonic = real(harmonic);

yy= ampDet.*harmonic;%.*modulation;
