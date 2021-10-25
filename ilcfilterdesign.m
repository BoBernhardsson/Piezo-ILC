function [l] = ilcfilterdesign(H,par)
% function [l] = ilcfilterdesign(H,par)
% Calculates timedomain ilcfilter l 
% length of filter = par.Lfilterlength
% Filter is noncausal with delay = par.Lskip

sigma = 0.1*abs(H(1));  % lower value gives a more high frequency inverse
Hinv = conj(H)./(sigma^2+abs(H.^2));
l = ifft(Hinv);
l = circshift(l,par.Lskip);
l = l(1:par.Lfilterlength);
end

