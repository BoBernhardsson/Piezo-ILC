function [h,H,Voffsethat] = piezoident(par)
%[h,H,Voffsethat] = function piezoident(simulationfrequency)
% Identify Piezo transfer function H(iw)
% Uses experiment with impulse response vpiezo(t)
% Alternative could be to use step or chirp experiments 
% For real operation, delta1 and delta2 need to be estimated using detuningestim.m

global u vpiezo Ib

% First  experiment
smallamp = 0.1;      % small cavity field to be able to measure detuning
u = smallamp*ones(par.nt,1);
vpiezo = 0*ones(par.nt,1);
Ib = 0*ones(par.nt,1);
xinit = zeros(2*par.cavity.nrmodes+1,1);
[x,delta1] = cavitysimulator(xinit,par);

% Second experiment, impulse response
ind1 = 2;  % dont use ind1=1 since interpolation problem in cavitysimulator
impamp = 100;
vpiezo(ind1) = impamp;
[x,delta2] = cavitysimulator(xinit,par);

% Calculate transfer function
diff = delta2 - delta1;
H = fft(diff)./fft(vpiezo);
h = ifft(H);

% Estimate of Voffset, i.e. Vpiezo value giving zero initial detuning error
Voffsethat = (par.deltaref-mean(delta1))/H(1);   % Assumes error-free detuning estimation
end

% To check the results
% 
% figure(1)
% loglog(par.fvec,abs(H),'-x')
% grid on
% figure(2)
% plot(par.tvec,h,'b',par.tvec,diff/sum(vpiezo),'r')
