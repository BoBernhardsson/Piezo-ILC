function [newvpiezo] = ilcupdate(par,l,deltahat,delta0hat,skip)
% function [newpiezopulse,newvoffset] = ilcupdate(par,V,l,deltahat,delta0hat,skip)

global vpiezo 
if nargin < 5
    skip = par.Lskip
end

e = par.deltaref - deltahat;
b = par.deltaref - delta0hat;
changev = filter(l,1,[e; zeros(skip,1)]);                
changev = changev(skip+1:skip+par.nt);
changev(1:round((par.tstart-0.2e-3)/par.dt)) = 0;    % only change Vpiezo  near pulse
changev(round((par.tend+par.tdecay/2)/par.dt):end) = 0;     % changed here   
   
[B3,A3] = butter(2,par.bw_ilc*par.dt*2);    % Q-filter in ILC
newvpiezo = filtfilt(B3,A3,vpiezo  + par.ilckappa1 * changev + par.ilckappa2 * b); 
end

