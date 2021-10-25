function [deltahat,delta0hat] = estimatedetuning(par,V,C1est,C2est,Vreg,B,A)
% function [deltahat,delta0hat] = estimatedetuning(V,u,Ib,C1,C2,Vreg)
% estimates the detuning offset deltahat(t), size(deltahat)=[par.nt,1]
% and the constant detuning offset delta0
global u Ib
if nargin <3
    C1est = par.cavity.C1;
    C2est = par.cavity.C2;
end
if  nargin<5
    Vreg = 0.05; % regularization parameter, larger value gives more robust estimates
end
if nargin<6
    bw = 2e3;    % bandwidth in Hz
    [B,A] = butter(1,bw*par.dt*2);
end

% filter V, dV/dt, u and Ib
Vf = filter(B,A,V);
dVf = filter(B,A,[0; V(3:end)-V(1:end-2); 0])/2/par.dt;
uf = filter(B,A,u);
Ibf = filter(B,A,Ib);

deltahatn = 1/2/pi*imag((conj(Vf).*dVf - C1est*conj(Vf).*uf - C2est*conj(Vf).*Ibf)) ./ (Vreg^2 + conj(Vf).*Vf);

% Only use deltahat when V is sufficiently large, else put deltahat=deltaref
% Filter to smooth out transients
largeV = (abs(V)>2*Vreg);
deltahat = filter(B,A,largeV.*deltahatn + (1-largeV)*par.deltaref);

indinitial = find(par.tstart < par.tvec & par.tvec < par.tstart+par.initialmeasuretime);
delta0hat = mean(deltahat(indinitial));
end

