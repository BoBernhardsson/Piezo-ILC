function [x,delta] = cavitysimulator(xinit,par)
% Cavitysimulator, noise free. Add noise outside
% xinit: initial statevector [V(0);delta1(0); delta1dot(0);...]
% par: Initialize parameters with initsystem.m
global u vpiezo Ib
[t,x] = ode45(@(t,x) cavityprocess(t,x,par) ,par.tvec, xinit, par.options);
delta = par.cavity.offsetHz + sum(x(:,2:2:end),2);
end


function dx = cavityprocess(t,x,par)
global u vpiezo Ib
% function dx = process(t,x,par)
% t: time (scalar)
% x: state vector for cavity + detuning dynamics
% x(1) = V 
% x(2k)= delta_k(t)
% x(2k+1) = deltadot_k(t)
% par.offset: Piezo voltage vpiezo = 0 gives detuning = offset in Hz in stationarity

% Dynamics
dx = zeros(par.cavity.nrmodes+1,1);    % a column vector
% Interpolate forward power u and piezo voltage vpiezo at time t
ut = lininterp1(par.tvec,u,t);
% Forward power, a.u.
vpiezot = lininterp1(par.tvec,vpiezo,t);  % Piezo voltage, a.u.
Ibt = lininterp1(par.tvec,Ib,t);          % Beam, a.u.
% 
delta = par.cavity.offsetHz + sum(x(2:2:end));   % Offset in Hertz
dx(1) = -(par.cavity.whalf-1i*2*pi*delta)*x(1) + par.cavity.C1*ut + par.cavity.C2*Ibt;
% Mechanical dynamics, here 2nd order dynamics with nrmodes
for k = 1:par.cavity.nrmodes    
    Omegak = par.cavity.Omegavector(k);
    Qk = par.cavity.Qvector(k);
    KLFDk = par.cavity.KLFDvector(k);
    Kpiezok = par.cavity.Kpiezovector(k);
    dx(2*k) = x(2*k+1);
    dx(2*k+1) = - Omegak/Qk*x(2*k+1) - Omegak^2*x(2*k) + Omegak^2*KLFDk*abs(x(1)^2) + Omegak^2*Kpiezok*vpiezot;
end
end

    
function v = lininterp1(X, V, x)
% linear interpolation, given set of X and V values, and an x query
% assumes X values are in strictly increasing order
%
% Differences from matlab built-in :
%       much, much faster
%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
%       extends values off the ends instead of giving NaN
%       
if length(X) ~= length(V), error('X and V sizes do not match'); end
pindex = find((x >= X), 1, 'last');
index = find((x <= X), 1, 'first');
if isempty(pindex)
    warning('interpolating before beginning');
    pindex = index;
    slope = 0;
elseif isempty(index)
    warning('interpolating after end');
    index = pindex;
    slope = 0;
elseif pindex == index
    slope = 0;
else
    Xp = X(pindex);
    slope = (x - Xp) / (X(index) - Xp);
end
v = V(pindex) * (1 - slope) + V(index) * slope;
end    
    
    