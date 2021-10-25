function dx = process(t,x,offset)
global tvec u vpiezo dt
% function dx = process(t,x,offset)
% t: time
% x: state vector for cavity + detuning dynamics
% offset: Piezo voltage vpiezo = 0 gives detuning = offset in stationarity

global whalf KLFD Kpiezo Omega1 Q1

% Interpolate forward power u and piezo voltage vpiezo at time t

ux = lininterp1(tvec,u,t);     % Forward power, a.u.
vpiezoint = lininterp1(tvec,vpiezo,t);             % Piezo voltage, a.u.

% Dynamics
n = length(x);      % state dimension, e.g. n = 3
dx = zeros(n,1);    % a column vector
V = x(1);           % cavity probe signal, a.u.
delta = x(2);       % detuning in Hertz
deltadot = x(3);    % derivative of x(2)

% Cavity dynamics
dx(1) = -(whalf-1i*2*pi*delta)*V + whalf*ux;
% Mechanical dynamics, here 2nd order dynamics
dx(2) = deltadot;
dx(3) = - Omega1/Q1*deltadot - Omega1^2*(delta-offset) + Omega1^2*KLFD*abs(V^2) + Omega1^2*Kpiezo*vpiezoint;
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

