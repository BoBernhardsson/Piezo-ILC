% Identifies Cavity parameters 
% dV/dt = -whalf * V + 2*pi*delta*V + C1*u + C2*Ib
% Where whalf is real and C1 and C2 are complex numbers
% Sufficiently exciting signals u and Ib are needed. 
% Constant signals give bad identifiability 

global u vpiezo Ib

% Initialize system and simulation parameters 
simulationfrequency = 1000e3;
par = initsystem(simulationfrequency);

% Initialize LLRF control , piezo control and beam 
u = par.ufill*(par.tstart < par.tvec & par.tvec <= par.tstart + par.tfill) + par.uduringbeam*(par.tstart+par.tfill < par.tvec & par.tvec < par.tend);
u = u.'; 
% Better identifiability with more varying inpout
%u = 0.3*sin(2*pi*100*par.tvec) + i*0.3*sin(2*pi*180*par.tvec) + 0.3*sin(2*pi*250*par.tvec); 
%u = u.'; 

vpiezo = 0*ones(size(par.tvec))';
Ib = par.beam*(par.tstart+par.tfill < par.tvec & par.tvec < par.tend);
Ib = Ib.';

% simulate cavity
xinit = zeros(2*par.cavity.nrmodes+1,1);
[x,delta] = cavitysimulator(xinit,par);
V = x(:,1);


bw = 100e3; % filter bandwidth in Hz
[b,a] = butter(1,bw*2*par.dt);      % Find LP-filter parameters
dVdt = (V(3:end)-V(1:end-2))/(2*par.dt); % Symmetric diff
dVdt = [0;dVdt;0];
dVdtfilt = filter(b,a,dVdt);    % 
% Could possibly filter V, u and Ib also

% Solving least squares problem
if max(abs(Ib))> 1e-8
    % Estimating all parameters by least squares if beam present
    y = real(conj(V).*dVdtfilt);
    A = [-conj(V).*V real(conj(V).*u) -imag(conj(V).*u) real(conj(V).*Ib) -imag(conj(V).*Ib)];
    theta = inv(A'*A)*(A'*y);      % Inverting 5x5 matrix
    whalfhat = theta(1)
    C1hat =  theta(2) + i*theta(3)
    C2hat = theta(4) + i*theta(5)
else
    % Estimating cavity parameters by least squares if beam not present
    y = real(conj(V).*dVdtfilt);
    A = [-conj(V).*V real(conj(V).*u) -imag(conj(V).*u) ];
    theta = inv(A'*A)*(A'*y);      % Inverting 3x3 matrix
    whalfhat = theta(1)
    C1hat =  theta(2) + i*theta(3)
end
