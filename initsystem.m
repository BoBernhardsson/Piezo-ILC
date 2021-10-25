function par = initsystem(simulationfrequency)
% Initialize parameters for simulation
%
% Cavity Parameters
% The KLFD, Kpiezo, Omega and Q vectors describe the Piezo mechanical modes

par.cavity.whalf = 2*pi*500;           % Cavity half bandwidth in rad/sec
par.cavity.cavphase = 0;               % Cavity phase in radians
par.cavity.C1 = par.cavity.whalf*exp(i*par.cavity.cavphase); % A.u.
par.cavity.C2 = par.cavity.C1;         % Assumed equal to C1 for simplicity
par.cavity.KLFDvector = [400];         % Lorentz force detuning constants per mode
par.cavity.Kpiezovector = [1];         % Piezo force constants Hz/Volt. A.u.
par.cavity.Omegavector = 2*pi*[350];   % Mechanical modes in rad/sec
par.cavity.Qvector = [3];              % Q-factors of modes
par.cavity.offsetHz = (-100);          % Piezo offset in Hz for Vpiezo=0
par.cavity.nrmodes = length(par.cavity.KLFDvector); % Nr of mechanical modes
if length(par.cavity.Kpiezovector) ~= par.cavity.nrmodes 
    error('wrong number of modes in Kpiezovector')
end
if length(par.cavity.Omegavector) ~= par.cavity.nrmodes 
    error('wrong number of modes in Omegavector')
end
if length(par.cavity.Qvector) ~= par.cavity.nrmodes 
    error('wrong number of modes in Qvector')
end

%
% Simulation Parameters
% We assume to fill the cavity to the (normalized) level V = 1
par.dt = 1/simulationfrequency;         % Piezo sample time, e.g. 20us
par.tmax = 20e-3 - par.dt;              % Piezo pulse length, here made shorter than 1/14 sec
par.tvec = 0 : par.dt : par.tmax;
par.nt = length(par.tvec);
par.fmax = 1/par.dt;                    % Piezo sample rate, e.g 50kHz
par.df = 1/(par.nt*par.dt);
par.fvec = 0 : par.df : par.fmax-par.df;
par.tstart = 4e-3;                      % Start of filling, e.g. 5ms in
par.tfill = 0.3e-3;                     % Fill time, e.g 0.3ms
par.tpulse = 2.86e-3;                   % Beam pulse length, e.g 2.86ms
par.tend = par.tstart + par.tfill + par.tpulse;  % time for end of pulse 
%par.ufill = 1.03*exp(-i*par.cavity.cavphase);
par.ufill = 1.638*exp(-i*par.cavity.cavphase);   % LLRF signal during filling of cavity, assumed constant
par.beam = 0*(-0.3-0.2*i);              % Complex valued beam,  assumed constant
par.uduringbeam = exp(-i*par.cavity.cavphase) - par.beam; % LLRF signal during beam, perfect control
par.options = odeset('RelTol',1e-6,'AbsTol',1e-5*ones(1,2*par.cavity.nrmodes+1)); % simu pars
par.noiselevel = 0*0.001;                 % Noiselevel on u and V, 0.001 = -60dB
%
% Algorithm Parameters
par.initialmeasuretime = 0.2e-3;        % Initial part of pulse used for offset estimation
par.tdecay = 1e-3;                      % Pulse decay used for estimation
par.Lskip = 150;                        % Noncausal part of L-filter, nr of samples
par.Lfilterlength = 301;                % Length of L-filter, nr of samples
par.bw_ilc = 3e3; 
par.ilckappa1 = 1;                      % ilc gain for detuning error e(t)
par.ilckappa2 = 0;                      % ilc gain for offset bias b
par.deltaref = 0;                       % setpoint for detuning [Hz]
