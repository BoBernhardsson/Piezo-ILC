%close all
clear
         
global whalf KLFD Kpiezo  Omega1  Q1 tvec u vpiezo dt

whalf = 2*pi*500;   % Cavity bandwidth 500Hz
KLFD = 400;         % Lorentz force detuning constant, 
Kpiezo = 1;         % Piezo force constant Hz/Volt, arbitrary
Omega1 = 2*pi*350;  % Main mechanical mode 350 Hz
Q1 = 3;             % its damping ratio zeta = 1/6 
offset = -300;         % Initial detuning offset in Hz, typically [-1000,1000];

dt = 2e-6;          % sample interval in simulation
tmax = 20e-3 - dt;  % total simulation time 
tstart = 5e-3;      % start of filling, e.g. 5ms
tfill = 0.3e-3;     % fill time, typically 0.3ms
tpulse = 2.86e-3;   % pulse length, typically 2.86ms
tdecay = 1e-3;      % pulse decay used for estimation
tend = tstart + tfill + tpulse;  % time for end of pulse 
tvec = 0 : dt : tmax;
nt = length(tvec);
fmax = 1/dt;
df = 1/(nt*dt);
f = 0 : df : fmax-df;

% Filling profile, no closed loop control assumed
% Find filling value by trial and error, depends on tfill and whalf
ufill = 1.638;
u = ufill*(tstart < tvec & tvec <= tstart + tfill) + 1*(tstart+tfill < tvec & tvec < tend);
u = u.';
vpiezo = 0*ones(size(tvec))';


%% initial simulation with v = 0, i.e. no Piezo control
options = odeset('RelTol',1e-6,'AbsTol',[1e-5 1e-5 1e-5]);     % setting accuracies
[t,x] = ode45(@(t,y) process2(t,y,offset) ,tvec,[0 offset 0],options);
t_ms = 1000*t;     % time vector in millisec
olddelta = (x(:,2));   % detuning
V = x(:,1);            % cavity probe signal 
%nextoffset = olddelta(end); % assume next pulse starts with end state of previous pulse

% Plotting
figure(1)
subplot(211)
plot(t_ms,abs(V),'b','Linewidth',2)
hold on; grid on
xlabel('time [ms]')
ylabel('abs(V)  [a.u.]')
subplot(212)
plot(t_ms,olddelta,'b','Linewidth',2)
hold on; grid on
xlabel('time [ms]')
ylabel('detuning [Hz]')
set(gca,'Fontsize',12)

%%%%%%%%%%%%%%%% System Identification *************
% Filtered derivative
% 
fnyq = 1/(2*dt); % Nyquist frequency
bw = 50e3; % filter bandwidth in Hz
[b,a] = butter(1,bw/fnyq);      % Find LP-filter parameters
dVdt = (V(3:end)-V(1:end-2))/(2*dt); % Symmetric diff
dVdt = [0;dVdt;0];
dVdtfilt = filtfilt(b,a,dVdt);    % Symmetric forward-backward filter

% Estimating parameters by least squares
y = real(conj(V).*dVdtfilt);
A = [-conj(V).*V real(conj(V).*u) -imag(conj(V).*u)];
x = inv(A'*A)*(A'*y);      % Inverting 3x3 matrix

fhat = x(1)/2/pi
a1hat = x(2)
b1hat = x(3)

%%%%%%%%%%%%%%%% . End Estimation Code . %%%%%%%%%%%%%%%%%


function dx = process2(t,x,offset)
global tvec u vpiezo dt
% function dx = process(t,x,offset)
% t: time
% x: state vector for cavity + detuning dynamics
% offset: Piezo voltage v = 0 gives detuning = offset in stationarity

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

a1 = 2000;          % true a and b parameters
b1 = 200;
% Cavity dynamics
dx(1) = -(whalf-1i*2*pi*delta)*V + (a1+i*b1)*ux;
% Mechanical dynamics, here 2nd order dynamics
dx(2) = deltadot;
dx(3) = - Omega1/Q1*deltadot - Omega1^2*(delta-offset) + Omega1^2*KLFD*abs(V^2) + Omega1^2*Kpiezo*vpiezoint;

end
