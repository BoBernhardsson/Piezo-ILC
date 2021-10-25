%close all
clear
         
global whalf KLFD Kpiezo  Omega1  Q1 tvec u vpiezo dt

whalf = 2*pi*500;   % Cavity bandwidth 500Hz
KLFD = 400;         % Lorentz force detuning constant, 
Kpiezo = 1;         % Piezo force constant Hz/Volt, arbitrary
Omega1 = 2*pi*350;  % Main mechanical mode 350 Hz
Q1 = 3;             % its damping ratio zeta = 1/6 
offset = (-100);      % Initial detuning offset in Hz, typically [-1000,1000];
noiselevel = 0*0.001;    % 0.001 corresponds to 60dB SNR (within fmax)

dt = 20e-6;         % sample interval in simulation
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


%% initial simulation with vpiezo = 0, i.e. no Piezo control
options = odeset('RelTol',1e-6,'AbsTol',[1e-5 1e-5 1e-5]);     % setting accuracies
[t,x] = ode45(@(t,y) process(t,y,offset) ,tvec,[0 offset 0],options);
t_ms = 1000*t;     % time vector in millisec
olddelta = (x(:,2));   % detuning
V = x(:,1);            % cavity probe signal 
nextoffset = olddelta(end); % assume next pulse starts with end state of previous pulse

% Plotting
% figure(1)
% subplot(211)
% plot(t_ms,abs(V),'b','Linewidth',2)
% hold on; grid on
% xlabel('time [ms]')
% ylabel('abs(V)  [a.u.]')
% subplot(212)
% plot(t_ms,olddelta,'b','Linewidth',2)
% hold on; grid on
% xlabel('time [ms]')
% ylabel('detuning [Hz]')
% set(gca,'Fontsize',12)


%% test pulse
% Uses input equal to half period of cosine
bliplength = 4e-4;    
nrb = round(bliplength/dt);
blipamp = 100;
blip = blipamp*0.5*(1-cos(2*pi/nrb*(0:nrb)'));


% Impulse response simulation
oldvpiezo = vpiezo;
ind1 = round((tstart + tfill + tpulse*0.5)/dt);

vpiezo(ind1-nrb/2:ind1+nrb/2) = vpiezo(ind1-nrb/2:ind1+nrb/2) + blip;
[t,x] = ode45(@(t,y) process(t,y,offset) ,tvec,[0 offset 0],options);
V = x(:,1);
delta = (x(:,2));
diff = delta - olddelta;
h1 = diff/sum(blip);

figure(2)
subplot(211)
plot(t_ms,vpiezo)
axis([0 0.001*tmax/dt -200 200])
subplot(212)
plot(t_ms,diff);
%axis([0 1000*tmax -5 5])
hold on
drawnow



%% Frequency domain approximate inverse (h)

H1 = fft(h1);
sigma = 0.05 * abs(H1(1));  % higher value gives more high frequency inverse
ufinv = conj(H1)./(sigma^(2)+abs(H1.^2));
uinv = real(ifft(ufinv));
[~,dd]=max(abs(uinv));
filterlength=301;
skip = 150;
l1 = uinv(dd-skip:dd-skip+filterlength-1);

res1 = filter(l1,1,h1); 
% [~, maxind1] = max(abs(vpiezo-oldvpiezo));
% [~, maxind2] = max(abs(res));
% skip1 = maxind2-maxind1;     % used for time alignment of inverse
% [~,skip] = max(abs(l));



%clear u1 uf ufinv uinv

%% Iterative Learning Control Loop

% Initializing v, delta, rmsvalue
%newv = v;
vpiezo = oldvpiezo;
delta = olddelta;
rmsvalue=[];
indpulse = find(tstart <= t & t <= tend);
rmsvalue(1) = rms(delta(indpulse));


% Calibrated gain and phase offsets 
% Gain parameters assumed known from previous experiment
alpha = whalf; 
beta = 0;
alphaest = alpha; 
betaest = beta;
delta0est = offset;

% Setting up some filters
% lowpass filter used for filtering u, V and deltahat 
cutoff1 = 2e3;    % bandwidth in Hz
[B1,A1] = butter(1,cutoff1*dt*2);

% Derivative filter 
cutoff2 = 2e3;   % bandwidth in Hz
[B2,A2]=butter(1,cutoff2*dt*2);

 
 cutoff3 = 3e3;      % Q-filter in ILC
 [B3,A3] = butter(2,cutoff3/(fmax/2));
    
% Old implementation of derivatefilter:
% ff = ifftshift(f)-fmax/2;
% ff1 = 1i*fmax*sin(2*pi*ff/fmax)./(1+(ff/cutoff2).^2);
% derivatefilter = fftshift(real(ifft(ff1)));
% cutoff2 = 3e3;   % bandwidth in Hz
% ff = ifftshift(f)-fmax/2;
% ff(end/2)=0;
% ff1 = 2*pi*1i*ff./(1+(ff/cutoff2).^2);
% derivatefilter = fftshift(real(ifft(ff1)));


%% Main ILC loop
for iter = 1:30

    iter

    % Adding noise to measurements
    
    un = u + noiselevel*max(abs(u))/sqrt(2)*(randn(nt,1)+1i*randn(nt,1));
    Vn = V + noiselevel*max(abs(V))*(randn(nt,1)+1i*randn(nt,1));

    % Estimate detuning
    % Filtering u, V, dV and deltahat 
    
    
    %uf = filtfilt(B,a,un);
    %dum = filter(B,A,un);
    %dum2 = filter(B,A,dum(end:-1:1));
    %uf = dum2(end:-1:1);
    
    unf = filter(B1,A1,un);
    
    %Vf = filtfilt(B,A,Vn);
    %dum = filter(B,A,Vn);
    %dum2 = filter(B,A,dum(end:-1:1));
    %Vf = dum2(end:-1:1);
    
    Vf = filter(B1,A1,Vn);
    
    
%     dVf = conv(derivatefilter,Vn);
%     dVf = dVf(nt/2+1:nt/2+nt);
%     
    
    dum = [0; Vn(3:end)-Vn(1:end-2); 0];
    dum1 = filter(B1,A1,dum);
    %dum2 = filter(BB1,AA1,dum1(end:-1:1));
    %dVf = 1/2/dt*dum2(end:-1:1);
    dVf = 1/2/dt*dum1;
    
    Vreg = 1e-1; % regularization parameter, larger value gives more robust estimates
    deltahatn = 1/2/pi*imag((conj(Vf).*dVf - (alphaest+1i*betaest)*conj(Vf).*unf))./(Vreg^2 + conj(Vf).*Vf);
    deltahatn(1:1+round((tstart+dt)/dt)) = 0;        % dont trust deltahat when V is low
    deltahatn(1+round((tend+tdecay)/dt):end) = 0;
    deltahat = filter(B2,A2,deltahatn);

    changev = conv(l1,deltahat);               % only change Vpiezo  near pulse 
    changev = changev(skip+1:skip+nt);
    changev(1:round((tstart-0.1e-3)/dt)) = 0;    
    changev(round((tend+tdecay/2)/dt):end) = 0;          
    
    % Estimating initial frequency offset 

    initialmeasuretime = 0.3e-3;
    indinitial = find(tstart < t & t < tstart+initialmeasuretime);
    mean1 = mean(deltahat(indinitial))
  
    % If you dont know the V_offset
%     if iter < 5
%         factor1 = 1;  factor2 = 0;
%         
%     else 
%         factor1 = 0.5; factor2 = 1;
%     end
%     vpiezo = filtfilt(B3,A3,vpiezo - factor1*mean1 - factor2*changev);
    
    % If you do know the V_offset
    
    if iter == 1
        vpiezo = -offset*ones(size(vpiezo));   % Guess of rest detuning offset
    else
        factor1 = 0;
        factor2 = 1;
        vpiezo = filtfilt(B3,A3,vpiezo - factor1*mean1 - factor2*changev);
    end
    
    % some plotting
    figure(3) 
    plot(t_ms,deltahatn,'g','Linewidth',1.5)
    hold on;grid on
    plot(t_ms,deltahat,'r','Linewidth',1.5)
    hold on;grid on
    plot(t_ms,delta,'b','Linewidth',1.5)
    xlabel('time [ms]','fontsize',14)
    ylabel('detune [Hz]','fontsize',14)
    set(gca,'Fontsize',14)
    hold off
      
    figure(4)
    subplot(211)
    plot(t_ms,vpiezo,'Linewidth',1.5)
    grid on;
    ylabel('Vpiezo')
    subplot(212)
    plot(t_ms,delta,'Linewidth',1.5)
    hold on
    plot(t_ms(indpulse),delta(indpulse),'r','Linewidth',2)
    grid on
    hold off
    ylabel('detune [Hz]')
    set(gca,'Fontsize',14)
    drawnow
    
    figure(5);
    plot(t_ms,filtfilt(B3,A3,changev))
    ylabel('Change in Vpiezo')
    grid on
    drawnow
    

    
    % simulate new pulse
    [t,x] = ode45(@(t,y) process(t,y,offset) ,tvec,[0 nextoffset 0],options);
    V = x(:,1);
    delta = (x(:,2));
    nextoffset = delta(end);
  
    rmsvalue(iter+1) = rms(delta(indpulse));
end

%% Some final plotting
figure(7)
subplot(211)
plot(t_ms,abs(V),'b','Linewidth',1.5)
grid on
ylabel('abs(V) [a.u]')
xlabel('Time [ms]')
subplot(212)
plot(t_ms,vpiezo,'b','Linewidth',1.5)
grid on
xlabel('Time [ms]')
ylabel('Vpiezo')
grid on

figure(8)
semilogy(rmsvalue,'Linewidth',1.5)
xlabel('iteration')
ylabel('detuning during pulse rms [Hz]')
set(gca,'Fontsize',14)
grid on
hold on