% Run simulation
global u Ib vpiezo 

% Initialize system and simulation parameters 
simulationfrequency = 50e3;
par = initsystem(simulationfrequency);
[h,H,Voffsethat] = piezoident(par);
l=ilcfilterdesign(H,par);

% Initialize signals
u = par.ufill*(par.tstart < par.tvec & par.tvec <= par.tstart + par.tfill) +...
    par.uduringbeam*(par.tstart+par.tfill < par.tvec & par.tvec < par.tend);
u = u.';
u + par.noiselevel*max(abs(u))/sqrt(2)*(randn(par.nt,1)+1i*randn(par.nt,1));
Ib = par.beam*(par.tstart+par.tfill < par.tvec & par.tvec < par.tend);
Ib = Ib.';
indpulse = find(par.tstart <= par.tvec & par.tvec <= par.tend); % index for pulse


% Choose alternative for estimation of initial offset:
% 1) Oracle, assuming H correct estimated
 Voffset = (par.deltaref-par.cavity.offsetHz) / sum(par.cavity.Kpiezovector);
% 2) From piezoident     
% Voffset = Voffsethat
vpiezo = Voffset * ones(size(par.tvec))';

% Initialize ILC simulation
xinit = zeros(2*par.cavity.nrmodes+1,1);
rmsvalue=[];

% ILC simulation
for iter = 1:30
    [x,delta] = cavitysimulator(xinit,par);
    rmsvalue(iter) = rms(delta(indpulse)-par.deltaref);
    V = x(:,1);
    V = V + par.noiselevel*max(abs(V))*(randn(par.nt,1)+1i*randn(par.nt,1)); 
    [deltahat,delta0hat] = estimatedetuning(par,V);
    
    % plot results
    figure(11);
    subplot(211)
    plot(par.tvec,abs(V));grid on; 
    axis([0 par.tmax 0 1.2])
    xlabel('Time [s]');ylabel('abs(V) [a.u.]')
    subplot(212)
    plot(par.tvec,vpiezo);grid on; 
    axis([0 par.tmax -400 400])
    xlabel('Time [s]');ylabel('Vpiezo')
    figure(12)
    plot(par.tvec,delta,'b',par.tvec(indpulse),delta(indpulse),'r',par.tvec,deltahat,'g')
    grid on; title('detuning [Hz]'); xlabel('Time [s]')
    yscale = max(20,max(abs(delta)));
    axis([0 par.tmax -yscale yscale]);  
    figure(13)
    plot(par.tvec,abs(u))
    grid on; title('LLRF u'); shg   
    
    vpiezo = ilcupdate(par,l, deltahat,delta0hat,par.Lskip+5);    % skip+5 even better !
    xinit = (x(end,:)).';      % for next pulse, needed for nonzero offsets
end
figure(20)
semilogy(rmsvalue)
grid on;  ylabel('rms error [Hz]'); hold on
xlabel('iteration')
axis([0 length(rmsvalue) 0.1 1000])
