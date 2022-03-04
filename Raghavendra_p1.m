% model fiber parameters
clear all; 
close all; 

CF    = [500.0 4.0e3]; % CF in Hz;   

cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 5e-3;   % rise/fall time in seconds

stimdb = 10; % stimulus intensity in dB SPL
intensities = -10:10:80;

% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

% ---------20 repeations---------
psth_both_freqs_20_iters = zeros(2, 20, 200);
for freq=1:2
        psth_20_iters = [];
        for iter=1:20
            pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*CF(freq)*t); % unramped stimulus
            pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
            pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
             % vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
             vihc = catmodel_IHC(pin,CF(freq),nrep,1/Fs,T*2,cohc,cihc); 
             [synout,psth] = catmodel_Synapse(vihc,CF(freq),nrep,1/Fs,fiberType,implnt); 
    
             timeout = (1:length(psth))*1/Fs;
             psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
             psthtime = timeout(1:psthbins:end); % time vector for psth
             pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
             Psth = pr/psthbinwidth; % psth in units of spikes/s
             psth_20_iters = [psth_20_iters; Psth];
        end

        psth_20_iters_reshaped = reshape(psth_20_iters, 1, 20, 200);
        psth_both_freqs_20_iters(freq, :, :) = psth_20_iters_reshaped;

end





% will write my own output - leave this
figure
subplot(4,1,1)
plot(timeout,[pin zeros(1,length(timeout)-length(pin))])
title('Input Stimulus')

subplot(4,1,2)
plot(timeout,vihc(1:length(timeout)))
title('IHC output')

subplot(4,1,3)
plot(timeout,synout); 
title('Synapse Output')
xlabel('Time (s)')

subplot(4,1,4)
plot(psthtime,Psth)
title('psth')
xlabel('Time (s)')