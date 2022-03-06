%%  ================== q1=====================================
% model fiber parameters
clear all; 
CF    = 10.0e3; % CF in Hz;   
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 5e-3;   % rise/fall time in seconds
stimdb = 10; % stimulus intensity in dB SPL
% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

best_freqs = [500 4000];
freqs = 125*2.^(0:1/8:8);
intensities = -10:10:80;

n_iters = 5;
rate_tensor_500hz = zeros(n_iters, length(intensities), length(freqs));
CF = 500;
% ----------- for 500Hz ------
for iter=1:n_iters
    fprintf("-----iter num %d \n", iter);
    for i=1:length(intensities)
        for f=1:length(freqs)
                pin = sqrt(2)*20e-6*10^(intensities(i)/20)*sin(2*pi*freqs(f)*t); % unramped stimulus
                pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
                pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
                
                vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
                [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
                
                timeout = (1:length(psth))*1/Fs;
                psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
                psthtime = timeout(1:psthbins:end); % time vector for psth
                pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
                Psth = pr/psthbinwidth; % psth in units of spikes/s
                avg_Psth = sum(Psth)/length(Psth);
                rate_tensor_500hz(iter, i, f) = avg_Psth;
        end
    end
end

rate_tensor_500hz_avg_iters = zeros(length(intensities), length(freqs));
for i=1:length(intensities)
    for f=1:length(freqs)
        rate_tensor_500hz_avg_iters(i,f) = sum(rate_tensor_500hz(:,i,f))/n_iters;
    end
end

figure
    hold on
        for i=1:length(intensities)
            y = rate_tensor_500hz_avg_iters(i,:);
            x=log10(freqs);
            plot(x, y);
        end
    hold off
    title('500 hz at 10 diff intensities')
grid


% -------- for 4 khz ----------
rate_tensor_4khz = zeros(n_iters, length(intensities), length(freqs));
CF = 4000;
for iter=1:n_iters
    fprintf("-----iter num %d \n", iter);
    for i=1:length(intensities)
        for f=1:length(freqs)
                pin = sqrt(2)*20e-6*10^(intensities(i)/20)*sin(2*pi*freqs(f)*t); % unramped stimulus
                pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
                pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
                
                vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
                [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
                
                timeout = (1:length(psth))*1/Fs;
                psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
                psthtime = timeout(1:psthbins:end); % time vector for psth
                pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
                Psth = pr/psthbinwidth; % psth in units of spikes/s
                avg_Psth = sum(Psth)/length(Psth);
                rate_tensor_4khz(iter, i, f) = avg_Psth;
        end
    end
end

rate_tensor_4khz_avg_iters = zeros(length(intensities), length(freqs));
for i=1:length(intensities)
    for f=1:length(freqs)
        rate_tensor_4khz_avg_iters(i,f) = sum(rate_tensor_4khz(:,i,f))/n_iters;
    end
end

figure
    hold on
        for i=1:length(intensities)
            y = rate_tensor_4khz_avg_iters(i,:);
            x=log10(freqs);
            plot(x, y);
        end
    hold off
    title('4k hz at 10 diff intensities')
grid

% ## rate vs intensity 
% --- for 500 hz -----
intensites_vs_rate_500hz = zeros(n_iters, length(intensities));
F0 = 500;
CF = 500;
for iter=1:n_iters
   for i=1:length(intensities)
                stimdb = intensities(i);
                
                pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
                pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
                pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
                
                vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
                [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
                
                timeout = (1:length(psth))*1/Fs;
                psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
                psthtime = timeout(1:psthbins:end); % time vector for psth
                pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
                Psth = pr/psthbinwidth; % psth in units of spikes/s
                avg_Psth = sum(Psth)/length(Psth);
                intensites_vs_rate_500hz(iter, i) = avg_Psth;
   end
end

intensites_vs_rate_500hz_avg = zeros(1, length(intensities));
for i=1:length(intensities)
    intensites_vs_rate_500hz_avg(1,i) = sum(intensites_vs_rate_500hz(:,i))/n_iters;
end

% for 4 khz 
intensites_vs_rate_4khz = zeros(n_iters, length(intensities));
F0 = 4000;
CF = 4000;
for iter=1:n_iters
   for i=1:length(intensities)
                stimdb = intensities(i);
                
                pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
                pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
                pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
                
                vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
                [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
                
                timeout = (1:length(psth))*1/Fs;
                psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
                psthtime = timeout(1:psthbins:end); % time vector for psth
                pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
                Psth = pr/psthbinwidth; % psth in units of spikes/s
                avg_Psth = sum(Psth)/length(Psth);
                intensites_vs_rate_4khz(iter, i) = avg_Psth;
   end
end

intensites_vs_rate_4khz_avg = zeros(1, length(intensities));
for i=1:length(intensities)
    intensites_vs_rate_4khz_avg(1,i) = sum(intensites_vs_rate_4khz(:,i))/n_iters;
end

figure
    hold on
        plot(intensities,intensites_vs_rate_500hz_avg);
        plot(intensities,intensites_vs_rate_4khz_avg);
        title('rate vs intensity')
        legend('500 hz','4 khz');
    hold off
grid

%% ================== q2 ===================
[fivewo, fs] = audioread('fivewo.wav');   %reading entire wav file
ah = fivewo(98250:113500);    % 'ah' sound
sound(ah, fs);
rms_ah = rms(ah);
I_0 = 20*10^(-6);
db_SPL = 20*log10(rms_ah/I_0);
% what is use of RMS??? - mulitply later to feed into ANFs
intensities = -20:10:80;
ah = ah';

CF    = 500; % CF in Hz;   
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = length(ah);  % stimulus duration in seconds
rt = 5e-3;   % rise/fall time in seconds
% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

rate_vs_intensity_500 = zeros(1, length(intensities));
for i=1:length(intensities)
    stimdb = intensities(i);
    pin = ah; % unramped stimulus - ah
    pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
    
    vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
    [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt); 
    
    timeout = (1:length(psth))*1/Fs;
    psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
    psthtime = timeout(1:psthbins:end); % time vector for psth
    pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
    Psth = pr/psthbinwidth; % psth in units of spikes/s
    Psth_avg = sum(Psth)/length(Psth);
    rate_vs_intensity_500(1, i) = Psth_avg;
end

figure
    plot(intensities,rate_vs_intensity_500)
    title('rate vs intensity for ah sound -500 hz BF')
grid
