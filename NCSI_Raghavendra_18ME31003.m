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
% what is use of RMS??? - mulitply later to feed into ANFs
intensities = -20:10:80;

[fivewo, fs] = audioread('fivewo.wav');   %reading entire wav file
ah = fivewo(98250:113499);    % 'ah' sound
sound(ah, fs);
rms_ah = rms(ah);
I_0 = 20*10^(-6);
db_SPL = 20*log10(rms_ah/I_0);
appropirate_factor = zeros(1,length(intensities));

for i=1:length(intensities)
   appropirate_factor(1,i) = (10^(intensities(i)/20)*(20 * (10^-6)))/rms_ah;
end



ah = ah';

CF    = 500; % CF in Hz;   
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = length(ah)/Fs;  % stimulus duration in seconds
% T = 50;
rt = 5e-3;   % rise/fall time in seconds
% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

n_iters = 50;
rate_vs_intensity_500 = zeros(n_iters, length(intensities));

for iter=1:n_iters
        for i=1:length(intensities)
            pin = ah*appropirate_factor(1,i); % unramped stimulus - ah
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
            rate_vs_intensity_500(iter, i) = Psth_avg;
    end
end

rate_vs_intensity_500_avg = zeros(1, length(intensities));
for i=1:length(intensities)
    rate_vs_intensity_500_avg(1, i) = sum(rate_vs_intensity_500(:,i))/n_iters;
end

figure
    plot(intensities,rate_vs_intensity_500_avg)
    title('rate vs intensity for ah sound -500 hz BF')
grid

%% q2-a intensties based on curve
i1 = 0; i2 = 40; i3 = 76; % 76 is also intensity of aah sound, it will be easier to compare 
bank_of_freqs = 125*2.^[0:1/8:6];
num_of_freqs_in_bank = length(bank_of_freqs);
levels_3 = [i1, i2, i3];

[fivewo, fs] = audioread('fivewo.wav');   %reading entire wav file
fivewo = fivewo';
rms5 = rms(fivewo);
I_0 = 20*10^(-6);
db_SPL = 20*log10(rms5/I_0);
appropirate_factor = zeros(1,3);

for i=1:3
   appropirate_factor(1,i) = (10^(levels_3(i)/20)*(20 * (10^-6)))/rms5;
end



cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = length(fivewo)/Fs;  % stimulus duration in seconds
% T = 50;
rt = 5e-3;   % rise/fall time in seconds
% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;




spike_trains = zeros(3, num_of_freqs_in_bank, 312500);

for i=1:length(levels_3)
    for b=1:num_of_freqs_in_bank
            fprintf("\n i=%d, b = %d\n",i,b);
            pin = fivewo*appropirate_factor(1,i); % unramped stimulus - ah
            CF = bank_of_freqs(b);
            vihc = catmodel_IHC(pin,CF,nrep,1/Fs,T*2,cohc,cihc); 
            [synout,psth] = catmodel_Synapse(vihc,CF,nrep,1/Fs,fiberType,implnt);
            psth_reshaped = reshape(psth, 1,1,312500);

            spike_trains(i,b,:) = psth_reshaped;
    end
 
end
%% -------- q2 c----
[fivewo, fs] = audioread('fivewo.wav');   %reading entire wav file
Fs = 100e3;
freqs = 125*2.^[0:1/8:6] ;
fivewo = fivewo';
figure
    spectrogram(fivewo, hann(25.6e-3*Fs), 12.8e-3*Fs, 1:8000, Fs, 'yaxis');
    title('spectrogram')
grid

s3 = spike_trains(3,:, :);
s3 = squeeze(s3);

window_sizes = [4,8,16,32,64,128];
window_sizes = window_sizes*1e-3*Fs;
half_window_sizes = floor(window_sizes)./2;

for w=1:length(window_sizes)
    interval = window_sizes(w)/2:half_window_sizes(w):length(fivewo)-window_sizes(w)/2;
    response_rates = zeros(49, length(interval));
    for f=1:length(freqs)
            for i=1:length(interval)
                response_rates(f, i) = sum(s3(f, interval(i)-half_window_sizes(w)+1:interval(i)+half_window_sizes(w)))/window_sizes(w);
            end
    end
    
    figure(60)
        subplot(2,3,w)
         [ t, f ] = meshgrid( interval, freqs);
         surf(t, f, response_rates,'edgecolor','none');
        colorbar;
        xlim([0,1.5e5]);
        view(2);


end

%% ------------ q3 -------------
[fivewo, fs] = audioread('fivewo.wav');   %reading entire wav file
Fs = 100e3;
bfs = [1, 5, 9,  13, 17, 21,   25, 29, 33,  37, 41];
cmap1 = hsv(11);
window_size = 12.8e-3*Fs;
half_window = floor(window_size/2);
interval = window_size/2:half_window:length(fivewo)-window_size/2;
max_pts = zeros(1, length(interval));
spike_train1 = spike_trains(3,:,:);
spike_train1 = squeeze(spike_train1);
figure(30)
    spectrogram(fivewo, hann(12.8e-3*Fs), 6.4e-3*Fs, 1:8000, Fs, 'yaxis');
    view(3);
    hold on;
for f=1:length(bfs)
    for i=1:length(interval)
        rate1 = spike_train1(bfs(f), (interval(i)-half_window+1) : (interval(i)+half_window));
        % subtracting mean
        mean_rate = mean(rate1);
        fft_rate = abs(fft(rate1 - mean_rate));
        [max_val, max_index] = max(squeeze(fft_rate(1:length(fft_rate)/2)));
        max_pts(i) = max_index*Fs/length(fft_rate);
    end

    scatter3(interval/Fs,max_pts/1000,zeros(1,length(interval))-50,[],cmap1(f,:),'filled', 'MarkerEdgeColor', 'k');
    ylim([0 3]);
end
view(2)

grid


%% part b 

% 250 Hz to 2 kHz
freqs = 125*2.^[0:1/8:4];

% digitize 10 khz
[fivewo, fs] = audioread('fivewo.wav');
audiowrite('fivewo_10k.wav', fivewo, 10e3);
[fivewo2, fs2] = audioread('fivewo_10k.wav');
cutoff = 500;
[b,a] = butter(4,0.5);
y = filtfilt(b,a,fivewo2);
plot(y);

