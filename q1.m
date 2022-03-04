% model fiber parameters
clear all; 
close all;

stimdb = 10; % stimulus intensity in dB SPL
CF    = 10.0e3; % CF in Hz; 
best_freqs = [500 4000];
freqs = 125*2.^(0:1/8:8);
intensities = -10:10:80;

cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 5e-3;   % rise/fall time in seconds
% PSTH parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;


n_iters = 2;
rate_tensor = zeros(length(best_freqs), n_iters, length(freqs), length(intensities), 200);
% 200 is the size of psth

for bf=1:length(best_freqs)
    for f=1:length(freqs)
        for intense=1:length(intensities)
                for i=1:n_iters
                    fprintf("iter num %d \n",i)
                    pin = sqrt(2)*20e-6*10^(intensities(intense)/20)*sin(2*pi*freqs(f)*t); % unramped stimulus
                    pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
                    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
                    
                    vihc = catmodel_IHC(pin,best_freqs(bf),nrep,1/Fs,T*2,cohc,cihc); 
                    [synout,psth] = catmodel_Synapse(vihc,best_freqs(bf),nrep,1/Fs,fiberType,implnt); 
                    
                    timeout = (1:length(psth))*1/Fs;
                    psthbins = round(psthbinwidth*Fs);  % number of psth bins per psth bin
                    psthtime = timeout(1:psthbins:end); % time vector for psth
                    pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
                    Psth = pr/psthbinwidth; % psth in units of spikes/s

                    Psth_reshaped = reshape(Psth, 1,1,1,1,200);
                    rate_tensor(bf,i,f,intense,:) = Psth_reshaped;
                end
        end
    end
end

bf_500_rates = rate_tensor(1,:, :, :, :);
bf_500_rates_psth_avg = zeros(1, n_iters, length(freqs), length(intensities), length(freqs));
for iter=1:n_iters
    for intense=1:length(intensities)
        for f=1:length(freqs)
            bf_500_rates_psth_avg(1, iter, f, intense, f) = sum(bf_500_rates(1,iter,f,intense,:))/200; 
        end
    end
end

bf_500_rates_iters_avg = zeros(1, 1, length(freqs), length(intensities), length(freqs));
for intense=1:length(intensities)
    for freq=1:length(freqs)
        for iter=1:n_iters
            bf_500_rates_iters_avg(1, 1,f,intense,f) = sum(bf_500_rates_psth_avg(1,:,f,intense,f))/n_iters;
        end
    end
end


% TODO figure 
figure
    hold on
        for intense=1:length(intensities)
            for f=1:length(freqs)
                y = bf_500_rates_iters_avg(1,1,f,intense,:);
                y = rehsape(y,1,65);
                y = bf_500_rates_iters_avg(1,1,)
            end
        end
    hold off
grid



% rate_tensor_avg = zeros(length(best_freqs), length(freqs), length(intensities), 200);
% for i=1:200
%     for bf=1:length(best_freqs)
%         for f=1:length(freqs)
%             for intense=1:length(intensities)
%                 rate_tensor_avg(bf, f, intense, i) = sum(rate_tensor(bf, :, f, intense, i))/n_iters;
%             end
%         end
%     end
% end
% 
% rate_tensor_psth_avg = zeros(length(best_freqs), length(freqs), length(intensities),length(freqs));
% for bf=1:length(best_freqs)
%     for intense=1:length(intensities)
%             for f=1:length(freqs)
%                 rate_tensor_psth_avg(bf,f,intense,f) = sum(rate_tensor_avg(bf,f,intense, :))/200;
%             end
%     end
% end

% for 500 Hz
figure
    for intense=1:length(intensities)
        hold on
%             x_axis = log(freqs);
            y_axis = reshape(rate_tensor_psth_avg(1,f,intense, :),1,length(freqs));
            semilogx(freqs,y_axis);
            title('500 Hz BF - Rate vs Log(freq)')
        hold off
    end
grid

% 4 Khz
figure
    for intense=1:length(intensities)
        hold on
%             x_axis = log(freqs);
            y_axis = reshape(rate_tensor_psth_avg(2,f,intense, :),1,length(freqs));
            semilogx(freqs,y_axis);
            title('4 KHz BF - Rate vs Log(freq)')
        hold off
    end
grid
% figure
% subplot(4,1,1)
% plot(timeout,[pin zeros(1,length(timeout)-length(pin))])
% title('Input Stimulus')
% 
% subplot(4,1,2)
% plot(timeout,vihc(1:length(timeout)))
% title('IHC output')
% 
% subplot(4,1,3)
% plot(timeout,synout); 
% title('Synapse Output')
% xlabel('Time (s)')
% 
% subplot(4,1,4)
% plot(psthtime,Psth)
% title('psth')
% xlabel('Time (s)')