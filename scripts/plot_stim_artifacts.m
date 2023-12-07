%% Import metadata

homeDir = ("D:\MATLAB\MPhil_scripts\stimulation_tasks\Mona");
cd(homeDir)
metadataSpreadsheet = "MEC.xlsx";

% Metadata
sheet = 2; %2-spikes, 1-voltage
voltageFlg = 0;
spikesFlg = 1;
baseline_recs = txt(:,2);
xlRange = 'A2:L17';
[num,txt,~] = xlsread(metadataSpreadsheet,sheet,xlRange);
samples = txt(:,1);
hub_stim_recs = txt(:,3);
per_stim_recs = txt(:,4);
hub_stim_nodes = num(:,6);
per_stim_nodes = num(:,8);

% Paths
% spikeDetectedDir = .\spikes;
% addpath(spikeDetectedDir)
addpath("voltage"); addpath("spikes"); addpath("4s_trials")
addpath("network_response\")
addpath("stim_artifacts\")
addpath("salpa")

% Parameters for segmenting recording into trials
n_trials = 60;
n_channels = 60; % TTX recordings = 12
fs = 25e3;
trial_duration = 5*fs;
total_samples = trial_duration * n_trials;
trial_onset_t = 0:trial_duration:total_samples;

stim_duration = 0;
new_trial_duration = 4*fs; % 4s trials
new_trial_onset_t = 0:trial_duration+stim_duration:...
    (n_trials-1)*(trial_duration+stim_duration);
new_trial_offset_t = new_trial_onset_t + new_trial_duration;

%% Plot
for channel = 1:60
    channelSpikesUnfiltd = mergedSpikeTimesUnfiltdFrames{channel};
    channelSpikesFiltd = mergedSpikeTimesFiltdFrames{channel};
    for trial = 1:trialsN
        start = trialOnT(trial) + 1;
        stop =  trialOffT(trial); %trial_offset_t(trial);
        trial_win_unfiltd = dat(start:stop,channel);
        trial_win_unfiltd = (trial_win_unfiltd-mean(trial_win_unfiltd));
        trial_win_filtd = filteredDat(start:stop,channel);
        trial_win_filtd = (trial_win_filtd-mean(trial_win_filtd));

        maxPoint = max(max([trial_win_unfiltd, trial_win_filtd]));

        trialSpikesTimesUnfiltd = channelSpikesUnfiltd((channelSpikesUnfiltd > start) & (channelSpikesUnfiltd < stop));
        trialSpikesUnfiltd = trialSpikesTimesUnfiltd - start;
        trialSpikesTimesFiltd = channelSpikesFiltd((channelSpikesFiltd > start) & (channelSpikesFiltd < stop));
        trialSpikesFiltd = trialSpikesTimesFiltd - start;
        
        figure("Position",[1, 49, 2560,1315])
        plot(trial_win_unfiltd)
        hold on
%             aesthetics
        xticks([625,1250,1875,2500])
        xticklabels([25,50,75,100])
        xlim([0,2500])
        xlabel("Time post-stimulus (ms)")
        ylabel("MEA output")
        aesthetics

        plot(trial_win_filtd)
        y1 = ones(1,length(trialSpikesUnfiltd))* maxPoint;
        scatter(trialSpikesUnfiltd,y1,"filled","v")
        y2 = ones(1,length(trialSpikesFiltd))*(maxPoint-5);
        scatter(trialSpikesFiltd,y2,"filled","v")
        hold off

        cd Stim_artifacts
        saveas(gcf, strcat("Channel_",num2str(channel), "_trial",num2str(trial)),'png')
        cd ..
    end
end

%     plot(trace)
%     hold on
%     yline([-7*thr -6*thr -5*thr -4*thr -3*thr -2*thr -thr thr 2*thr 3*thr 4*thr 5*thr 6*thr 7*thr])
    
%     end
%     
%     cd(homeDir)
%     clear dat
% end