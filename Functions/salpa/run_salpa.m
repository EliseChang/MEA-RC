function [filteredTrace] = runSALPA(trace, lostTime)

%% Import metadata

homeDir = ("D:\MATLAB\MPhil_scripts\stimulation_tasks\Mona");
cd(homeDir)
metadataSpreadsheet = "MEC.xlsx";

% .xlsx spreadsheet
sheet = 1; %2-spikes, 1-voltage
xlRange = 'A2:H79';
[num,txt,~] = xlsread(metadataSpreadsheet,sheet,xlRange);
voltage_filenames = txt(:,1);
% spike_filenames = txt(:,2);

stim_nodes = num(:,5);

% Parameters for segmenting recording into trials
n_trials = 60;
fs = 25e3;
lost_time = 1e-3*fs; % 1 ms

% Paths
cd(homeDir)
addpath('salpa')
addpath("4s_trials")
addpath("filtered_voltage")
addpath("stim_artifacts")
cd(homeDir)

% Get indices for electrodes
[channel_IDs,channel_pos,~] = get_MEA_coords;
channel_idx = channel_pos(~isnan(channel_pos));
n_channels = length(channel_idx);

for n = 1:length(voltage_filenames)

    disp(voltage_filenames{n})
    cd 4s_trials
    try
        load(strcat(voltage_filenames{n},'.mat'),'dat')
    catch
        cd(homeDir)
        continue
    end
    cd(homeDir)

%     % Clean
%     stim_node_idx = channel_idx(channel_IDs == stim_nodes(n));
%     dat(:,[stim_node_idx, 15]) = nan;

    % Segment recording into trials
    total_samples = length(dat);
    trial_duration = total_samples / n_trials;
    trial_onset_t = [0:trial_duration:total_samples];

    filtered_data = zeros(trial_duration-lost_time,n_channels,n_trials);

    cd stim_artifacts
    mkdir(voltage_filenames{n}); 
    cd(homeDir)

    for trial = 1:n_trials

        start = trial_onset_t(trial) + 1;
        stop =  trial_onset_t(trial+1);
        trial_win = dat(start:stop,:);
%         plotting_win_duration = trial_duration/50;

        for channel = 1:n_channels
            
            if ~any(isnan(trial_win(:,channel)))
%                 subplot_idx = find(channel_pos == channel);
%                 subplot(8, 8, subplot_idx)

                % Get time before depegging
                channel_data = mean_centre(trial_win(:,channel));
                thr = 6*std(channel_data);
                [~,peg_delay] = findpeaks(channel_data,"NPeaks",1, "MinPeakHeight",thr);

                if isempty(peg_delay)
                    warning("Trial %d, channel %d: Stimulation saturation not detected. Check artifact removal.",...
                        trial, channel)
                    peg_delay = 5;
                end

                filtered_trace = salpa(channel_data,'t_blankdepeg',peg_delay);

                % Remove lost time

                filtered_data(:,channel,trial) = filtered_trace(lost_time+1:end);

%                 cd stim_artifacts
%                 cd(voltage_filenames{n})
%                 figure('units','normalized','outerposition',[0 0 1 1])
%                 plot(channel_data)
%                 hold on
%                 plot(filtered_trace)
%                 hold off
%                 aesthetics
%                 xticks([0:trial_duration/8:trial_duration])
%                 xticklabels([0:0.5:4])
%                 xlabel("Time (s)")
%                 ylabel("MEA readout (microvolt)")
%                 saveas(gcf, strcat("Trial_",num2str(trial),"_channel_",num2str(channel),"_filtered.png"))
%                 close all
%                 clear filtered_trace
%                 cd(homeDir)
            else
                continue
            end
        end
        
        clear trial_win

    end
    
    filtered_data = reshape(permute(filtered_data, [1 3 2]), [], n_channels, 1);
    cd filtered_voltage
    save(strcat(voltage_filenames{n},'_filtd.mat'),"filtered_data")
    cd(homeDir)
    clear dat filtered_data

end