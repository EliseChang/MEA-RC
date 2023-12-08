% Process spike-detected MEA data output from MEA-NAP for reservoir
% computing analysis
%   - Get and/or merge spike time variables
%   - Get evoked activity (spike counts/spike rates) within defined windows
%       for reservoir computing analysis
%   - Visualise evoked activity: PSTHs, trial-stacked raster plots,
%       heatmaps

% Created: Elise Chang, February 2023

%% Set parameters and import metadata

% Activity state variable
stateVar = 'spike_count'; % 'spike_rates' | 'spike_prob' | 'latency'
varLabel = "Spike counts"; % 'Spike rates' | 'Spike probabilities' | 'Spike latencies'

% Directories
homeDir = ("D:\MATLAB\MEA-RC");
cd(homeDir)
metadataSpreadsheet = "mecp2RecordingsListNew.xlsx";
preSpikeDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData06Dec2023\1_SpikeDetection\Pre';
postSpikeDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData06Dec2023\1_SpikeDetection\Post';
spikeDir = postSpikeDir;
baselineDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData01Nov2023\1_SpikeDetection\1A_SpikeDetectedData';
voltageDir = 'D:\MATLAB\MEA-NAP\organoids\Nov2023DIV150Stim';
figDir = fullfile(homeDir,'figs');
outputDir = fullfile('C:\Users\elise\Python\ReservoirComputing\data\MEA\organoids',stateVar);

% Metadata
xlSheet = 'Stim';
xlRange = 'A2:M11';
spreadsheetDir = "D:\MATLAB\MEA-NAP";
[num,txt,~] = xlsread(fullfile(spreadsheetDir,metadataSpreadsheet),xlSheet,xlRange);
samples = txt(:,1); % name of sample
ages = num(:,1);
genotypes = txt(:,3);
patternAStimID = txt(:,4);
patternBStimID = txt(:,5); % idx
groundElecID = txt(:,6);
startFrame = num(:,6);
circShift = num(:,7);
startTrial = num(:,8);
endTrial = num(:,9);
nTrials = num(:,10);
preStimFile = txt(:,12);
postStimFile = txt(:,13);

% Paths
addpath(preSpikeDir)
addpath(postSpikeDir)
% addpath(outputDir)
addpath spikes
addpath Functions
funcSubFolders = genpath('/Functions');
addpath(funcSubFolders)
addpath(figDir)
figSubFolders = genpath(figDir);
addpath(figSubFolders);
addpath("D:\MATLAB\MEA-NAP\Functions")

% Parameters for segmenting recording into trials
channelsN = 60;
fs = 25e3; % sampling frequency
ISI = fs; % inter-stimulus interval in frames
psthBin = 5e-3*fs; % in frames
preStimRecLength = 600; % in secs
postStimRecLength = 300; % in secs
minFiringRate = 0.02;

windows = {'0_100ms'}; % {'window0_20','window20_40','window40_60','window60_80','window80_100','window100_120','window120_140','window140_160','window160_180'};
% window_start_frames = 1; % [3*1e-3*fs, 20*1e-3*fs]; % :20*1e-3*fs:180*1e-3*fs;
windowLengths = 0.1*fs; % [17*1e-3*fs, 80*1e-3*fs];
lostTime = 0.003*fs; 

% Load in stimulation protocol
stimProt = double(~readmatrix("MeCP2OrgStimProt1.csv"));
stimProt(1) = []; % remove extra entry in first pos because of zero indexing

% Readout electrodes
outputElecs = [7 5 2 59 56 54]; % indices

figPos = [1 49 1920 955];

%% Optional: if spikes only available in spike matrix and not spike times format, convert spike matrix to spike times
    
for n = 1:length(samples)

    cd(spikeDir)
    
    try
        spikeMatrix = full(load(strcat(samples{n},'.mat'),'cSpikes')).cSpikes;
    catch
        try
            spikeMatrix = full(load(strcat(samples{n},'.mat'),'spikes')).spikes;
        catch
            spikeMatrix = full(load(strcat(samples{n},'.mat'),'spikeMatrix')).spikeMatrix;
        end
    end

    disp(samples{n})

    cd(homeDir)

    mergedSpikeTimes = spikeMatrixToSpikeTimes(spikeMatrix,fs); % produces spikeTimes variable in miliseconds
    
    cd(spikeDir)
    save(strcat(samples{n},'.mat'),'spikeMatrix','mergedSpikeTimes')
    cd(homeDir)

    clear spikeMatrix mergedSpikeTimes


end

%% Plot spike detection results pre- and post-SALPA

saveFigDir = fullfile(figDir, 'stimArtifactRemoval', 'recoveredSpikesWaveforms');
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

methods = {'thr4', 'thr5', 'bior1p5', 'bior1p3', 'db2'}; % spike detection methods used in MEA-NAP
methodsN = length(methods);

for w = 1:length(windows)
    window = windows{w};
%     spikeRatios = zeros(length(samples),methodsN);
    
    disp(['Comparing pre- and post-artifact removal spike detection results: ', window])
    for n = 1:length(samples)

        disp(samples{n})

        % Load voltage data
        cd(voltageDir)
        stimDat = load(strcat(samples{n},'.mat'), 'stimDat').stimDat;
        preFieldName = strcat('preSALPA',window);
        cd(homeDir)

        % Load spike time data
        cd(preSpikeDir)
        preSpikeTimes = load(strcat(samples{n},'_spikes.mat'), 'spikeTimes').spikeTimes;
        preSpikeWaveforms = load(strcat(samples{n},'_spikes.mat'), 'spikeWaveforms').spikeWaveforms;
        cd(homeDir)

        cd(postSpikeDir)
        postSpikeTimes = load(strcat(samples{n},'_spikes.mat'), 'spikeTimes').spikeTimes;
        postSpikeWaveforms = load(strcat(samples{n},'_spikes.mat'), 'spikeWaveforms').spikeWaveforms;
        cd(homeDir)
        clear spikeTimes

        % Adjust stim prot as needed
        startTrialIdx = startTrial(n);
        endTrialIdx = endTrial(n);
        shiftTrials = circShift(n);
        if ~isnan(shiftTrials)
            patternSeq = circshift(stimProt, trialsN - shiftTrials);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end
        firstPattern = patternSeq(1) + 1; % zero-indexing

        % Find stimulation electrodes
        patternAStimIdx = getElectrodeIdx(str2num(patternAStimID{n}));
        patternBStimIdx = getElectrodeIdx(str2num(patternBStimID{n}));
        patterns = [patternAStimIdx;patternBStimIdx];
        stimElecs = patterns(firstPattern,:);
        
        % Get start and end time of first trial
        trialsN = nTrials(n);
        windowLength = windowLengths(w); % in frames -- length of window
        framesN = windowLength*trialsN;
        trialT = 1:windowLength:framesN+windowLength;
        times = [trialT(1),trialT(2)];
        
        % Find additional spikes postSALPA
        incSpikes = cell(1, channelsN); % initialise cell array to filter the post-SALPA spike time struct
        recovSpikes = cell(1, channelsN); % initialise cell array to store times of recovered spikes
        preSpikeCounts = zeros(1, methodsN); % initialise vector to store times channel spike counts
        postSpikeCounts = zeros(1, methodsN); % initialise vector to store times channel spike counts
       
        for ch = 1:channelsN
            for m = 1:methodsN
                method = methods{m};
                preSpikeTrain = preSpikeTimes{ch}.(method);
                postSpikeTrain = postSpikeTimes{ch}.(method);
                preSpikeCounts(m) = preSpikeCounts(m) + numel(preSpikeTrain);
                postSpikeCounts(m) = postSpikeCounts(m) + numel(postSpikeTrain);
                incSpikes{1, ch}.(method) = ~ismember(postSpikeTrain,preSpikeTrain); % logical variable for each of original spike times post-SALPA
                recovSpikes{1, ch}.(method) = setdiff(postSpikeTrain,preSpikeTrain); % actual times of additional spikes
                clear preSpikeTrain postSpikeTrain
            end
        end

%         spikeRatios(n,:) = postSpikeCounts ./ preSpikeCounts;
        clear postSpikeCounts preSpikeCounts

        % Plot additional spikes on output electrodes
        
        cd(saveFigDir)
        if ~isfolder(samples{n})
            mkdir(samples{n})
        end
        cd(samples{n})

        for ch = 1:numel(outputElecs)
            channel = outputElecs(ch);
            t = plotSpikeWaveforms(channel, incSpikes, recovSpikes, postSpikeWaveforms, stimDat.(preFieldName), fs, methods, figPos);
            saveas(t, strcat("Electrode ", num2str(channel), ".png"))
            close all
        end
        cd(homeDir)
        
        clear incSpikes recovSpikes stimDat

%         f = plotSpikeDetectionChecksStim(stimDat, window, preSpikeTimes, postSpikeTimes, stimElecs, outputElecs,... 
%             times, fs, figPos);

    end

end

clear saveFigDir

disp("Finished running section.\n")

%% Merge spike times across spike detection methods
    
cd(spikeDir)
disp("MERGING SPIKES")

for n = 1:length(samples)

    load(strcat(samples{n},'_spikes.mat'));
    disp(samples{n})

    mergedSpikeTimes = cell(1,channelsN);

    for ch = 1:channelsN
        [mergedSpikes,~, ~] = mergeSpikes(spikeTimes{ch}, 'all'); % in seconds
        mergedSpikeTimes{1,ch} = mergedSpikes;
        clear mergedSpikes
    end

    save(strcat(samples{n},'_spikes.mat'), "mergedSpikeTimes", "-append")
    clear spikeDetectionResult spikeTimes mergedSpikeTimes spikeWaveforms thresholds

end

cd(homeDir)

disp("Finished running section.\n")

%% % Get baseline values

disp("GETTING PRE- AND POST-STIM FIRING RATES")

prePostStimFiringRates = zeros(length(samples),3);

for n = 1:length(samples)

    disp(samples{n})

    cd(baselineDir)

    % Load pre-stim baseline spike data and merge spikes
    load(strcat(preStimFile{n},'_spikes.mat'));
    mergedSpikeTimes = cell(1,channelsN);
    for ch = 1:channelsN
        [mergedSpikes,~, ~] = mergeSpikes(spikeTimes{ch}, 'all'); % in seconds
        mergedSpikeTimes{1,ch} = mergedSpikes;
        clear mergedSpikes
    end
    
    % Calculate spiking probabilities for whole network and each channel
    totalSpikes = cell2mat(mergedSpikeTimes);
    preStimNetSpikeRate = length(totalSpikes) / preStimRecLength; % spikes across all channels
%     [binnedNetSpikeCounts,~] = histcounts(totalSpikes, 0:1e-3:preStimRecLength);
%     preStimNetSpikeProb = sum(binnedNetSpikeCounts > 0) / (preStimRecLength * 1e3 + 1);
%     channelSpikeCounts = cell2mat(cellfun(@(a) length (a), mergedSpikeTimes, 'UniformOutput', false)); 
%     channelSpikeRates = channelSpikeCounts / preStimRecLength;
%     channelSpikeRates(channelSpikeRates < minFiringRate) = [];
%     preStimNetSpikeRate = mean(channelSpikeRates);
    prePostStimFiringRates(n,1) = preStimNetSpikeRate;
%     preStimChSpikeProb = channelSpikeCounts / preStimRecLength * 1e-3;
    clear spikeTimes mergedSpikeTimes totalSpikes totalSpikeCount uniqueSpikeCount channelSpikeCounts channelSpikeRates preStimNetSpikeRate

    % Load post-stim baseline spike data and merge spikes
    load(strcat(postStimFile{n},'_spikes.mat'));
    mergedSpikeTimes = cell(1,channelsN);
    for ch = 1:channelsN
        [mergedSpikes,~, ~] = mergeSpikes(spikeTimes{ch}, 'all'); % in seconds
        mergedSpikeTimes{1,ch} = mergedSpikes;
        clear mergedSpikes
    end

    % Calculate spiking probabilities for each channel
    totalSpikes = cell2mat(mergedSpikeTimes);
    postStimNetSpikeRate = length(totalSpikes) / postStimRecLength; % spikes across all channels
%     channelSpikeCounts = cell2mat(cellfun(@(a) length (a), mergedSpikeTimes, 'UniformOutput', false)); 
%     channelSpikeRates = channelSpikeCounts / preStimRecLength;
%     channelSpikeRates(channelSpikeRates < minFiringRate) = [];
%     postStimNetSpikeRate = mean(channelSpikeRates);
    prePostStimFiringRates(n,2) = postStimNetSpikeRate;
%     postStimChSpikeProb = channelSpikeCounts / postStimRecLength * 1e-3;
    clear spikeTimes mergedSpikeTimes channelSpikeCounts channelSpikeRates postStimNetSpikeRate
    
%     cd(spikeDir)
%     try
%         save(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeRate','preStimNetSpikeProb','preStimChSpikeProb',...
%             'postStimChSpikeProb', '-append'); % 'thresholds'
%     catch 'MATLAB:save:couldNotWriteFile' % if file for stimulation recording does not yet exist
%         save(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeRate','preStimNetSpikeProb','preStimChSpikeProb',...
%             'postStimChSpikeProb')
%     end
%     clear mergedSpikeTimes totalSpikeCount preStimSpikeRate preStimSpikeProbNet postStimSpikeRate postStimSpikeProbNet
%     cd(homeDir)

end

cd(homeDir)

disp("Finished running section.\n")

%% % Create stimulus-locked trials and get trial activity

% Make folder to save fig
saveFigDir = fullfile(figDir, 'networkResponse', 'rasterPlot');
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

for w = 1:length(windows)
    window = windows{w};
    windowLength = windowLengths(w); % in frames -- length of window
        
    disp(['GETTING TRIAL ACTIVITY: ', window, '\n'])
    
    for n = 1:length(samples)

        trialsN = nTrials(n);
        framesN = windowLength*trialsN;
        trialOnT = lostTime:windowLength:framesN;
        trialOffT = windowLength:windowLength:framesN;
            
        cd(spikeDir)
        load(strcat(samples{n},'_spikes.mat'), 'mergedSpikeTimes', 'allX')
        if ~exist('allX','var')
            allX = struct();
        end

        cd(homeDir)
       
        disp(samples{n})
        
        % If needed, rearrange trial order in stimProt
        startTrialIdx = startTrial(n);
        endTrialIdx = endTrial(n);
        if startTrialIdx ~= 1
            k = length(stimProt) - startTrialIdx + 1; % ensures that original kth trial becomes last trial
            patternSeq = circshift(stimProt,k);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end

        % Get stim electrodes indices
        patternAStimIdx = getElectrodeIdx(str2num(patternAStimID{n}));
        patternBStimIdx = getElectrodeIdx(str2num(patternBStimID{n}));
        
        % Initiate variable x for storing trial data and convert spike
        % time units
        x = zeros(trialsN,channelsN);
        mergedSpikeTimesFrames = cellfun(@(a) a*fs, mergedSpikeTimes, 'UniformOutput', false); % convert from s to frames
        
        % for PSTH, get counts across all non-stim channels to exclude
        % channels as necessary
        if ~strcmp(stateVar, 'latency')
            spikeCounts = cell2mat(cellfun(@(a) length (a), mergedSpikeTimes, 'UniformOutput', false));
            excludeUnits = find((spikeCounts > framesN/fs*100) |...
            (spikeCounts < 10)); % exclude units with mean firing rate >100 Hz across recording
            % and units which record fewer than 10 spikes across the whole
            % recording -- assumed to be noise
            incChMergedSpikeTimesFrames = mergedSpikeTimesFrames;
            incChMergedSpikeTimesFrames([excludeUnits,patternAStimIdx,patternBStimIdx]) = [];
            allSpikeTimes = cell2mat(incChMergedSpikeTimesFrames);
            unqiueSpikeTimes = unique(allSpikeTimes);
            nBins = windowLength/psthBin - 1;
            psthData = zeros(trialsN,nBins); % store mean binned spike counts across electrodes per trial
            
        end
    
        for trial = 1:trialsN
    
            start = trialOnT(trial);
            stop = trialOffT(trial);
            if ~strcmp(stateVar, 'latency')
                % define PSTH bins
                binEdges = start:psthBin:stop;
            end
                
            for ch = 1:channelsN
    
                channelSpikes = mergedSpikeTimesFrames{ch}; % all spikes for channel
                trialMask = (channelSpikes <= stop) & (channelSpikes > start); % which spikes fall within this trial
                spikeCount = sum(trialMask);
                if strcmp(stateVar, "spike_count")
                    x(trial,ch) = spikeCount;
                elseif strcmp(stateVar, "spike_rates")
                    spikeRate = spikeCount / (windowLength/fs);
                    x(trial,ch) = spikeRate;
                elseif strcmp(stateVar, "spike_prob")
                    x(trial,ch) = double(any(trialMask));
                elseif strcmp(stateVar, "latency")
                    if any(trialMask)
                        trialSpikeTimes = channelSpikes(trialMask);
                        firstSpikeFrame = trialSpikeTimes(1) - start + lostTime;
                    else
                        firstSpikeFrame = NaN;
                    end
                    x(trial,ch) = firstSpikeFrame;
                end
                clear t
        
            end

                % for PSTH, bin spike counts across all channels
                
                if ~strcmp(stateVar, 'latency')
                    if strcmp(stateVar, 'spike_prob')
                        % for spike probabilities, do not count spikes across electrodes arriving in the same 1 ms bin
                        trialSpikes = unqiueSpikeTimes((unqiueSpikeTimes <= stop) & (unqiueSpikeTimes > start));
                    else
                        trialSpikes = allSpikeTimes((allSpikeTimes <= stop) & (allSpikeTimes > start));
                    end
                    [counts,~] = histcounts(trialSpikes,binEdges);
                    psthData(trial, :) = counts;
                end
        end
        
        allX.(strcat(stateVar, "_", window)) = x;

        % Save variables
        saveDataDir = fullfile(outputDir,window);
        if ~isfolder(saveDataDir)
            mkdir(saveDataDir)
        end
        cd(saveDataDir)
        
        writematrix(x, strcat(samples{n},".csv"))
        cd(homeDir)
        
        cd(spikeDir)
        save(strcat(samples{n},'_spikes.mat'),'allX','-append');
        cd(homeDir)

        % Plot basic raster
        % If needed, rearrange trial order in stimProt
        startTrialIdx = startTrial(n);
        endTrialIdx = endTrial(n);
        if startTrialIdx ~= 1
            k = length(stimProt) - startTrialIdx + 1; % ensures that original kth trial becomes last trial
            patternSeq = circshift(stimProt,k);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end
        mask = patternSeq == 0;
        p = plotTrialRaster(x',mask, ["Pattern A", "Pattern B"], varLabel,figPos);
        cd(saveFigDir)
        saveas(p,samples{n},'_trial_raster_',window,'.png')
        close all
        cd(homeDir)

%         % Plot PSTH
%         % get baseline
%         cd(baselineDir)
%         if strcmp(stateVar, 'spike_prob')
%             baseline = load(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeProb').preStimNetSpikeProb;
%         else
%             baseline = load(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeRate').preStimNetSpikeRate;
%         end
%         cd(homeDir)
% 
%         % Plot PSTH
%         p = psth(psthData, stateVar, stimProt, psthBin/fs*1e3, baseline);
%         cd .\network_response\PTSH
%         saveas(p,samples{n},'png')
%         cd(homeDir)
%         close all
%         prePostSALPACounts(n,2) = sum(psthData, "all");
        
%         % Plot trial heatmap
%         t = plotTrialHeatmap(figPos,patternSeq,psthData,psthBin,fs);
%         cd .\network_response\trialHeatmaps
%         saveas(t,strcat(samples{n},'_trial_heatmap.png'))
%         cd(homeDir)
%         close all
       
    end

end
        
clear x allX mergedSpikeTimes mergedSpikeTimesFrames psthData saveFigDir

disp("Finished running section.\n")

%% Plot trial-stacked raster plot

for n = 1:length(samples)

    if types{n} ~= "BASELINE"

        cd(spikeDir)
        spikeTimes = load(strcat(samples{n},'.mat')).mergedSpikeTimes;
        cd(homeDir)
    
%         % Plot spike counts at each electrode
%         cd network_response/Raster_plots
%         mkdir(samples{n}); 
%         cd(homeDir)
%         for ch = 1:channelsN
%             t = seconds(spikeTimes{1,ch})*1e-3; % convert from ms to s
%             if ~isempty(t)
%                 trialAssignment = assignSpikeTimeTrials(t, seconds(stimOnT));
%                 activeTrials = unique(trialAssignment);
%                 alignmentTimes = seconds(activeTrials-1)*5;
%                 h = spikeRasterPlot(t,trialAssignment,'AlignmentTimes',alignmentTimes,...
%                     'XLimits',seconds([lostTime/fs,0.1]),'XLabelText', "Time (s)", "YLabelText", "Trial",...
%                     'TitleText',strcat(samples{n},"_electrode_",num2str(ch)));
%                 saveas(h, strcat("network_response/Raster_plots/",samples{n},"/",...
%                     samples{n},"_electrode_",num2str(ch),".png"))
%                 close all
%             end
%         end

        % Plot spike counts summed across all electrodes
        cd network_response/Raster_plots
        mkdir(samples{n});
        cd(homeDir)
        activeChannelSpikeTimes = spikeTimes(~cellfun('isempty',spikeTimes));
        allSpikeTimes = cell2mat(activeChannelSpikeTimes);
        t = seconds(allSpikeTimes)*1e-3; % convert from ms to s
        if ~isempty(t)
            trialAssignment = assignSpikeTimeTrials(t, seconds(stimOnT));
            activeTrials = unique(trialAssignment);
            alignmentTimes = seconds(activeTrials-1)*5;
            h = spikeRasterPlot(t,trialAssignment,'AlignmentTimes',alignmentTimes,...
                'XLimits',seconds([startFrame/fs,0.1]),'XLabelText', "Time (s)", "YLabelText", "Trial",...
                'TitleText',samples{n});
            saveas(h, strcat("network_response/Raster_plots/",samples{n},".png"))
            close all
        end
    end
end

%% Electrophysiological characterisation
% 
% [~,Params.GrpNm] = findgroups(genotypes);
% [~,Params.DivNm] = findgroups(ages);
% Params.bursting = 1;
% Params.fs = fs;
% Params.ephysMetricsToPlot = {'numActiveElec','FRmean','FRstd',...
%     'meanNBstLengthS','NBurstRate','meanNumChansInvolvedInNbursts'};
% Info.duration_s = durationS;
% 
% for n = 1:length(samples)
% 
%     cd(spikeDir)
%     try
%         load(strcat(samples{n},'_spikes.mat'),'mergedSpikeTimes')
%     catch
%         load(strcat(samples{n},'.mat'),'mergedSpikeTimes')
%     end
%     cd(homeDir)
%    
%     disp(samples{n})
%     
%     mergedSpikeTimesFrames = cellfun(@(a) a*fs*1e-3, mergedSpikeTimes, 'UniformOutput', false);
%     spikeMatrix = spikeTimesToSpikeMatrix(mergedSpikeTimesFrames, framesN);
%     Ephys = firingRatesBursts(spikeMatrix,Params,Info);
% 
%     cd(spikeDir)
%     try
%         save(strcat(samples{n},'.mat'),'Ephys','-append');
%     catch
%         save(strcat(samples{n},'_spikes.mat'),'Ephys','-append');
%     end
%     cd(homeDir)
%     
% end
% 
% saveEphysStats(samples, Params, genotypes, ages, homeDir)

%% Plot heatmaps

% Make folder to save fig
saveFigDir = fullfile(figDir, 'heatmaps', stateVar);
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

for w = 1:length(windows)

    windowLength = windowLengths(w);

    for n = 1:length(samples)
        
        disp(samples{n})
                
        % Get activity
        cd(spikeDir)
        allX = load(strcat(samples{n},'_spikes.mat'), 'allX').allX;
        cd(homeDir)
        
        % Get stim nodes
        patternAStimIdx = getElectrodeIdx(str2num(patternAStimID{n}));
        patternBStimIdx = getElectrodeIdx(str2num(patternBStimID{n}));
        groundElecIdx = getElectrodeIdx(str2num(groundElecID{n})); %#ok<ST2NM>
        removeElecs = [patternAStimIdx,patternBStimIdx,groundElecIdx];

        % If needed, rearrange trial order in stimProt
        startTrialIdx = startTrial(n);
        endTrialIdx = endTrial(n);
        if startTrialIdx ~= 1
            k = length(stimProt) - startTrialIdx + 1; % ensures that original kth trial becomes last trial
            patternSeq = circshift(stimProt,k);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end

        x = allX.(strcat(stateVar, "_", windows{w}));
        if strcmp(stateVar, 'latency')
            x = x / fs * 1000; % convert from frames to milliseconds
             minVal = [];
             maxVal = [];
%             minVal = lostTime / fs * 1000; % min. latency in milliseconds
%             maxVal = windowLength / fs * 1000; % max. latency in milliseconds
        else
            minVal = [];
            maxVal = [];
        end
        maskA = patternSeq == 0;
        maskB = patternSeq == 1;
        datA = mean(x(maskA,:),1, "omitnan");
        datB = mean(x(maskB,:),1, "omitnan");
        
%         [channels,coords] = getCoordsFromLayout('MCS60');
        % Plot heatmap
%         f = electrodeHeatMapsEC(samples{n},NaN,datA,datB,removeElecs,channels,coords);
        f = plotMEAHeatMap(datA,datB,removeElecs,stateVar,["Pattern A", "Pattern B"],minVal,maxVal);
        cd(saveFigDir)
        saveas(f,strcat(samples{n},'_',stateVar,'_',window,'pattern_comparison.png'))
        cd(homeDir)
        close all
    end        
    clear allX x
end
clear saveFigDir

disp("Finished running section.\n")

%% Plot pre- and post-stim heatmaps

% Make folder to save fig
saveFigDir = fullfile(figDir, 'heatmaps', stateVar);
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

for n = 1:length(samples)
    
    disp(samples{n})
    
    % Get activity
    cd(spikeDir)
    load(strcat(samples{n},'_spikes.mat'),'preStimChSpikeProb', 'postStimChSpikeProb');
    cd(homeDir)
    groundElecIdx = groundElecID(n); % getElectrodeIdx(str2num(groundElecID{n})); %#ok<ST2NM> 

    % Plot heatmap
    %TODO: change to highlight/remove electrodes
    f = plotMEAHeatMap(preStimChSpikeProb,postStimChSpikeProb,groundElecIdx,samples{n},'p');
    cd(saveFigDir)
    saveas(f,strcat(samples{n},'_',stateVar,'_',window,'_pre_vs_post-stim.png'))
    cd(homeDir)
    close all

    clear groundElecIdx preStimChSpikeRate postStimChSpikeRate

end

clear saveFigDir

disp("Finished running section.\n")


%% Plot distances between activity state matrices
% 
for w=1:length(windows)
    distances = zeros(length(pairName),1);
    for n = 1:length(pairName)
        
        disp(pairName{n})
        
        cd(spikeDir)
        
        % Get activity for pattern 0
        cd(spikeDir)
        try
            allX1 = load(strcat(pattern0Rec{n},'_spikes.mat')).allX;
        catch
            allX1 = load(strcat(pattern0Rec{n},'.mat')).allX;
        end
        cd(homeDir)
    
        try
            stimNode1 = pattern0StimElec{n};
        catch
            stimNode1 = pattern0StimElec(n);
        end
    
        if isa(stimNode1,'char')
            stimNode1 = str2num(stimNode1);
        end

        x1 = allX1.(windows{w});
    
        % Get activity for pattern 1
        cd(spikeDir)
        try
            allX2 = load(strcat(pattern1Rec{n},'_spikes.mat')).allX;
        catch
            allX2 = load(strcat(pattern1Rec{n},'.mat')).allX;
        end
        cd(homeDir)
    
        try
            stimNode2 = pattern1StimElec{n};
        catch
            stimNode2 = pattern1StimElec(n);
        end
        
        if isa(stimNode2,'char')
            stimNode2 = str2num(stimNode2);
        end
        x2 = allX2.(windows{w});
    
        dist = vecnorm(x1'-x2');
        distances(n) = mean(dist);
    end

%     dat1 = mean(x1)/baselineSpikeRate; % fold change in firing rate relative to baseline
%     dat2 = mean(x2)/baselineSpikeRate; % fold change in firing rate relative to baseline
%     
%     f = plotMEAHeatMap(dat1,dat2,stimNode1,stimNode2,pairName{n},'foldChange');
%     cd .\network_response\heatmaps
%     saveas(f,pairName{n},'png')
%     cd(homeDir)
%     close all

end
% for n = 1:length(pairName)
%     
%     disp(pairName{n})
%     
%     % Get spike times for pattern 0
%     cd(spikeDir)
%     try
%         spikeTimes = load(strcat(pattern0Rec{n},'_spikes.mat')).mergedSpikeTimes;
%     catch
%         spikeTimes = load(strcat(pattern0Rec{n},'.mat')).mergedSpikeTimes;
%     end
% 
%     cd(homeDir)
% 
%     spikeMatrix0 = spikeTimesToSpikeMatrix(spikeTimes, durationS, fs); % spike times given in ms 
%     stimNode = stimElecIdx(1,:);
%     spikeMatrix0(:,stimNode) = 0;
%     clear spikeTimes stimNode
% 
%     % Get spike times for pattern 1
%     cd(spikeDir)
%     try
%         spikeTimes = load(strcat(pattern1Rec{n},'_spikes.mat')).mergedSpikeTimes;
%     catch
%         spikeTimes = load(strcat(pattern1Rec{n},'.mat')).mergedSpikeTimes;
%     end
%         
%     cd(homeDir)
% 
%     spikeMatrix1 = spikeTimesToSpikeMatrix(spikeTimes, durationS, fs);
%     stimNode = stimElecIdx(1,:);
%     spikeMatrix1(:,stimNode) = 0;
%     clear spikeTimes stimNode
%     
%     dsF = 25; % downsample factor of 25 means that spike rates given per ms
%     spikeMatrix0 = downSampleSum(spikeMatrix0,framesN/dsF);
%     spikeMatrix1 = downSampleSum(spikeMatrix1,framesN/dsF);
%     plotWin = [1, 0.1*fs/dsF]; % up to 100 ms
%     f = plotStateDist(spikeMatrix0, spikeMatrix1, trialsN, trialLength, dsF, plotWin);
%     cd .\network_response\State_distances
%     saveas(f,pairName{n},'png')
%     cd(homeDir)
%     close all
% 
%     clear spikeMatrix0 spikeMatrix1
% 
% end