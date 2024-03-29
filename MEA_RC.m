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
stateVar = 'spike_count'; % 'spike_rates' | 'spike_count' | 'spike_prob' | 'latency'
varLabel = "Spike counts"; % 'Spike rates' | 'Spike counts' | 'Spike probabilities' | 'Spike latencies'

% Directories
homeDir = ("D:\MATLAB\MEA-RC");
cd(homeDir)
metadataSpreadsheet = 'mecp2RecordingsListNew.xlsx'; % file name
spreadsheetDir = "D:\MATLAB\MEA-NAP\metadata";
preSpikeDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData08Jan2024\1_SpikeDetection\1A_SpikeDetectedData\Pre';
postSpikeDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData08Jan2024\1_SpikeDetection\1A_SpikeDetectedData\Post';
saveSpikeDir = 'D:\MATLAB\MEA-RC\spikes\organoids'; % save spikes after combinging and merging
baselineDir = 'D:\MATLAB\MEA-NAP\outputs\OutputData20Nov2023\ExperimentMatFiles'; % ExperimentMatFolder from pipeline run of baseline data
baselineDirDate = '20Nov2023';
voltageDir = 'D:\MATLAB\MEA-NAP\organoids\Nov2023DIV150Stim';
figDir = fullfile(homeDir,'figs');
activityOutputDir = fullfile(homeDir,'activityData');
RCOutputDir = fullfile('C:\Users\elise\Python\ReservoirComputing\data\MEA\organoids',stateVar);

% Metadata
xlSheet = 'Stim';
xlRange = 'A2:M11';
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
% addpath(preSpikeDir)
addpath(postSpikeDir)
% addpath(outputDir)
addpath spikes
addpath Functions
funcSubFolders = genpath('/Functions');
addpath(funcSubFolders)
addpath(figDir)
addpath("D:\MATLAB\MEA-NAP\Functions")

% Parameters for segmenting recording into trials
channelsN = 60;
fs = 25e3; % sampling frequency
ISI = fs; % inter-stimulus interval in frames
psthBin = 5e-3*fs; % in frames
preStimRecLength = 600; % in secs
postStimRecLength = 300; % in secs
minFiringRate = 1; % minimum firing rate for baseline activity, Hz

windows = {'0_40ms', '0_80ms','0_120ms','0_160ms','0_200ms'};  % {'window0_20','window20_40','window40_60','window60_80','window80_100','window100_120','window120_140','window140_160','window160_180'};
window_start_frames = 0; % [3*1e-3*fs, 20*1e-3*fs]; % :20*1e-3*fs:180*1e-3*fs;
windowLengths = [40e-3*fs, 80e-3*fs, 120e-3*fs, 160e-3*fs, 200e-3*fs]; % 200e-3*fs [17*1e-3*fs, 80*1e-3*fs];
lostTime = 3e-3*fs;
wholeWindow = '0_200ms';
wholeWindowLength = 200e-3*fs;

% Load in stimulation protocol
stimProt = double(~readmatrix(fullfile('stimuli', 'MeCP2OrgStimProt1.csv')));
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

%% Validating stimulus artifact removal and spike detection

saveFigDir = fullfile(figDir, 'stimArtifactRemoval', 'spikesWaveforms');
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

methods = {'thr4', 'thr5', 'bior1p5', 'bior1p3', 'db2'}; % spike detection methods used in MEA-NAP
methodsN = length(methods);

%     spikeRatios = zeros(length(samples),methodsN); % stores ratio of post:pre-SALPA spike rates for each method
spikeRatios = zeros(channelsN,length(samples));
disp(['Comparing pre- and post-artifact removal spike detection results: ', wholeWindow])
for n = 1:length(samples)

    disp(samples{n})

    % Load voltage data
    cd(voltageDir)
    stimDat = load(strcat(samples{n},'.mat'), 'stimDat').stimDat;
    preFieldName = strcat('preSALPA',wholeWindow);
    postFieldName = strcat('postSALPA',wholeWindow);
    cd(homeDir)

    % Load spike time data
    cd(preSpikeDir)
    preSpikeTimes = load(strcat(samples{n},'_spikes_',wholeWindow,'.mat'), 'spikeTimes').spikeTimes;
    preSpikeWaveforms = load(strcat(samples{n},'_spikes_',wholeWindow,'.mat'), 'spikeWaveforms').spikeWaveforms;
    cd(homeDir)

    cd(postSpikeDir)
    postSpikeTimes = load(strcat(samples{n},'_spikes_',wholeWindow,'.mat'), 'spikeTimes').spikeTimes;
    postSpikeWaveforms = load(strcat(samples{n},'_spikes_',wholeWindow,'.mat'), 'spikeWaveforms').spikeWaveforms;
    cd(homeDir)

    % Adjust stim prot as needed
    trialsN = nTrials(n);
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
    times = [1,wholeWindowLength];
    
    % Find additional spikes postSALPA
    recovSpikeMask = cell(1, channelsN); % initialise cell array to store masks of the unique post-SALPA spikes
    recovSpikeTimes = cell(1, channelsN); % initialise cell array to store times of recovered spikes
    lostSpikeMask = cell(1, channelsN); % initialise cell array to store masks of the unique pre-SALPA spikes
    lostSpikeTimes = cell(1, channelsN); % initialise cell array to store times of the unique pre-SALPA spikes

%     preSpikeCounts = zeros(1, methodsN); % initialise vector to store spike counts by method
%     postSpikeCounts = zeros(1, methodsN); % initialise vector to store spike counts by method
   
    combinedSpikeTimes = cell(1, channelsN);
    for ch = 1:channelsN
        [preMergedSpikes,~, ~] = mergeSpikes(preSpikeTimes{ch}, 'all');
        preSpikeChCount = numel(preMergedSpikes);
        [postMergedSpikes,~, ~] = mergeSpikes(postSpikeTimes{ch}, 'all');
        postSpikeChCount = numel(postMergedSpikes);
        spikeRatios(ch,n) = postSpikeChCount / preSpikeChCount;
        for m = 1:methodsN
            method = methods{m};
            preSpikes = preSpikeTimes{ch}.(method);
            postSpikes = postSpikeTimes{ch}.(method);
             
%             % Running total of spikes counts across electrode for each method
%             preSpikeCounts(m) = preSpikeCounts(m) + numel(preSpikes);
%             postSpikeCounts(m) = postSpikeCounts(m) + numel(postSpikes);

            % Check if any 'recovered' spikes have been jittered by +/- 0.5 ms
            checkSpikes = ~ismember(postSpikes,preSpikes);
            for s = 1:length(postSpikes)
                if checkSpikes(s)
                    offsets = abs(preSpikes - postSpikes(s));
                    if any(offsets < 1e-3)
                        checkSpikes(s) = 0;
                    end
                end
            end
            recovSpikeMask{1, ch}.(method) = checkSpikes;
            recovSpikeTimes{1, ch}.(method) = postSpikes(checkSpikes);
            clear checkSpikes

            % Check if any 'lost' spikes have been jittered by +/- 0.5 ms
            checkSpikes = ~ismember(preSpikes,postSpikes);
            for s = 1:length(preSpikes)
                if checkSpikes(s)
                    offsets = abs(postSpikes - preSpikes(s));
                    if any(offsets < 1e-3)
                        checkSpikes(s) = 0;
                    end
                end
            end
            lostSpikeMask{1, ch}.(method) = checkSpikes;
            lostSpikeTimes{1, ch}.(method) = preSpikes(checkSpikes); % actual times of additional spikes
            % Store in combinedSpikeTimes
            if any(checkSpikes)
                combinedSpikeTimes{ch}.(method) = sort([postSpikes; preSpikes(checkSpikes)]);
            else
                combinedSpikeTimes{ch}.(method) = postSpikes;
            end
            clear checkSpikes preSpikes postSpikes
        end
    end

%     spikeRatios(:,n) = postSpikeCounts ./ preSpikeCounts;
%     clear postSpikeCounts preSpikeCounts

    % Plot recovered spikes on output electrodes
    cd(saveFigDir)
    if ~isfolder(samples{n})
        mkdir(samples{n})
    end
    cd(samples{n})
    for ch = 1:numel(outputElecs)
        channel = outputElecs(ch);
        origTrace = stimDat.(preFieldName);
        t = plotSpikeWaveforms(channel, recovSpikeMask, recovSpikeTimes, postSpikeWaveforms, origTrace, fs, methods, figPos); % 1 flag for recovered spikes
        if ~isempty(t)
            saveas(t, strcat("Electrode ", num2str(channel), " Recovered Spikes.png"))
        end
        close all
        clear origTrace
    end
    cd(homeDir)

    % Plot lost spikes on output electrodes
    cd(saveFigDir)
    if ~isfolder(samples{n})
        mkdir(samples{n})
    end
    cd(samples{n})
    for ch = 1:numel(outputElecs)
        channel = outputElecs(ch);
        origTrace = stimDat.(postFieldName);
        t = plotSpikeWaveforms(channel, lostSpikeMask, lostSpikeTimes, preSpikeWaveforms, origTrace, fs, methods, figPos); % 0 flag for lost pikes
        if ~isempty(t) % empty plots without spikes will return nan
            saveas(t, strcat("Electrode ", num2str(channel), " Lost Spikes.png"))
        end
        close all
        clear origTrace
    end
    cd(homeDir)
    
    cd(saveSpikeDir)
    save(strcat(samples{n},'_spikes_',wholeWindow,'.mat'), 'combinedSpikeTimes')
    cd(homeDir)
    
    clear preSpikeTimes postSpikeTimes recovSpikeMask lostSpikeMask recovSpikeTimes lostSpikeTimes combinedSpikeTimes stimDat preSpikeWaveforms postSpikeWaveforms postSpikeTimes

%         f = plotSpikeDetectionChecksStim(stimDat, window, preSpikeTimes, postSpikeTimes, stimElecs, outputElecs,... 
%             times, fs, figPos);

end

%     for n = 1:length(samples)
%         % still rough -- ignore RHS plot, need to define removeElecs
%         cd figs/stimArtifactRemoval/firingRateHeatmaps
%         f = plotMEAHeatMap(spikeRatios(:,n),ones(60,1),removeElecs, 'foldChangePostSALPA');
%         saveas(f, samples{n}, 'png')
%         cd(homeDir)
%         close all
%     end

clear saveFigDir

disp("Finished running section", newline)

%% Merge spike times across spike detection methods

disp("MERGING SPIKES")

cd(saveSpikeDir)

% for w = 1:length(windows)
%     window = windows{w};
for n = 1:length(samples)

    spikeTimes = load(strcat(samples{n},'_spikes_',wholeWindow,'.mat'), 'combinedSpikeTimes').combinedSpikeTimes;

    disp(samples{n})

    mergedSpikeTimes = cell(1,channelsN);

    for ch = 1:channelsN
        [mergedSpikes,~, ~] = mergeSpikes(spikeTimes{ch}, 'all'); % in seconds
        mergedSpikeTimes{1,ch} = mergedSpikes;
        clear mergedSpikes
    end

    save(strcat(samples{n},'_spikes_', wholeWindow,'.mat'), "mergedSpikeTimes", "-append")
    clear spikeDetectionResult spikeTimes mergedSpikeTimes spikeWaveforms thresholds

end
% end

cd(homeDir)

disp("Finished running section", newline)

%% % Get baseline values and compare pre- and post-stimulation

disp("GETTING PRE- AND POST-STIM FIRING RATES")

% Make folder to save fig
saveFigDir = fullfile(figDir, 'prePostStim');
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

ephysVars = {'FRmean','meanNBstLengthS',...%'numActiveElecs',...
'numNbursts',...
'meanNumChansInvolvedInNbursts'};
prePostStimEphysStats = struct;
for v = 1:length(ephysVars)
    var = ephysVars{v};
    prePostStimEphysStats.(var) = zeros(length(samples),2);
end

for n = 1:length(samples)

    disp(samples{n})

    cd(baselineDir)
    
    % Extract ephys vars
    preStimEphys = load(strcat(preStimFile{n},'_', baselineDirDate, '.mat'), 'Ephys').Ephys;
    postStimEphys = load(strcat(postStimFile{n},'_', baselineDirDate, '.mat'), 'Ephys').Ephys;
    for v = 1:length(ephysVars)
        var = ephysVars{v};
        prePostStimEphysStats.(var)(n,1) = preStimEphys.(var);
        prePostStimEphysStats.(var)(n,2) = postStimEphys.(var);
    end

    preStimNetSpikeRate = preStimEphys.FRmean;
    preStimNetSpikeCount = preStimNetSpikeRate * preStimRecLength;
    postStimNetSpikeRate = preStimEphys.FRmean;
    postStimNetSpikeCount = postStimNetSpikeRate * postStimRecLength;

    cd(activityOutputDir)  
    try
        save(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeCount','preStimNetSpikeRate', ...
            'postStimNetSpikeCount', 'postStimNetSpikeRate','-append');
    catch 'MATLAB:save:couldNotWriteFile' % if file for stimulation recording does not yet exist
        save(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeCount','preStimNetSpikeRate',...
            'postStimNetSpikeCount', 'postStimNetSpikeRate')
    end
    clear preStimNetSpikeCount preStimNetSpikeRate postStimNetSpikeCount postStimNetSpikeRate
    cd(homeDir)

end

ephysVarNames = {'Mean firing rate (Hz)', 'N. active electrodes',...
    'N. bursts', 'Burst size (electrodes)'};
f = figure('Position', figPos);
t = tiledlayout("flow");
for v = 1:length(ephysVars)
    var = ephysVars{v};
    nexttile
    boxplot(prePostStimEphysStats.(var))
    xticklabels({'Pre-stim', 'Post-stim'})
    ylabel(ephysVarNames{v})
end

cd(saveFigDir)
saveas(t,'all_samples_pre_post_stim_comparison.png')
cd(homeDir)
close all
disp("Finished running section", newline)

%% % Create stimulus-locked trials and get trial activity

% Make folder to save fig
saveFigDir = fullfile(figDir, 'networkResponse', 'rasterPlot');
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

for w = 1:length(windows)
    window = windows{w};
    windowLength = windowLengths(w); % in frames -- length of window
        
    disp(['GETTING TRIAL ACTIVITY: ', window, newline])
    
    for n = 1:length(samples)

        trialsN = nTrials(n);
        framesN = wholeWindowLength*trialsN;
        trialOnT = lostTime:windowLength:framesN;
        trialOffT = trialOnT + windowLength - lostTime - 1;
            
        cd(saveSpikeDir)
        load(strcat(samples{n},'_spikes_', wholeWindow, '.mat'), 'mergedSpikeTimes')
%         load(strcat(samples{n},'_spikes_',window,'.mat'), 'mergedSpikeTimes')
        cd(homeDir)

        cd(activityOutputDir)
        load(strcat(samples{n},'_spikes.mat'));
        if ~exist('allX','var')
            allX = struct();
        end
        cd(homeDir)
       
        disp(samples{n})
        
        % If needed, rearrange trial order in stimProt
        startTrialIdx = startTrial(n);
        endTrialIdx = endTrial(n);
        shiftTrials = circShift(n);
        if ~isnan(shiftTrials)
            patternSeq = circshift(stimProt, trialsN - shiftTrials);
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
            psthData = zeros(trialsN, nBins, channelsN); % store mean binned spike counts across electrodes per trial
            
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
                trialSpikes = channelSpikes(trialMask);
                [psthCounts,~] = histcounts(trialSpikes,binEdges);
                psthData(trial, :, ch) = psthCounts;
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

%             % for PSTH, bin spike counts across all channels
%             
%             if ~strcmp(stateVar, 'latency')
%                 if strcmp(stateVar, 'spike_prob')
%                     % for spike probabilities, do not count spikes across electrodes arriving in the same 1 ms bin
%                     trialSpikes = unqiueSpikeTimes((unqiueSpikeTimes <= stop) & (unqiueSpikeTimes > start));
%                 else
%                     trialSpikes = allSpikeTimes((allSpikeTimes <= stop) & (allSpikeTimes > start));
%                 end
%                 [counts,~] = histcounts(trialSpikes,binEdges);
%                 psthData(trial, :) = counts;
%             end

        end
        
        allX.(strcat(stateVar, "_", window)) = x;

        % Save variables
        saveDataDir = fullfile(RCOutputDir,window);
        if ~isfolder(saveDataDir)
            mkdir(saveDataDir)
        end
        cd(saveDataDir)
        
        writematrix(x, strcat(samples{n},".csv"))
        cd(homeDir)
        
        cd(activityOutputDir)
        try
            save(strcat(samples{n},'_spikes.mat'),'allX','-append');
        catch 'MATLAB:save:couldNotWriteFile' % if file for stimulation recording does not yet exist
            save(strcat(samples{n},'_spikes.mat'),'allX');
        end
        cd(homeDir)

        % ACTIVITY PLOTS
        % sum across all or subset of electrodes/exclude trials
%         sumPsthData = sum(psthData, 3);

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
        mask = patternSeq == 0; % 1 is flag for trials always plotted on LHS

%         % Plot trial-by-trial raster
%         % Sort channels in order of stimulation electrodes then distance from stimulation electrodes
        [~, allChannelsOrder, ~] = getMEACoords;
        allChannelsOrder = allChannelsOrder(~isnan(allChannelsOrder));
        nonStimElecs = setdiff(allChannelsOrder, [patternAStimIdx, patternBStimIdx], 'stable');
        channelsOrder = [patternAStimIdx, patternBStimIdx, nonStimElecs];
        ylines = [numel(patternAStimIdx), numel(patternAStimIdx) + numel(patternBStimIdx)];
        labels = {'Pattern A', 'Pattern B'};
%         p = plotTrialRaster(x', mask, channelsOrder, ylines, labels, ["Pattern A", "Pattern B"], varLabel,figPos);
%         cd(saveFigDir)
%         saveas(p,strcat(samples{n},'_trial_raster_',window,'.png'))
%         close all
%         cd(homeDir)

%         % Plot PSTH
%         % get baseline
%         if ~strcmp(stateVar,'latency')
%             cd(activityOutputDir)
%             if strcmp(stateVar,'spike_count')
%                 preStimBaseline = load(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeCount').preStimNetSpikeCount;
%                 postStimBaseline = load(strcat(samples{n},'_spikes.mat'),'postStimNetSpikeCount').postStimNetSpikeCount;
%             elseif strcmp(stateVar,'spike_rates')
%                 preStimBaseline = load(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeRate').preStimNetSpikeRate;
%                 postStimBaseline = load(strcat(samples{n},'_spikes.mat'),'postStimNetSpikeRate').postStimNetSpikeRate;
%             elseif strcmp(stateVar,'spike_prob')
%                 preStimBaseline = load(strcat(samples{n},'_spikes.mat'),'preStimNetSpikeProb').preStimNetSpikeProb;
%                 postStimBaseline = load(strcat(samples{n},'_spikes.mat'),'postStimNetSpikeProb').postStimNetSpikeProb;
%             end
%             cd(homeDir)
%             preStimBaselinePerWindow = preStimBaseline*((psthBin/fs) / preStimRecLength);
%             postStimBaselinePerWindow = postStimBaseline*((psthBin/fs) / postStimRecLength);
%             p = psth(sumPsthData, stateVar, mask, psthBin/fs*1e3, [preStimBaselinePerWindow, postStimBaselinePerWindow]);
%             cd figs\PTSH
%             saveas(p,strcat(samples{n}, ' ', wholeWindow), 'png')
%             cd(homeDir)
%             close all
%         end

        % Plot trial heatmap
        t = plotTrialHeatmap(figPos,patternSeq,sumPsthData,psthBin,fs);
        cd figs\PTSH
        saveas(t,strcat(samples{n},'_trial_heatmap.png'))
        cd(homeDir)

        close all
    end

end
        
clear x allX mergedSpikeTimes mergedSpikeTimesFrames psthData saveFigDir

disp("Finished running section", newline)

%% Plot trial-stacked raster plot

for n = 1:length(samples)

    if types{n} ~= "BASELINE"

        cd(activityOutputDir)
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
        cd(activityOutputDir)
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
        shiftTrials = circShift(n);
        if ~isnan(shiftTrials)
            patternSeq = circshift(stimProt, trialsN - shiftTrials);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end

        x = allX.(strcat(stateVar, "_", windows{w}));
        maskA = patternSeq == 0;
        maskB = patternSeq == 1;
        
        if strcmp(stateVar, 'latency')
            x = x / fs * 1000; % convert from frames to milliseconds
            trialsA = x(maskA,:);
            trialsB = x(maskB,:);
            minVal = 3;
            maxVal = 40;
            activeTrialsA = trialsA(min(trialsA, [], 2) <= maxVal, :);
            activeTrialsB = trialsB(min(trialsB, [], 2) <= maxVal, :);
            activeTrialsA(activeTrialsA > maxVal) = NaN;
            activeTrialsB(activeTrialsB > maxVal) = NaN;
            datA = mean(activeTrialsA, 1, "omitnan");
            datB = mean(activeTrialsB, 1, "omitnan");
%             minVal = lostTime / fs * 1000; % min. latency in milliseconds
%             maxVal = windowLength / fs * 1000; % max. latency in milliseconds
        else
            minVal = [];
            maxVal = [];
            trialsA = x(maskA,:);
            trialsB = x(maskB,:);
            datA = mean(trialsA,1, "omitnan");
            datB = mean(trialsB,1, "omitnan");
        end
        
        
        
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

disp("Finished running section", newline)

%% Compare windows

% Make folder to save fig
saveFigDir = fullfile(figDir, 'networkResponse', 'windowComparisons', stateVar);
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

% for w = 1:length(windows)
data = zeros(length(samples),length(windows)*2);
colsB = [1:length(windows)]*2;
colsA = colsB - 1;

for w = 1:length(windows)

    windowLength = windowLengths(w);

    for n = 1:length(samples)
        
        disp(samples{n})
                
        % Get activity
        cd(activityOutputDir)
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
        shiftTrials = circShift(n);
        if ~isnan(shiftTrials)
            patternSeq = circshift(stimProt, trialsN - shiftTrials);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end

        x = allX.(strcat('spike_count', "_", windows{w}));

        maskA = patternSeq == 0;
        maskB = patternSeq == 1;
        datA = mean(x(maskA,:),"all");
        datB = mean(x(maskB,:),"all");
        data(n,colsA(w)) = datA;
        data(n,colsB(w)) = datB;
    end        
    clear allX x
end
clear saveFigDir

disp("Finished running section", newline)

%% Compare windows

% Make folder to save fig
saveFigDir = fullfile(figDir, 'networkResponse', 'windowComparisons', stateVar);
if ~isfolder(saveFigDir)
    mkdir(saveFigDir)
end

% for w = 1:length(windows)
data = zeros(length(samples),trialsN);
colsB = [1:length(windows)]*2;
colsA = colsB - 1;

for w = 1:length(windows)

    windowLength = windowLengths(w);

    for n = 1:length(samples)
        
        disp(samples{n})
                
        % Get activity
        cd(activityOutputDir)
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
        shiftTrials = circShift(n);
        if ~isnan(shiftTrials)
            patternSeq = circshift(stimProt, trialsN - shiftTrials);
            patternSeq = patternSeq(startTrialIdx:endTrialIdx);
        else
            patternSeq = stimProt(startTrialIdx:endTrialIdx);
        end

        x = allX.(strcat('spike_count', "_", windows{w}));

        maskA = patternSeq == 0;
        maskB = patternSeq == 1;
        datA = mean(x(maskA,:),"all");
        datB = mean(x(maskB,:),"all");
        data(n,colsA(w)) = datA;
        data(n,colsB(w)) = datB;
    end        
    clear allX x
end
clear saveFigDir

disp("Finished running section", newline)

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