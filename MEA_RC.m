% Process spike-detected MEA data output from MEA-NAP for reservoir
% computing analysis
%   - Get and/or merge spike time variables
%   - Get evoked activity (spike counts/spike rates) within defined window
%       for reservoir computing analysis
%   - Visualise evoked activity: PSTHs, trial-stacked raster plots,
%       heatmaps

% Created: Elise Chang, February 2023

%% Set parameters and import metadata

homeDir = ("D:\MATLAB\MEA-RC");
cd(homeDir)
metadataSpreadsheet = "mecp2timepoints.csv";
spikeDir = fullfile(homeDir,"spikes"); % 'E:\spikes';
stateVar = 'spike_rates'; % spike_counts
% window = ''; % '3-20';
% outputDir = fullfile('C:\Users\elise\conn2res\data\MEA\MEA-Mecp2_Mona_stim_dataset_conn2res',stateVar, window);

% Metadata
xlSheet = 'spikes';
xlRange = 'A2:G39';
[num,txt,~] = xlsread(metadataSpreadsheet,xlSheet,xlRange);
samples = txt(:,1); % name of sample
ages = num(:,1);
genotypes = txt(:,3);
types = txt(:,4); % baseline, hub or peripheral stimulation
% stimNodes = num(:,5); % single node
stimNodes = txt(:,6); % multiple nodes

% Paths
addpath(spikeDir)
% addpath(outputDir)
addpath spikes
addpath Functions
addpath network_response\
addpath("D:\MATLAB\MPhil_scripts\MEA-NAP\Functions")

% Parameters for segmenting recording into trials
trialsN = 60; % number of stimuli delivered; TTX recordings = 12
channelsN = 60;
fs = 25e3; % sampling frequency
ISI = 5*fs; % inter-stimulus interval in frames

% windows={'window0_20','window20_40','window40_60','window60_80','window80_100','window100_120','window120_140','window140_160','window160_180'};
% window_start_times=0:20*1e-3*fs:160*1e-3*fs;
stimLength = 5; % number of frames
% for w=1:length(windows)

%     lostTime = window_start_times(w);
    lostTime = 0e-3*fs; %3e-3*fs; % in frames -- time removed due to stimulus artifacts and early response component
    trialLength = 100e-3*fs; % in frames -- length of trial
    
    % % For samples which have already been segmented into trials of the correct
    % % length
    % trialOnT = 1:trialLength:framesN; % is actually cut off at the end
    % trialOffT = trialLength:trialLength:framesN;
    % stimOnT = trialOnT/fs;
    % framesN = trialsN*trialLength;
    
    % For samples which have NOT been segmented into trials
    trialOnT = lostTime:ISI+stimLength:...
        lostTime+(trialsN-1)*(ISI+stimLength);
    trialOffT = trialOnT + trialLength;
    stimOnT = [1,ISI:ISI:(trialsN-1)*ISI]/fs;
    framesN = trialsN*ISI;
    
    durationS = framesN/fs;
    
    %% Optional: if spikes only available in spike matrix and not spike times format, convert spike matrix to spike times
    % 
    % for n = 1:length(samples)
    % 
    %     cd spikes
    %     
    %     try
    %         spikeMatrix = full(load(strcat(samples{n},'.mat'),'cSpikes')).cSpikes;
    %     catch
    %         try
    %             spikeMatrix = full(load(strcat(samples{n},'.mat'),'spikes')).spikes;
    %         catch
    %             spikeMatrix = full(load(strcat(samples{n},'.mat'),'spikeMatrix')).spikeMatrix;
    %         end
    %     end
    % 
    %     disp(samples{n})
    % 
    %     cd(homeDir)
    % 
    %     mergedSpikeTimes = spikeMatrixToSpikeTimes(spikeMatrix,fs); % produces spikeTimes variable in miliseconds
    %     
    %     cd(spikeDir)
    %     save(strcat(samples{n},'.mat'),'spikeMatrix','mergedSpikeTimes')
    %     cd(homeDir)
    % 
    %     clear spikeMatrix mergedSpikeTimes
    % 
    % 
    % end
    % % Merge spike times across spike detection methods
    % 
    % cd(spikeDir)
    % disp("MERGING SPIKES")
    % 
    % for n = 1:length(samples)
    % 
    %     try
    %         load(strcat(samples{n},'_spikes.mat'));
    %     catch
    %         load(strcat(samples{n},'.mat'));
    %     end
    % 
    %     disp(samples{n})
    % 
        mergedSpikeTimesFiltd = cell(1,channelsN);
    
        for ch = 1:channelsN
            [mergedSpikes,~, ~] = mergeSpikes(spikeTimesFiltd{ch}, 'all');
            mergedSpikeTimesFiltd{1,ch} = mergedSpikes;
        end
    % 
    %     save(strcat(samples{n},'.mat'),'spikeDetectionResult','spikeTimes','mergedSpikeTimes','spikeWaveforms'); % 'thresholds'
    % 
    %     clear spikeDetectionResult spikeTimes mergedSpikeTimes spikeWaveforms thresholds
    % 
    % end
    % 
    % cd(homeDir)
    
    %% Create stimulus-locked trials and get trial activity for reservoir computing analysis
    
    disp("GETTING TRIAL ACTIVITY")
    nActiveTrials = zeros(1,length(samples));
    % networkX = zeros(length(samples),trialsN);
    cvX = zeros(length(samples),1);
    
    for n = 1:length(samples)
        
        if ~strcmp(types{n}, "BASELINE")
    
        cd(spikeDir)
%         try
%             load(strcat(samples{n},'_spikes.mat'),'mergedSpikeTimes','spikeMatrix')
%         catch
%             load(strcat(samples{n},'.mat'),'mergedSpikeTimes','spikeMatrix')
%         end
        try
            load(strcat(samples{n},'_spikes.mat'))
        catch
            load(strcat(samples{n},'.mat'))
        end

        cd(homeDir)
       
        disp(samples{n})
%         if w == 1
% %             baselineX = struct();
%             allX = struct(); %'window 0_20',{},'window 20_40',{},'window 40_60',{},'window 60_80',{},'window 80_100',{},'window 100_120',{},'window 120_140',{},'window 140_160',{});
%         end
        if strcmp(stateVar, "spike_rates")
    
            x = zeros(trialsN,channelsN);
            mergedSpikeTimesFiltdFrames = cellfun(@(a) a*fs*1e-3, mergedSpikeTimesFiltd, 'UniformOutput', false); % convert from ms to frames
    
        elseif strcmp(stateVar, "spike_counts")
    
            x = zeros(trialsN,floor(trialLength/fs*1e3)); % in ms
            spikeMatrix = spikeTimesToSpikeMatrix(mergedSpikeTimes, durationS*1e3); % downsample to 1 ms bins
            stimNode = stimNodes(n);
            if isnan(stimNode)
                spikeMatrix(:,15) = []; % remove ground electrode
            else
                spikeMatrix(:,[15,stimNode]) = [];
            end
            spikeCounts = sum(spikeMatrix,2);
    
        end
    
        for trial = 1:trialsN
    
            start = trialOnT(trial);
            stop = trialOffT(trial);
    
            if strcmp(stateVar, "spike_rates")
    
%                 % TEMP: get network-summed rates + CV
%                 
                    activeChannelSpikeTimes = mergedSpikeTimesFrames(~cellfun('isempty',mergedSpikeTimesFrames));
                    allSpikeTimesFrames = cell2mat(activeChannelSpikeTimes);
                    spikeCount = sum((allSpikeTimesFrames <= stop) & (allSpikeTimesFrames > start));
                    spikeRate = spikeCount / (trialLength/fs);
                    networkX(trial) = spikeRate;
    
                for ch = 1:channelsN
        
                    t = mergedSpikeTimesFrames{ch};
                    spikeCount = sum((t <= stop) & (t > start));
                    spikeRate = spikeCount / (trialLength/fs);
                    x(trial,ch) = spikeRate;
        
                    clear t
        
                end
    
            elseif strcmp(stateVar, "spike_counts")
    
                msStart = floor(start/fs*1e3) + 1;
                msStop = floor(stop/fs*1e3);
                trialSpikeCounts = spikeCounts(msStart:msStop);
                x(trial,:) = trialSpikeCounts;
    
            end
    
        end
%         baselineX.(windows{w}) = x;
%         allX.(windows{w}) = x;
        cvX(n) = std(networkX) / mean(networkX);
    %     nActiveTrials(n) = length(find(find(~all(x == 0,2))));
        
    %     % Save variables
    %     if strcmp(types{n}, "BASELINE")
    % 
    %         % Get average baseline spiking rate
    %         activeChannelSpikeTimes = mergedSpikeTimes(~cellfun('isempty',mergedSpikeTimes));
    %         spikesN = length(cell2mat(activeChannelSpikeTimes));
    %         baselineSpikeRate = spikesN / (durationS*1000); % get network spike count per ms
    % 
    %         cd(spikeDir)
    %         try
    %             save(strcat(samples{n},'.mat'),'mergedSpikeTimes','spikeMatrix','baselineSpikeRate','x')
    %         catch
    %             save(strcat(samples{n},'.mat'),'mergedSpikeTimes','baselineSpikeRate','x')
    %         end
    %         cd(homeDir)
    % 
    %         x = [x;x];
    %     else
%             cd(spikeDir)
%             try
%                 save(strcat(samples{n},'.mat'),'allX','-append');
%             catch
%                 save(strcat(samples{n},'_spikes.mat'),'allX','-append');
%             end

    %         try
    %             save(strcat(samples{n},'.mat'),'mergedSpikeTimes','spikeMatrix','x');
    %         catch
    %             save(strcat(samples{n},'.mat'),'mergedSpikeTimes','x');
    %         end
            cd(homeDir)
        end
    end
    % 
    %     cd(outputDir)
    %     writematrix(x, strcat(samples{n},".csv"))
    %     cd(homeDir)
        
        clear x  mergedSpikeTimes mergedSpikeTimesFrames spikeMatrix

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
                'XLimits',seconds([lostTime/fs,0.1]),'XLabelText', "Time (s)", "YLabelText", "Trial",...
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

%% Plot x over trials

% for w=1:length(windows)
    window='window0_180';
    stimRecs = samples(~strcmp(types,'BASELINE'));
    stimRecsN = length(stimRecs);
    trialX = zeros(trialsN,stimRecN);

    baselineRecs = samples(strcmp(types,'BASELINE'));
    baselineRecsN = length(baselineRecs);
    trialX = zeros(trialsN,baselineRecsN);

    for n = 1:baselineRecsN
        
        disp(baselineRecs{n})
        
        cd(spikeDir)
        
        % Get activity for pattern 0
        cd(spikeDir)
        try
            allX = load(strcat(baselineRecs{n},'_spikes.mat')).baselineX;
        catch
            allX = load(strcat(baselineRecs{n},'.mat')).baselineX;
        end
        cd(homeDir)
    
%         try
%             stimNode = pattern0StimElec{n};
%         catch
%             stimNode = pattern0StimElec(n);
%         end
    
%         if isa(stimNode1,'char')
%             stimNode1 = str2num(stimNode1);
%         end

        x = baselineX.(window); % s{w}
        trialX(:,n) = mean(x,2);
    
    end

    colData = reshape(trialX, [], 1);

% end
%% For side-by-side plotting, load in metadata in stimulation pair format

% Metadata -- see example spreadsheet
xlSheet = 'paired_crit'; % paired
xlRange = 'A2:M15';
[num,txt,~] = xlsread(metadataSpreadsheet,xlSheet,xlRange);
pairName = txt(:,2);
baselineRec = txt(:,3);
pattern0Rec = txt(:,4);
pattern1Rec = txt(:,5);

pattern0StimElec = num(:,6); % for a single stimulation electrode per pattern
pattern1StimElec = num(:,8); % for a single stimulation electrode per pattern

% pattern0StimElec = txt(:,10); % for a single stimulation electrode per pattern
% pattern1StimElec = txt(:,12); % for a single stimulation electrode per pattern

%% Plot PSTH

% PSTH parameters
binsize = 1; % bin size for plotting histogram in miliseconds


for n = 1:length(pairName)
    
    disp(pairName{n})
    
    cd(spikeDir)
    
    % Get recording-averaged baseline firing rate
    load(strcat(baselineRec{n},'.mat'),'baselineSpikeRate')
    
    % Get spike times for pattern 0
    cd(spikeDir)
    try
        spikeTimes = load(strcat(pattern0Rec{n},'_spikes.mat')).mergedSpikeTimes;
    catch
        spikeTimes = load(strcat(pattern0Rec{n},'.mat')).mergedSpikeTimes;
    end

    cd(homeDir)
   
    % Get stimulation electrode(s)
    try
        stimNode = pattern0StimElec{n};
    catch
        stimNode = pattern0StimElec(n);
    end
    if isa(stimNode,'char')
        stimNode = str2num(stimNode);
    end
        
    spikeCounts = cell2mat(cellfun(@(a) length (a), spikeTimes, 'UniformOutput', false));
    excludeUnits = find((spikeCounts > framesN/fs*100) |...
    (spikeCounts < 10)); % exclude units with mean firing rate >100 Hz across recording
    % and units which record fewer than 10 spikes across the whole
    % recording -- assumed to be noise
    spikeTimes([excludeUnits,stimNode]) = [];
    spikeTimes0 = cell2mat(spikeTimes)*fs/1000; % convert from ms to frames
    clear spikeCounts spikeTimes stimNode

    % Get spike times for pattern 1
    cd(spikeDir)
    try
        spikeTimes = load(strcat(pattern1Rec{n},'_spikes.mat')).mergedSpikeTimes;
    catch
        spikeTimes = load(strcat(pattern1Rec{n},'.mat')).mergedSpikeTimes;
    end
    cd(homeDir)

%     stimNode = stimElecIdx(2,:);
    try
        stimNode = pattern1StimElec{n};
    catch
        stimNode = pattern1StimElec(n);
    end

    if isa(stimNode,'char')
        stimNode = str2num(stimNode);
    end

    spikeCounts = cell2mat(cellfun(@(a) length (a), spikeTimes, 'UniformOutput', false));
    excludeUnits = find((spikeCounts > framesN/fs*100) |...
    (spikeCounts < 10)); % exclude units with mean firing rate >100 Hz across recording
    % and units which record fewer than 10 spikes across the whole recording
    spikeTimes([excludeUnits,stimNode]) = [];
    spikeTimes1 = cell2mat(spikeTimes)*fs/1000; % convert from ms to frames

    t = plotPSTH(spikeTimes0, spikeTimes1, binsize, fs, trialsN, ISI, pairName{n}, baselineSpikeRate);
    cd .\network_response\PTSH
    saveas(t,pairName{n},'svg')
    cd(homeDir)
    close all

end

%% Plot heatmaps

distances = zeros(length(samples),1);
for n = 1:length(pairName)
    
    disp(pairName{n})
    
    cd(spikeDir)
    
    % Get recording-averaged baseline firing rate
    load(strcat(baselineRec{n},'.mat'),'baselineSpikeRate')
    
    % Get activity for pattern 0
    cd(spikeDir)
    try
        x1 = load(strcat(pattern0Rec{n},'_spikes.mat')).x;
    catch
        x1 = load(strcat(pattern0Rec{n},'.mat')).x;
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

    % Get activity for pattern 1
    cd(spikeDir)
    try
        x2 = load(strcat(pattern1Rec{n},'_spikes.mat')).x;
    catch
        x2 = load(strcat(pattern1Rec{n},'.mat')).x;
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
    

    dist = vecnorm(x1'-x2');
    distances(n) = mean(dist);

%     dat1 = mean(x1)/baselineSpikeRate; % fold change in firing rate relative to baseline
%     dat2 = mean(x2)/baselineSpikeRate; % fold change in firing rate relative to baseline
%     
%     f = plotMEAHeatMap(dat1,dat2,stimNode1,stimNode2,pairName{n},'foldChange');
%     cd .\network_response\heatmaps
%     saveas(f,pairName{n},'png')
%     cd(homeDir)
%     close all

end
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