%% Avalanche detection

%% Metadata and paths

addpath NCCToolboxV1

homeDir = ("D:\MATLAB\MPhil_scripts\MEA-RC");
cd(homeDir)
metadataSpreadsheet = "D:\MATLAB\MPhil_scripts\MEA-NAP\Mecp2_2022_dataset_MEA-NAP.xlsx";
spikeDir = "D:\MATLAB\MPhil_scripts\MEA-NAP\OutputData09Nov2022\ExperimentMatFiles"; % fullfile(homeDir,"spikes"); % 'E:\spikes';
outputDir = fullfile(homeDir, "spikes/"); % for savig asdf2 files
stateVar = 'spike_rates'; % spike_rates
% window = ''; % '3-20';

% Metadata
xlSheet = 'fully_connected_mature';
xlRange = 'A2:C52';
[num,txt,~] = xlsread(metadataSpreadsheet,xlSheet,xlRange);
samples = txt(:,1);
ages = num(:,1);
genotypes = txt(:,3);
sampleN = length(samples);
Params.date = '09Nov2022'; % if using spikeTimes from a pipeline output folder
% types = txt(:,4);
% stimNodes = txt(:,6); % for multiple stimulation electrodes per pattern
% stimNodes = num(:,5); % for a single stimulation electrode per pattern

% Paths
addpath(spikeDir)
addpath network_response\avalancheAnalysis
addpath Functions\

framesN = 25e3*60*10; % length of whole recording in frames

%% Hyperparameters

binSizes = [1]; % ms
channelsN = 60;
% % Parameters for segmenting recording into trials
% trialsN = 60; % TTX recordings = 12
% channelsN = 60;
fs = 25e3;
% ISI = 5*fs; % inter-stimulus interval
% 
% stimLength = 5; % number of frames
% lostTime = 3e-3*fs; % time removed due to stimulus artifacts and early response component
% trialLength = 17e-3*fs; % CHECK THIS IS CORRECT

% % For samples which have already been segmented into trials of the correct
% % length
% trialOnT = 1:trialLength:framesN; % is actually cut off at the end
% trialOffT = trialLength:trialLength:framesN;
% stimOnT = trialOnT/fs;
% framesN = trialsN*trialLength;

% % For samples which have NOT been segmented into trials
% trialOnT = lostTime:ISI+stimLength:...
%     lostTime+(trialsN-1)*(ISI+stimLength);
% trialOffT = trialOnT + trialLength;
% stimOnT = [1,ISI:ISI:(trialsN-1)*ISI]/fs;
% framesN = trialsN*ISI;
% 
% durationS = framesN/fs;

%% Get avalanche properties
branchingRatio = zeros(sampleN,1);
for n = 1:sampleN

    disp(samples{n})

    cd(spikeDir)
    try
        mergedSpikeTimes = load(strcat(samples{n},'.mat'),'mergedSpikeTimes').mergedSpikeTimes;
    catch
        spikeTimes = load(strcat(samples{n},'_',Params.date,'.mat')).spikeTimes;    
        mergedSpikeTimes = cell(1,channelsN);
        for ch = 1:channelsN
            [mergedSpikes,~, ~] = mergeSpikes(spikeTimes{ch}, 'all');
            mergedSpikeTimes{1,ch} = mergedSpikes;
        end
    end

    cd(homeDir)
    % Note that Mecp2 2022 spikeTimes are in seconds
    mergedSpikeTimesFrames = cellfun(@(a) a*fs, mergedSpikeTimes, 'UniformOutput', false); % convert from seconds to frames
    spikeMatrix = spikeTimesToSpikeMatrix(mergedSpikeTimesFrames, framesN);

    expSys = "Mecp2_2022Baseline";
    dataType = "spikes";
    dataID = samples{n};
    binSize = 1/25; % unit of time represented by each bin in ms
    asdf2 = rastertoasdf2(spikeMatrix', binSize, expSys, dataType, dataID);
    % TODO: add physical electrode coordinates?
    asdf2Binned = rebin(asdf2, 1);
    rasterBinned = asdf2toraster(asdf2Binned);
    
    cd(outputDir)
    save(samples{n},'asdf2Binned','-mat')
    cd(homeDir)

%     %% Plot binary raster plot at 1 ms resolution
%     h = imagesc(rasterBinned);
%     savefig(fullfile(h,figDir,strcat(dataID,"_raster_1ms")),'compact')
%     close all

% end

% %     %% Plot IEI distribution
% % 
% %     allSpikeTimes = sort(cell2mat(asdf2Binned.raster'));
% %     intervals = diff(allSpikeTimes);
% % %     intervals(intervals > 100) = [];
% %     avalancheData.meanIEI = mean(intervals);
% %     h = histogram(intervals);
% % %     xlim([0,100])
% %     hold on
% %     xline(mean(intervals))
% %     hold off
% %     savefig(fullfile(h,figDir,strcat(dataID,"_IEI_distribution")),'compact')
% %     close all

    % Branching ratio estimate
    [br,~,~] = brestimate(asdf2Binned);
    branchingRatio(n) = br;

%     Avalanche = avprops(asdf2);
%     maxSize = max(Avalanche.size);
%     for s = 1:maxSize 
%         if sum(Avalanche.size == s) == 1
%             idx = find(Avalanche.size == s);
%             Avalanche.size(idx) = [];
%             Avalanche.duration(idx) = [];
%             Avalanche.shape(idx) = [];
%         end
%     end
% 
%     [tau, xmin, xmax, sigmaTau, p, pCrit, ks, PlotData] = avpropvals(Avalanche.duration, 'duration','plot');
%     [tau, xmin, xmax, sigmaTau, p, pCrit, ks, PlotData] = avpropvals(Avalanche.size, 'size','plot');
 
end