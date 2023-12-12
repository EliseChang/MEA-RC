% Process .mat voltage data from MEA stimulation experiments ready for spike detection and analysis with
% MEA-NAP (https://github.com/SAND-Lab/MEA-NAP)

    % - trim recordings to specified number of trials
    % - apply SALPA to each stimulation trial (Wagenaar & Potter, 2002) for stimulus artifact removal
    % - recompile stimulation trials into one recording (note that some
        % sections of the recording may have been removed as part of stimulus
        % artifact suppression)

% Created: Elise Chang, May 2023

%% Set parameters and import metadata

homeDir = "D:\MATLAB\MEA-RC";
cd(homeDir)
metadataSpreadsheet = "mecp2RecordingsListNew.xlsx";
dataDir = 'D:\MATLAB\MEA-NAP\organoids\Nov2023DIV150Stim'; % "E:\231102_Mecp2_organoids"

% Metadata -- see example spreadsheet
xlSheet = 'Stim';
xlRange = 'A2:M11';
spreadsheetDir = "D:\MATLAB\MEA-NAP\metadata";
[num,txt,~] = xlsread(fullfile(spreadsheetDir,metadataSpreadsheet),xlSheet,xlRange);
samples = txt(:,1);

% Paths
addpath(dataDir)
addpath(genpath("Functions"))
figDir = fullfile('figs', "stimArtifactRemoval");
addpath(figDir)

% Parameters for segmenting recording into trials
patternsN = 2;
channelsN = 60;
trialsN = 300; % TTX recordings = 12
fs = 25e3;
ISI = 1*fs; % inter-stimulus interval
framesN = trialsN*ISI;
lostTime = 3e-3*fs; % for SALPA algorithm -- but note this is still included in trial
trialStartFrame = 160e-3*fs;
trialLength = 40e-3*fs; % in frames
window = '160_200ms';

% idealised stimulation times -- will be adjusted
stimLength = 5; % number of frames
trialOnT = 0:ISI+stimLength:...
    trialsN*(ISI+stimLength);

% % Set stimulation electrode indices
patternAStimID = txt(:,4);
patternBStimID = txt(:,5); % idx
groundElecID = txt(:,6);
circShift = num(:,7);
startTrial = num(:,8);
endTrial = num(:,9);
% nTrials = num(:,10);
preStimFile = txt(:,12);
postStimFile = txt(:,13);

% NOTE: currently 0 = A and 1 = B
stimProt = double(~readmatrix(fullfile('stimuli', "MeCP2OrgStimProt1.csv")));
stimProt(1) = []; % remove first element due to zero-indexing

figPos = [1 49 1920 955];

%% Find stimulation trial onset times

load("channels.mat")
fs = 25e3;

for n = 1:length(samples)
    
    cd(dataDir)
    load(strcat(samples{n},'.mat'))
    cd(homeDir)
    disp(samples{n})
    
%     % Actually sometimes useful not to truncate
%     if length(dat) > framesN
%         dat(framesN+1:end,:) = [];
%     end

    % mean-centre voltage traces at 0 for easier interpretation of voltage deflections
    centredDat = bsxfun(@minus, dat, sum(dat)/size(dat, 1));
    
    % get stimulation electrodes used for each pattern and ground electrodes
    patternAStimIdx = getElectrodeIdx(str2num(patternAStimID{n})); %#ok<ST2NM> 
    patternBStimIdx = getElectrodeIdx(str2num(patternBStimID{n})); %#ok<ST2NM>
    stimElecs = [patternAStimIdx,patternBStimIdx];
    groundElecIdx = getElectrodeIdx(str2num(groundElecID{n})); %#ok<ST2NM> 
    excludeElecs = [patternAStimIdx;patternBStimIdx;groundElecIdx];

%     % check trial sequence
%     for e = 1:length(stimElecs)
%         plot(dat(:,e))
%         TODO: save these
%         shiftTrials = input("Enter the actual index of the final trial: ")
%         startTrialIdx = input("Enter start trial: ")
%         endTrialIdx = input("Enter end trial: ")
%         close all
%     end

    % Adjust stim prot as needed
    shiftTrials = circShift(n);

    if ~isnan(shiftTrials)
        patternSeq = circshift(stimProt, trialsN - shiftTrials);
%         patternSeq = patternSeq(startTrialIdx:endTrialIdx);
    else
        patternSeq = stimProt; %(startTrialIdx:endTrialIdx);
    end
    firstPattern = patternSeq(1) + 1; % zero-indexing
    checkStimElecs = stimElecs(:,firstPattern);
    stimElecsN = length(checkStimElecs);

    % visualising voltage traces at stimulation electrodes for first trial
    f = figure("Position",figPos);
    t = tiledlayout(1,stimElecsN);
    artifactT = zeros(stimElecsN,1); % store start frames for each electrode
            
    % plot and find stimulation electrode which registers artifact
    % earliest -- we will lock stimulation trials to this electrode
    for i = 1:length(checkStimElecs)
        
        e = checkStimElecs(i);
        nexttile
        firstTrial = centredDat(1:ISI,e);
        plot(firstTrial)
        [~,maxIdx] = max(firstTrial); % find artifact peak
        [~,minIdx] = min(firstTrial); % or sometimes trough
        artifactT(i) = maxIdx; % but we assume peak, can be corrected if necessary
        xline(maxIdx, '-red', num2str(maxIdx)) % plot and label time frame of peak
        xline(minIdx, '-blue', num2str(minIdx))
        aesthetics
        title(["Pattern ", firstPattern, " Electrode ", e])
        clear firstTrial maxIdx

    end    
    
    [startFrame,earlyElecIdx] = min(artifactT); % find electrode with earliest artifact
    lockToElec = checkStimElecs(earlyElecIdx); % store electrode index for each pattern
    override = input("Check artifact polarity. Specify start frame and electrode index to lock to as vector if required, otherwise no input required. ");
    if ~isempty(override)
        startFrame = override(1);
        lockToElec = override(2);
    end    

    % save fig
    xlabel(t,"Sampling frame")
    ylabel(t, "Voltage")
    aesthetics
    cd figs\stimArtifactRemoval\firstArtifact
    saveas(t, samples{n}, 'png')
    cd(homeDir)
    close all
    
%     % For selected electrode in each pattern, find artifact times
%     stimT = zeros(trialsN/patternsN,patternsN); % store actual trial times for each pattern
%     trialGuides = trialOnT + startFrames(p); % idealised trial times
%     for p = 1:patternsN
%         % get trials
%         mask = patternSeq == p-1;
%         fprintf('Pattern %d: Electrode %d \n',p,lockToElecs(p));
%         chDat = centredDat(:,lockToElecs(p));
%         adjust = 'y';
%         while ~isempty(adjust)
%             plot(chDat)
%             xline(trialGuides(mask))
%             thr = input("Enter minimum peak amplitude for stimulus artifacts: ");
%             [peaks,locs] = findpeaks(chDat,"NPeaks",trialsN/patternsN,"MinPeakHeight",thr,"MinPeakDistance",ISI*0.75);
%             hold on
%             scatter(locs,peaks,'v')
%             hold off
%             adjust = input("Check stimulus artifact identification. Enter 1 to adjust threshold, otherwise no input required. ");
%             close all
%         end
%         clear chDat
%         stimT(:,p) = locs;
%     end
%     stimTimes = sort(reshape(stimT,[trialsN,1]));
%     fprintf("ISIs: %d", diff(stimTimes))
%     pause

    stimTimes = trialOnT + startFrame;
    stimTimes = stimTimes(startTrialIdx:endTrialIdx);
    cd(dataDir)
    save(strcat(samples{n},'.mat'), 'stimTimes', 'centredDat','-append')
    cd(homeDir)
    clear dat stimTimes

end

%% Get stimulation cut-outs across all electrodes and trials and run SALPA

for n = 1:length(samples)
    
    cd(dataDir)
    centredDat = load(strcat(samples{n},'.mat'), 'centredDat').centredDat;
    stimTimes = load(strcat(samples{n},'.mat'), 'stimTimes').stimTimes;
    try
        stimDat = load(strcat(samples{n},'.mat'), 'stimDat').stimDat;
    catch
        stimDat = struct;
    end
    cd(homeDir)
    disp(samples{n})
    
    if stimTimes(end) + trialStartFrame + trialLength > framesN
        disp("Final trial exceeds recording length. Remove trials.")
        continue
    end

    mask = zeros(length(centredDat),1);
    for s = 1:length(stimTimes)
        trialStart = stimTimes(s) + trialStartFrame;
        trialEnd = trialStart + trialLength - 1;
        mask(trialStart:trialEnd) = 1;
    end
    cutoutDat = centredDat(logical(mask),:);
    fieldName = strcat('preSALPA',window);
    stimDat.(fieldName) = cutoutDat;
    clear fieldName

    % get stimulation electrodes used for each pattern and ground electrodes
    patternAStimIdx = getElectrodeIdx(str2num(patternAStimID{n})); %#ok<ST2NM> 
    patternBStimIdx = getElectrodeIdx(str2num(patternBStimID{n})); %#ok<ST2NM>
    stimElecs = [patternAStimIdx,patternBStimIdx];
    groundElecIdx = getElectrodeIdx(str2num(groundElecID{n})); %#ok<ST2NM> 
    excludeElecs = [patternAStimIdx;patternBStimIdx;groundElecIdx];

    % Run SALPA
    % set new trial times after 'cutting out' data
    newFramesN = length(cutoutDat);
    newTrialTOn = 1:trialLength:newFramesN;
    newTrialTOff = [newTrialTOn(2:end)-1, newFramesN];
    tic
    filteredDat = runSALPA(cutoutDat, excludeElecs, lostTime, trialLength, newTrialTOn, newTrialTOff);
    toc
    fieldName = strcat('postSALPA',window);
    stimDat.(fieldName) = filteredDat;

    cd(dataDir)
    save(strcat(samples{n},'.mat'), 'stimDat', '-append')
    cd(homeDir)

    clear dat filteredDat stimTimes stimDat

end
