% Process .mat voltage data from MEA stimulation experiments ready for spike detection and analysis with
% MEA-NAP (https://github.com/SAND-Lab/MEA-NAP)

    % - trim recordings to specified number of trials
    % - apply SALPA to each stimulation trial (Wagenaar & Potter, 2002) for stimulus artifact removal
    % - recompile stimulation trials into one recording (note that some
        % sections of the recording may have been removed as part of stimulus
        % artifact suppression)

% Created: Elise Chang, May 2023

%% Set parameters and import metadata

homeDir = ("D:\MATLAB\MPhil_scripts\MEA-NAP");
cd(homeDir)
metadataSpreadsheet = "Mona_MEC_MEA-NAP.xlsx";
dataDir = fullfile(homeDir, 'voltage');

% Metadata -- see example spreadsheet
xlSheet = 'Sheet1';
xlRange = '';
[num,txt,~] = xlsread(metadataSpreadsheet,xlSheet,xlRange);
samples = txt(:,1);
patterns = txt(:,3);

% Paths
addpath(dataDir)
addpath Functions\

% Parameters for segmenting recording into trials
trialsN = 60; % TTX recordings = 12
channelsN = 60;
fs = 25e3;
ISI = 5*fs; % inter-stimulus interval

stimLength = 5; % number of frames
lostTime = 3e-3*fs; % for SALPA algorithm
trialLength = ISI;
trialOnT = 1:ISI+stimLength:...
    trialsN*(ISI+stimLength); % align trial start times
trialOffT = trialOnT + trialLength - 1;
framesN = trialOffT(end);

% % Set stimulation electrode indices
% stimElecIdx = getElectrodeIdx([[12, 21, 22]; [17, 27, 28]]); % Pattern 1 on first row, etc.
groundElecIdx = []; %getElectrodeIdx([15, 16, 78, 87]);

%% Run SALPA

for n = 1:length(samples)

    if patterns{n} ~= "BASELINE" % only stimulation, not baseline recordings
    
        cd(dataDir)
        load(strcat(samples{n},'.mat'))
        cd(homeDir)

%         if length(dat) < framesN
%             dat(end:framesN,:) = 0;
%         end
        
        stimElecs = patterns(n); % if input node positions are not set

        excludeElecs = [stimElecs,groundElecIdx];
        disp("Running SALPA")
        filteredDat = runSALPA(dat, excludeElecs, lostTime, trialLength, trialOnT, trialOffT);
        cd(dataDir)
        save(strcat(samples{n},'.mat'),"dat","filteredDat")
        cd(homeDir)

    end

    clear dat filteredDat

end