function f = plotSpikeDetectionChecksStim(stimDat, window, preSpikeTimes,postSpikeTimes, stimElecs, outputElecs,... 
    times, fs, figPos)
%
% Plot voltage traces pre- and post-stimulus artifact removal
% within specific windows
% 
% Parameters 
% ----------
% spikeTimes (both pre and post) : cell array of structures
%    1 X N cell array, where N is the number of electrodes / units 
%    each item in a cell should be a 1 x 1 struct
%    which contain fields corresponding to the spike detection methods used,
%    eg. bior1p5 [numSpikes x 1 double] and timing data given in
%    SECONDS
% stimDat : (struct)
%         structure output from stimPrePipelineProcessing containing
%         cut-out data
% window : (str)
%         name of post-stimulation cutout window e.g. '0-100ms'
% stimElecs : (array)
%           indices of stimulation electrodes (will be plotted in first
%           row)
% outputElecs : (array)
%           indices of output electrodes (will be plotted in second row)
% times: (1 x 2 array)
%           start and end frame number for post-stim window
% figFolder : str
%     specify where the full path to the folder to save the spike detection
%     check plots
% Returns 
% -------

%}

%% Define plotting colours

spikeMethodColors = ...
    [  0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ... 
    0.4660    0.6740    0.1880; ... 
    0.3010    0.7450    0.9330; ... 
    0.6350    0.0780    0.1840];

%% Get voltage data pre- and post-SALPA

preFieldName = strcat('preSALPA',window);
postFieldName = strcat('postSALPA',window);
preDat = stimDat.(preFieldName);
postDat = stimDat.(postFieldName);
 
%% Filter voltage data -- should do this before?

% filtered_data = zeros(size(dat));
% num_chan = size(dat,2);
% 
% for ch = 1:num_chan
%     lowpass = Params.filterLowPass;
%     highpass = Params.filterHighPass;
%     wn = [lowpass highpass] / (fs / 2); % seems like the fs here comes from workspace???
%     filterOrder = 3;
%     [b, a] = butter(filterOrder, wn);
%     filtered_data(:,ch) = filtfilt(b, a, dat(:,ch));
% end

%% Get spike detection methods

preMethods = fieldnames(preSpikeTimes{1});
postMethods = fieldnames(postSpikeTimes{1});

% Plot spikes for only methods used for both data
methods = intersect(preMethods,postMethods);
methodsN = numel(methods);
methodsList = strrep(methods, 'p','.'); % for legend

%% Plot voltage traces with spike times

stimElecsN = numel(stimElecs);
outputElecsN = numel(outputElecs);
colN = max(stimElecsN,outputElecsN);
start = times(1);
stop = times(2) - 1;

f = figure('Position',figPos);

% Stimulation electrodes
for p = 1:stimElecsN

    subplot(2,colN,p)
    channel = stimElecs(p);
    preTrace = preDat(start:stop, channel);
    plot(preTrace, 'k-')
    hold on
    postTrace = postDat(start:stop, channel);
    plot(postTrace, 'r--') % , 'Color', [0.5020, 0.5020, 0.5020]) % grey
    
    % add markers to show spike times 
    maxVal = max([preTrace,postTrace],[],'all');
    minVal = min([preTrace,postTrace],[],'all');

    for m = 1:length(methods)
        method = methods{m};
        preSpikeTrain = preSpikeTimes{channel}.(method);
        preSpikeTrainFrames = preSpikeTrain*fs;
        postSpikeTrain = postSpikeTimes{channel}.(method);
        postSpikeTrainFrames = postSpikeTrain*fs;
    end
   
    if ~isempty(preSpikeTrain)
        scatter(preSpikeTrainFrames, repmat(maxVal,numel(preSpikeTrain),1), ...
           15, 'v', 'filled', ... % 15 is marker size
            'markerfacecolor',spikeMethodColors(m, :), ... 
            'markeredgecolor', spikeMethodColors(m, :), 'linewidth',0.1);
        maxVal = 0.9*maxVal; % ensures that spike time markers are not plotted over each other
    end

    if ~isempty(postSpikeTrain)
        scatter(postSpikeTrainFrames, repmat(minVal,numel(postSpikeTrain),1), ...
           15, '^', 'filled', ... % 15 is marker size
            'markerfacecolor',spikeMethodColors(m, :), ... 
            'markeredgecolor', spikeMethodColors(m, :), 'linewidth',0.1);
        minVal = 0.9*minVal;
    end
    
    % format
    hold off
    aesthetics
    xlim([start stop])
    ylabel('Amplitude (\muV)');
    xlabel("Time post-stimulus (ms)")
    xticklabels(xticks/fs*1e3) % convert from frames to ms
    title(["Electrode ", num2str(channel)])
    set(gca,'TickDir','out');

end
% TODO: add title for first row

% Output electrodes
for p = 1:outputElecsN
    
    channel = outputElecs(p);
    subplot(2,colN,p+colN)
    preTrace = preDat(start:stop, channel);
    plot(preTrace, 'k-')
    hold on
    postTrace = postDat(start:stop, channel);
    plot(postTrace, 'r--') % , 'Color', [0.5020, 0.5020, 0.5020]) % grey
    
    % add markers to show spike times 
    maxVal = max([preTrace,postTrace],[],'all');
    minVal = min([preTrace,postTrace],[],'all');
    
    for m = 1:length(methods)
        method = methods{m};
        preSpikeTrain = preSpikeTimes{channel}.(method);
        preSpikeTrainFrames = preSpikeTrain*fs;
        postSpikeTrain = postSpikeTimes{channel}.(method);
        postSpikeTrainFrames = postSpikeTrain*fs;
        if ~isempty(preSpikeTrain)
            scatter(preSpikeTrainFrames, repmat(maxVal,numel(preSpikeTrain),1), ...
               15, 'v', 'filled', ... % 15 is marker size
                'markerfacecolor',spikeMethodColors(m, :), ... 
                'markeredgecolor', spikeMethodColors(m, :), 'linewidth',0.1);
            maxVal = 0.9*maxVal;
        end
        if ~isempty(postSpikeTrain)
            scatter(postSpikeTrainFrames, repmat(minVal,numel(postSpikeTrain),1), ...
               15, '^', 'filled', ... % 15 is marker size
                'markerfacecolor',spikeMethodColors(m, :), ... 
                'markeredgecolor', spikeMethodColors(m, :), 'linewidth',0.1);
            minVal = 0.9*minVal;
        end
    end

    % format
    hold off
    aesthetics
    xlim([start stop])
    ylabel('Amplitude (\muV)');
    xlabel("Time post-stimulus (ms)")
    xticklabels(xticks/fs*1e3) % convert from frames to ms
    title(["Electrode ", num2str(channel)])
    set(gca,'TickDir','out');

end

% TODO: add title for second row

% add legend
hL = legend('pre-SALPA', 'postSALPA', methodsList{:});
% newPosition = [0.6 0.12 0.1 0.1];
% set(hL,'Position', newPosition,'Units', newUnits,'Box','off');
% title({strcat(regexprep(Info.FN{1},'_','','emptymatch')),' '});

