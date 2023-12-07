function [filteredDat] = runSALPA(dat, excludeElecs, lostTime, trialLength, trialOnT, trialOffT)

trialN = length(trialOnT);
channelsN = size(dat,2);
filteredDat = zeros(trialLength,channelsN,trialN);

% For plotting voltage signals pre- and post-SALPA -- not currently used as
% separate function for plotting after spike detection

% % plotTrialsN = 5;
% plotWin = 0.1*fs; % Plot 100 ms
% plottingFiltDat = zeros(plotWin,plotTrialsN,length(outputElecs));
% plottingRawDat = zeros(plotWin,plotTrialsN,length(outputElecs));

% Split into trials, find and remove stimulus artifact from each trial
for trial = 1:trialN

    start = trialOnT(trial);
    stop =  trialOffT(trial);
    trialWin = dat(start:stop,:);

    for channel = 1:channelsN
        channelDat = trialWin(:,channel);
        if ~ismember(channel, excludeElecs)
            
            filteredTrace = salpa(channelDat,'tau',lostTime,'t_blankdepeg',0); % ,'t_blankdepeg',pegDelay
            filteredTrace(isnan(filteredTrace)) = 0;
            filteredDat(:,channel,trial) = filteredTrace;
                        
%             if trial <= plotTrialsN && ismember(channel, outputElecs)
%                 outputElecN = find(outputElecs == channel);
%                 plottingFiltDat(:,trial,outputElecN) = filteredTrace(1:plotWin);
%                 plottingRawDat(:,trial,outputElecN) = channelDat(1:plotWin);
%             end

        else
%             if ismember(channel, outputElecs)
%                 warning("Channel %d is a stimulation/grounded electrode. Voltage traces will not be plotted.", channel)
%             end
            filteredDat(:,channel,trial) = channelDat; % don't filter data from excluded electrodes
        end

        clear filteredTrace channelDat

    end
    clear trialWin
end

filteredDat = reshape(permute(filteredDat, [1 3 2]), [], channelsN, 1);

% % Plot
% t = tiledlayout(2,3); % Position=[1, 49, 1920, 955]
% for p = 1:length(outputElecs)
%     z1 = plottingFiltDat(:,:,p);
%     z2 = plottingRawDat(:,:,p);
%     nexttile
%     stackedPlot(z1, 4, 1, 4, "Color", "red"); % 4 = style, 1 = spacing, 4 = labels
%     hold on
%     stackedPlot(z2, 4, 1, 4, "Color", "black");
%     hold off
%     title(strcat("Electrode ", num2str(outputElecs(p))))
%     clear z1 z2
% end
% 
% % Create one pair of legend labels
% legendLabels = strings(plotTrialsN*2,1);
% legendLabels(1) = "Filtered";
% legendLabels(plotTrialsN+1) = "Raw";
% l = legend(legendLabels);
% 
% l.Layout.Tile = 'east';
% saveas(t,fullfile(figDir, strcat(fileName,".png")))
% close all
% clear plottingFiltDat plottingRawDat