function [filteredDat] = runSALPA(dat, excludeElecs, lostTime, trialLength, trialOnT, trialOffT)

trialN = length(trialOnT);
channelsN = size(dat,2);
filteredDat = zeros(trialLength,channelsN,trialN);

for trial = 1:trialN

    start = trialOnT(trial);
    stop =  trialOffT(trial);
    trialWin = dat(start:stop,:);

    for channel = 1:channelsN
        
        if ~ismember(channel, excludeElecs)

            % Get time before depegging
            channelDat = meanCentre(trialWin(:,channel));
            thr = 6*std(channelDat); % threshold for finding time of saturation
            [~,pegDelay] = findpeaks(channelDat,"NPeaks",1, "MinPeakHeight",thr);

            if isempty(pegDelay)
                warning("Trial %d, channel %d: Stimulation saturation not detected. Check artifact removal.",...
                    trial, channel)
                pegDelay = 5;
            end

            filteredTrace = salpa(channelDat,'tau',lostTime,'t_blankdepeg',pegDelay);
            filteredDat(:,channel,trial) = filteredTrace(lostTime+1:end);

        else
            filteredDat(:,channel,trial) = channelDat(lostTime+1:end);
        end
        
        clear filteredTrace
    end
    
    clear trialWin

end

filteredDat = reshape(permute(filteredDat, [1 3 2]), [], channelsN, 1);
clear dat

end