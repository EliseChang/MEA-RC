function [minAmp, minSlope, maxSlope] = getSpikeWaveform(spikeWaveforms, method)
% TODO: exclude grounded electrodes
% Create arrays to store values for each channel
chMinAmp = zeros(1,60);
chMinSlope = zeros(1,60);
chMaxSlope = zeros(1,60);

for ch = 1:60
    spikes = spikeWaveforms{1, ch}.(method);
    [~, allSpikesMinIdx] = min(spikes, [], 2);
    inclSpikes = spikes(allSpikesMinIdx == 25, :); % the 'clean' spike waveforms are centred on the midpoint of the cutout window, 
                                                    %  so just consider these for now
    if isempty(inclSpikes)
        continue
    end
    
    [spikeMinAmp,spikeMinIdx] = min(inclSpikes, [], 2);
    
    preSpikes = inclSpikes(:,1:25); % find the max. value and index before the spike minimum
    [spikeMaxAmp,spikeMaxIdx] = max(preSpikes, [], 2);

    amps = spikeMaxAmp - spikeMinAmp;
    chMinAmp(1, ch) = prctile(amps, 10); % set the minimum amplitude as the 10th percentile
    slopeDur = (spikeMinIdx - spikeMaxIdx) * 40; % time between spike max. and min. in microseconds
    slopes = amps ./ slopeDur; % slope in microvolts / microseconds
    chMinSlope(1, ch) = prctile(slopes, 10);
    chMaxSlope(1, ch) = prctile(slopes, 90);

    clear spikes allSpikesMinIdx inclSpikes spikeMinAmp spikeMinIdx preSpikes spikeMaxAmp spikeMaxIdx amps slopeDur slopes

end

minAmp = min(chMinAmp(chMinAmp ~= 0));
minSlope = min(chMinSlope(chMinSlope ~= 0));
maxSlope = max(chMaxSlope(chMaxSlope ~= 0));

clear chMinAmp chMinSlope chMaxSlope

end