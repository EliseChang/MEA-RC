function [trialAssignment] = assignSpikeTimeTrials(spikeTimes,stims)
% Assigns spike times to trial number

% Output
    % trials: categorical vector trial assignment of each spike in spikeTimes

trialAssignment = zeros(length(spikeTimes),1);

for trial = 2:length(stims)
    first_spike_idx = find((spikeTimes>=stims(trial-1)),1,'first');
    last_spike_idx = find((spikeTimes<stims(trial)),1,'last');
    trialAssignment(first_spike_idx:last_spike_idx,:) = trial-1;
end

end