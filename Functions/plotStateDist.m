function f = plotStateDist(spikeMatrix0, spikeMatrix1, trialsN, trialLength, dsF, plotWin)
    
    % TODO: check spike matrices have consistent dimensions
    % TODO: option to pass subset of electrodes
    elecDist = sqrt(sum((spikeMatrix0 - spikeMatrix1).^2,2));
    trialDist = reshape(elecDist, [trialsN,trialLength/dsF]);
    
    f = figure;
    stdshade(trialDist)
    aesthetics
    xlim(plotWin) % given in downsampled frames
    
    xlabel("Time (ms)")
    ylabel("Distance between output states")
    

end