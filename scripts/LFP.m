
%% Import metadata

homeDir = "D:\MATLAB\MPhil_scripts\stimulation_tasks\Mona";
cd(homeDir)
metadataSpreadsheet = "OWT.xlsx";

% Metadata
sheet = 3; %2-spikes, 1-voltage=
xlRange = 'A2:M15';
[num,txt,~] = xlsread(metadataSpreadsheet,sheet,xlRange);
samples = txt(:,2);
baseline_recs = txt(:,3);


% Paths
addpath("spikes")

% Parameters for binning recording 
fs = 25e3;

%% Get spike rates

for n = 1:length(samples)

    cd spikes

    try
        load(strcat(baseline_recs{n},'.mat'),'cSpikes')
        cSpikes = full(cSpikes);
    catch
        load(strcat(baseline_recs{n},'.mat'),'spikes')
        cSpikes = spikes;
        clear spikes
    end    
    cd(homeDir)

    disp(baseline_recs{n})
    total_samples = length(cSpikes);
    bin_duration = 1 * fs;
    bin_onset_t = 0:bin_duration:total_samples;
    n_bins = total_samples/bin_duration;
    n_channels = 60;

    x = zeros(n_channels, n_bins);
    
    for bin = 1:n_bins
    
        start = bin_onset_t(bin) + 1;
        stop = bin_onset_t(bin + 1);
    
        bin_win = cSpikes(start:stop,:);
        spike_rates = sum(bin_win, 1);
        x(:,bin) = spike_rates;
        clear bin
    
    end

    % Plot raster
    Params.fs = 1;
    Params.rasterPlotUpperPercentile = 100;
    Params.figExt = {'.png'};
    Params.fullSVG = 0;
    spikeFreqMax = max(x,[],'all');
    figFolder = 'D:\MATLAB\MPhil_scripts\stimulation_tasks\Mona\network_response\LFPs';
    rasterPlot(baseline_recs{n},x',Params,spikeFreqMax,figFolder)

end

