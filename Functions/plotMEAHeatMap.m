function t = plotMEAHeatMap(dat1,dat2,removeElecs,metric,subtitles,minVal,maxVal) 
%Input = any metric e.g. latency, fr , AC for 60 electrodes from the MEA
%recordings 1x60
%options: lat for latency (ms post stim) ;counts (mean spike count),
%p (ppt firing ), AC (average controllablity)
%adapted from makeHeatMap_AD by Alex Dunn
arguments
    dat1
    dat2
    removeElecs
    metric
    subtitles = [];
    minVal = [];
    maxVal = [];
end


figure("Position",[100, 100, 1500, 750])
t = tiledlayout(1,2,'TileSpacing','Compact');

pltIndex = [21,19,16,15,12,10,...
    24,22,20,17,14,11,9,7,26,...
    25,23,18,13,8,6,5,...
    29,30,28,27,4,3,1,2,...
    32,31,33,34,57,58,60,59,...
    35,36,38,43,48,53,55,56,...
    37,39,41,44,47,50,52,54,...
    40,42,45,46,49,51];

template = zeros(8, 8); 
template(1, 1) = NaN; 
template(1, 8) = NaN; 
template(8, 1) = NaN;
template(8, 8) = NaN;

dat1(removeElecs) = NaN;
dat1 = dat1(pltIndex);
heatMatrix1 = template;
heatMatrix1(2:7) = dat1(1:6);
heatMatrix1(9:56) = dat1(7:54); 
heatMatrix1(58:63) = dat1(55:60);

dat2(removeElecs) = NaN;
dat2 = dat2(pltIndex);
heatMatrix2 = template;
heatMatrix2(2:7) = dat2(1:6);
heatMatrix2(9:56) = dat2(7:54); 
heatMatrix2(58:63) = dat2(55:60);

h(1) = nexttile(t);
% plot grid
hm = imagesc(heatMatrix1); 
set(hm, 'AlphaData', ~isnan(heatMatrix1))
aesthetics
if ~isempty(subtitles) && numel(subtitles) == 2
    title(subtitles(1), FontSize=16)
end
h(1).XAxis.Visible = 'off';
h(1).YAxis.Visible = 'off';
% title("Pre-stim")

h(2) = nexttile(t);
% plot grid
hm = imagesc(heatMatrix2); 
set(hm, 'AlphaData', ~isnan(heatMatrix2))
aesthetics
% title("Post-stim")
if ~isempty(subtitles) && numel(subtitles) == 2
    title(subtitles(2), FontSize=16)
end
h(2).XAxis.Visible = 'off';
h(2).YAxis.Visible = 'off';

% title(t, strcat(pairName,"_heatmap"),'FontSize',28)

if ~isempty(minVal) && ~isempty(maxVal)
    cmin = minVal;
    cmax = maxVal;
else
    cmin = min([dat1,dat2],[],'all');
    cmax = max([dat1,dat2],[],'all');
end
if strcmp(metric, 'latency')
    set(h, 'Colormap', gray, 'CLim', [cmin cmax])
else
    set(h, 'Colormap', parula, 'CLim', [cmin cmax])
end

% Set colour bar properties
cb = colorbar;
if strcmp(metric, 'AC')
     ylabel(cb, 'Average Controllability')
elseif strcmp(metric, 'latency') 
    ylabel(cb, 'Latency post-stimulation (ms)')
elseif strcmp(metric, 'FR') 
    ylabel(cb, 'Spike rate (Hz)')
elseif strcmp(metric, 'p')
    %ylabel(cb, 'Log10 spike count')   
    ylabel(cb, 'Probablility of firing')
%     ylimit_cbar = 1;
%     caxis([0,ylimit_cbar])
elseif strcmp(metric,'foldChangeStim')
    ylabel(cb,'Fold change in firing rate vs baseline')
elseif strcmp(metric,'foldChangePostStim')
    ylabel(cb,'Fold change in firing rate post- vs pre-stim')
elseif strcmp(metric,'foldChangePostSALPA')
    ylabel(cb,'Fold change in spike count post- vs pre-artifact removal')
elseif strcmp(metric, 'spike_count')
    %ylabel(cb, 'Log10 spike count')   
    ylabel(cb, 'Mean spike counts')
end
    cb.FontSize = 16;
   cb.TickDirection = 'out';
%     cb.Location = 'Southoutside';
    cb.Box = 'off';
    cb.Layout.Tile = 'east'; 

end 