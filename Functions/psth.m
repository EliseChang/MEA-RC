function t = psth(data, var, stimProt, binSize, baseline)
% binSize: ms
% baselineRate: Hz

figure('Position', [100, 100, 1500, 750])
t = tiledlayout(1, 2, 'TileSpacing','Compact');

maskA = logical(stimProt);
maskB = logical(~stimProt);

if strcmp(var, 'spike_prob')
    plotBaseline = baseline * binSize;
else
    plotBaseline = baseline*binSize*1e-3;
end


%% Plot Pattern A
ax1 = nexttile;
trials = data(maskA,:);
if strcmp(var, 'spike_prob')
    trials = trials > 0;
end
histData = mean(trials,1);
cvData = std(trials,1) ./ histData;

yyaxis left
psth = bar(histData, 1); % bars touching
set(psth,'edgecolor','black','facecolor',"#0072BD");
xticklabels(xticks*binSize)
if strcmp(var, 'spike_prob')
    ylabel("Spike probability")
else
    ylabel("Mean network spike count across trials")
end

yyaxis right
plot(cvData, "--", "Color", "#D95319", "LineWidth", 1)
if strcmp(var, 'spike_prob')
    ylabel("Fano factor")
else
    ylabel("CV of network spike count across trials")
end

yline(plotBaseline, "black", "LineWidth", 2, "Label", "Baseline")

title("Pattern A")
aesthetics

clear mask trials
%% Plot Pattern B
ax2 = nexttile;
trials = data(maskB,:);

if strcmp(var, 'spike_prob')
    trials = trials > 0;
end
histData = mean(trials,1);
cvData = std(trials,1) ./ histData;

yyaxis left
psth = bar(histData, 1); % bars touching
set(psth,'edgecolor','black','facecolor',"#0072BD");
xticklabels(xticks*binSize)
if strcmp(var, 'spike_prob')
    ylabel("Spike probability")
else
    ylabel("Mean network spike count across trials")
end

yyaxis right
plot(cvData, "--", "Color", "#D95319", "LineWidth", 1)
if strcmp(var, 'spike_prob')
    ylabel("Fano factor")
else
    ylabel("CV of network spike count across trials")
end

yline(plotBaseline, "black", "LineWidth", 2, "Label", "Baseline")

title("Pattern B")
aesthetics

clear mask trials

%% Format whole figure
yl1 = ax1.YAxis(1); yr1 = ax1.YAxis(2);
yl2 = ax2.YAxis(1); yr2 = ax2.YAxis(2);
linkprop([yl1 yl2], 'Limits');
linkprop([yr1 yr2], 'Limits');
xticklabels(xticks*binSize)

xlabel(t,"Time post-stimulus (ms)")
title(t,"Peri-stimulus Time Histogram")


end
