function t = plotPSTH(spikeTimes0, spikeTimes1, trialsN, trialLength, pairName, options)

arguments
    spikeTimes0
    spikeTimes1
    trialsN
    trialLength % frames
    pairName
    options.binsize = 1; % miliseconds
    options.limits = [0,100]; % x-limits
    options.baselineSpikeCount = nan % baseline spike count
end

earlyTimepoints = 0:20;
lateTimepoints = 21:500;

figure("Position",[1, 49, 2560,1315])
t = tiledlayout(1,2,'TileSpacing','Compact');

% Pattern 0
nexttile
hold on
[~,counts] = psth(spikeTimes0, options.binsize, options.fs, trialsN, trialLength, gca);
xlim(options.limits) % look at first 100 ms

% Plot lines showing average firing rates
if ~isnan(baselineSpikeCount)
    a = yline(baselineSpikeRate,'LineStyle','--','LineWidth',5,'Color','red','DisplayName','Average baseline');
%     earlySpikeRate = sum(counts(1:earlyTimepoints(end)))/length(earlyTimepoints);
%     b = line(earlyTimepoints,repmat(earlySpikeRate,length(earlyTimepoints)),'LineStyle','--','LineWidth',5,'Color','blue',...
%         'DisplayName','Average early component');
%     lateSpikeRate = sum(counts(lateTimepoints))/length(lateTimepoints);
%     c = line(lateTimepoints,repmat(lateSpikeRate,length(lateTimepoints)),'LineStyle','--','LineWidth',5,'Color','magenta',...
%         'DisplayName','Average late component');
end

title("Pattern A")
aesthetics
hold off

% Pattern 1
nexttile
hold on
[~,counts] = psth(spikeTimes1, binsize, fs, trialsN, trialLength, gca);
xlim([limits]) % look at first 100 ms

% Plot lines showing average firing rates
if ~isnan(baselineSpikeCount)
    a = yline(baselineSpikeRate,'LineStyle','--','LineWidth',5,'Color','red','DisplayName','Average baseline');earlySpikeRate = sum(counts(1:earlyTimepoints(end)))/length(earlyTimepoints);
%     b = line(earlyTimepoints,repmat(earlySpikeRate,length(earlyTimepoints)),'LineStyle','--','LineWidth',5,'Color','blue',...
%         'DisplayName','Average early component');
%     lateSpikeRate = sum(counts(lateTimepoints))/length(lateTimepoints);
%     c = line(lateTimepoints,repmat(lateSpikeRate,length(lateTimepoints)),'LineStyle','--','LineWidth',5,'Color','magenta',...
%         'DisplayName','Average late component');
end

title("Pattern B")
aesthetics
hold off

xlabel(t, "Time post-stimulus (ms)",'FontSize',28)
ylabel(t, "Trial-averaged spike count",'FontSize',28)
title(t, strcat(pairName,"_PTSH"),'FontSize',28)
legend([a]) %b(1),c(1)