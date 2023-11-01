function t = plotPSTH(spikeTimes0, spikeTimes1, binsize, fs, trialsN, trialLength, pairName, baselineSpikeRate)

earlyTimepoints = 0:20;
lateTimepoints = 21:100;

figure("Position",[1, 49, 2560,1315])
t = tiledlayout(1,2,'TileSpacing','Compact');

% Pattern 0
nexttile
hold on
[~,counts] = psth(spikeTimes0, binsize, fs, trialsN, trialLength, gca);
xlim([0 180]) % look at first 100 ms

% % Plot lines showing average firing rates
% if ~isnan(baselineSpikeRate)
%     a = yline(baselineSpikeRate,'LineStyle','--','LineWidth',5,'Color','red','DisplayName','Average baseline');
%     earlySpikeRate = sum(counts(1:earlyTimepoints(end)))/length(earlyTimepoints);
%     b = line(earlyTimepoints,repmat(earlySpikeRate,length(earlyTimepoints)),'LineStyle','--','LineWidth',5,'Color','blue',...
%         'DisplayName','Average early component');
%     lateSpikeRate = sum(counts(lateTimepoints))/length(lateTimepoints);
%     c = line(lateTimepoints,repmat(lateSpikeRate,length(lateTimepoints)),'LineStyle','--','LineWidth',5,'Color','magenta',...
%         'DisplayName','Average late component');
% end

% title("Pattern A")
aesthetics(40)
hold off

% Pattern 1
% nexttile
% hold on
% [~,counts] = psth(spikeTimes1, binsize, fs, trialsN, trialLength, gca);
% xlim([0 180]) % look at first 100 ms

% % Plot lines showing average firing rates
% if ~isnan(baselineSpikeRate)
%     a = yline(baselineSpikeRate,'LineStyle','--','LineWidth',5,'Color','red','DisplayName','Average baseline');earlySpikeRate = sum(counts(1:earlyTimepoints(end)))/length(earlyTimepoints);
%     b = line(earlyTimepoints,repmat(earlySpikeRate,length(earlyTimepoints)),'LineStyle','--','LineWidth',5,'Color','blue',...
%         'DisplayName','Average early component');
%     lateSpikeRate = sum(counts(lateTimepoints))/length(lateTimepoints);
%     c = line(lateTimepoints,repmat(lateSpikeRate,length(lateTimepoints)),'LineStyle','--','LineWidth',5,'Color','magenta',...
%         'DisplayName','Average late component');
% end

% title("Pattern B")
% aesthetics(40)
% hold off

xlabel(t, "Time post-stimulus (ms)",'FontSize',28)
ylabel(t, "Trial-averaged spike count",'FontSize',28)
% title(t, strcat(pairName,"_PTSH"),'FontSize',28)
% legend([a b(1),c(1)])