function p = plotTrialRaster(x, mask, channelOrder, ylines, labels, subtitles, varLabel, figPos)
% x: [channel X trial]
% mask: trials to plot in LHS plot

p = figure("Position",figPos);
t = tiledlayout(1,2, 'TileSpacing','compact');
trialsN = sum(mask);
channelsN = size(x,1);
trials1 = x(channelOrder,mask);
trials2 = x(channelOrder,~mask);

nexttile
imagesc(trials1)
xticks([1,50:50:trialsN])
yticks([1,10:10:channelsN])
hold on
yline(ylines, '--', labels, 'Color', [0.7 0.7 0.7])
hold off
aesthetics
title(subtitles(1))
set(gca, 'FontSize',24)

nexttile
imagesc(trials2)
xticks([1,50:50:trialsN])
yticks([1,10:10:channelsN])
hold on
yline(ylines, '--', labels, 'Color', [0.7 0.7 0.7])
hold off
aesthetics
title(subtitles(2))
set(gca, 'FontSize',24)

cb = colorbar;
cb.Label.String = varLabel;
cb.FontSize = 24;
cb.TickDirection = 'out';
ylabel(t, 'Electrode', 'FontSize',28)
xlabel(t, 'Trial', 'FontSize',28)
end

