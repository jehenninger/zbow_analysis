function [barH] = clusterBarGraph(clusterInfo, clusterColor)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

sortTab = flipud(sortrows(clusterInfo,3));
barH = bar([transpose(sortTab(:,3));nan(size(transpose(sortTab(:,3))))],'Stacked');
hold on
for ii = 1:numel(sortTab(:,3))
    set(barH(ii),'FaceColor',clusterColor(sortTab(ii,1),:));
end
ylim([0 100]);
xlim([0 2]);
set(gca,'ygrid','on');
set(gca,'xtick',1);
set(gca,'FontSize',20,'FontName','Arial');
set(gca,'XTickLabel',' ');
title([num2str(size(sortTab,1)), ' clusters']);

end

