
[fileNames,pathName] = uigetfile('E:\zon_lab\FACS\*.fcs','MultiSelect','on');

if iscell(fileNames) == 1
    numFiles = numel(fileNames);
else
    numFiles = 1;
    fileNames = {fileNames};
end
multiWaitbar('Analyzing flow cytometry data...',0);
rho = cell(numFiles,1);
delta = cell(numFiles,1);
nneigh = cell(numFiles,1);
for kk = 1:numFiles
    file = fullfile(pathName,fileNames{kk});
    %Step 1 use zbow_logicle to transform FACS data.
    [data,normData,customData,ternColor,ternCoords,sampleName] = zbow_logicle(file,[],20000);
    
    [N{kk}, Xedges, Yedges] = histcounts2(ternCoords(:,1),ternCoords(:,2),15,'Normalization','pdf');
%     figure, hist3(ternCoords);
%     set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%     colormap('parula');
    

    
%     [rho{kk}, delta{kk},nneigh{kk}] = deltarho(normData,1);
%     
%     figure, scatter3(ternCoords(:,1),ternCoords(:,2),rho{kk},15,customData,'filled');
%     hh = title(sampleName);
%     hh.Interpreter = 'none';
    
%     %Step 2 Get delta, rho, and nneigh from deltarho function
%     
%     [rho, delta, nneigh] = deltarho(ternCoords,1);
%     
%     clusterMetric = rho.*delta';
%     
%     metricStd = std(clusterMetric);
%     metricMean = mean(clusterMetric);
%     
%     threshold = metricMean + 1.25*metricStd;
%     
%     figure,
%     subplot(2,2,3),
%     plotSpread(clusterMetric);
%     
%     hold on
%     plot([0 1 2], [metricMean metricMean metricMean],'--r');
%     plot([0 1 2], [threshold threshold threshold],'--g');
%     hold off
%     
%     set(gca,'YScale','log');
%     sampleName = strrep(sampleName,'_','.');
%     t = suptitle(sampleName);
%     % set(t,'Interpreter','none');
%     
%     clusterCenterIdx = find(clusterMetric > threshold);
%     
%     clusterIdx = assign_cluster(rho,nneigh,clusterCenterIdx);
%     
%     clusterInfo = tabulate(clusterIdx);
%     
%     clusterCenterColor = customData(clusterCenterIdx,:);
%     
%     subplot(2,2,4),
%     clusterBarGraph(clusterInfo,clusterCenterColor);
%     
%     subplot(2,2,1)
figure,
ternPlot(ternCoords,customData,'false');
t2 = title(sampleName);
%     
%     groupColorMap = distinguishable_colors(max(clusterInfo(:,1)));
%     groupColor = zeros(size(clusterIdx,1),3);
%     
%     for ii = 1:size(clusterIdx,1)
%         groupColor(ii,:) = groupColorMap(clusterIdx(ii),:);
%     end
%     subplot(2,2,2),
%     ternPlot(ternCoords, groupColor,'false');
% %     t3 = title({sampleName, 'ternary plot with cluster color'});
    drawnow;
    multiWaitbar('Analyzing flow cytometry data...',kk/numFiles);
end

for jj = 1:numFiles
    maxProb(jj) = max(N{jj}(:));
end

maxProb = maxProb';
figure, scatter(1:numFiles,maxProb);
multiWaitbar('CloseAll');