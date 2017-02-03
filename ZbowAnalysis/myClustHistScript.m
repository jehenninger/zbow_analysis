% myclear all
%TO ADD: DISTINGUISH XLS AND CSV INPUTS
clear all


folder_name = uigetdir('E:\Jon\Google Drive\');
dataType = questdlg('Excel or csv?','Excel or csv','Excel','csv','Excel');

switch dataType
    case 'Excel'
        file_list = dir(fullfile(folder_name,'*_Summary*.xlsx'));
    case 'csv'
        file_list = dir(fullfile(folder_name,'*_Summary*.csv'));
end

label_length = 12;

for k = 1:length(file_list)
    file_names{k} = file_list(k).name;
    
    switch dataType
        case 'Excel'
            data{k} = xlsread(fullfile(folder_name,file_list(k).name));
        case 'csv'
            data{k} = csvread(fullfile(folder_name,file_list(k).name));
    end
    
end

labels = {file_list(:).name};

for kk = 1:length(labels)
    cellContents = labels{kk};
    labels{kk} = cellContents(1:label_length);
end
% hold off

for k = 1:length(data)
    currentData = data{k};
    
    red = currentData(:,4);
    green = currentData(:,5);
    blue = currentData(:,6);
    
    redTern = red./(red+green+blue);
    redOnly{k} = currentData(redTern == max(redTern(:)),:);
    currentData(redTern == max(redTern(:)),:) = [];
    
    currentData(currentData(:,2)<5,:) = [];
    
    dataPro{k} = currentData;
    clusterPercentages{k} = dataPro{k}(:,3);
    
    medianCluster(k) = median(dataPro{k}(:,3));
    
    
    
    %For Quartiles +/- IQR*1.5
    percentiles{k} = prctile(zscore(clusterPercentages{k}), [25 75]);
    outliers{k} = clusterPercentages{k} < (percentiles{k}(1)-1.5*iqr(clusterPercentages{k}))...
        | clusterPercentages{k} > (percentiles{k}(2)+1.5*iqr(clusterPercentages{k}));
    noOutliers{k} = clusterPercentages{k}(~outliers{k});
    
    % %     %For MAD
    %     scaleFactor(k) = 1/prctile(zscore(clusterPercentages{k}),75);
    % %     scaleFactor = 1.4826;
    %     madStat(k) = scaleFactor(k).*mad(clusterPercentages{k},1);
    % %     outliers{k} = clusterPercentages{k} < medianCluster(k)-2.5*medianAbs{k}...
    % %                   | clusterPercentages{k} > medianCluster(k)+2.5*medianAbs{k};
    %     testStat{k} = abs(clusterPercentages{k}-median(clusterPercentages{k}))./madStat(k);
    %     outliers{k} = clusterPercentages{k}(testStat{k}>=3);
    %     noOutliers{k} = clusterPercentages{k}(~outliers{k});
    
%     %For Sn MAD
%     
%     distance{k} = squareform(pdist(clusterPercentages{k}));
%     
%     for ii = 1:numel(clusterPercentages{k})
%         distance{k}(numel(clusterPercentages{k})+1,ii) = median(distance{k}(:,ii));
%     end
%     sN(k) = median(distance{k}(numel(clusterPercentages{k})+1,:));
%     
%     if mod(numel(clusterPercentages{k}),2) == 0
%         cN(k) = numel(clusterPercentages{k})/(numel(clusterPercentages{k})+3.8);
%     else
%         cN(k) = numel(clusterPercentages{k})/(numel(clusterPercentages{k})-0.9);
%     end
%     
%     sN(k) = cN(k).*1.1926.*sN(k);
%     
%     testStat{k} = abs(clusterPercentages{k}-median(clusterPercentages{k}))./sN(k);
%     outliers{k} = testStat{k}>=2;
%     noOutliers{k} = clusterPercentages{k}(~outliers{k});
%     
    meanClusterSize(k) = mean(noOutliers{k});
    stdClusterSize(k) = std(noOutliers{k})./sqrt(numel(noOutliers{k}));
    predictedClones(k) = 100/meanClusterSize(k);
    individualClones{k} = 100./noOutliers{k};
    individualClonesStd(k) = std(individualClones{k})./sqrt(numel(noOutliers{k}));
    
end

[~,medianIndex] = sort(medianCluster);

clear dataProCombined
dataProCombined(:,1) = dataPro{1}(:,3);
for k = 2:length(dataPro)
    dataProCombined = padadd(dataProCombined,dataPro{k}(:,3));
end

% dataNoOutliersCombined(:,1) = noOutliers{1};
% for k = 2:length(noOutliers)
%     dataNoOutliersCombined = padadd(dataNoOutliersCombined,noOutliers{k});
% end


figure,
subplot(2,2,2)
boxPlotHandleSorted = boxplot(dataProCombined(:,medianIndex));
title1 = title([file_names{1,1}(1:6),' Sorted Cluster Distribution']);
set(title1,'interpreter','none');

subplot(2,2,1)
boxPlotHandle = boxplot(dataProCombined);
title2 = title([file_names{1,1}(1:6),' Cluster Distribution']);
set(title2,'interpreter','none');

subplot(2,2,3)
clusterPlot = errorbar(1:numel(predictedClones),predictedClones,individualClonesStd,'.','MarkerSize',25);
title3 = title([file_names{1,1}(1:6),' Predicted Clones']);
set(title3,'interpreter','none');

[~, sortIndex] = sort(predictedClones);

subplot(2,2,4)
clusterPlotSorted = errorbar(1:numel(predictedClones),predictedClones(sortIndex), individualClonesStd(sortIndex),'.',...
    'MarkerSize',25);
title4 = title([file_names{1,1}(1:6),' Sorted Predicted Clones']);
set(title4,'interpreter','none');




for ll = 1:length(noOutliers)
    medianTest(ll) = median(noOutliers{ll});
    medianPredictedTest(ll) = 100/medianTest(ll);
end

disp('Predicted clones from median');
averagePredictedCloneFromMedian = mean(medianPredictedTest)
disp('Predicted clones from average');
averagePredictedCloneFromMean = mean(predictedClones)

myClusterPlot(dataProCombined,dataPro,labels,1);
