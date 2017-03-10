
[fileNames,pathName] = uigetfile('E:\zon_lab\FACS\*.fcs','MultiSelect','on');

if iscell(fileNames) == 1
    numFiles = numel(fileNames);
else
    numFiles = 1;
    fileNames = {fileNames};
end

plotQuest = questdlg('Plot ternary graphs?','Plot question','Yes','No','Yes');
multiWaitbar('Analyzing flow cytometry data...',0);
rho = cell(numFiles,1);
delta = cell(numFiles,1);
nneigh = cell(numFiles,1);
recombination = cell(numFiles,1);
clusterMetric = zeros(numFiles,1);
clusterCenterColor = cell(numFiles,1);
meanBigClustColor = zeros(numFiles,3);
sampleName = cell(numFiles,1);
normRho = cell(numFiles,1);
normDelta = cell(numFiles,1);
clusterSize = cell(numFiles,1);
meanClusterColor = cell(numFiles,1);

for kk = 1:numFiles
    file = fullfile(pathName,fileNames{kk});
    %Step 1 use zbow_logicle to transform FACS data.
    [data,normData,customData,ternColor,ternCoords,sampleName{kk}] = zbow_logicle(file,[],20000);
    sampleName{kk} = strrep(sampleName{kk},'_','.');
    sampleSize = 2500;
    
    totalCells = size(data,1);
    
    
    %exclude red cells - probably could be tweaked a bit
    red_cells= ternColor(:,1)>=0.8;
    red_count = sum(red_cells);
    nonred_count = totalCells - red_count;
    recombination{kk} = nonred_count/totalCells;
    
    if nonred_count > 10
        
        ternCoords = ternCoords(~red_cells,:);
        customData = customData(~red_cells,:);
        
        if nonred_count < sampleSize
            dataSamp = ternCoords;
            dataSampIdx = 1:nonred_count;
        else
            [dataSamp, dataSampIdx] = datasample(ternCoords,sampleSize,'Replace',false);
        end
        
        [rho{kk}, delta{kk}, nneigh{kk}] = deltarho(dataSamp,1);
        
        normRho{kk} = myNorm(rho{kk});
        normDelta{kk} = myNorm(delta{kk});
        
        rhoThresh = find(normRho{kk} > 0.1);
        deltaThresh = find(normDelta{kk} > 0.2);
        clusterCentersSampleIdx = intersect(rhoThresh,deltaThresh);
        
        clusterCentersIdx = dataSampIdx(clusterCentersSampleIdx);
        
        clusterCenters = ternCoords(clusterCentersIdx,:);
        
        %For top density
        %         [maxRho, maxRhoIdx] = max(rho{kk});
        %
        %         maxPointIdx = dataSampIdx(maxRhoIdx);
        %         maxPoint3D = customData(maxPointIdx,:);
        %         maxPointTern = ternCoords(maxPointIdx,:);
        
        
        
        %calculates pairwise distance of cluster centers with all cells
        switch plotQuest
            case 'Yes'
                figure,
            case 'No'
        end
        clusterSize{kk} = zeros(size(clusterCenters,1),1);
        subPlotNum = numSubplots(size(clusterCenters,1));
        
        meanClusterColor{kk} = zeros(size(clusterCenters,1),3);
        for nn = 1:size(clusterCenters,1)
            %             D = pdist2(maxPointTern,ternCoords);
            %         D = pdist2(maxPoint3D,customData);
            D = pdist2(clusterCenters(nn,:),ternCoords);
            D = D';
            %     set cutoff max radius, everything within this sphere gets counted
            %     towards willmetric, 0.2 works better than 0.3, also could be tweaked
            Dclose = D<=0.08;
            Dsums = sum(Dclose);
            clusterSize{kk}(nn) = 100.*Dsums./totalCells;
            clusterColor = 0.8*ones(nonred_count,3);
            clusterColor(Dclose,1) = customData(Dclose,1);
            clusterColor(Dclose,2) = customData(Dclose,2);
            clusterColor(Dclose,3) = customData(Dclose,3);
            
            meanClusterColor{kk}(nn,1) = mean(customData(Dclose,1));
            meanClusterColor{kk}(nn,2) = mean(customData(Dclose,2));
            meanClusterColor{kk}(nn,3) = mean(customData(Dclose,3));
            switch plotQuest
                case 'Yes'
                    subplot(subPlotNum(1),subPlotNum(2),nn),
                    ternPlot(ternCoords,clusterColor,'false');
                    drawnow;
                    text(0,0.75,sampleName{kk});
                case 'No'
            end
            
        end
        
        
        
        %     [N{kk}, Xedges, Yedges] = histcounts2(ternCoords(:,1),ternCoords(:,2));
        %
        %
        %     [maxN(kk), maxNIdx(kk)] = max(N{kk}(:));
        %     bigClust(kk) = 100.*maxN(kk)./(size(data,1));
        
        %     [xx, yy] = ind2sub(maxNIdx(kk));
        
        %     figure, hist3(ternCoords);
        %     set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        %     colormap('parula');
        
        
        
        
        %     %Step 2 Get delta, rho, and nneigh from deltarho function
        %
        %     [rho, delta, nneigh] = deltarho(ternCoords,1);
        %
        %     clusterMetric = rho.*delta';
        %
        %     metricStd = nanstd(clusterMetric);
        %     metricMean = nanmean(clusterMetric);
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
        %     set(t,'Interpreter','none');
        %
        %     clusterCenterIdx = find(clusterMetric > threshold);
        %
        %     clusterIdx = assign_cluster(rho,nneigh,clusterCenterIdx);
        %
        %     clusterInfo = tabulate(clusterIdx);
        %
        %     clusterCenterColor = ternCoords(clusterCenterIdx,:);
        %     clusterCenterColor3D = customData(clusterCenterIdx,:);
        %
        %     subplot(2,2,4),
        %     clusterBarGraph(clusterInfo,clusterCenterColor3D);
        %
        %     subplot(2,2,1)
        %     ternPlot(ternCoords,customData,'false');
        %     t2 = title(sampleName);
        %
        %     %calculates pairwise distance of cluster centers with all cells
        %     %     D = pdist2(clusterCenterColor{kk},customData);
        %     D = pdist2(clusterCenterColor3D,customData);
        %     D = D';
        %     % set cutoff max radius, everything within this sphere gets counted
        %     % towards willmetric, 0.2 works better than 0.3, also could be tweaked
        %     Dclose = D<=0.2;
        %     DIdx = find(sum(Dclose,2) > 0);
        %     Dsums = sum(Dclose);
        %     willmetric{kk} = Dsums./nonred_count;
        %
        %     groupColorMap = distinguishable_colors(max(clusterInfo(:,1)));
        %
        %
        %     DColor = 0.8*ones(size(ternCoords,1),3);
        %     DColor(DIdx,1) = 0;
        %     DColor(DIdx,2) = 0.2;
        %     DColor(DIdx,3) = 0.8;
        %
        %
        %     groupColor = zeros(size(clusterIdx,1),3);
        %
        %     for ii = 1:size(clusterIdx,1)
        %         groupColor(ii,:) = groupColorMap(clusterIdx(ii),:);
        %     end
        %     subplot(2,2,2),
        %     ternPlot(ternCoords, groupColor,'false');
        %     t3 = title({sampleName, 'ternary plot with cluster color'});
        %
        %     figure,
        %     ternPlot(ternCoords,DColor,'false');
        %
        drawnow;
        multiWaitbar('Analyzing flow cytometry data...',kk/numFiles);
    else
    end
end

if numFiles > 1
    switch plotQuest
        case 'Yes'
            cascade;
        case 'No'
    end
end

%
% Plotting the cluster sizes
figure,
hold on
for jj = 1:numFiles
    plot([jj jj], [0 max(clusterSize{jj})],'--k');
    xx = jj*ones(size(clusterSize{jj},1),1);
    scatter(xx,clusterSize{jj},200,meanClusterColor{jj},'filled');
end

ylim([0,70]);
hold off

set(gca,'XTick',1:numFiles);
xticklabels(sampleName);
xtickangle(90);



multiWaitbar('CloseAll');