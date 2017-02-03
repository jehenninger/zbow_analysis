function [h] = myClusterPlot(data,myColor,labels, transform,sorted)
%Takes input data and and plots individual fish with dots for each clusters
%on linear or log scale
% data is matrix with rows corresponding to clusters and columns
% corresponding to individual fish. Use NaN to fill in elements.
%
% color is cell array of individual fish with columns 4 to 6 equal to the
% RGB values

% transform is 1 for linear and 2 for log2

% sorted is a cell array specifying controls (1 cell, indices of control samples)
% and experimental (1 cell, indices of experimental samples), and will sort
% the clusters by median percentage OR it is a vector specifying a sort
% index
if exist('sorted','var')
    if iscell(sorted) == 1 %Only if sorted cell array is provided
        
        sortType = questdlg('Mean or median?','Mean or median','Median','Mean','Mean');
        groupNum = size(sorted,2);
        
        switch sortType
            case 'Median'
                
                for kk = 1:groupNum
                    
                    groupData{kk} = data(:,sorted{kk});
                    groupMedian{kk} = nanmedian(groupData{kk});
                    [~, groupMedianIndex{kk}] = sort(groupMedian{kk});
                    groupColor{kk} = myColor(sorted{kk});
                    groupColorSorted{kk} = groupColor{kk}(groupMedianIndex{kk});
                end
                data = cell2mat(groupData);
                myColor = cell2mat(groupColorSorted);
                new_labels = labels(cell2mat(groupMedianIndex));
                labels = new_labels;
                
%                 OLD WAY
%                 controlData = data(:,sorted{1}); %Parse control data from data matrix
%                 conMedian = nanmedian(controlData); %Find median of the control population without NaN
%                 [~, conMedianIndex] = sort(conMedian); %Sort by median and get indices
%                 
%                 expData = data(:,sorted{2}); %Same as above but with experimental
%                 expMedian = nanmedian(expData);
%                 [~, expMedianIndex] = sort(expMedian);
%                 
%                 controlColor = myColor(sorted{1}); %Parse color of control samples from color matrix
%                 controlColorSorted = controlColor(conMedianIndex); %Sort color matrix by the median index
%                 
%                 expColor = myColor(sorted{2}); %Same as control above
%                 expColorSorted = expColor(expMedianIndex);
%                 
%                 data = [controlData(:,conMedianIndex), expData(:,expMedianIndex)]; %Concatenate into one matrix for graphing below
%                 myColor = [controlColorSorted, expColorSorted]; %Concatenate into one matrix for graphing below
%                 
                title1 = 'Sorted Linear Cluster Percentages Median';
                title2 = 'Sorted Log2 Cluster Percentages Median';
                
            case 'Mean'
                
                for kk = 1:groupNum
                    
                    groupData{kk} = data(:,sorted{kk});
                    for ii = 1:size(groupData{kk},2)
                        if ii == 1
                            groupDataVector{kk} = log2(groupData{kk}(:,ii));
                        else
                            groupDataVector{kk} = [groupDataVector{kk};log2(groupData{kk}(:,ii))];
                        end
                    end
                    groupDataVector{kk}(isnan(groupDataVector{kk})) = [];
                    groupMean{kk} = nanmean(groupData{kk});
                    [~, groupMeanIndex{kk}] = sort(groupMean{kk});
                    
                    adjustedLabel(kk) = numel(groupMeanIndex{kk});
                    if kk == 1
                        adjustedLabelIndex{kk} = groupMeanIndex{kk};
                    else
                        adjustedLabelIndex{kk} = groupMeanIndex{kk}+adjustedLabel(kk-1);
                    end
                    
                    
                    groupDataSorted{kk} = groupData{kk}(:,groupMeanIndex{kk});
                    groupColor{kk} = myColor(sorted{kk});
                    groupColorSorted{kk} = groupColor{kk}(groupMeanIndex{kk});
                    new_labels{kk} = labels(sorted{kk});
                    new_labels{kk} = new_labels{kk}(groupMeanIndex{kk});
                end
                
                
                
                myColor = {[groupColorSorted{:}]};
                myColor = [myColor{:}];
                
                data = cell2mat(groupDataSorted);
                
                new_labels = [new_labels{:}];
%                 adjustedLabelIndex = [adjustedLabelIndex{:}];
%                 new_labels = new_labels(adjustedLabelIndex);
                labels = new_labels;
                
%                 OLD WAY
%                 controlData = data(:,sorted{1}); %Parse control data from data matrix
%                 conMean = nanmean(controlData); %Find median of the control population without NaN
%                 [~, conMeanIndex] = sort(conMean); %Sort by median and get indices
%                 
%                 expData = data(:,sorted{2}); %Same as above but with experimental
%                 expMean = nanmean(expData);
%                 [~, expMeanIndex] = sort(expMean);
%                 
%                 controlColor = myColor(sorted{1}); %Parse color of control samples from color matrix
%                 controlColorSorted = controlColor(conMeanIndex); %Sort color matrix by the median index
%                 
%                 expColor = myColor(sorted{2}); %Same as control above
%                 expColorSorted = expColor(expMeanIndex);
%                 
%                 data = [controlData(:,conMeanIndex), expData(:,expMeanIndex)]; %Concatenate into one matrix for graphing below
%                 myColor = [controlColorSorted, expColorSorted]; %Concatenate into one matrix for graphing below
%                 
                title1 = 'Sorted Linear Cluster Percentages Mean';
                title2 = 'Sorted Log2 Cluster Percentages Mean';
        end
        
    else
        data = data(:,sorted);
        myColor = myColor(sorted);
        labels = labels(sorted);
        title1 = 'Sorted by Index Linear Cluster Percentages';
        title2 = 'Sorted by Index Log2 Cluster Percentages';
    end
    
else
    title1 = 'Unsorted Linear Cluster Percentages';
    title2 = 'Unsorted Log2 Cluster Percentages';
    
end

boxPlotQuestion = questdlg('Also plot boxplots?','Boxplot Option','Yes','No','Yes');

if transform == 1
    figure,
    switch boxPlotQuestion
        case 'Yes'
            boxplot(data,'colors','k','whisker',1000);
    end
    
    hold on
    for jj = 1:numel(data(1,:))
        h = scatter(jj*ones(numel(data(:,jj))-sum(isnan(data(:,jj))),1),...
            data(~isnan(data(:,jj)),jj),110,myColor{jj}(:,4:6),'filled');
        hold on
    end
    
    hold off
    title(title1);
    ax = gca;
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 45;

    
    
elseif transform == 2
    
    figure,
    switch boxPlotQuestion
        case 'Yes'
            boxplot(log2(data),'colors','k','whisker',1000);
    end
    hold on
    for jj = 1:numel(data(1,:))
        h = scatter(jj*ones(numel(data(:,jj))-sum(isnan(data(:,jj))),1),...
            log2(data(~isnan(data(:,jj)),jj)),110,myColor{jj}(:,4:6),'filled');
        hold on
    end
    hold off
    title(title2);
    ax = gca;    
    ax.XTickLabel = labels;
    ax.XTickLabelRotation = 45;
    

    
    yLabels = str2double(ax.YTickLabel);
    yLabels = 2.^yLabels;
    for jj = 1:numel(yLabels)
        yLabelsNew{jj} = num2str(yLabels(jj));
    end
    ax.YTickLabel = yLabelsNew;
end

if exist('groupNum','var') == 1
    %To make histograms of grouped data
    figure,
    for kk = 1:groupNum
        x_values = 0:0.1:100;
        columnData = reshape(data(:,sorted{kk}),numel(data(:,sorted{kk})),1);
        pd{kk} = fitdist(columnData,'Kernel','BandWidth',1);
        pdfunc{kk} = pdf(pd{kk},x_values);
        legendInfo{kk} = ['Group = ' num2str(kk)];
        plot(x_values,pdfunc{kk},'LineWidth',2.5);
        hold on
    end
    hold off
    legend(legendInfo);
    
    %To make box plots of distributions
    figure,
    for kk = 1:groupNum
        means{kk} = nanmean(data(:,sorted{kk}));
        plot(kk*ones(1,numel(means{kk})),means{kk},'.');
        hold on
    end
    
    combinedMeans(:,1) = means{1};
    for kk = 2:groupNum
        combinedMeans = padadd(combinedMeans, means{kk}');
    end
    
    boxplot(combinedMeans,'colors','k');
    hold off
    xlim([0, groupNum+1]);
end

if exist('groupNum','var') == 1
    distColors = distinguishable_colors(groupNum);
    figure,
    distributionPlot(groupDataVector,'color',mat2cell(distColors,ones(groupNum,1),3),...
        'histOpt',1,'addSpread',1);
end

end

