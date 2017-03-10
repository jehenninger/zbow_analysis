myColor = distinguishable_colors(10);
count = 1;
h1 = figure;
h2 = figure;

for jj = 20
    
    theorHSC  = jj;
    
    autoDiversity = zeros(numFiles,3);
    for i = 1:numFiles
        clusters = zeros(theorHSC,1);
        nonRedClusters = clusterSize{i}./100;
        nonClusterPercent = 1 - sum(nonRedClusters);
        if sum(nonRedClusters) >1
            nonClusterPercent = 0;
        end
        numClusters = size(nonRedClusters(:,1),1);
        percentPerClust = nonClusterPercent./(theorHSC-numClusters);


        if percentPerClust < 0
            clusters(1:numClusters) = nonRedClusters;
            autoDiversity(i,3) = ginicoeff(clusters);
            clustersAdded = 0;
        else
        clusters(1:numClusters) = nonRedClusters;
        clusters((numClusters+1):theorHSC) = ones((theorHSC - numClusters),1).*percentPerClust;
        clustersCombined{i} = clusters;
        
        
        %SHannon Entropy
        %     autoDiversity(i,1) = shannon_entro(clusters);
        %Shannon ENS
        %     autoDiversity(i,2) = exp(autoDiversity(i,1));
        
        %gini index:
        autoDiversity(i,3) = ginicoeff(clusters);
        end


    end
    figure(h1);
    plotSpread(clustersCombined);
%     hold on
%     for nn = 1:numFiles
%         plot([nn nn], [0 max(clusterSize{nn})],'--k');
%         xx = nn*ones(size(clusterSize{nn},1),1);
%         scatter(xx,clusterSize{nn},200,meanClusterColor{nn},'filled');
%     end
%     
%     ylim([0,70]);
%     hold off
    
    set(gca,'XTick',1:numFiles);
    xticklabels(sampleName);
    xtickangle(90);
    
    figure(h2),
    scatter(1:19,autoDiversity(:,3),50,myColor(count,:),'filled');
    title([num2str(jj),' ','theoretical HSCs']);
    hold on
    
    count = count + 1;
    
    % plotSpread(autoDiversity(1:4,1),autoDiversity(5:19,1));
    % title('Auto')
    % figure;
    % plotSpread(Diversity_2015_12_2(:,1));
    % title('Manual')
end

hold off