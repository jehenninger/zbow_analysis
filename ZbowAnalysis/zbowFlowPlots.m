function [h] = zbowFlowPlots(fileName)
%UNTITLED Summary of this function goes here
%   This function will read fcs data and make color plots that you would
%   see from FACS Diva.

%Load and read data
% try
    [fcsdat, fcshdr] = fca_readfcs(fileName);
    
    
    channelNames = {fcshdr.par.name};
    
    %For Aria
    red = fcsdat(:,find(strcmp('FJComp-PE-A',channelNames)));
    green = fcsdat(:,find(strcmp('FJComp-FITC-A',channelNames)));
    blue = fcsdat(:,find(strcmp('FJComp-mCFP-A',channelNames)));
    
    
    sessionData = [red green blue];
    linData = [red green blue];
    
    
    %Default logicle values
    T = 2^18;
    M = 4.5;
    A = 0;
    
    %Default width parameters
    clustRedW = (M - log10(T/abs(min(sessionData(:,1)))))/2;
    clustGreenW = (M - log10(T/abs(min(sessionData(:,2)))))/2;
    clustBlueW = (M - log10(T/abs(min(sessionData(:,2)))))/2;
    
    if clustRedW < 0
        clustRedW = 0;
    end
    
    if clustGreenW < 0
        clustGreenW = 0;
    end
    
    if clustBlueW < 0
        clustBlueW = 0;
    end
    
    %Custom width parameters for Zebrabow
    ternRedW = 1.5;
    ternGreenW = 1.75;
    ternBlueW = 1.75;
    
    %Logicle transformed data with default widths
    
    logDataDefault = [logicleTransform(sessionData(:,1),T,clustRedW,M,A),...
        logicleTransform(sessionData(:,2),T,clustGreenW,M,A),...
        logicleTransform(sessionData(:,3),T,clustBlueW,M,A)];
    
    
    
    
    %Logicle transformed data with custom widths
    logDataTern = [logicleTransform(sessionData(:,1),T,ternRedW,M,A),...
        logicleTransform(sessionData(:,2),T,ternGreenW,M,A),...
        logicleTransform(sessionData(:,3),T,ternBlueW,M,A)];
    
    
    [normDataTern, ternCoords] = myDataProcessing(logDataTern);
    
    [normDataDefault, ternCoordsLog] = myDataProcessing(logDataDefault);
    
    h = figure;
    bins = round(sqrt(length(normDataDefault(:,1))));
    
    %Red vs Green Plot
    subplot(1,3,1)
    scatter(normDataDefault(:,1),normDataDefault(:,2),10,normDataTern,'filled');
    hold on
    dscatter(normDataDefault(:,1),normDataDefault(:,2),'BINS',[bins bins],'PLOTTYPE','contour');
    xlim([0 1]), ylim([0 1]);
    xlabel('Red'),ylabel('Green');
    axis('equal')
    hold off
    
    %Red vs Blue Plot
    subplot(1,3,2)
    scatter(normDataDefault(:,1),normDataDefault(:,3),10,normDataTern,'filled');
    hold on
    dscatter(normDataDefault(:,1),normDataDefault(:,3),'BINS',[bins bins],'PLOTTYPE','contour');
    xlim([0 1]), ylim([0 1]);
    xlabel('Red'),ylabel('Blue');
    axis('equal')
    hold off
    
    %Green vs Blue Plot
    subplot(1,3,3)
    scatter(normDataDefault(:,2),normDataDefault(:,3),10,normDataTern,'filled');
    hold on
    dscatter(normDataDefault(:,2),normDataDefault(:,3),'BINS',[bins bins],'PLOTTYPE','contour');
    xlim([0 1]), ylim([0 1]);
    xlabel('Green'),ylabel('Blue');
    axis('equal')
    hold off
    
% catch
%     msgbox('Failed to load file');
% end


    function [normDataTern,ternCoords] = myDataProcessing(logDataTern)
    %Normalized data
    normDataTern(:,1) = myNorm(logDataTern(:,1));
    normDataTern(:,2) = myNorm(logDataTern(:,2));
    normDataTern(:,3) = myNorm(logDataTern(:,3));
    
    
    %Relative contributions of each color
    relData = myColorRatio(normDataTern(:,1),normDataTern(:,2),normDataTern(:,3));
    R = relData(:,1);
    G = relData(:,2);
    
    %Generate 2D Ternary Graph Coordinates
    y = G*sin(pi/3);
    x = R + y*cot(pi/3);
    
    ternCoords = [x, y];
    end
end

