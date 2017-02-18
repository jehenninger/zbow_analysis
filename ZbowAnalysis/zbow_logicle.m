function [sessionData,normDataTern,colorm,ternColor,ternCoords,sampleName] = zbow_logicle(file,cellType,sampleSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('file','var') || isempty(file)
    [file,path] = uigetfile('E:\zon_lab\FACS\*.fcs');
    file = fullfile(path,file);
end

if ~exist('cellType','var') || isempty(cellType)
    
    cellType = 'Minimal';
    
end

if ~exist('sampleSize','var')
    sampleSize = 20000;
end

[~,sampleName,~] = fileparts(file);
[fcsdat, fcshdr] = fca_readfcs(file);

if fcsdat == 0
    sessionData = NaN;
    normDataTern = NaN;
    colorm = NaN;
    ternColor = NaN;
    ternCoords = NaN;
else
    channelNames = {fcshdr.par.name};
    
    
    
    
    %     %For Aria
    %     red = fcsdat(:,find(strcmp('PE-A',channelNames)));
    %     green = fcsdat(:,find(strcmp('FITC-A',channelNames)));
    %     blue = fcsdat(:,find(strcmp('CFP-A',channelNames)));
    %     sessionData = [red green blue];
    
    cambridgeTest = strfind(channelNames,'FJComp-FITC-A');
    cambridgeTest = [cambridgeTest{:}];
    childrensTest = strfind(channelNames,'FJComp-GFP-A');
    childrensTest = [childrensTest{:}];
    
    if cambridgeTest == 1
        %For Aria
        red = fcsdat(:,find(strcmp('FJComp-PE-A',channelNames)));
        green = fcsdat(:,find(strcmp('FJComp-FITC-A',channelNames)));
        blue = fcsdat(:,find(strcmp('FJComp-CFP-A',channelNames)));
        sessionData = [red green blue];
        
    elseif childrensTest == 1
        %For Children's Aria (self-sorter)
        red = fcsdat(:,find(strcmp('FJComp-DsRed-A',channelNames)));
        green = fcsdat(:,find(strcmp('FJComp-GFP-A',channelNames)));
        blue = fcsdat(:,find(strcmp('FJComp-DAPI-A',channelNames)));
        sessionData = [red green blue];
    else
        error('Could not identify fluorescence parameters...');
    end
    
    % %For Fortessa un-compensated
    % red = fcsdat(:,find(strcmp('PE-A',channelNames)));
    % green = fcsdat(:,find(strcmp('FITC-A',channelNames)));
    % blue = fcsdat(:,find(strcmp('mCFP-A',channelNames)));
    % sessionData = [red green blue];
    
    % %For Fortessa
    % red = fcsdat(:,find(strcmp('FJComp-PE-A',channelNames)));
    % green = fcsdat(:,find(strcmp('FJComp-FITC-A',channelNames)));
    % blue = fcsdat(:,find(strcmp('FJComp-mCFP-A',channelNames)));
    % sessionData = [red green blue];
    
    
    % Rough process outliers by removing above 10^5 and removing below 5000
    % lowIndex = sessionData(:,1)<3000 & sessionData(:,2)<2000 & sessionData(:,3)<2000;
    % lowIndex = sessionData(:,1)<500 & sessionData(:,2)<500 & sessionData(:,3)<500;
    
    switch cellType
        % For myeloid
        case 'Myeloid'
            sessionData(sessionData(:,1)<5000 & sessionData(:,2)<5000 & sessionData(:,3)<5000,:) = [];
            sessionData(sessionData(:,1)>1.5*10e5 |...
                sessionData(:,2) > 1.5*10e5 |...
                sessionData(:,3) > 1.5*10e5, :) = [];
            
        case 'Erythroid'
            sessionData(sessionData(:,1)<500 & sessionData(:,2)<500 & sessionData(:,3)<500,:) = [];
            sessionData(sessionData(:,1)>1.5*10e3 |...
                sessionData(:,2) > 1.5*10e3 |...
                sessionData(:,3) > 1.5*10e3, :) = [];
            
        case 'Minimal'
            sessionData(sessionData(:,1)<1000 & sessionData(:,2)<1000 & sessionData(:,3)<1000,:) = [];
            sessionData(sessionData(:,1)>2.0*10e5 |...
                sessionData(:,2) > 2.0*10e5 |...
                sessionData(:,3) > 2.0*10e5, :) = [];
            
            % sessionData(lowIndex==1,:) = [];
            % sessionData(highIndex,:) = [];
    end
    
    if sampleSize < length(sessionData(:,1))
        sessionData = datasample(sessionData,sampleSize,'Replace',false);
    end
    
    T = 2^18;
    M = 4.5;
    A = 0;
    
    %%Zbow parameters
    redW = 1.5;
    greenW = 1.75;
    blueW = 1.75;
    logColor = [logicleTransform(sessionData(:,1), T, redW, M, A),...
        logicleTransform(sessionData(:,2),T,greenW,M,A),...
        logicleTransform(sessionData(:,3),T,blueW,M,A)];
    
    colorm(:,1) = (logColor(:,1)-min(logColor(:,1)))./(max(logColor(:,1))-min(logColor(:,1)));
    colorm(:,2) = (logColor(:,2)-min(logColor(:,2)))./(max(logColor(:,2))-min(logColor(:,2)));
    colorm(:,3) = (logColor(:,3)-min(logColor(:,3)))./(max(logColor(:,3))-min(logColor(:,3)));
    
    % Auto width parameters
    redW = (M - log10(T/abs(min(sessionData(:,1)))))/2;
    greenW = (M - log10(T/abs(min(sessionData(:,2)))))/2;
    blueW = (M - log10(T/abs(min(sessionData(:,2)))))/2;
    
    if redW < 0
        redW = 0;
    end
    
    if greenW < 0
        greenW = 0;
    end
    if blueW < 0
        blueW = 0;
    end
    
    logData = [logicleTransform(sessionData(:,1), T, redW, M, A),...
        logicleTransform(sessionData(:,2),T,greenW,M,A),...
        logicleTransform(sessionData(:,3),T,blueW,M,A)];
    
    normDataTern(:,1) = (logData(:,1)-min(logData(:,1)))./(max(logData(:,1))-min(logData(:,1)));
    normDataTern(:,2) = (logData(:,2)-min(logData(:,2)))./(max(logData(:,2))-min(logData(:,2)));
    normDataTern(:,3) = (logData(:,3)-min(logData(:,3)))./(max(logData(:,3))-min(logData(:,3)));
    
    % R = normDataTern(:,1)./(normDataTern(:,1)+ normDataTern(:,2)+ normDataTern(:,3));
    % G = normDataTern(:,2)./(normDataTern(:,1)+ normDataTern(:,2)+ normDataTern(:,3));
    
    totalColor = colorm(:,1) + colorm(:,2) + colorm(:,3);
    R = colorm(:,1)./(totalColor(:));
    G = colorm(:,2)./(totalColor(:));
    B = colorm(:,3)./(totalColor(:));
    
    %Generate 2D ternary graph coordinates
    y = G*sin(pi/3);
    x = R + y*cot(pi/3);
    
    ternCoords = [x,y];
    ternColor = [R G B];
    
end

end

