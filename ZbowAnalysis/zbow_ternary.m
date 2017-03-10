function [handle3D, handleTern, densityHandle] = zbow_ternary(fileName,sampleSize, cellType, plotTypeChoice,combinedChoice)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[pathName, sampleName, ext] = fileparts(fileName);

if isequal(ext,char('.fcs'))
    %Load and read data
    [fcsdat, fcshdr] = fca_readfcs(fileName);
    channelNames = {fcshdr.par.name};
    
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
    case 'None'
        
      
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

R = colorm(:,1)./(colorm(:,1)+ colorm(:,2)+ colorm(:,3));
G = colorm(:,2)./(colorm(:,1)+ colorm(:,2)+ colorm(:,3));

%Generate 2D ternary graph coordinates
y = G*sin(pi/3);
x = R + y*cot(pi/3);

% switch combinedChoice
%     case 'Yes'
%         switch plotTypeChoice
%             case 'Both'
%                 figure,
%                 scrollsubplot(
%             case 'Ternary'
%             case '3D Scatter'
%         end
%
%     case 'No'
%         figure,
% end

figure,

switch plotTypeChoice
    case 'Both'
        
        set(gcf,'Units','normalized','Position',[0.5, 0.5, 0.25, 0.16]);
        a = subplot(1,2,1);
        handle3D = scatter3(normDataTern(:,1),...
            normDataTern(:,2),...
            normDataTern(:,3),...
            20,colorm,'filled');
        xlabel('Red'),ylabel('Green'),zlabel('Blue');
        view(135,30);
        
        subplot(1,2,2);
        handleTern = ternPlot([x, y], colorm,'false','true',10,'false');
        b = title(a,[sampleName, ' ', 'n = ',num2str(length(normDataTern(:,1)))]);
        set(b,'interpreter','none');
        
        %         binNumber = 500;
        %
        %         hold on
        %         if length(normDataTern(:,1))>500
        %             densityHandle = dscatter(x,y,...
        %                 'PLOTTYPE','contour','BINS',[binNumber binNumber]);
        %         else
        %             densityHandle = 'Sample size too low';
        %         end
        %
        %         hold off
        
    case 'Ternary'
        
        set(gcf,'Units','normalized','Position',[0.5, 0.5, 0.25, 0.25]);
        handleTern = ternPlot([x, y], colorm,'true','true',10);
        b = title([sampleName, ' ', 'n = ',num2str(length(normDataTern(:,1)))]);
        set(b,'interpreter','none');
        set(findall(gca,'type','text'),'visible','on');
        
        handle3D = 0;
                
        %         hold on
        %         if length(normDataTern(:,1))>500
        %             densityHandle = dscatter(x,y,...
        %                 'PLOTTYPE','contour','BINS',[binNumber binNumber]);
        %         else
        %             densityHandle = 'Sample size too low';
        %         end
        %
        %         hold off
        
    case '3D Scatter'
        
        set(gcf,'Units','normalized','Position',[0.5, 0.5, 0.25, 0.25]);
        
        handle3D = scatter3(normDataTern(:,1),...
            normDataTern(:,2),...
            normDataTern(:,3),...
            20,colorm,'filled');
        xlabel('Red'),ylabel('Green'),zlabel('Blue');
        view(135,30);
        b = title([sampleName,' ', 'n = ', num2str(length(normDataTern(:,1)))]);
        set(b, 'interpreter','none');
        set(findall(gca,'type','text'),'visible','on');
        
        handleTern = 0;
        densityHandle = 0;
        
        hold off
end


end

