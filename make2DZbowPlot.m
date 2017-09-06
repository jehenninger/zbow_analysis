function [hFig, hPlot, file] = make2DZbowPlot(fileInput,sampleSize, markerSize, alpha, visible)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~exist('fileInput','var') || isempty(fileInput)
    error('Could not find or load file...');
end

if ~exist('sampleSize','var')|| isempty(sampleSize)
    sampleSize = 20000;
end

if ~exist('markerSize','var') || isempty(markerSize)
    markerSize = 15;
end

if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end

file = fileInput;

[cellColor,~,sampleName,idx] = flowTransformCustom(file);
sampleName = strrep(sampleName,'_',' ');


cellColor = cellColor(:,idx);
cellColor = normalize_var(cellColor,0,1);

if sampleSize < size(cellColor,1)
    [~, sampleIdx] = datasample(cellColor,sampleSize,'Replace', false);
    cellColorSample = cellColor(sampleIdx,:);
else
    cellColorSample = cellColor;
end

[~, ternCoords] = ternaryProcessing(cellColorSample);

switch visible
    case 'on'
        hFig = figure('Units','normalized','Position',[0.5 0 0.5 1], 'Visible','on');
    case 'off'
        hFig = figure('Units','normalized','Position',[0.5 0 0.5 1], 'Visible','off');
end

hPlot = ternPlot(ternCoords, cellColorSample, 'false', 'true', markerSize, 'true', alpha);

title(sampleName);

    function [relData,ternCoords] = ternaryProcessing(normData)
        
        %Relative contributions of each color
        relData = myColorRatio(normData(:,1),normData(:,2),normData(:,3));
        R = relData(:,1);
        G = relData(:,2);
        
        %Generate 2D Ternary Graph Coordinates
        b = G*sin(pi/3);
        a = R + b*cot(pi/3);
        
        ternCoords = [a, b];
    end

end
