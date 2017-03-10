
function varargout = zbow_analysisV2(varargin)
% ZBOW_ANALYSISV2 MATLAB code for zbow_analysisV2.fig
%      ZBOW_ANALYSISV2, by itself, creates a new ZBOW_ANALYSISV2 or raises the existing
%      singleton*.
%
%      H = ZBOW_ANALYSISV2 returns the handle to a new ZBOW_ANALYSISV2 or the handle to
%      the existing singleton*.
%
%      ZBOW_ANALYSISV2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZBOW_ANALYSISV2.M with the given input arguments.
%
%      ZBOW_ANALYSISV2('Property','Value',...) creates a new ZBOW_ANALYSISV2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before zbow_analysisV2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to zbow_analysisV2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help zbow_analysisV2

% Last Modified by GUIDE v2.5 06-Dec-2016 22:23:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @zbow_analysisV2_OpeningFcn, ...
    'gui_OutputFcn',  @zbow_analysisV2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before zbow_analysisV2 is made visible.
function zbow_analysisV2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to zbow_analysisV2 (see VARARGIN)

% Choose default command line output for zbow_analysisV2
handles.output = hObject;

handles.T = 2^18;
handles.M = 4.5;
handles.A = 0;

handles.scatterFigHandle = [];
handles.ternFigHandle = [];
handles.scatterPlotHandle = [];
handles.ternPlotHandle = [];
handles.densityOverlay = [];

handles.silCount = 0;
handles.clusterCount = 0;

handles.clusterEvalCount = 0;

annotation('textbox','FontSize',15,'LineStyle','none','Units','normalized',...
    'Position',[0.72 0.162 0.04 0.032],'String','\rho',...
    'FontUnits','normalized','FontWeight','bold');
annotation('textbox','FontSize',15,'LineStyle','none','Units','normalized',...
    'Position',[0.419 0.649 0.033 0.045],'String','\delta',...
    'FontUnits','normalized','FontWeight','bold');

set(gcf,'units','normalized','OuterPosition', [0.01, 0.5, 0.5, 0.6]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes zbow_analysisV2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = zbow_analysisV2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
% hObject    handle to loadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.sortTab = [1 1 1];
handles.meanValues = [0.8 0.8 0.8];


% try
%% Load FSC file
[fileName, pathName]  = uigetfile('E:\zon_lab\FACS\*.fcs','Select FCS file');
handles.pathName = pathName;
[~, handles.sampleName, ~] = fileparts(fileName);

[handles.dataLin, handles.dataDef, handles.dataCus,handles.fcsdat,handles.channelNames] = loadFCSFile(...
    fullfile(pathName,fileName), handles);

%Show sample size
numOfCells = size(handles.dataLin,1);
set(handles.originalSampleSize,'String',num2str(numOfCells));


sampleSize = str2double(handles.sampleSizeCluster.String);
if size(handles.dataLin,1)> sampleSize
    [~, sampleIdx] = datasample(handles.dataLin,sampleSize,'Replace',false);
    handles.dataLin = handles.dataLin(sampleIdx,:);
    handles.dataDef = handles.dataDef(sampleIdx,:);
    handles.dataCus = handles.dataCus(sampleIdx,:);
    handles.fcsdat = handles.fcsdat(sampleIdx,:);
end

handles.clusterIndex = ones(size(handles.dataLin,1),1);

handles.normDataLin = myDataProcessing(handles.dataLin);
handles.normDataDef = myDataProcessing(handles.dataDef);
handles.normDataCus = myDataProcessing(handles.dataCus);


%Update file name text box
set(handles.fileNameText,'String',fileName);

%% Make new figures/axes for scatter and ternary plots
handles.scatterFigHandle = figure('Units','normalized','OuterPosition',[0.51 0.405 0.35 0.57]);
handles.ternFigHandle = figure('Units','normalized','OuterPosition',[0.47 0.046 0.42 0.36]);

handles = updatePlots(handles);

guidata(hObject,handles);


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% try

saveFolder = uigetdir(handles.pathName);
sampleName = handles.sampleName;
% [saveFileName,savePathName] = uiputfile([handles.pathName,'\*.xlsx'],...
%     'Choose save file...',[handles.sampleName,'_Summary.xlsx']);

% [~, sampleName,~] = fileparts(fullfile(savePathName,saveFileName));

hWait = waitbar(0,'Saving... Please wait...');
%Decision Graph
if isfield(handles,'decisionGraph')
    % ff = figure('PaperUnits','inches','PaperPositionMode','manual',...
    %     'PaperPosition',[0 0 10 7],'Visible','off');
    ff = figure('units','normalized','outerposition',[0 0 0.5 1],'Visible','off');
    new_handle = copyobj(handles.decisionGraph,ff);
    set(new_handle,'Units','normalized','outerposition',[0.01 0.01 0.98 0.98]);
    print(ff,fullfile(saveFolder, [sampleName, '_decision_graph.eps']),'-painters','-depsc','-r300');
    print(ff,fullfile(saveFolder, [sampleName, '_decision_graph.jpg']),'-opengl','-djpeg','-r300');
    close gcf
end
waitbar(1/4);

%Ternary Graph
% ff = figure('PaperUnits','inches','PaperPositionMode','manual',...
%     'PaperPosition',[0 0 10 7],'Visible','off');
ff = figure('units','normalized','outerposition',[0 0 0.5 1],'Visible','off');
new_handle = copyobj(handles.sub1,ff);
set(new_handle,'Units','normalized','outerposition',[0.01 0.01 0.9 0.9]);
print(ff,fullfile(saveFolder, [sampleName, '_ternary_graph.eps']),'-painters','-depsc','-r300');
print(ff,fullfile(saveFolder, [sampleName, '_ternary_graph.jpg']),'-opengl','-djpeg','-r300');
close gcf
waitbar(2/4);

%Bar Graph
% ff = figure('PaperUnits','inches','PaperPositionMode','manual',...
%     'PaperPosition',[0 0 7 10],'Visible','off');
ff = figure('units','normalized','outerposition',[0 0 0.5 1],'Visible','off');
new_handle = copyobj(handles.sub2,ff);
set(new_handle,'Units','normalized','outerposition',[0.01 0.01 0.98 0.98]);
print(ff,fullfile(saveFolder, [sampleName, '_bar_graph.eps']),'-painters','-depsc','-r300');
print(ff,fullfile(saveFolder, [sampleName, '_bar_graph.jpg']),'-opengl','-djpeg','-r300');
close gcf
waitbar(3/4);

%Silhouette plot (if it exists)
if isfield(handles,'silFigHandle')
    if isvalid(handles.silFigHandle)
        ff = figure('units','normalized','outerposition',[0 0 0.5 1],'Visible','off');
        new_handle = copyobj(get(handles.silPlotHandle,'Children'),ff);
        set(new_handle,'Units','normalized','outerposition',[0.01 0.01 0.98 0.98]);
        print(ff,fullfile(saveFolder, [sampleName, '_sil_graph.eps']),'-painters','-depsc','-r300');
        print(ff,fullfile(saveFolder, [sampleName, '_sil_graph.jpg']),'-opengl','-djpeg','-r300');
        close gcf
    end
end


%Save excel file of tabulated clusters with color
if isfield(handles,'excelOutputSummary')
    excelOutputSummary = handles.excelOutputSummary;
    excelOutputSummary = num2cell(excelOutputSummary);
    excelOutputSummaryHeader = {'Clone ID','Num of Cells','Percentage of Total','Red','Green','Blue'};
    excelOutputSummary = [excelOutputSummaryHeader; excelOutputSummary];
    xlswrite(fullfile(saveFolder,[sampleName,'_Summary.xlsx']),excelOutputSummary);
    waitbar(4/4);
    excelOutputRaw = handles.excelOutputRaw;
    excelOutputRaw = num2cell(excelOutputRaw);
    % excelOutputRawHeader = {'RAW RED', 'RAW GREEN', 'RAW BLUE', 'TRANSFORM RED', 'TRANSFORM GREEN', 'TRANSFORM BLUE',...
    %     'NORMAL RED','NORMAL GREEN','NORMAL BLUE','COLOR RED','COLOR GREEN','COLOR BLUE',...
    %     'X-COORD','Y-COORD','CLUSTER ID'};
    excelOutputRawHeader = {'Norm Linear Red', 'Norm Linear Green', 'Norm Linear Blue',...
        'Norm Default Red', 'Norm Default Green', 'Norm Default Blue',...
        'Norm Custom Red','Norm Custom Green','Norm Custom Blue',...
        'Ternary X','Ternary Y','Cluster ID'};
    excelOutputRaw = [excelOutputRawHeader; excelOutputRaw];
    xlswrite(fullfile(saveFolder,[sampleName,'_RawValues.xlsx']),excelOutputRaw);
end
close(hWait);


% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% arrayfun(@cla,findall(0,'type','axes'))
% set(handles.clusterList,'Value',[],'String',[]);
% clear

if ishandle(handles.scatterFigHandle)
    if isvalid(handles.scatterFigHandle)
        close(handles.scatterFigHandle);
    end
end

if ishandle(handles.ternFigHandle)
    if isvalid(handles.ternFigHandle)
        close(handles.ternFigHandle)
    end
end

close(gcbf);
zbow_analysisV2;



function sampleSizeCluster_Callback(hObject, eventdata, handles)
% hObject    handle to sampleSizeCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function sampleSizeCluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleSizeCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','20000');

% --- Executes on button press in removeOutliers.
function removeOutliers_Callback(hObject, eventdata, handles)
% hObject    handle to removeOutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.scatterDisplayType.Value
    case 1
        plotDataScatter = handles.normDataCus;
    case 2
        plotDataScatter = handles.normDataDef;
    case 3
        plotDataScatter = handles.normDataLin;
end


idx = removeOutliers(plotDataScatter,handles.scatterPlotHandle);

handles.dataDef = handles.dataDef(idx,:);
handles.dataCus = handles.dataCus(idx,:);
handles.dataLin = handles.dataLin(idx,:);
handles.fcsdat = handles.fcsdat(idx,:);

handles.normDataLin = myDataProcessing(handles.dataLin);
handles.normDataDef = myDataProcessing(handles.dataDef);
handles.normDataCus = myDataProcessing(handles.dataCus);

newSampleSize = size(handles.normDataLin,1);
handles.outlierSampleSize.String = num2str(newSampleSize);

handles.clusterIndex = ones(newSampleSize,1);

handles = updatePlots(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in makeDecisionGraph.
function makeDecisionGraph_Callback(hObject, eventdata, handles)
% hObject    handle to makeDecisionGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.clusterCount = 0;
hWait = waitbar(0,'Please wait...');

set(handles.clusterList,'String', []);

cla(handles.decisionGraph,'reset');

% Choose transform of data

switch handles.transformType.Value
    case 1
        transformData = handles.normDataCus;
    case 2
        transformData = handles.normDataDef;
    case 3
        transformData = handles.normDataLin;
end

waitbar(.2);

switch handles.clusteringMethod.Value
    case 1
        inputData = transformData;
    case 2
        [~,~,inputData] = myDataProcessing(transformData);
end

[rho, delta, nneigh] = deltarho(inputData, 1);
handles.rho = rho;
handles.delta = delta;
handles.nneigh = nneigh;

waitbar(.8);

%% Normalize and log2 scale delta and rho for plotting

modrho = rho;
modrho(modrho>2) = log2(modrho(modrho>2));
normModRho = modrho./max(modrho(:));

moddelta = 10*delta;
moddelta(moddelta>2) = log2(moddelta(moddelta>2));
normModDelta = moddelta./max(moddelta(:));

scale = 600;

markerSize = normModRho.*normModDelta'*scale;
markerSize(markerSize<25)=30;

handles.modrho = normModRho;
handles.moddelta = normModDelta;

scatter(handles.decisionGraph,normModRho(:),normModDelta(:),...
    markerSize,handles.normDataCus,'filled');
title ('Decision Graph','FontName','Arial','FontSize',15.0)
xlim(handles.decisionGraph,[-0.1 1.1]);
ylim(handles.decisionGraph,[-0.1 1.1]);
xlabel ('\rho')
ylabel ('\delta')

drawnow;
close(hWait);

guidata(hObject,handles);

% --- Executes on selection change in clusterList.
function clusterList_Callback(hObject, eventdata, handles)
% hObject    handle to clusterList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clusterList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusterList


% --- Executes during object creation, after setting all properties.
function clusterList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusterList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addCenter.
function addCenter_Callback(hObject, eventdata, handles)
% hObject    handle to addCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load data that we will need
modrho = handles.modrho;
moddelta = handles.moddelta;

%Draw rectangle on plot and get coordinates
clusterRect = imrect(handles.decisionGraph);

clusterRectCoord = getPosition(clusterRect);
clusterRectX = [clusterRectCoord(1);...
    clusterRectCoord(1);...
    clusterRectCoord(1)+clusterRectCoord(3);...
    clusterRectCoord(1)+clusterRectCoord(3)];
clusterRectY = [clusterRectCoord(2);...
    clusterRectCoord(2)+clusterRectCoord(4);...
    clusterRectCoord(2)+clusterRectCoord(4);...
    clusterRectCoord(2)];

% Identify cluster centers as points within rectangle
cluster_center = inpoly([modrho,moddelta.'],[clusterRectX,clusterRectY]);

% Test if cluster centers already exist, and if so, add to the list
if handles.clusterCount == 0
    clusterCenterIndex = unique(find(cluster_center == 1));
    handles.clusterCenterIndex = clusterCenterIndex;
    handles.clusterCount = numel(handles.clusterCenterIndex);
else
    oldClusterCenterIndex = handles.clusterCenterIndex;
    newClusterCenterIndex = find(cluster_center == 1);
    clusterCenterIndex = [oldClusterCenterIndex; newClusterCenterIndex];
    clusterCenterIndex = unique(clusterCenterIndex);
    handles.clusterCenterIndex = clusterCenterIndex;
    handles.clusterCount = numel(clusterCenterIndex);
    clusterText = handles.clusterText;
    for k = 1:numel(clusterText)
        delete(clusterText{k});
    end
    
end

%Assign clusters by nearest neighbors
clusterIndex = assign_cluster(1, handles.nneigh, handles.clusterCenterIndex);
handles.clusterIndex = clusterIndex;

% Place text of cluster number to the right of the cluster center
clusterString = cell(handles.clusterCount,1);
for k = 1:handles.clusterCount
    clusterText{k} = text(modrho(clusterCenterIndex(k))+.035,moddelta(clusterCenterIndex(k)),num2str(k),...
        'FontUnits','normalized',...
        'FontWeight','bold','HorizontalAlignment','center');
    clusterString{k} = [num2str(k),' ','(',num2str(sum(clusterIndex == k)),')'];
end

set(handles.clusterList,'String', clusterString);
handles.clusterText = clusterText;
delete(clusterRect);

handles = updatePlots(handles);

guidata(hObject,handles);


function [figure_handle, plot_handle, density_handle, bar_handle,subhandle1,subhandle2] = initiateTernPlot(...
    data,color,sortTab,meanValues,previous_handle)

if ~exist('sortTab','var')
    sortTab = [1 1 1];
elseif isempty(sortTab)
    sortTab = [1 1 1];
end

if ~exist('meanValues','var')
    meanValues = [0.8 0.8 0.8];
elseif isempty(meanValues)
    meanValues = [0.8 0.8 0.8];
end

if ~exist('previous_handle','var')
    figure_handle = figure('Units','normalized','OuterPosition',[0.47 0.046 0.42 0.36]);
elseif ~isvalid(previous_handle)
    figure_handle = figure('Units','normalized','OuterPosition',[0.47 0.046 0.42 0.36]);
else
    figure_handle = previous_handle;
    figure(figure_handle);
end


[~,~,ternCoords] = myDataProcessing(data);
subhandle1 = subplot(1,5,1:4);
[plot_handle, density_handle] = ternPlot(ternCoords, color,'true','true');


subhandle2 = subplot(1,5,5);
bar_handle = bar([transpose(sortTab(:,3));nan(size(transpose(sortTab(:,3))))],'Stacked');
hold on
for ii = 1:numel(sortTab(:,3))
    set(bar_handle(ii),'FaceColor',meanValues(sortTab(ii,1),:));
end
ylim([0 100]);
xlim([0 2]);
set(gca,'ygrid','on');
set(gca,'xtick',1);
set(gca,'FontSize',20,'FontName','Arial');
set(gca,'XTickLabel',' ');
title([num2str(size(sortTab,1)), ' clusters']);
hold off


% --- Executes on button press in removeCenter.
function removeCenter_Callback(hObject, eventdata, handles)
% hObject    handle to removeCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% try

modrho = handles.modrho;
moddelta = handles.moddelta;
nneigh = handles.nneigh;
clusterCenterIndex = handles.clusterCenterIndex;


clusterToDelete = get(handles.clusterList,'Value');
set(handles.clusterList,'Value',1);


clusterCenterIndex(clusterToDelete) = [];
handles.clusterCenterIndex = clusterCenterIndex;

clusterCount = numel(clusterCenterIndex);
handles.clusterCount = clusterCount;

clusterText = handles.clusterText;
for k = 1:numel(clusterText)
    delete(clusterText{k});
end

clusterIndex = assign_cluster(1, nneigh, clusterCenterIndex);
handles.clusterIndex = clusterIndex;

clusterString = cell(clusterCount);
for k = 1:clusterCount
    clusterText{k} = text(modrho(clusterCenterIndex(k))+.035,moddelta(clusterCenterIndex(k)),num2str(k),...
        'FontUnits','normalized',...
        'FontWeight','bold','HorizontalAlignment','center','Parent',handles.decisionGraph);
    clusterString{k} = [num2str(k),' ','(',num2str(sum(clusterIndex == k)),')'];
end

set(handles.clusterList,'String', clusterString);
handles.clusterText = clusterText;

handles = updatePlots(handles);

guidata(hObject,handles);


% --- Executes on button press in highlightCluster.
function highlightCluster_Callback(hObject, eventdata, handles)
% hObject    handle to highlightCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clusterIndex = handles.clusterIndex;
normDataTern = handles.normDataCus;
normData = handles.normDataDef;
normDataLin = handles.normDataLin;

clusterSelection = get(handles.clusterList,'Value');

colorControl = get(handles.ternColorType,'Value');

highColorTern = 0.9*ones(length(clusterIndex),3);
switch colorControl
    case 1
        highColorTern(clusterIndex == clusterSelection,:)...
            = normDataTern(clusterIndex == clusterSelection,:);
    case 2
        highColorTern(clusterIndex == clusterSelection,:)...
            = normData(clusterIndex == clusterSelection,:);
    case 3
        highColorTern(clusterIndex == clusterSelection,:)...
            = normDataLin(clusterIndex == clusterSelection,:);
    case 4
        highColorTern(clusterIndex == clusterSelection,:)...
            = normDataTern(clusterIndex == clusterSelection,:);
        
end

color3DControl = get(handles.scatterColorType,'Value');

highColorScatter = 0.9*ones(length(clusterIndex),3);
switch color3DControl
    case 1
        highColorScatter(clusterIndex == clusterSelection,:)...
            = normDataTern(clusterIndex == clusterSelection,:);
    case 2
        highColorScatter(clusterIndex == clusterSelection,:)...
            = normData(clusterIndex == clusterSelection,:);
    case 3
        highColorScatter(clusterIndex == clusterSelection,:)...
            = normDataLin(clusterIndex == clusterSelection,:);
    case 4
        highColorScatter(clusterIndex == clusterSelection,:)...
            = normDataTern(clusterIndex == clusterSelection,:);
end

handles = updatePlots(handles,highColorTern,highColorScatter);

guidata(hObject,handles);



% --- Executes on selection change in ternColorType.
function ternColorType_Callback(hObject, eventdata, handles)
% hObject    handle to ternColorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ternColorType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ternColorType

handles = updatePlots(handles);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ternColorType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ternColorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in transformType.
function transformType_Callback(hObject, eventdata, handles)
% hObject    handle to transformType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns transformType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from transformType


% --- Executes during object creation, after setting all properties.
function transformType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transformType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in clusteringMethod.
function clusteringMethod_Callback(hObject, eventdata, handles)
% hObject    handle to clusteringMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clusteringMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusteringMethod


% --- Executes during object creation, after setting all properties.
function clusteringMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusteringMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',2);



% --- Executes on selection change in scatterColorType.
function scatterColorType_Callback(hObject, eventdata, handles)
% hObject    handle to scatterColorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns scatterColorType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scatterColorType

handles = updatePlots(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function scatterColorType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scatterColorType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in scatterDisplayType.
function scatterDisplayType_Callback(hObject, eventdata, handles)
% hObject    handle to scatterDisplayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns scatterDisplayType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scatterDisplayType

handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function scatterDisplayType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scatterDisplayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ternDisplayType.
function ternDisplayType_Callback(hObject, eventdata, handles)
% hObject    handle to ternDisplayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ternDisplayType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ternDisplayType


handles = updatePlots(handles);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ternDisplayType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ternDisplayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in clusterDataCopy.
function clusterDataCopy_Callback(hObject, eventdata, handles)
% hObject    handle to clusterDataCopy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clusterIndex = handles.clusterIndex;
clusterNumber = get(handles.clusterList,'Value');

numOfCells = sum(clusterIndex == clusterNumber);

percentOfBlood = 100*numOfCells./numel(clusterIndex);
output = [clusterNumber, numOfCells,percentOfBlood];
num2clip(output);

function [linData, logDataDefault, logDataTern, fcsdat, channelNames] = loadFCSFile(file, handles)
%Load and read data
[fcsdat, fcshdr] = fca_readfcs(file);
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

linData = [red green blue];

%Default logicle values
T = handles.T;
M = handles.M;
A = handles.A;

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

function [figure_handle, plot_handle] = initiateScatter(data, color, previous_handle,currentView)


% color = myColorProcessing(color,1);

if  ~exist('previous_handle','var')
    figure_handle = figure('Units','normalized','OuterPosition', [0.51 0.405 0.35 0.57]);
elseif ~isvalid(previous_handle)
    figure_handle = figure('Units','normalized','OuterPosition', [0.51 0.405 0.35 0.57]);
else
    figure_handle = previous_handle;
    figure(figure_handle);
end


plot_handle = scatter3(data(:,1),...
    data(:,2),...
    data(:,3),...
    8,color,'filled');

xlabel('Red'),ylabel('Green'),zlabel('Blue');
if ~exist('currentView','var')
    view(135,30);
else
    view(currentView(1),currentView(2));
end
axis('equal');


title('3D Scatter');


function [data] = myColorProcessing(data,gamma)
%Regular normalization
data(:,1) = myNorm(data(:,1));
data(:,2) = myNorm(data(:,2));
data(:,3) = myNorm(data(:,3));

%Gamma correction
if ~exist('gamma','var')
    gamma = 1;
end

data(:,1) = data(:,1).^(1/gamma);
data(:,2) = data(:,2).^(1/gamma);
data(:,3) = data(:,3).^(1/gamma);



function [normData,relData,ternCoords] = myDataProcessing(data)
%Normalized data
normData(:,1) = myNorm(data(:,1));
normData(:,2) = myNorm(data(:,2));
normData(:,3) = myNorm(data(:,3));

%Relative contributions of each color
relData = myColorRatio(normData(:,1),normData(:,2),normData(:,3));
R = relData(:,1);
G = relData(:,2);

%Generate 2D Ternary Graph Coordinates
y = G*sin(pi/3);
x = R + y*cot(pi/3);

ternCoords = [x, y];

function [idx] = removeOutliers(data, plot_handle)

newData = [plot_handle.XData',plot_handle.YData',plot_handle.ZData'];

idx = find(ismember(data,newData,'rows'));

function handles = updatePlots(handles,highColorTern,highColorScatter)

if ishandle(handles.scatterFigHandle)
    if isvalid(handles.scatterFigHandle)
        figure(handles.scatterFigHandle)
        [currentAz,currentEl] = view;
    end
end
currentView = [135,30];

if ishandle(handles.scatterPlotHandle)
    if isvalid(handles.scatterPlotHandle)
        currentView = [currentAz,currentEl];
    end
end

if isvalid(handles.ternFigHandle)
    if ~isempty(handles.ternFigHandle)
        cla(handles.ternPlotHandle);
        cla(handles.densityOverlay);
    end
end


%Get ternary user options
switch handles.ternDisplayType.Value
    case 1
        plotDataTern = handles.normDataCus;
    case 2
        plotDataTern = handles.normDataDef;
    case 3
        plotDataTern = handles.normDataLin;
end


if ~exist('highColorTern','var')
    switch handles.ternColorType.Value
        case 1
            newColorTern = handles.normDataCus;
        case 2
            newColorTern = handles.normDataDef;
        case 3
            newColorTern = handles.normDataLin;
        case 4
            clusterColor = distinguishable_colors(max(handles.clusterIndex));
            newColorTern = clusterColor(handles.clusterIndex(:),:);
    end
else
    newColorTern = highColorTern;
end

%Update outputs
meanValues = getMeanColorValues(newColorTern,handles.clusterIndex);
handles.meanValues = meanValues;
[~,~,ternCoords] = myDataProcessing(plotDataTern);
tab = tabulate(handles.clusterIndex);
sortTab = flipud(sortrows(tab,3));
handles.sortTab = sortTab;

meanValuesExcel = zeros(max(sortTab(:,1)),3);
for ii = 1:max(sortTab(:,1))
    meanValuesExcel(ii,:) = meanValues(sortTab(ii,1),:);
end

handles.excelOutputSummary = [sortTab, meanValuesExcel];

handles.excelOutputRaw = [handles.normDataLin, handles.normDataDef,handles.normDataCus,...
    ternCoords, handles.clusterIndex];


%Update ternary graph
[handles.ternFigHandle,handles.ternPlotHandle,handles.densityOverlay,...
    handles.barPlotHandle,handles.sub1,handles.sub2]...
    = initiateTernPlot(plotDataTern,newColorTern, handles.sortTab,...
    handles.meanValues, handles.ternFigHandle);

%Get scatter user options
switch handles.scatterDisplayType.Value
    case 1
        plotDataScatter = handles.normDataCus;
    case 2
        plotDataScatter = handles.normDataDef;
    case 3
        plotDataScatter = handles.normDataLin;
end

if ~exist('highColorScatter','var')
    switch handles.scatterColorType.Value
        case 1
            newColorScatter = handles.normDataCus;
        case 2
            newColorScatter = handles.normDataDef;
        case 3
            newColorScatter = handles.normDataLin;
        case 4
            clusterColor = distinguishable_colors(max(handles.clusterIndex));
            newColorScatter = clusterColor(handles.clusterIndex(:),:);
    end
else
    newColorScatter = highColorScatter;
end


%Update scatter graph
[handles.scatterFigHandle, handles.scatterPlotHandle] = initiateScatter(...
    plotDataScatter,newColorScatter,handles.scatterFigHandle,currentView);


function meanValues = getMeanColorValues(transformData,cluster_idx)
meanValues = zeros(max(cluster_idx),3);

for kk = 1:(max(cluster_idx))
    meanValues(kk,1) = mean(transformData(cluster_idx == kk,1));
    meanValues(kk,2) = mean(transformData(cluster_idx == kk,2));
    meanValues(kk,3) = mean(transformData(cluster_idx == kk,3));
end


% --- Executes on button press in makeSilPlot.
function makeSilPlot_Callback(hObject, eventdata, handles)
% hObject    handle to makeSilPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Load data that we need
clusterIndex = handles.clusterIndex;
meanValues = handles.meanValues;
sampleName = handles.sampleName;

numOfClusters = numel(unique(clusterIndex));
switch handles.transformType.Value
    case 1
        transformData = handles.normDataCus;
    case 2
        transformData = handles.normDataDef;
    case 3
        transformData = handles.normDataLin;
end

if numOfClusters > 1
    [handles.silFigHandle,handles.silPlotHandle] = generateSilPlot(transformData,clusterIndex,meanValues,sampleName);
else
    errordlg('There must be more than one cluster','Error')
end

guidata(hObject,handles);




function [hFig, hPlot] = generateSilPlot(transformData, clusterIndex, meanValues,sampleName)
numOfClusters = numel(unique(clusterIndex));

hWait = waitbar(0,'Evaluating clustering... please wait...');
% Sort by cluster size
clusterTable = tabulate(clusterIndex);
sortedClusterTable = sortrows(clusterTable,3);
sortedClusters = sortedClusterTable(:,1);


sortedData = transformData(clusterIndex == sortedClusters(1),:);
sortedClusterIndex = clusterIndex(clusterIndex == sortedClusters(1));

for jj = 1:numOfClusters
    if jj == 1
        sortedData = transformData(clusterIndex == sortedClusters(1),:);
        sortedClusterIndex = clusterIndex(clusterIndex == sortedClusters(1));
    else
        sortedData = [sortedData; transformData(clusterIndex == sortedClusters(jj),:)];
        sortedClusterIndex = [sortedClusterIndex; clusterIndex(clusterIndex == sortedClusters(jj))];
    end
    
end

sortedMeanColor = meanValues(sortedClusters,:);

waitbar(0.5,hWait,'Generating graph...');
hFig = figure('Visible','off','Units','Normalized','Position',[0.08 0.08 0.38 0.62]);
[silData, hPlot,ticks] = mySilhouette(sortedData, sortedClusterIndex,...
    'sqEuclidean',sortedClusters);

set(get(get(hPlot,'Children'),'Children'),...
    'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.5);
hold on
title(['Sample ', sampleName],'interpreter','none');

silMeans = zeros(1,numOfClusters);
for kk = 1:numOfClusters
    numberOfBars = clusterTable(sortedClusters(kk),2);
    silMeans(kk) = mean(silData(sortedClusterIndex == sortedClusters(kk)));
    rectangle('Position',[0 ticks(kk)-0.5*numberOfBars 0.09 numberOfBars],...
        'Curvature', [0.2 0.2],...
        'FaceColor',sortedMeanColor(kk,:),'EdgeColor','none','LineStyle','none');
    rectangle('Position',[silMeans(kk) ticks(kk) 0.005 500],...
        'FaceColor','k');
end
hold off

close(hWait);
hFig.Visible = 'on';


% --- Executes on button press in drawGate.
function drawGate_Callback(hObject, eventdata, handles)
% hObject    handle to drawGate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'ternGate')
    if isvalid(handles.ternGate)
        delete(handles.ternGate);
    end
end


[handles.ternGate,index] = drawGate(get(handles.ternPlotHandle,'Parent'),handles.ternPlotHandle);

if isfield(handles,'backGateHandle')
    if isvalid(handles.backGateHandle)
        delete(handles.backGateHandle);
    end
end

handles.backGateHandle = makeBackGatePlots(handles, index);

% Update handles structure
guidata(hObject, handles);


function [h,index] = drawGate(ax,plot_handle)
axes(ax);
h = impoly;
api = iptgetapi(h);
api.setColor([0.7 0.5 0]);

set(h,'Tag','oldTernGate');

gateCoord = getPosition(h);
% xData = get(get(ax,'Children'),'XData');
% yData = get(get(ax,'Children'),'YData');
xData = plot_handle.XData';
yData = plot_handle.YData';
dataInGate = inpoly([xData,yData],gateCoord);
index = find(dataInGate == 1);

function [index] = updateGate(ax,plot_handle, poly_handle)
axes(ax);
% xData = get(get(ax,'Children'),'XData');
% yData = get(get(ax,'Children'),'YData');
xData = plot_handle.XData';
yData = plot_handle.YData';

gateCoord = getPosition(poly_handle);
dataInGate = inpoly([xData,yData],gateCoord);
index = find(dataInGate == 1);

function [h] = makeBackGatePlots(handles, gateIndex)
%Find red, green, and blue FACS parameters

channelNames = handles.channelNames;
fcsdat = handles.fcsdat;

%Default logicle values
T = handles.T;
M = handles.M;
A = handles.A;

cambridgeTest = strfind(channelNames,'FJComp-FITC-A');
    cambridgeTest = [cambridgeTest{:}];
    childrensTest = strfind(channelNames,'FJComp-GFP-A');
    childrensTest = [childrensTest{:}];
    
    if cambridgeTest == 1
        %For Aria
        fscA = fcsdat(:,find(strcmp('FSC-A',channelNames)));
        fscH = fcsdat(:,find(strcmp('FSC-H',channelNames)));
        fscW = fcsdat(:,find(strcmp('FSC-W',channelNames)));
        
        sscA = fcsdat(:,find(strcmp('SSC-A',channelNames)));
        sscAW = calculateLogicleWidth(handles,sscA);
        sscA = logicleTransform(sscA,T,sscAW,M,A);
        
        sscH = fcsdat(:,find(strcmp('SSC-H',channelNames)));
        sscHW = calculateLogicleWidth(handles,sscH);
        sscH = logicleTransform(sscH,T,sscHW,M,A);
        
        
        sscW = fcsdat(:,find(strcmp('SSC-W',channelNames)));
        sscWW = calculateLogicleWidth(handles,sscW);
        sscW = logicleTransform(sscW,T,sscWW,M,A);
        
        liveDye = fcsdat(:,find(strcmp('FJComp-APC-A',channelNames)));
        liveDyeW = calculateLogicleWidth(handles,liveDye);
        liveDye = logicleTransform(liveDye,T,liveDyeW,M,A);
        
        red = fcsdat(:,find(strcmp('FJComp-PE-A',channelNames)));
        redW = calculateLogicleWidth(handles,red);
        red = logicleTransform(red,T,redW,M,A);
        
        redAuto = fcsdat(:,find(strcmp('PE-Cy7-A',channelNames)));
        redAutoW = calculateLogicleWidth(handles,redAuto);
        redAuto = logicleTransform(redAuto,T,redAutoW,M,A);
        
        green = fcsdat(:,find(strcmp('FJComp-FITC-A',channelNames)));
        greenW = calculateLogicleWidth(handles,green);
        green = logicleTransform(green,T,greenW,M,A);
        
        greenAuto = fcsdat(:,find(strcmp('PerCP-A',channelNames)));
        greenAutoW = calculateLogicleWidth(handles,greenAuto);
        greenAuto = logicleTransform(greenAuto,T,greenAutoW,M,A);
        
        blue = fcsdat(:,find(strcmp('FJComp-CFP-A',channelNames)));
        blueW = calculateLogicleWidth(handles,blue);
        blue = logicleTransform(blue,T,blueW,M,A);
        
        blueAuto = fcsdat(:,find(strcmp('Chromomycin-A',channelNames)));
        blueAutoW = calculateLogicleWidth(handles,blueAuto);
        blueAuto = logicleTransform(blueAuto,T,blueAutoW,M,A);
        
    elseif childrensTest == 1
        %For Children's Aria (self-sorter)
        fscA = fcsdat(:,find(strcmp('FSC-A',channelNames)));
        fscH = fcsdat(:,find(strcmp('FSC-H',channelNames)));
        fscW = fcsdat(:,find(strcmp('FSC-W',channelNames)));
        
        sscA = fcsdat(:,find(strcmp('SSC-A',channelNames)));
        sscAW = calculateLogicleWidth(handles,sscA);
        sscA = logicleTransform(sscA,T,sscAW,M,A);
        
        sscH = fcsdat(:,find(strcmp('SSC-H',channelNames)));
        sscHW = calculateLogicleWidth(handles,sscH);
        sscH = logicleTransform(sscH,T,sscHW,M,A);
        
        
        sscW = fcsdat(:,find(strcmp('SSC-W',channelNames)));
        sscWW = calculateLogicleWidth(handles,sscW);
        sscW = logicleTransform(sscW,T,sscWW,M,A);
        
        liveDye = fcsdat(:,find(strcmp('FJComp-APC-Cy7-A',channelNames)));
        liveDyeW = calculateLogicleWidth(handles,liveDye);
        liveDye = logicleTransform(liveDye,T,liveDyeW,M,A);
        
        red = fcsdat(:,find(strcmp('FJComp-DsRed-A',channelNames)));
        redW = calculateLogicleWidth(handles,red);
        red = logicleTransform(red,T,redW,M,A);
        
        redAuto = fcsdat(:,find(strcmp('PE-Cy7-A',channelNames)));
        redAutoW = calculateLogicleWidth(handles,redAuto);
        redAuto = logicleTransform(redAuto,T,redAutoW,M,A);
        
        green = fcsdat(:,find(strcmp('FJComp-GFP-A',channelNames)));
        greenW = calculateLogicleWidth(handles,green);
        green = logicleTransform(green,T,greenW,M,A);
        
        greenAuto = fcsdat(:,find(strcmp('PerCP-A',channelNames)));
        greenAutoW = calculateLogicleWidth(handles,greenAuto);
        greenAuto = logicleTransform(greenAuto,T,greenAutoW,M,A);
        
        blue = fcsdat(:,find(strcmp('FJComp-DAPI-A',channelNames)));
        blueW = calculateLogicleWidth(handles,blue);
        blue = logicleTransform(blue,T,blueW,M,A);
        
        blueAuto = fcsdat(:,find(strcmp('Hoechst Red-A',channelNames)));
        blueAutoW = calculateLogicleWidth(handles,blueAuto);
        blueAuto = logicleTransform(blueAuto,T,blueAutoW,M,A);
    else
        msgbox('Could not identify fluorescence parameters...');
    end



totalData = [fscA, fscH, fscW, sscA, sscH, sscW, liveDye,...
    red, redAuto, green, greenAuto, blue, blueAuto];

gatedData = totalData(gateIndex,:);

h = figure;

subplot(2,5,1),
scatter(totalData(:,1),totalData(:,4),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,1),gatedData(:,4),15,'b','filled');
hold off
xlabel('FSC-A'), ylabel('SSC-A');
xlim([0, 2^18]), ylim([0 1]);

subplot(2,5,2),
scatter(totalData(:,1),totalData(:,2),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,1),gatedData(:,2),15,'b','filled');
hold off
xlabel('FSC-A'), ylabel('FSC-H');
xlim([0, 2^18]), ylim([0 2^18]);

subplot(2,5,3),
scatter(totalData(:,6),totalData(:,5),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,6),gatedData(:,5),15,'b','filled');
hold off
xlabel('SSC-W'), ylabel('SSC-H');
xlim([0, 1]), ylim([0 1]);

subplot(2,5,4),
scatter(totalData(:,1),totalData(:,7),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,1),gatedData(:,7),15,'b','filled');
hold off
xlabel('FSC-A'), ylabel('Live dye');
xlim([0, 2^18]), ylim([0 1]);

subplot(2,5,6),
scatter(totalData(:,12),totalData(:,13),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,12),gatedData(:,13),15,'b','filled');
hold off
xlabel('CFP-A'), ylabel('Auto CFP-A');
xlim([0, 1]), ylim([0 1]);

subplot(2,5,7),
scatter(totalData(:,10),totalData(:,11),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,10),gatedData(:,11),15,'b','filled');
hold off
xlabel('YFP-A'), ylabel('Auto YFP-A');
xlim([0, 1]), ylim([0 1]);

subplot(2,5,8),
scatter(totalData(:,8),totalData(:,9),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,8),gatedData(:,9),15,'b','filled');
hold off
xlabel('RFP-A'), ylabel('Auto RFP-A');
xlim([0, 1]), ylim([0 1]);

subplot(2,5,9),
scatter(totalData(:,10),totalData(:,12),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,10),gatedData(:,12),15,'b','filled');
hold off
xlabel('YFP-A'), ylabel('CFP-A');
xlim([0, 1]), ylim([0 1]);

subplot(2,5,10),
scatter(totalData(:,8),totalData(:,10),15,[0.8 0.8 0.8],'filled');
hold on
scatter(gatedData(:,8),gatedData(:,10),15,'b','filled');
hold off
xlabel('RFP-A'), ylabel('YFP-A');
xlim([0, 1]), ylim([0 1]);






function w = calculateLogicleWidth(handles,data)
%Default logicle values
T = handles.T;
M = handles.M;

w = (M - log10(T/abs(min(data))))/2;

if w < 0
    w = 0;
end