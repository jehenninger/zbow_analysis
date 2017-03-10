function varargout = myeloidLymphoidComparison(varargin)
% MYELOIDLYMPHOIDCOMPARISON MATLAB code for myeloidLymphoidComparison.fig
%      MYELOIDLYMPHOIDCOMPARISON, by itself, creates a new MYELOIDLYMPHOIDCOMPARISON or raises the existing
%      singleton*.
%
%      H = MYELOIDLYMPHOIDCOMPARISON returns the handle to a new MYELOIDLYMPHOIDCOMPARISON or the handle to
%      the existing singleton*.
%
%      MYELOIDLYMPHOIDCOMPARISON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYELOIDLYMPHOIDCOMPARISON.M with the given input arguments.
%
%      MYELOIDLYMPHOIDCOMPARISON('Property','Value',...) creates a new MYELOIDLYMPHOIDCOMPARISON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before myeloidLymphoidComparison_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to myeloidLymphoidComparison_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help myeloidLymphoidComparison

% Last Modified by GUIDE v2.5 15-Jun-2016 19:05:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @myeloidLymphoidComparison_OpeningFcn, ...
                   'gui_OutputFcn',  @myeloidLymphoidComparison_OutputFcn, ...
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


% --- Executes just before myeloidLymphoidComparison is made visible.
function myeloidLymphoidComparison_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to myeloidLymphoidComparison (see VARARGIN)

% Choose default command line output for myeloidLymphoidComparison
handles.output = hObject;

% To turn on/off erythroid, change this value to 1/0
handles.includeErythroid = 1;

handles.T = 2^18;
handles.M = 4.5;
handles.A = 0;

handles.backGateCount = 0;
handles.scatter3DFig = [];
handles.ternPlotFig = [];

handles.clusterCount = 0;

handles.clusterEvalCount = 0;

set(gcf,'units','normalized','OuterPosition', [0.45, 0.3, 0.55, 0.7]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes myeloidLymphoidComparison wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = myeloidLymphoidComparison_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadFile.
function loadFile_Callback(hObject, eventdata, handles)
% hObject    handle to loadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Load myeloid file
[fileNameMyeloid, pathNameMyeloid] = uigetfile('*myeloid*.fcs','Select myeloid data','E:\Zon Lab\FACS\');
set(handles.fileNameMyeloid,'String',fileNameMyeloid);
handles.fullFileNameMyeloid = fullfile(pathNameMyeloid,fileNameMyeloid);

[linDataMyeloid, logDataDefaultMyeloid, logDataTernMyeloid] = loadFSCFile(...
    fullfile(pathNameMyeloid,fileNameMyeloid), handles);

handles.linDataMyeloid = linDataMyeloid;

handles.logDataDefaultMyeloid = logDataDefaultMyeloid;
handles.logDataDefaultMyeloidCopy = handles.logDataDefaultMyeloid;

handles.logDataTernMyeloid = logDataTernMyeloid;
handles.logDataTernMyeloidCopy = handles.logDataTernMyeloid;



%% Load lymphoid file
[fileNameLymphoid, pathNameLymphoid] = uigetfile('*lymphoid*.fcs','Select lymphoid data',pathNameMyeloid);
set(handles.fileNameLymphoid,'String',fileNameLymphoid);
handles.fullFileNameLymphoid = fullfile(pathNameLymphoid,fileNameLymphoid);

[linDataLymphoid, logDataDefaultLymphoid, logDataTernLymphoid] = loadFSCFile(...
    fullfile(pathNameLymphoid,fileNameLymphoid), handles);

handles.linDataLymphoid = linDataLymphoid;

handles.logDataDefaultLymphoid = logDataDefaultLymphoid;
handles.logDataDefaultLymphoidCopy = handles.logDataDefaultLymphoid;

handles.logDataTernLymphoid = logDataTernLymphoid;
handles.logDataTernLymphoidCopy = handles.logDataTernLymphoid;



%% Load erythroid file
if handles.includeErythroid == 1
    [fileNameErythroid, pathNameErythroid] = uigetfile([pathNameMyeloid,'*precursor.fcs'],'Select erythroid data');
    set(handles.fileNameErythroid,'String',fileNameErythroid);
    
    [linDataErythroid, logDataDefaultErythroid, logDataTernErythroid] = loadFSCFile(...
        fullfile(pathNameErythroid,fileNameErythroid), handles);
    
    handles.linDataErythroid = linDataErythroid;
    
    handles.logDataDefaultErythroid = logDataDefaultErythroid;
    handles.logDataDefaultErythroidCopy = handles.logDataDefaultErythroid;
    
    handles.logDataTernErythroid = logDataTernErythroid;
    handles.logDataTernErythroidCopy = handles.logDataTernErythroid;
end



%%Make 3D scatter plots for observation and outlier removal
[handles.scatterFigureHandle, handles.scatterPlotHandleMyeloid] = initiate3DScatter(...
    handles, logDataDefaultMyeloid,logDataTernMyeloid,'Myeloid');
[~, handles.scatterPlotHandleLymphoid] = initiate3DScatter(...
    handles, logDataDefaultLymphoid,logDataTernLymphoid,'Lymphoid',handles.scatterFigureHandle,2);

if handles.includeErythroid == 1
    [~, handles.scatterPlotHandleErythroid] = initiate3DScatter(...
        handles, logDataDefaultErythroid,logDataTernErythroid,'Erythroid',handles.scatterFigureHandle,3);
end

% Update handles structure
guidata(hObject, handles);

updateTernPlots_Callback(handles.updateTernPlots,eventdata,handles);


% --- Executes on button press in drawMyeloidGate.
function drawMyeloidGate_Callback(hObject, eventdata, handles)
% hObject    handle to drawMyeloidGate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles,'myeloidGate');
    delete(handles.myeloidGate);
end

[handles.myeloidGate,index] = drawGate(handles.myeloidTern,handles.ternPlotHandleMyeloid, 'Myeloid');

handles.outputTernCoordXMyeloid = mean(handles.ternPlotHandleMyeloid.XData(index));
handles.outputTernCoordYMyeloid = mean(handles.ternPlotHandleMyeloid.YData(index));
handles.outputColorMyeloid = [...
    mean(handles.ternPlotHandleMyeloid.CData(index,1)),...
    mean(handles.ternPlotHandleMyeloid.CData(index,2)),...
    mean(handles.ternPlotHandleMyeloid.CData(index,3))];


totalCellCount = size(handles.logDataTernMyeloid,1);
myeloidCount = size(index,1);
gatePercentage = 100*myeloidCount/totalCellCount;

handles.numOfCellsInGateMyeloid.String = num2str(myeloidCount);
handles.totalNumOfCellsMyeloid.String = num2str(totalCellCount);
handles.gatePercentageMyeloid.String = num2str(gatePercentage);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in updateMyeloidGate.
function updateMyeloidGate_Callback(hObject, eventdata, handles)
% hObject    handle to updateMyeloidGate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[index] = updateGate(handles.myeloidTern, handles.ternPlotHandleMyeloid, handles.myeloidGate);

handles.outputTernCoordXMyeloid = mean(handles.ternPlotHandleMyeloid.XData(index));
handles.outputTernCoordYMyeloid = mean(handles.ternPlotHandleMyeloid.YData(index));
handles.outputColorMyeloid = [...
    mean(handles.ternPlotHandleMyeloid.CData(index,1)),...
    mean(handles.ternPlotHandleMyeloid.CData(index,2)),...
    mean(handles.ternPlotHandleMyeloid.CData(index,3))];


totalCellCount = size(handles.logDataTernMyeloid,1);
myeloidCount = size(index,1);
gatePercentage = 100*myeloidCount/totalCellCount;

handles.numOfCellsInGateMyeloid.String = num2str(myeloidCount);
handles.totalNumOfCellsMyeloid.String = num2str(totalCellCount);
handles.gatePercentageMyeloid.String = num2str(gatePercentage);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in drawLymphoidGate.
function drawLymphoidGate_Callback(hObject, eventdata, handles)
% hObject    handle to drawLymphoidGate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'lymphoidGate');
    delete(handles.lymphoidGate);
end

[handles.lymphoidGate,index] = drawGate(handles.lymphoidTern,handles.ternPlotHandleLymphoid, 'Lymphoid');

handles.outputTernCoordXLymphoid = mean(handles.ternPlotHandleLymphoid.XData(index));
handles.outputTernCoordYLymphoid = mean(handles.ternPlotHandleLymphoid.YData(index));
handles.outputColorLymphoid = [...
    mean(handles.ternPlotHandleLymphoid.CData(index,1)),...
    mean(handles.ternPlotHandleLymphoid.CData(index,2)),...
    mean(handles.ternPlotHandleLymphoid.CData(index,3))];

totalCellCount = size(handles.logDataTernLymphoid,1);
lymphoidCount = size(index,1);
gatePercentage = 100*lymphoidCount/totalCellCount;

handles.numOfCellsInGateLymphoid.String = num2str(lymphoidCount);
handles.totalNumOfCellsLymphoid.String = num2str(totalCellCount);
handles.gatePercentageLymphoid.String = num2str(gatePercentage);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in updateLymphoidGate.
function updateLymphoidGate_Callback(hObject, eventdata, handles)
% hObject    handle to updateLymphoidGate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[index] = updateGate(handles.lymphoidTern,handles.ternPlotHandleLymphoid, handles.lymphoidGate);

handles.outputTernCoordXLymphoid = mean(handles.ternPlotHandleLymphoid.XData(index));
handles.outputTernCoordYLymphoid = mean(handles.ternPlotHandleLymphoid.YData(index));
handles.outputColorLymphoid = [...
    mean(handles.ternPlotHandleLymphoid.CData(index,1)),...
    mean(handles.ternPlotHandleLymphoid.CData(index,2)),...
    mean(handles.ternPlotHandleLymphoid.CData(index,3))];

totalCellCount = size(handles.logDataTernLymphoid,1);
lymphoidCount = size(index,1);
gatePercentage = 100*lymphoidCount/totalCellCount;

handles.numOfCellsInGateLymphoid.String = num2str(lymphoidCount);
handles.totalNumOfCellsLymphoid.String = num2str(totalCellCount);
handles.gatePercentageLymphoid.String = num2str(gatePercentage);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in drawGateErythroid.
function drawGateErythroid_Callback(hObject, eventdata, handles)
% hObject    handle to drawGateErythroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.includeErythroid == 1
    if isfield(handles,'erythroidGate');
        delete(handles.erythroidGate);
    end
    
    [handles.erythroidGate,index] = drawGate(handles.erythroidTern,handles.ternPlotHandleErythroid, 'Erythroid');
    
    totalCellCount = size(handles.logDataTernErythroid,1);
    erythroidCount = size(index,1);
    gatePercentage = 100*erythroidCount/totalCellCount;
    
    handles.numOfCellsInGateErythroid.String = num2str(erythroidCount);
    handles.totalNumOfCellsErythroid.String = num2str(totalCellCount);
    handles.gatePercentageErythroid.String = num2str(gatePercentage);
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in updateGateErythroid.
function updateGateErythroid_Callback(hObject, eventdata, handles)
% hObject    handle to updateGateErythroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.includeErythroid == 1
    [index] = updateGate(handles.erythroidTern,handles.ternPlotHandleErythroid, handles.erythroidGate);
    
    totalCellCount = size(handles.logDataTernErythroid,1);
    erythroidCount = size(index,1);
    gatePercentage = 100*erythroidCount/totalCellCount;
    
    handles.numOfCellsInGateErythroid.String = num2str(erythroidCount);
    handles.totalNumOfCellsErythroid.String = num2str(totalCellCount);
    handles.gatePercentageErythroid.String = num2str(gatePercentage);
    
    % Update handles structure
    guidata(hObject, handles);
end



% --- Executes on button press in updateTernPlots.
function updateTernPlots_Callback(hObject, eventdata, handles)
% hObject    handle to updateTernPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.myeloidTern);
cla(handles.lymphoidTern);
cla(handles.erythroidTern);

logicleChoiceMyeloid = handles.logicleChoiceMyeloid.Value;
logicleChoiceLymphoid = handles.logicleChoiceLymphoid.Value;

markerSize = 5;

%%  Myeloid
% Process data for making ternary graphs
[normDataTernMyeloid, ternCoordsMyeloid] = myDataProcessing(handles.logDataTernMyeloid);
[normDataDefaultMyeloid, ternCoordsDefaultMyeloid] = myDataProcessing(handles.logDataDefaultMyeloid);
[normDataLinMyeloid, ternCoordsLinMyeloid] = myDataProcessing(handles.linDataMyeloid);
[newColorMyeloid] = myColorProcessing(normDataTernMyeloid);


%Plot ternary graph
axes(handles.myeloidTern);
% ternPlot(ternCoordsMyeloid,normDataTernMyeloid,'false');

switch logicleChoiceMyeloid
    case 1
        plotDataMyeloid = ternCoordsDefaultMyeloid;
        %plotDataMyeloid3D = normDataDefaultMyeloid;
        plotDataMyeloid3D = handles.logDataDefaultMyeloid;
    case 2
        plotDataMyeloid = ternCoordsMyeloid;
%         plotDataMyeloid3D = normDataTernMyeloid;
        plotDataMyeloid3D = handles.logDataTernMyeloid;
    case 3
        plotDataMyeloid = ternCoordsLinMyeloid;
%         plotDataMyeloid3D = normDataLinMyeloid;
        plotDataMyeloid3D = handles.linDataMyeloid;
end

sampleSize = 20000;
if size(plotDataMyeloid,1) < sampleSize
    handles.ternPlotHandleMyeloid = ternPlot(plotDataMyeloid,newColorMyeloid,'false','true',markerSize);
else
    [~, dataSample] = datasample(plotDataMyeloid,sampleSize,'Replace',false);
    handles.ternPlotHandleMyeloid = ternPlot(plotDataMyeloid(dataSample,:),...
        newColorMyeloid(dataSample,:),'false','true',markerSize);
end

% Update 3D scatter

cla(get(handles.scatterPlotHandleMyeloid,'Parent'))
cla(get(handles.scatterPlotHandleLymphoid,'Parent'))
if isfield(handles,'scatterPlotErythroid')
    cla(get(handles.scatterPlotHandleErythroid,'Parent'));
end

[handles.scatterFigureHandle, handles.scatterPlotHandleMyeloid] = initiate3DScatter(...
    handles, plotDataMyeloid3D,handles.logDataTernMyeloid,'Myeloid',handles.scatterFigureHandle);



%%  Lymphoid
% Process data for making ternary graphs
[normDataTernLymphoid, ternCoordsLymphoid] = myDataProcessing(handles.logDataTernLymphoid);
[normDataDefaultLymphoid, ternCoordsDefaultLymphoid] = myDataProcessing(handles.logDataDefaultLymphoid);
[normDataLinLymphoid, ternCoordsLinLymphoid] = myDataProcessing(handles.linDataLymphoid);
[newColorLymphoid] = myColorProcessing(normDataTernLymphoid,1);


%Plot ternary graph

axes(handles.lymphoidTern);
% ternPlot(ternCoordsLymphoid,normDataTernLymphoid,'false');

switch logicleChoiceLymphoid
    case 1
        plotDataLymphoid = ternCoordsDefaultLymphoid;
        %         plotDataLymphoid3D = normDataDefaultLymphoid;
        plotDataLymphoid3D = handles.logDataDefaultLymphoid;
    case 2
        plotDataLymphoid = ternCoordsLymphoid;
        %         plotDataLymphoid3D = normDataTernLymphoid;
        plotDataLymphoid3D = handles.logDataTernLymphoid;
    case 3
        plotDataLymphoid = ternCoordsLinLymphoid;
        %         plotDataLymphoid3D = normDataLinLymphoid;
        plotDataLymphoid3D = handles.linDataLymphoid;
end

if size(plotDataLymphoid,1) < sampleSize
    handles.ternPlotHandleLymphoid = ternPlot(plotDataLymphoid,...
        newColorLymphoid,'false','true',markerSize);
else
    [~,dataSample] = datasample(plotDataLymphoid,sampleSize,'Replace',false);
    handles.ternPlotHandleLymphoid = ternPlot(plotDataLymphoid(dataSample,:),...
        newColorLymphoid(dataSample,:),'false','true',markerSize);
end

% Update 3D scatter

[~, handles.scatterPlotHandleLymphoid] = initiate3DScatter(...
    handles, plotDataLymphoid3D,handles.logDataTernLymphoid,'Lymphoid',handles.scatterFigureHandle,2);

%%  Erythroid

if handles.includeErythroid == 1
    % Process data for making ternary graphs
    [normDataTernErythroid, ternCoordsErythroid] = myDataProcessing(handles.logDataTernErythroid);
    [normDataDefaultErythroid, ternCoordsDefaultErythroid] = myDataProcessing(handles.logDataDefaultErythroid);
    [normDataLinErythroid, ternCoordsLinErythroid] = myDataProcessing(handles.linDataErythroid);
    [newColorErythroid] = myColorProcessing(normDataTernErythroid,1.2);
    
    
    %Plot ternary graph
    axes(handles.erythroidTern);
    % ternPlot(ternCoordsErythroid,normDataTernErythroid,'false');
    
    if size(ternCoordsDefaultErythroid,1) < sampleSize
        handles.ternPlotHandleErythroid = ternPlot(ternCoordsLinErythroid,...
            newColorErythroid,'false','true',markerSize);
    else
        [~,dataSample] = datasample(ternCoordsLinErythroid,sampleSize,'Replace',false);
        handles.ternPlotHandleErythroid = ternPlot(ternCoordsLinErythroid(dataSample,:),...
            newColorErythroid(dataSample,:),'false','true',markerSize);
    end
    
    % Update 3D scatter
    
    [handles.scatterFigureHandle, handles.scatterPlotHandleErythroid] = initiate3DScatter(...
        handles, normDataLinErythroid,handles.logDataTernErythroid,'Erythroid',handles.scatterFigureHandle,3);
    
end


% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in removeOutliersMyeloid.
function removeOutliersMyeloid_Callback(hObject, eventdata, handles)
% hObject    handle to removeOutliersMyeloid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx = removeOutliers(handles.logDataDefaultMyeloid,handles.scatterPlotHandleMyeloid);

handles.logDataDefaultMyeloid = handles.logDataDefaultMyeloid(idx,:);
handles.logDataTernMyeloid = handles.logDataTernMyeloid(idx,:);
handles.linDataMyeloid = handles.linDataMyeloid(idx,:);


% Update handles structure
guidata(hObject, handles);

updateTernPlots_Callback(handles.updateTernPlots,eventdata,handles);


% --- Executes on button press in removeOutliersLymphoid.
function removeOutliersLymphoid_Callback(hObject, eventdata, handles)
% hObject    handle to removeOutliersLymphoid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

idx = removeOutliers(handles.logDataDefaultLymphoid,handles.scatterPlotHandleLymphoid);

handles.logDataDefaultLymphoid = handles.logDataDefaultLymphoid(idx,:);
handles.logDataTernLymphoid = handles.logDataTernLymphoid(idx,:);
handles.linDataLymphoid = handles.linDataLymphoid(idx,:);

% Update handles structure
guidata(hObject, handles);

updateTernPlots_Callback(handles.updateTernPlots,eventdata,handles);


% --- Executes on button press in removeOutliersErythroid.
function removeOutliersErythroid_Callback(hObject, eventdata, handles)
% hObject    handle to removeOutliersErythroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.includeErythroid == 1
    idx = removeOutliers(handles.logDataDefaultErythroid,handles.scatterPlotHandleErythroid);
    
    handles.logDataDefaultErythroid = handles.logDataDefaultErythroid(idx,:);
    handles.logDataTernErythroid = handles.logDataTernErythroid(idx,:);
    handles.linDataErythroid = handles.linDataErythroid(idx,:);
    
    % Update handles structure
    guidata(hObject, handles);
    
    updateTernPlots_Callback(handles.updateTernPlots,eventdata,handles);
    
end



function [data] = myColorProcessing(data,gamma)


%% For normalization to a certain percentage of top and bottom cells
% numOfCells = size(data,1);
% % max_threshold = round(0.05*numOfCells);
% max_threshold = 1;
% min_threshold = 1;
% minColor = 0;
% 
% data(:,1) = myNorm(data(:,1));
% data(:,2) = myNorm(data(:,2));
% data(:,3) = myNorm(data(:,3));
% 
% %Red
% s = sort(data(:,1),'descend');
% redMax = s(max_threshold);
% s = sort(data(:,1),'ascend');
% redMin = s(min_threshold);
% 
% %Green
% s = sort(data(:,2),'descend');
% greenMax = s(max_threshold);
% s = sort(data(:,2),'ascend');
% greenMin = s(min_threshold);
% 
% %Blue
% s = sort(data(:,3),'descend');
% blueMax = s(max_threshold);
% s = sort(data(:,3),'ascend');
% blueMin = s(min_threshold);
% 
% %Red
% data(data(:,1)>redMin | data(:,1)<redMax,1) = myNorm(data(data(:,1)>redMin | data(:,1)<redMax,1)...
%     ,1,minColor);
% data(data(:,1)>redMax,1) = 1;
% data(data(:,1)<redMin,1) = minColor;
% 
% %Green
% data(data(:,2)>greenMin | data(:,2)<greenMax,2) = myNorm(data(data(:,2)>greenMin | data(:,2)<greenMax,2)...
%     ,1,minColor);
% data(data(:,2)>greenMax,2) = 1;
% data(data(:,2)<greenMin,2) = minColor;
% 
% %Blue
% data(data(:,3)>blueMin | data(:,3)<blueMax,3) = myNorm(data(data(:,3)>blueMin | data(:,3)<blueMax,3)...
%     ,1,minColor);
% data(data(:,3)>blueMax,3) = 1;
% data(data(:,3)<blueMin,3) = minColor;

%% Regular normalization
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



function [normData,ternCoords] = myDataProcessing(data)
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



function [h,index] = drawGate(ax,plot_handle, gateType)
axes(ax);
h = impoly;
api = iptgetapi(h);
api.setColor([0.7 0.5 0]);

switch gateType
    case 'Myeloid'
        set(h,'Tag','oldMyeloidGate');
    case 'Lymphoid'
        set(h,'Tag','oldLymphoidGate');
    case 'Erythroid'
        set(h,'Tag','oldErythroidGate');
end

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

function [idx] = removeOutliers(data, plot_handle)

newData = [plot_handle.XData',plot_handle.YData',plot_handle.ZData'];

idx = find(ismember(data,newData,'rows'));

function [figure_handle, plot_handle] = initiate3DScatter(handles, data,color,type,previous_handle,subplot_num)
% %Normalize data
% data(:,1) = myNorm(data(:,1));
% data(:,2) = myNorm(data(:,2));
% data(:,3) = myNorm(data(:,3));

%Normalize color
switch type
    case 'Myeloid'
        color = myColorProcessing(color);
    case 'Lymphoid'
        color = myColorProcessing(color,2.2);
    case 'Erythroid'
        color = myColorProcessing(color,2.2);
end

color = myColorProcessing(color);

if  nargin == 5
    figure_handle = previous_handle;
    figure(figure_handle);
    subplot_num = 1;
elseif nargin > 5
    figure_handle = previous_handle;
    figure(figure_handle);
else
    figure_handle = figure;
    set(figure_handle,'Units','normalized','Position', [0, 0.3, 0.4, 0.6]);
    subplot_num = 1;
end

if handles.includeErythroid == 1
    subplot(1,3,subplot_num)
else
    subplot(1,2,subplot_num)
end

plot_handle = scatter3(data(:,1),...
    data(:,2),...
    data(:,3),...
    8,color,'filled');

xlabel('Red'),ylabel('Green'),zlabel('Blue');
view(135,30);
axis('equal');

switch type
    case 'Myeloid'
        title('Myeloid');
    case 'Lymphoid'
        title('Lymphoid');
    case 'Erythroid'
        title('Erythroid');
end

function [linData, logDataDefault, logDataTern] = loadFSCFile(file, handles)
%Load and read data
[fcsdat, fcshdr] = fca_readfcs(file);
channelNames = {fcshdr.par.name};
%For Aria
red = fcsdat(:,find(strcmp('FJComp-PE-A',channelNames)));
green = fcsdat(:,find(strcmp('FJComp-FITC-A',channelNames)));
blue = fcsdat(:,find(strcmp('FJComp-CFP-A',channelNames)));

sessionData = [red green blue];
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


% --- Executes on selection change in logicleChoiceMyeloid.
function logicleChoiceMyeloid_Callback(hObject, eventdata, handles)
% hObject    handle to logicleChoiceMyeloid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns logicleChoiceMyeloid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from logicleChoiceMyeloid


% --- Executes during object creation, after setting all properties.
function logicleChoiceMyeloid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to logicleChoiceMyeloid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in logicleChoiceLymphoid.
function logicleChoiceLymphoid_Callback(hObject, eventdata, handles)
% hObject    handle to logicleChoiceLymphoid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns logicleChoiceLymphoid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from logicleChoiceLymphoid


% --- Executes during object creation, after setting all properties.
function logicleChoiceLymphoid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to logicleChoiceLymphoid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in copyToClip.
function copyToClip_Callback(hObject, eventdata, handles)
% hObject    handle to copyToClip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


numOfCellsInGateMyeloid = str2double(handles.numOfCellsInGateMyeloid.String);
totalNumOfCellsMyeloid = str2double(handles.totalNumOfCellsMyeloid.String);
gatePercentageMyeloid = str2double(handles.gatePercentageMyeloid.String);

numOfCellsInGateLymphoid = str2double(handles.numOfCellsInGateLymphoid.String);
totalNumOfCellsLymphoid = str2double(handles.totalNumOfCellsLymphoid.String);
gatePercentageLymphoid = str2double(handles.gatePercentageLymphoid.String);

outputTernCoordXMyeloid = handles.outputTernCoordXMyeloid;
outputTernCoordYMyeloid = handles.outputTernCoordYMyeloid;
outputColorMyeloid = handles.outputColorMyeloid;

outputTernCoordXLymphoid = handles.outputTernCoordXLymphoid;
outputTernCoordYLymphoid = handles.outputTernCoordYLymphoid;
outputColorLymphoid = handles.outputColorLymphoid;

output = [numOfCellsInGateMyeloid,totalNumOfCellsMyeloid,gatePercentageMyeloid,...
    numOfCellsInGateLymphoid, totalNumOfCellsLymphoid,gatePercentageLymphoid,...
    outputTernCoordXMyeloid,outputTernCoordYMyeloid,outputColorMyeloid,...
    outputTernCoordXLymphoid,outputTernCoordYLymphoid,outputColorLymphoid];

num2clip(output);


% --- Executes on button press in resetButton.
function resetButton_Callback(hObject, eventdata, handles)
% hObject    handle to resetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ishandle(handles.scatterFigureHandle) == 1
    close(handles.scatterFigureHandle);
end
close(gcbf);
myeloidLymphoidComparison;