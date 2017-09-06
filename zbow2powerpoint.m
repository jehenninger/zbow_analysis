function [] = zbow2powerpoint(files, sampleSize)
% TO DO: Make file size smaller for powerpoint.


% close any open pptx instances. This is in case there was an error in this
% function where a pptx instance was still open

pptxTest = exportToPPTX('query');

if ~isempty(pptxTest)
    exportToPPTX('close');
end

% setting parameters

markerSize = 50;
alpha = 0.5;

% if file is not provided, then allow user to select multiple files
% code below generates full file names + path
if ~exist('files','var') || isempty(files)
    
    [fileNames, pathName] = uigetfile('*.fcs', 'Multiselect','on');
    
    if ~iscell(fileNames) && fileNames == 0
        error('Could not load files');
    end
    
end

if iscell(fileNames) == 1
    numFiles = numel(fileNames);
    for ii = 1:numFiles
        files{ii} = fullfile(pathName, fileNames{ii});
    end
else
    numFiles = 1;
    files = fullfile(pathName, fileNames);
    files = {files};
end

% set default for sample size
if ~exist('sampleSize','var') || isempty(sampleSize)
    sampleSize = 20000;
end

%initiate wait bar
multiWaitbar('Generating output...', 0);

%initialize powerpoint
exportToPPTX('new');

for kk = 1:numFiles
    
    [~, sampleName, ~] = fileparts(files{kk});
    [hFig3D, ~, ~, ~] = make3DZbowPlot(files{kk},sampleSize, markerSize, alpha,[], 'off');
    
    exportToPPTX('addslide');
    exportToPPTX('addpicture', get(hFig3D,'Children')); %should this be figure or axis?
    exportToPPTX('addtext', sampleName, 'FontWeight','bold');
    
    [hFig2D, ~, ~] = make2DZbowPlot(files{kk},sampleSize, markerSize, alpha,'off');
    
    exportToPPTX('addslide');
    exportToPPTX('addpicture', get(hFig2D,'Children')); %should this be figure or axis?
    exportToPPTX('addtext', sampleName, 'FontWeight','bold');
    
    close(hFig3D);
    close(hFig2D);
    
    multiWaitbar('Generating output...', kk/numFiles);
end

exportToPPTX('save',fullfile(pathName,'exported_plots.pptx'));
exportToPPTX('close');

fprintf(['Successfully generated and saved the file to\n ', fullfile(pathName, 'exported_plots.pptx'),'\n']);

multiWaitbar('CloseAll');


end

