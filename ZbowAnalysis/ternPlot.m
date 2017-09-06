function [h,hDen] = ternPlot(ternCoords, colorm,labels,outline, markerSize, density, alpha)
%UNTITLED3 Summary of this function goes here
%   ternCoords is the x and y coordinates of the RGB relative values
%   colorm is a mx3 color matrix specifying color for each cell in
%   labels is 'true' to have labels or 'false' to show no labels
%   outline is 'true' to show outline and 'false' to have no outline
%   markerSize is the desired marker size
%   density is 'true' to show contour and 'false' to not show contour


if nargin < 3
    labels = 'true';
    outline = 'true';
    markerSize = 15;
    density = 'true';
    
elseif nargin < 4
    outline = 'true';
    markerSize = 15;
    density = 'true';
    
elseif nargin < 5
    markerSize = 15;
    density = 'true';
    
elseif nargin < 6
    density = 'true';
end

if ~exist('alpha','var') || isempty('alpha')
    alpha = 1;
end

x = ternCoords(:,1);
y = ternCoords(:,2);

majors=5;%number of major ticks
%plot axis lines

switch outline
    case 'true'
        patch('xdata', [0 1 0.5 0], 'ydata', [0 0 sin(pi/3) 0], 'edgecolor','k','linewidth',2,'facecolor','w','handlevisibility','off');
        hold on
        plot ([0 1],[-0.001 -0.001], 'color', 'k', 'linewidth',2,'handlevisibility','off');
    case 'false'
end



switch labels
    case 'true'
        % Generate labels
        majorticks = linspace(0, 1, majors + 1);
        majorticksLess = [0,0.5,1];
        % majorticks = majorticks(1:end);
        labels = num2str(majorticks'*100);
        labelsLess = num2str(majorticksLess'*100);
        labelsleft = [num2str(majorticks'*100) repmat('   ', length(labels), 1)];
        labelsleftLess = [num2str(majorticksLess'*100) repmat('   ',length(labelsLess),1)];
        labelsright = [repmat('  ', length(labels), 1) num2str(majorticks'*100)];
        labelsrightLess = [repmat('  ', length(labelsLess), 1) num2str(majorticksLess'*100)];
        zerocomp = zeros(size(majorticks)); % represents zero composition
        zerocompLess = zeros(size(majorticksLess));
        
        % Plot right labels
        lyc = majorticks*sin(pi/3);
        lxc = 1-majorticks + lyc*cot(pi/3);
        lycLess = majorticksLess*sin(pi/3);
        lxcLess = 1-majorticksLess + lycLess*cot(pi/3);
        %text(lxc, lyc, labelsright,'FontWeight','bold','FontSize',20,'color',[0 0.7 0]);
        text(lxcLess, lycLess, labelsrightLess,'FontName','Arial','FontWeight','bold','FontSize',20,'color',[0 0.7 0]);
        
        % Plot bottom labels
        lyb = zerocomp;
        lybLess = zerocompLess;
        lxb = majorticks + lyb*cot(pi/3);
        lxbLess = majorticksLess + lybLess*cot(pi/3);
        %text(lxb, lyb-0.01, labels,'FontWeight','bold','FontSize',20,'color',[0.7 0 0],'VerticalAlignment','top','HorizontalAlignment','center');
        text(lxbLess, lybLess-0.01, labelsLess,'FontName','Arial',...
            'FontWeight','bold','FontSize',20,'color',[0.7 0 0],...
            'VerticalAlignment','top','HorizontalAlignment','center');
        % Plot left labels
        lya = (1-majorticks)*sin(pi/3);
        lxa = lya*cot(pi/3);
        lyaLess = (1-majorticksLess)*sin(pi/3);
        lxaLess = lyaLess*cot(pi/3);
        % text(lxa, lya, labelsleft,'FontWeight','bold','FontSize',20,'color',[0 0 0.7],'HorizontalAlignment','right');
        text(lxaLess, lyaLess, labelsleftLess,'FontName','Arial',...
            'FontWeight','bold','FontSize',20,'color',[0 0 0.7],...
            'HorizontalAlignment','right');
        nlabels = length(labels)-2;
        
        for i = 1:nlabels
            plot([lxa(i+1) lxb(nlabels - i + 2)], [lya(i+1) lyb(nlabels - i + 2)],'LineStyle',':','Color',[0.8 0.8 0.8],'linewidth',1,'handlevisibility','off');
            plot([lxb(i+1) lxc(nlabels - i + 2)], [lyb(i+1) lyc(nlabels - i + 2)],'LineStyle',':','Color',[0.8 0.8 0.8],'linewidth',1,'handlevisibility','off');
            plot([lxc(i+1) lxa(nlabels - i + 2)], [lyc(i+1) lya(nlabels - i + 2)],'LineStyle',':','Color',[0.8 0.8 0.8],'linewidth',1,'handlevisibility','off');
        end
    case 'false'
end

if nargin<5
    h = scatter(x, y, 15, colorm,'filled','Linewidth',1, 'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha',alpha);
end

if nargin >= 5
    h = scatter(x, y, markerSize, colorm, 'filled', 'Linewidth', 1, 'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha',alpha);
end

switch density
    case 'true'
        hold on
        binSize = round(sqrt(size(ternCoords,1)));
        
       if binSize > 5
           hDen = dscatter(ternCoords(:,1),ternCoords(:,2),'PLOTTYPE','contour','BINS',[binSize,binSize]);
           %     dscatter(ternCoords(:,1),ternCoords(:,2),'PLOTTYPE','contour','BINS',[500 500]);
           %dscatter(ternCoords(:,1),ternCoords(:,2),'MARKER','.','BINS',[1000
           %1000])
       end
        hold off
    case 'false'
end


hold off;

%text(0.5, -0.075,'% Red', 'horizontalalignment', 'center','FontSize',60,'FontWeight','bold','color',[0.7 0 0]);
%text(1-0.45*sin(pi/6)+0.07, 0.5,'% Green', 'rotation', -60, 'horizontalalignment', 'center','FontSize',60,'FontWeight','bold','color',[0 0.7 0]);
%text(0.45*sin(pi/6)-0.07, 0.5, '% Blue', 'rotation', 60, 'horizontalalignment', 'center','FontSize',60,'FontWeight','bold','color',[0 0 0.7]);
set(gca,'dataaspectratio',[1 1 1]);
set(gca, 'visible', 'off');



end

