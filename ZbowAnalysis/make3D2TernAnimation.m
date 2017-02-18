function [M] = make3D2TernAnimation(normDataTern,ternColor,ternCoords)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ternCoords(:,3) = 0;



for ii = 1:size(ternCoords,1)
    C{ii} = getInterpolatingPoints(normDataTern(ii,:),ternCoords(ii,:));
end

viewAngle = getInterpolatingPoints([135 50], [0 90]);

timePoint = cell(101,1);
for nn = 1:101
    timePoint{nn} = zeros(size(ternCoords,1),3);
    for jj = 1:size(ternCoords,1)
        timePoint{nn}(jj,:) = C{jj}(nn,:);
    end
end

figure,
axis manual
xlim([0 1]), ylim([0 1]), zlim([0 1]);
view(135,50);
for kk = 1:101
    
    scatter3(timePoint{kk}(:,1),timePoint{kk}(:,2),timePoint{kk}(:,3),20,ternColor,'filled');
    xlim([0 1]), ylim([0 1]), zlim([0 1]);
    view(viewAngle(kk,1), viewAngle(kk,2));
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    M(kk) = getframe;
    
end

figure,
set(gca,'Visible','off');
movie(M,1,24);


    function C = getInterpolatingPoints(point1, point2)
        t = 0:0.01:1;
        C = repmat(point1,length(t),1)'+(point2-point1)'*t;
        C = C';
    end

end

