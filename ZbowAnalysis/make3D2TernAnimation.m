function [M] = make3D2TernAnimation(normDataTern,ternColor,ternCoords)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ternCoords(:,3) = 0;


for ii = 1:size(ternCoords,1)
    C{ii} = getInterpolatingPoints(normDataTern(ii,:),ternCoords(ii,:));
end

viewAngle = getInterpolatingPoints([45 50], [0 90]);

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
view(45,50);
count = 1;

%rotation of 3D plot
scatter3(timePoint{1}(:,1),timePoint{1}(:,2),timePoint{1}(:,3),20,ternColor,'filled');
for nn = 1:90
    [az, el] = view;
    view(az+5,el);
    M(count) = getframe;
    count = count + 1;
end


%Pause after rotation
for mm = 1:15
    scatter3(timePoint{1}(:,1),timePoint{1}(:,2),timePoint{1}(:,3),20,ternColor,'filled');
    M(count) = getframe;
    count = count + 1;
end

%transform into 2D plot
for kk = 1:101
    scatter3(timePoint{kk}(:,1),timePoint{kk}(:,2),timePoint{kk}(:,3),20,ternColor,'filled');
    xlim([0 1]), ylim([0 1]), zlim([0 1]);
    view(viewAngle(kk,1), viewAngle(kk,2));
    set(gca,'XTick',[],'YTick',[],'ZTick',[]);
    M(count) = getframe;
    count = count + 1;
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

