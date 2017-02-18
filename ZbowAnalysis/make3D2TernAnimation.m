function [h] = make3D2TernAnimation(normDataTern,ternColor,ternCoords)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ternCoords(:,3) = 0;


figure,
for ii = 1:size(ternCoords,1)
    C{ii} = getInterpolatingPoints(normDataTern(ii,:),ternCoords(ii,:));
end

timePointCom = cat(2,C{:});

for nn = 1:101
    timePoint{nn} = reshape(timePointCom(nn,:),[size(ternCoords,1),3]);
end


for kk = 1:101
    
    scatter3(timePoint{kk}(:,1),timePoint{kk}(:,2),timePoint{kk}(:,3),20,'b','filled');
    
    M(kk) = getframe;
end

figure,
movie(M,5);

    function C = getInterpolatingPoints(point1, point2)
        t = 0:0.01:1;
        C = repmat(point1,length(t),1)'+(point2-point1)'*t;
        C = C';
    end

end

