function [rho, delta, nneigh] = deltarho(data, percent)
%DELTARHO Finds the density (rho) and distance from high densities (delta)
%   "Data" is a distance matrix of the normalized intensities
%   "Percent" is relative to the average number of neighbors
xx = pdist(data);
sda = sort(xx);

N = size(xx,2);

position=round(N*percent/100);
clear N

dc=sda(position);

clear sda
clear position

dist = squareform(xx);
%Zxx = squareform(xx);

clear xx




%% NEW WAY
% %Gaussian kernel to find rho
rhoElement = exp(-(dist/dc).*(dist/dc));
rho = sum(rhoElement,1)'-1;
clear rhoElement

% Dc method
% rho = dist - dc;
% rho(rho>=0) = 0;
% rho(rho<0) = 1;
% rho = sum(rho,2);




%Calculate delta for each point

delta1 = repmat(rho,[1 length(rho)]); %Make rho vector into rho x rho matrix
delta2 = bsxfun(@minus,delta1,rho'); %Subtract each column by rho(i) of that column
clear delta1
delta3 = delta2<=0; %Take all values less than zero, which are the points with lower density than rho(i)
clear delta2
%delta4 = dist; %make a copy of dist
dist(delta3) = NaN; %Set points with densities less than rho(i) to NaN
clear delta3
dist(dist == 0) = NaN; %Set diagonals, which are 0, to NaN to avoid issues with min()
[delta, nneigh] = nanmin(dist); %Find min of each column in the dist matrix of points with higher densities than rho(i)
                                   %and get indices of the min values,
                                   %which are the nearest neighbors

% %Define neighbor of max density point
[~,maxIndex] = max(rho);
nneigh(maxIndex) = 0;
delta(maxIndex) = max(delta(:));


%% OLD WAY
% Gaussian kernel
% ND = length(data(:,1));
% rho = zeros(ND,1); %no need to preallocate anymore since we aren't using a
% % %for loop
% 
% for i=1:ND-1
%   for j=i+1:ND
%      rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
%      rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
%   end
% end
% [~,ordrho]=sort(rho,'descend');
% delta(ordrho(1))=-1.;
% nneigh(ordrho(1))=0;
% 
% maxd=max(dist(:)); 
% 
% for ii=2:ND
%    delta(ordrho(ii))=maxd;
%    for jj=1:ii-1
%      if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
%         delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
%         nneigh(ordrho(ii))=ordrho(jj);
%      end
%    end
% end
% delta(ordrho(1))=max(delta(:));



end

