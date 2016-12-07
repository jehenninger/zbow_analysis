function [cl] = assign_cluster(rho, nneigh, index)
% ASSIGN_CLUSTER will assign cells to clusters based on deltarho()
%   Rho, delta, and nneigh are from deltarho(). index is a logical
%   matrix from determining what cells are inside polygon in Decision Graph

clust = numel(index);


%% NEW WAY
nneigh(index) = index;
nneighcopy = nneigh;

uniqueN = numel(unique(nneighcopy));

while uniqueN > clust
    nneighcopy = nneighcopy(nneighcopy);
    uniqueN = numel(unique(nneighcopy));
end

cl = zeros(numel(nneigh),1);

cl(index) = 1:clust;
cl = cl(nneighcopy);



%% OLD WAY
% % maxrho = find(max(rho(index)));
% % rho(maxrho) = rho(maxrho)+5;
% ND = length(rho);
% [rho_sorted,ordrho]=sort(rho,'descend');
% NCLUST=0;
% 
% cl = -1.*ones(ND,1);
% 
% 
% 
% 
% for ii = 1:clust
%     NCLUST = NCLUST+1;
%     cl(index(ii)) = NCLUST;
% end
% 
% %assignation
% for ii=1:ND
%   if (cl(ordrho(ii))==-1)
%     cl(ordrho(ii))=cl(nneigh(ordrho(ii)));
%   end
% end
% 
% clusterIndex = cl;
end



