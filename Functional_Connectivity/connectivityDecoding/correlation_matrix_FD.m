function [CM, CMdist]=correlation_matrix(TCS,FD,s)
% compute a pairwise pearson linear correlation matrix on voxel time courses
%
% IN -  TCS:    a nRegions x nTimepoints matrix of voxel time courses.
% OUT - CM:     a nRegions x nRegions (symmetric) matrix of pearson
%                   correlations
%       CMdist: a nPairs (nchoosek(nRegions,2)) x 1 vector of pairwise
%                   correlation coefficients
% v1.0 Dimitri van de Ville
% - initial release

for i = 1:90
    for k = 1:150
        if FD{s,1}(k,1) == 0
        TCS(k,i)= 0;
        end 
    end 
end

TCS_new = TCS(any(TCS,2),:);

CMdist=[]; %bCMdist=[];
%v2.0 SW
for iter=1:size(TCS_new,2) 
    TCS_new(:,iter)=TCS_new(:,iter)-mean(TCS_new(:,iter));
end


% compute crossproducts divided by products of variances
for iter1=1:size(TCS_new,2)
    CM(iter1,iter1)=1;
    for iter2=iter1+1:size(TCS_new,2)
        tmp=sum(TCS_new(:,iter1).*TCS_new(:,iter2))/((sum(TCS_new(:,iter1).^2).*sum(TCS_new(:,iter2).^2))^(0.5));
        if isnan(tmp)
            tmp=0;
        end
        CM(iter1,iter2)=tmp;
        CM(iter2,iter1)=tmp;
        CMdist=[CMdist tmp];
    end
% % compute crossproducts divided by products of variances
% for iter1=1:size(TCS,2)
%     CM(iter1,iter1)=1;
%     for iter2=iter1+1:size(TCS,2)
%         tmp=sum(TCS(iter1,:).*TCS(iter2,:))/((sum(TCS(iter1,:).^2).*sum(TCS(iter2,:).^2))^(0.5));
%         if isnan(tmp)
%             tmp=0;
%         end
%         CM(iter1,iter2)=tmp;
%         CM(iter2,iter1)=tmp;
%         CMdist=[CMdist tmp];
%     end

% % v1.0 DvdV
% for iter=1:size(TCS,1), 
%     TCS(iter,:)=TCS(iter,:)-mean(TCS(iter,:),2);
% end;
% 
% 
% % compute crossproducts divided by products of variances
% for iter1=1:size(TCS,1),
%     CM(iter1,iter1)=1;
%     for iter2=iter1+1:size(TCS,1),
%         tmp=sum(TCS(iter1,:).*TCS(iter2,:),2)/((sum(TCS(iter1,:).^2,2).*sum(TCS(iter2,:).^2,2))^(0.5));
%         if isnan(tmp),
%             tmp=0;,
%         end;
%         CM(iter1,iter2)=tmp;
%         CM(iter2,iter1)=tmp;
%         CMdist=[CMdist tmp];
%     end;
end
