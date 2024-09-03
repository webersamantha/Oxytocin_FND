function [pVals, sig, ci, stats]=jCorrmatSignificanceMulti2_SW(CMs1,CMs2,qVal, varargin)
% test significance of difference of correlation coefficients between
% two corrmats, controlling for multiple comparisons via FDR
% 
% IN:
%   CMs1: a nVars x n Vars x nSubjects array of correlation matrices for
%       condition 1
%   CMs2: likewise for condition 2
%   qVal: desired error-controlled p-value
%
% OUT:
%   pVals: a nVars x nVars matrix of p-values for the two-sample t-test
%   sig:  a nVars x nVars binary matrix indicating which elements in the
%           input correlation matrix are significantly different
%
% DEPENDENCIES:
% - needs the statistics toolbox
% - needs process_options, jUpperTriMatToVec
% 
% REFERENCES:
% - Logan, Rowe, An evaluation of thresholding techniques in fMRI analysis,
% NeuroImage, 2004 [DL OK]
% - Benjamini and Yekutieli, THE CONTROL OF THE FALSE DISCOVERY RATE IN
% MULTIPLE TESTING UNDER DEPENDENCY, Annals of Statistics, 2001 [DL OK]
%
% VERSION:
% 1.0 Jonas Richiardi
% - initial release
% 1.1 Samantha Weber
% - new FDR function

[plotStuff]=process_options(varargin,'plotStuff',false);



nCMs1=size(CMs1,3);
nCMs2=size(CMs2,3);

nVars1=size(CMs1,2);
nVars2=size(CMs2,2);


if (nVars1~=nVars2)
    error('Input CMs must have the same dimension in all conditions');
end

disp(['Doing ' num2str(nchoosek(nVars1,2)) ' two-sample t-tests...']);
pVals=zeros(nVars1,nVars1);
for r=1:nVars1
    for c=r:nVars1
        %h=1: reject null hypothesis of no difference -> sample HAVE a different mean 
        [h,p,ci,stats]=ttest2(CMs1(r,c,:),CMs2(r,c,:),qVal,'both','unequal');
        pVals(r,c)=p;
        pVals(c,r)=p;
    end
end
disp('...done');

%% do FDR control
% pVals_v=jUpperTriMatToVec(pVals,1);
% nComparisons=numel(pVals_v);
% %c=log(nComparisons)+0.5; 
% c=1; % FDR control
% slope=qVal/(c*nComparisons);
% pVals_v_sort=sort(pVals_v,'ascend');
% 
% % sometimes can have numerical problems with first value > 1*slope while
% % second is < 2*slope...
% curLev=1;
% for i=1:numel(pVals_v_sort)
%     if (pVals_v_sort(i)-eps) <= i*slope
%         curLev=i;
%     else
%         if (curLev==1)
%             %warning('cmMulti2:passThrough','Allowing first value to pass through despite failed test');
%             curLev=1;
%         else
%             break
%         end
%     end
% end
% qVal=pVals_v_sort(curLev);
% disp(['Controlling multiple-comparison error rate at q=' ...
%     num2str(qVal,'%0.5f') ', keeping first ' num2str(curLev) ' vals.' ]);
% 
% 
% if (plotStuff==true)
%     % plot pVals vs. FDR-slope
%     figure; plot(pVals_v_sort); hold on;
%     plot(1:nComparisons,(1:nComparisons).*slope,'r:');
%     scatter(curLev,qVal,'o');
%     legend('p - values','FDR slope','qVal');
% end
% 
% if (curLev<2)
%     %warning('Insufficient number of values kept');
% end
% 
% %sig=pVals<qVal; 
% sig=pVals<qVal;

%for d = 1:size(pVals,2)
    [pthr,pcor,padj] = fdr(pVals(:,d));
    %                    fdr(pval,q)
    %                    fdr(pval,q,cV)

  
    sig=pVals<qVal;

%end
