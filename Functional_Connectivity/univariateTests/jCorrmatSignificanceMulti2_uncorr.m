function [pVals, sig]=jCorrmatSignificanceMulti2_uncorr(CMs1,CMs2,qVal, varargin)
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

[plotStuff]=process_options(varargin,'plotStuff',true);



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
        [h,p]=ttest2(CMs1(r,c,:),CMs2(r,c,:),qVal,'both','unequal');
        pVals(r,c)=p;
        pVals(c,r)=p;
    end
end
disp('...done');


sig=pVals<qVal;

% kill insignificant vals in corrmats
%CMout=zeros(size(CMs1));
%for c=1:nCMs
%    CMout(:,:,c)=CMs(:,:,c).*sig;
%end
