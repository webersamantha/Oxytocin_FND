function [qVal,sig]=jFDRcorrection(pVals_v,qValTheo,varargin)
% compute FDR correction for multiple comparisons on a set of p values
%
%   IN
%       pVals_v: a vector of p-values from multiple comparisons
%       qValTheo: theoretical desired FDR-corrected p-value (typically
%       0.05)
%
%   OUT
%       qVal: the FDR-corrected significance level that could be achieved
%       sig: a locical mask vector of the same dimension as pVals_v,
%       indicating which pVals pass the FDR-corrected significance level
% 
%   REFERENCES
%       - Benjamini and Yekutieli, THE CONTROL OF THE FALSE DISCOVERY RATE
%       IN MULTIPLE TESTING UNDER DEPENDENCY, Annals of Statistics, 2001
%
%   VERSION
%   v1.0 Mar 2011 Jonas Richiardi
%   v1.0.1 Mar 2013 JR
%       - pvals <= qval (not only <) are significant


[plotStuff]=process_options(varargin,'plotStuff',false);

%% sanity check
pVals_v=pVals_v(:);
if size(pVals_v,2)>1
    error('Function can only deal with vector inputs');
end

%% do FDR control

nComparisons=numel(pVals_v);
%c=log(nComparisons)+0.5; 
c=1; % FDR control
slope=qValTheo/(c*nComparisons);
pVals_v_sort=sort(pVals_v,'ascend');

% sometimes can have numerical problems with first value > 1*slope while
% second is < 2*slope...
curLev=1;
for i=1:numel(pVals_v_sort)
    if (pVals_v_sort(i)-eps) <= i*slope
        curLev=i;
    else
        if (curLev==1)
            warning('jFDRcorrection:passthrough','Allowing first value to pass through despite failed test');
            curLev=1;
        else
            break
        end
    end
end
qVal=pVals_v_sort(curLev);
disp(['Controlling multiple-comparison error rate at q=' ...
    num2str(qVal,'%0.5f') ', keeping first ' num2str(curLev) ' vals.' ]);

if (plotStuff==true)
    % plot pVals vs. FDR-slope
    figure; plot(pVals_v_sort); hold on;
    plot(1:nComparisons,(1:nComparisons).*slope,'r:');
    scatter(curLev,qVal,'o');
    legend('p - values','FDR slope','qVal');
end

sig=pVals_v<=qVal;