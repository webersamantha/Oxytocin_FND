function W=jWinsorizeMatrix(M,pct)
% Winsorize a matrix
% IN
%   M: N x T, a matrix of N observations / subjects with T points
%   (typically in time)
%   pct: a struct with fields
%       .up: upper percentile (1-99)
%       .low: lower percentile (1-99)
% OUT
%   W: N x T, the same matrix winsorized to within specified percentile
%
% v1.0 Nov 2010 Jonas Richiardi
%

% sanity check
if pct.low < 0.01 || pct.up > 99.99 || pct.low >= pct.up
    error('Lower and upper percentile bounds must be between 0.1 and 99.9, and low < up');
end

[N T]=size(M);
W=M;

% compute percentile thresholds for all rows in matrix
allThresh=prctile(M,[pct.low pct.up],2);


for n=1:N
    W(n,W(n,:)>allThresh(n,2))=max(W(n,W(n,:)<=allThresh(n,2)));
    W(n,W(n,:)<allThresh(n,1))=min(W(n,W(n,:)>=allThresh(n,1)));
end