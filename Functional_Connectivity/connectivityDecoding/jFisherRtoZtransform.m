function zprime=jFisherRtoZtransform(r)
% convert Pearson correlations (whose sampling distrubtion is not normal)
% to normally distributed variable zprime by applying Fisher's 
% z' transformation.
%
% IN: r: a matrix of Pearson r values
% OUT: z: a matrix of fisher z values, approximately normally distributed
%       and with a standard error of 1/sqrt(N-3)
%
% v1.0 Oct 2009 Jonas Richiardi

if (any(r>1 | r<-1))
    error('Please input a correlation matrix');
end

% limit range to avoid infs in mapping
% for practical purposes anything above 0.999999 is equal to that ceiling.
r(r>0.999999)=0.999999;
r(r<-0.999999)=-0.999999;
% do transform
%zprime=0.5*(log(1+r)-log(1-r));
zprime=atanh(r);