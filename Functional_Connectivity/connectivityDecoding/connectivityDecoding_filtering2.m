function TCS_filt=connectivityDecoding_filtering(TCS, C_MAX, C_S)
% ** Connectivity decoding toolkit **
% 
% Compute Wavelet decomposition of region-averaged time series
%
% SYNOPSIS
% This filters time series into orthogonal subbands, which allows to
% analyse activity/connectivity at different frequencies.
% It is a thin wrapper around Dimitri Van De Ville's my_wavelet.m script,
% itself a wrapper around Michael Unser and Thierry Blu's fractional
% wavelet transform package
% (http://bigwww.epfl.ch/demo/fractsplines/matlab.html)
%
% IN
%   TCS: nRegions x T matrix of real numbers, region-averaged time series
%       (see connectivityDecoding_preproc.m)
%   C_MAX: scalar, maximum number of subbands 
%   C_S:  scalar,  choice of subband to use
%
% REFERENCE
% If you use this code please cite
%   - M. Unser, T. Blu, "Fractional Splines and Wavelets," SIAM Review, vol.
%   42, no. 1, pp. 43-67, March 2000
%	- Jonas Richiardi, Hamdi Eryilmaz, Sophie Schwartz, Patrik Vuilleumier,
%	and Dimitri Van De Ville, Decoding Brain States from fMRI Connectivity Graphs,
%	NeuroImage 56: 616-626, 2011
%
%
% REQUIREMENTS
%   - Signal Processing Toolbox
%   - connectivity decoding toolkit helper functions (cdtkhelper)
%
% VERSION
% 1.0, oct 2010
% Jonas Richiardi + Dimitri Van De Ville, Medical Image Processing Laboratory
% - initial public release. 

addpath(genpath(fullfile('.','cdtkhelper'))); % add access to necessary tools


TCS_l=size(TCS,2);
TCS_l_p2=2^ceil(log(TCS_l)/log(2)); % power of 2 above current length

% zero-mean
TCSmu=mean(TCS,2);
TCS=TCS-repmat(TCSmu,1,size(TCS,2));
%pad if needed
padLength=TCS_l_p2-TCS_l;
if(padLength~=0)
    TCS=[TCS zeros(size(TCS,1),padLength)];
end;
% filter
TCS_filt=my_wavelet2(TCS,3,0,C_MAX,C_S);
% unpad
TCS_filt=TCS_filt(:,1:end-padLength); % remove padding