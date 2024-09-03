function CM=connectivityDecoding_dependency(TCS)
% ** Connectivity decoding toolkit **
% 
% Compute correlation of region-averaged, potentially filtered time series
%
% SYNOPSIS
% This computes pairwise Pearson correlation coefficients between time
% series. The corr function could be used instead but requires the stats
% toolbox.
%
% IN
%   TCS: nRegions x T matrix of real numbers, region-averaged time series
%       (see connectivityDecoding_preproc.m,
%       connectivityDecoding_filtering.m)
%
% REFERENCE
% If you use this code please cite
% Jonas Richiardi, Hamdi Eryilmaz, Sophie Schwartz, Patrik Vuilleumier,
% and Dimitri Van De Ville, Decoding Brain States from fMRI Connectivity Graphs,
% NeuroImage 56: 616-626, 2011
%
% Where the preprocessing and connectivity matrix computation is based on
%
% S. Achard, R. Salvador, B. Whitcher, J. Suckling, E. Bullmore, A resilient,
% low frequency, small-world human brain functional network with highly
% connected association cortical hubs. J.of Neuroscience 26(1):63–72
%
% REQUIREMENTS
%   - connectivity decoding toolkit helper functions (cdtkhelper)
%
% VERSION
% 1.0, oct 2010
% Jonas Richiardi + Dimitri Van De Ville, Medical Image Processing Laboratory
% - initial public release, only pairwise Pearson correlation supported. 

[CM CMdist_ignore]=correlation_matrix(TCS);