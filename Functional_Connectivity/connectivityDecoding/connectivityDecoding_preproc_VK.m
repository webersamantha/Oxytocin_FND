function TCS=connectivityDecoding_preproc_VK(fMean_filename, sa_filename, ...
    fa_filename,maxRegionNidx, f_direc,prefix, atlas)
% ** Connectivity decoding toolkit **
% 
% Structural and functional preprocessing for connectivity-based decoding
% 
% SYNOPSIS
% This implements the necessary preprocessing code to obtain
% region-averaged time series. After running this, feature extraction
% typically proceeds by filtering the averaged time series into frequency
% subbands of interest (wavelet transform or bandpass filter), then
% computing a correlation or other dependency matrix.
%
% IN
%   fMean_filename: string, mean realigned functional volume
%   sa_filename: string, structural atlas volume
%   fa_filename: string, functional atlas volume (will be created)
%   maxRegionNidx: scalar, maximum region index to keep (e.g. 90: keep up to the
%       90th atlas region)
%   f_direc: string, path to the functional directory containing realigned
%   images
%
% OUT
%   TCS: nRegions x nTimepoints matrix, region-averaged time courses.
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
% connected association cortical hubs. J.of Neuroscience 26(1):63???72
% 
% REQUIREMENTS
% - SPM 5 or above (http://www.fil.ion.ucl.ac.uk/spm/)
% - IBASPM toolbox
% - connectivity decoding toolkit helper functions (cdtkhelper)
%
% VERSION
% Canonical script 1.0, oct 2010
% Jonas Richiardi, Medical Image Processing Laboratory
% - initial public release. Thanks to William Pettersson-Yeo at KCL for
% beta-testing this.
% 1.0.1, 21 oct 2010
% - path corrections and file renaming
% 1.0.2, 21 oct 2010
% - function
% 1.0.3, April 2012
% - fix for variable number of regions in atlas conversion

% The algorithm is as follows, for each subject:
% 0. Inputs:
%   A series of functional volumes:              f1.nii ... fN.nii
%   A structural volume (typically T1 MPRAGE):   s.nii
% 1. Realign and reslice functional data, obtain the mean functional
%   image
%   -> fMean.nii
% 2. Co-register structural to functional data using fMean.nii as
%   the Ref. and s.nii as the Source.
% 3. Compute the structural atlas. This is done by segmenting the
%   structural image, normalising it, and using the IBASPM toolbox or any
%   other that can compute subject-specific atlases
%   -> sa.nii
% 4. Map the structural atlas into a functional atlas
%   -> fa.nii
% 5. Load all functional data in memory as a 4D matrix. 
% 6. Loop over all your atlas labels, finding voxel indices in the
%   functional atlas that match this atlas label, and average by atlas
%   region.
%
% Steps 1-3 are standard SPM operations. Code for steps 4-6 are provided
% below.

switch atlas
    case 'AAL'
        nRois=90;
    case 'Shirer'
        nRois=90;
    case 'Hammers'
        nRois=83;
    case 'Brainnetome'
        nRois=246;
    case 'Yeo'
        nRois = 17;
end

nAtlasRegions=nRois; %VK

% 4. Map the structural atlas into a functional atlas. The functional mean
% image is only used to create a functional atlas volume with the same
% characteristics.

disp('* Converting atlas...');
%convert_atlas(fMean_filename,sa_filename,fa_filename,maxRegionNidx);
convert_atlas_VK(fMean_filename,sa_filename,fa_filename,maxRegionNidx); %VK

% 5. Create a 4D matrix of all the functional data - nPixels x nPixels x nSlices x nTimePoints


%[myListofFiles,dontcare]=spm_select('List', f_direc, '^s5wclean_*.*'); % make list of all realigned functional volumes - replace this by a programmatic call if preferred
[myListofFiles,dontcare]=spm_select('List', f_direc, ['^' prefix '*.*']); % make list of all realigned functional volumes - replace this by a programmatic call if preferred
disp('* Reading headers...');
tmp=pwd; 
cd(f_direc) 
V0i=spm_vol(myListofFiles);            % load all functional volumes headers
disp(['Selected ' num2str(length(V0i)) ' functional volumes.']);
V0=zeros(V0i(1).dim(1),V0i(1).dim(2),V0i(1).dim(3),length(V0i),'single'); % initialise storage
% store all volumes as 4D matrix
disp(['* Reading functional data...']);
for iter=1:length(V0i)
    V0(:,:,:,iter)=spm_read_vols(V0i(iter));
end    
nTimePoints=length(V0i);     %  save length of time series (number of volumes) in a variable
TCS=zeros(nAtlasRegions,nTimePoints,'single');   %  we will save each averaged timecourse in this 2D matrix

% 6. Loop over all atlas labels, finding voxel indices in the functional atlas that match this atlas label.
disp('* Doing averaging...');
cd(tmp);
FAi=spm_vol(fa_filename); % read functional atlas header
FA=spm_read_vols(FAi); % read in functional atlas
[x1,x2,x3]=ndgrid(1:FAi.dim(1),1:FAi.dim(2),1:FAi.dim(3));  %  FAi is the header of the functional atlas
for al=1:nAtlasRegions  %  al: atlas label
    idx=find(FA==al);            %  find voxels that have this label
    tc=zeros(1,nTimePoints);
    % accumulate all time series of voxels that have this atlas label
    for iter2=1:length(idx)
        tmpTC=squeeze(V0(x1(idx(iter2)),x2(idx(iter2)),x3(idx(iter2)),:)); %  V0 is a matrix of timecourses (4D of nPixels x nPixels x nSlices x nTimePoint volumes)
        tc= tc+tmpTC(:)';
    end
    tc=tc/length(idx);  %  divide by its length to get average time course
    TCS(al,:)=tc;    %  copy the averaged timecourse for this region into the TCS matrix.
end

%   plot averaged time courses of first 10 regions.
%   figure;plot(TCS(1:50,:)'); xlabel('volume number'); %VK
%   ylabel('region-averaged BOLD response'); %VK

