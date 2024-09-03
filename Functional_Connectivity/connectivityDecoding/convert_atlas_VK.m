function oobList=convert_atlas_VK(funcFN,atlasStructFN,outFN,maxRegionNidx)
% CONVERT_ATLAS Convert structural image atlas into functional image atlas,
% all in native space
%
% IN
%   funcFN: filename of realigned functional image (typically mean image),
%       whose header will be used to write out functional atlas.
%   atlasStructFN: filename of atlased structural image, as created by
%       IBASPM. Should be in register with the functional atlas.
%   outFN: filename of output functional atlas file
%   maxRegionNidx: maximum region index to keep (e.g. 90: keep up to the
%       90th atlas region)
% OUT
%   oobList: list of out-of-boundary structural voxels that were mapped
%       to edge of functional volume
% 
% EXAMPLE USE
%   myOOB=convert_atlas('fImage.nii','sAtlas.nii','fAtlas.nii');
%
% HISTORY
% v1.0 Dimitri Van de Ville (DVDV)
% - initial release
% v1.1 Jonas Richiardi (JR)
% - function
% - more robust to out-of-boundary voxels etc
% - stats on oob voxels for quality control
% v1.1.1 DVDV + Isik Karahanoglu + JR
% - fix intensity scaling
% v1.1.2 JR
% - don't assume contiguous and/or start-at-1 atlas values
% v1.1.3 2010Dec02 JR
% - bugfix release: regions were clipped at maxRegionNidx-1 in most cases.
% v1.1.4 2011Jun13 JR
% - more informative error message in case too many regions
% v1.1.5 2011Jun21 JR
% - allow for atlases with fewer than specified max number of regions
%
% REFERENCE
% If you use this code, please cite
% Jonas Richiardi, Hamdi Eryilmaz, Sophie Schwartz, Patrik Vuilleumier,
% and Dimitri Van De Ville, Decoding Brain States from fMRI Connectivity Graphs,
% NeuroImage 56: 616-626, 2011

DEBUGMODE=false;

% Read and load functional image - just for saving time, don't care about contents
Fi=spm_vol(funcFN);
F=spm_read_vols(Fi);
if DEBUGMODE==true
    figure;hist(F(:),unique(F(:)));
    title(['Unique voxel values in source functional image ' ...
        num2str(numel(unique(F)))]);
end
% read and load atlased image in native space
Ai=spm_vol(atlasStructFN);
A=spm_read_vols(Ai);

% number of unique regions: with 116-regions atlas should be 117 if all
% regions are visible in the volume
atlasCodes=unique(A); 
atlasCodes(any(isnan(atlasCodes),2),:)=[]; % VK
% don't assume atlas values start at 1 and / or are contiguous
if numel(atlasCodes) <= maxRegionNidx
    warning(['There are fewer than ' num2str(maxRegionNidx) ... 
        ' regions in the atlas, keeping at most ' ...
        num2str(numel(atlasCodes))]);
    maxRegionNidx=numel(atlasCodes);
end
maxRegionN=atlasCodes(maxRegionNidx);
if atlasCodes(1)==0 % take into account 0
    maxRegionN=maxRegionN+1;
end
nAregions=numel(atlasCodes);
disp(['There are ' num2str(nAregions) ' unique atlas region in structural atlas']);

[x1,x2,x3]=ndgrid(1:Fi.dim(1),1:Fi.dim(2),1:Fi.dim(3));
idx=1:numel(F); % let the atlasing set the non-intracranial voxels to 0.

% take every voxel in the volume spanned by the functional images,
% compute its real-world position in mm
oobList=zeros(0,4); % list of out-of-bound structural voxels
for iter=1:length(idx)
    oob=false;
    % recover world-space position of this voxel in mm from affine
    % transform matrix
    mm=Fi.mat*[x1(idx(iter)) x2(idx(iter)) x3(idx(iter)) 1]';
    % convert this position into index of the closest structural voxel
    vx=round(Ai.mat\[mm(1) mm(2) mm(3) 1]');
    vx(vx<=0)=1;
    vxOri=vx;
    % remap out-of-bounds voxels to last  
    if vx(1)>Ai.dim(1), vx(1)=Ai.dim(1); oob=true; end
    if vx(2)>Ai.dim(2), vx(2)=Ai.dim(2); oob=true; end
    if vx(3)>Ai.dim(3), vx(3)=Ai.dim(3); oob=true; end
    if (oob==true), oobList(end+1,:)=vxOri; end
    % idx(iter): current voxel
    F(idx(iter))=A(vx(1),vx(2),vx(3));
    if any(F(idx(iter))<0)
        disp(['Negative voxel values!']);
    end
end
u=unique(F); % VK
u(any(isnan(u),2),:)=[]; % VK
nFregions=numel(unique(u)); % VK
%nFregions=numel(unique(F));

if (size(oobList,1)>0)
    disp([num2str(size(oobList,1)) ' structural voxels were out of bounds']);
    %disp(['You can inspect them with figure; '... 
    %    'scatter3(myOOB(:,1),myOOB(:,2),myOOB(:,3))']);
end

if (nFregions > nAregions)
    warning('convertAtlas:tooManyFregions',['Something is wrong, there are '...
        num2str(nFregions) ' functional regions.']);
end

if DEBUGMODE==true
figure; subplot(211); hist(F(:),unique(F)); title('Full functional atlasing');
end

disp(['There are ' num2str(nFregions) ' unique atlas region in functional atlas']);
fprintf(['Removing regions higher than ' num2str(maxRegionN) ' \n']);
F(F>maxRegionN)=0;

uu=unique(F); % VK
uu(any(isnan(uu),2),:)=[]; % VK
nUregions=numel(unique(uu)); % VK
%nUregions=numel(unique(F));

disp(['There are now ' num2str(nUregions) ' unique atlas region in functional atlas']);
if (nUregions >maxRegionNidx+1)
    warning('convert_atlas:TooManyFunctionalRegions',...
        ['There should be at most ' num2str(maxRegionNidx) ' unique' ...
        'functional atlas regions, but there are ' num2str(nUregions) ...
        ': ' mat2str(unique(F))]);
end

if DEBUGMODE==true
subplot(212); hist(F(:),unique(F)); title('Restricted functional atlasing');
end

% save functional atlas
Fi.fname=outFN;
% using the mean functional image to save volume header creation time has
% one perverse effect: scaling can be non-integer. Fix that to avoid having
% non-integer atlas intensity values
if Fi.pinfo(1)~=1
    warning('convert_atlas:changingScaling','Changing plane intensity scaling to 1');
end
Fi.pinfo(1)=1; Fi.pinfo(2)=0; % scale planes by 1, zero intensity offset, don't touch "offset into image in bytes"
%Fi.dim(4)=4;
spm_write_vol(Fi,F);

