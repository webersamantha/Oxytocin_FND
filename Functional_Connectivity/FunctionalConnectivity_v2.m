%% FUNCTIONAL CONNECTIVITY DECODING 

% 1. Set Parameters
% 2. Preprocessing (already done)
% 3. Extraction of average WM and CSF signal. Extraction of grey matter signal
% 4. Build correlation matrices
% 5. Fisher-Z transform matrix

% Newest version implemented November 2023 by Samantha Weber
% This code uses functions from Dimitri Van De Ville and Jonas Richiardi 

% v1.0 Samantha Weber (March 2021, GCB Symposium)
% v2.0 Samantha Weber (March 2022, CAPs Paper)
% v3.0 Samantha Weber (November 2023, BioGen Teaching)
 
clear all
% Necessary: Download SPM12 tool on https://www.fil.ion.ucl.ac.uk/spm/software/download/
% Download cbrewer2
% Download BrainNetViewer

%% Load paths
%Scripts path
scriptsPath = '/Users/samanthaweber/Documents/GitHub/BioGen_Oxytocin_SW/02_Analysis/Functional_Connectivity';
addpath (genpath(scriptsPath));
%spm Path
spmPath = '/Users/samanthaweber/spm12';
addpath (spmPath);
%RootPath
RootPath= '/Volumes/T7_Samantha/Analyses/Oxytocin/00_Data/fMRI Data';
addpath (genpath ('RootPath'));%(3)
%Atlas path
atlasPath= fullfile(scriptsPath,'Atlases');
% % Mask Path
% masksPath=fullfile(scriptsPath,'preprocessing','Masks');
% addpath (genpath(masksPath));

% Correlation Matrices path
CMpath='/Volumes/T7_Samantha/Analyses/Oxytocin/02_ConnectivityAnalyses/01_CorrelationMatrices';
if ~exist(CMpath),mkdir(CMpath); end

%% 1. Set parameters
% Load subjects
Labelling = 'MNI'; % 'native' or 'MNI'
myGroups={'HC','FND'};
% load P-codes HC;
    filelist=dir(fullfile(RootPath, myGroups{1},'P*'));
    N_HC=size(filelist,1);
    myHC=cell(1,N_HC);
    for f=1:N_HC
        myHC{f}=filelist(f).name;
    end
    clear filelist f

%load P-codes FND
    filelist=dir(fullfile(RootPath, myGroups{2},'P*'));
    N_FND=size(filelist,1);
    myFND=cell(1,N_FND);
    for f=1:N_FND
        myFND{f}=filelist(f).name;
    end
    clear filelist f

mySubjects = [myHC,myFND];
N_HC = length(myHC); N_FND = length(myFND); nSubj=N_HC+N_FND;
       
%% Define all parameters

% I have the full codebook for AAL, that's why it's easiest to use the AAL
% if you want to visualize it with the BrainNetViewer. Also, brainnetviewer
% has AAL coordinates included already!

atlas = 'AAL' ; %choose between:  'AAL', 'Shirer', 'Hammers', 'Yeo'

% Set the atlas
switch atlas
    case 'AAL'
        nRois=90; %ROI = region of interest
        atlasFile='AAL116_correctLR.nii';
    case 'Shirer'
        nRois = 90;
        atlasFile='ShirerAtlas.nii';
    case 'Hammers'
        nRois=83;
        atlasFile='Hammers_mith_atlas_n30r83_SPM5.nii';
    case 'Brainnetome'
        nRois=246;
        atlasFile='BN_Atlas_246_1mm.nii';
    case 'HCP-ICA'
        nRois=32;
        atlasFile='networks.nii';
    case 'Yeo'
        nRois = 17;
        atlasFile = 'Yeo2011_17Networks_MNI152_1mm_LiberalMask.nii';
end

thisAtlas=fullfile(atlasPath,atlas,atlasFile);

%Define the folders --> this only works, if you have ALWAYS the same
%structure. You can of course have different group names, but you have to
%define those in the myGroups = {} above. 

%Data
    % FND
        %P...
        %P...
    % HC
        %P...
        %P...

for i = 1:N_HC
    folders{i} = fullfile(RootPath,myGroups{1},mySubjects{i});
end

for i = N_HC+1:nSubj
    folders{i} = fullfile(RootPath,myGroups{2},mySubjects{i});
end


%% 2. Preprocessing

% Already done using same pipeline as for CAPs analyes or as discussed
% previously. 
% https://github.com/FND-ResearchGroup/CAPs_SW/tree/main/01_CAPs_Preparation

%% 3. Extract average signal of WM and CSF, motion, and region-averaged timecourses

% The region-averaged time-courses are extracted (see INPUT/OUTPUT description below)

%INPUT: 
%        One previously cleaned functional file (s5wclean_Det_rf*) -->
%        native space would be something like clean_Det_rf*
%        If you work in native space, it will warp the atlas into native
%        space using the inverse deformation map, else it will just use the
%        MNI template atlas.
%        Atlas file
%
%OUTPUT: x = 3- dimensional matrix
%            [timepoints x 6 movement parameter x Number of Subjects]
%        y = Matrix containing the extracted BOLD signal for all 90 Regions of Interest (ROIs)for each Subject
%            [timepoints x ROI x Number of Subjects]

for s=1:nSubj
    disp(['***Extraction of average Signal from Subject ' num2str(s) ' out of ' num2str(nSubj) ' Subjects***']);
    % Define folders
    structPath=fullfile(folders{s},'struct');
    functPath=fullfile(folders{s},'funct'); 

    %%  Extract region-averaged time series (actual BOLD signal)%%
    %   The BOLD signal gets extracted from the whole brain (motion and WM/CSF is not yet regressed out!)
    %   The function 'connectivityDecoding_preproc_VK' is used:
    %       1. Map the structural atlas into a functional atlas -> fa*.nii
    %       2. Load all functional data in memory as a 4D matrix. 
    %       3. Loop over all your atlas labels, finding voxel indices in the functional atlas that match this atlas label, and average by atlas region.
    %   For further information, open the function!
    %
    %   INPUT
    %       fMean_filename: string, mean realigned functional volume
    %       sa_filename: string, structural atlas volume
    %       fa_filename: string, functional atlas volume (will be created)
    %       maxRegionNidx: scalar, maximum region index to keep (e.g. 90: keep up to the 90th atlas region)
    %       f_direc: string, path to the functional directory containing realigned images
    %
    %   OUTPUT
    %       TCS: nTimepoints matrix x nRegions, region-averaged time courses. (Time course series)
    %       y: saves TCS in a 3-dimensional matrix for all the subjects [ntimepoints x nROI x nSubjects]
    
    switch Labelling
        case 'MNI'
            prefix = 's5wclean';
            fMean_file=dir(fullfile(functPath,'s5wclean_Det_rf*.nii')); %Files have previously been cleaned using Preprocessing_RS_long_SW_v1
            fMean_filename=fullfile(functPath,fMean_file(1).name);
            sa_filename=fullfile(atlasPath,atlas,atlasFile);
            fa_filename=fullfile(folders{s},['FAtlas_' atlas '.nii']); %fa = functional atlased
            TCS=connectivityDecoding_preproc_VK(fMean_filename,sa_filename,fa_filename,nRois,functPath,prefix,atlas); %TCS = region-averaged time-series
            y(:,:,s)=TCS(:,:)';
            
        case 'native'
            prefix = 'clean';
            fMean_file=dir(fullfile(functPath,'clean_Det_rf*.nii')); %Files have previously been cleaned using Preprocessing_RS_long_SW_v1
            fMean_filename=fullfile(functPath,fMean_file(1).name);
            sa_file=dir(fullfile(structPath,['c1*' atlas '*.nii'])); %sa = structural atlased
            sa_filename=fullfile(structPath,sa_file(1).name);
            fa_filename=fullfile(folders{s},['FAtlas_native_' atlas '.nii']); %fa = functional atlased
            TCS=connectivityDecoding_preproc_VK(fMean_filename,sa_filename,fa_filename,nRois,functPath,prefix, atlas); %TCS = region-averaged time-series
            y(:,:,s)=TCS(:,:)';
    end
end
y(isnan(y))=0; %for problems with Shirers to exclude regions that overlap

clear pathRS mot RP TCS s fMean_file fMean_filename sa_dir sa_file sa_filename fa_filename masksPath i 

% Save preprocessed time series (y). You can uncomment the save or the load
% functions if you want to save intermediate steps. I normally do that, as
% then you can just start the code at a later stage and don't need to run
% the full connectivity decoding again. 

switch Labelling
    case 'MNI'
        %save(fullfile(CMpath,['preprocessedTimeSeries_' atlas '.mat']))
        %load((fullfile(CMpath,['preprocessedTimeSeries_' atlas '.mat'])))
    case'native'
        %save(fullfile(CMpath,['preprocessedTimeSeries_native_' atlas '.mat']))
        %load((fullfile(CMpath,['preprocessedTimeSeries_native_' atlas '.mat'])))
end

%% Grouping into RSNs
%not possible for Hammer, AAL, but already done for Yeo

% This part will only be done if you work with Shirer atlas. As I think the
% Shirer atlas is crap compared to the Yeo, I would anyway not use it. 
y_new = [];
switch atlas 
    case 'Shirer'
    % ROIs are grouped into 14 resting-state networks (RSNs) based on Shirer atlas
    headers = {'Auditory', 'Basal Ganglia', 'dDMN','high Visual','Language','LECN','Sensorimotor','postSalience','Precuneus'...
    ,'Primary Visual','RECN','Salience','vDMN','Visuospatial'}; % just to have it saved somewhere
        for i = 1:nSubj
            %Auditory
            y_new(:,1,i) = mean(y(:,1:3,i),2);
            % Basal Ganglia
            y_new(:,2,i) = mean(y(:,4:8,i),2);
            % dDMN
            y_new(:,3,i) = mean(y(:,9:17,i),2);
            %high Visual
            y_new(:,4,i) = mean(y(:,18:19,i),2);
            %Language
            y_new(:,5,i) = mean(y(:,20:26,i),2);
            %LECN
            y_new(:,6,i) = mean(y(:,27:32,i),2);
            %Sensorimotor
            y_new(:,7,i) = mean(y(:,33:38,i),2);
            %postSalience
            y_new(:,8,i) = mean(y(:,39:50,i),2);
            %Precuneus
            y_new(:,9,i) = mean(y(:,51:54,i),2);
            %Primary Visual
            y_new(:,10,i) = mean(y(:,55:56,i),2);
            %RECN
            y_new(:,11,i) = mean(y(:,57:62,i),2);
            % Salience
            y_new(:,12,i) = mean(y(:,63:69,i),2);
            %vDMN
            y_new(:,13,i) = mean(y(:,70:79,i),2);
            %Visuospatial
            y_new(:,14,i) = mean(y(:,80:90,i),2);
                 
        end
    nRois = 14;
y_net = y_new;    
end
y_net = y;

%% 4. Build correlation matrices
%   Correlation Matrix (based on the function 'connectivityDecoding_dependency' by Dimitri van de Ville)
%   The function computes pairwise Pearson correlation coefficients between time series
%   In case 'wavelets' filter is used, pairwise Pearson correlation coefficients will be calculated for each subband
%   INPUT: 
%        CorrTCS:       nRegions x nTimepoints matrix of voxel time courses
%   OUTPUT
%        CM_filt:       CORRELATION MATRIX (nRegions x nRegions (symmetric) matrix of pearson correlations for each subband)
%
%   Correlation Matrix will be saved as .mat file ('CorrelationMatrices_' filtering '_' atlas '_' mySubjects{s} '.mat'])

for s=1:N_HC
    % Put timeseries regressed out for motion parameters and WM/CSF signal into CorrTCS
    CorrTCS=y_net(:,:,s);   
    % Build correlation matrices for all frequency subbands
    
    disp('*** Computing dependencies...');
    CM_filt = connectivityDecoding_dependency(CorrTCS);
    %Save CM of all Subjects in a 3D Matrix
    my3DCMs(:,:,s)= CM_filt;
    %Save Correlation matrix
    save(fullfile(CMpath,['CorrelationMatrices_' atlas '_' mySubjects{s} '_.mat']),'CM_filt')
end

for s=N_HC+1:nSubj
    % Put timeseries regressed out for motion parameters and WM/CSF signal into CorrTCS
    CorrTCS=y_net(:,:,s);   
    % Build correlation matrices for all frequency subbands
    
    disp('*** Computing dependencies...');
    CM_filt = connectivityDecoding_dependency(CorrTCS);
    %Save CM of all Subjects in a 3D Matrix
    my3DCMs(:,:,s)= CM_filt;
    %Save Correlation matrix
    save(fullfile(CMpath,['CorrelationMatrices_' atlas '_' mySubjects{s} '_.mat']),'CM_filt')
end

clear structPath functPath folders i mvtParams 
disp('****Functional Connectivity Calculation: DONE****');

%% If exist you can load here --> uncomment to load
% In this case, you just run the very first parts of the code, where you
% define your paths etc (until line 120) and then you run the parts here
% below manually! 

for s=1:N_HC
    load(fullfile(CMpath,['CorrelationMatrices_' atlas '_' mySubjects{s} '_.mat']),'CM_filt')
    my3DCMs(:,:,s)= CM_filt;
end

for s=N_HC+1:nSubj
    load(fullfile(CMpath,['CorrelationMatrices_' atlas '_' mySubjects{s} '_.mat']),'CM_filt')
    my3DCMs(:,:,s)= CM_filt;
end

%% 5. Fisher-Z transform connectivity matrices & vectors (for bandpass filtered data)
%  make them gaussian with fisher Z transformation

% Convert Pearson correlations (whose sampling distrubtion is not normal)to
% normally distributed variable zprime by applying Fisher's z' transformation.
%
% INPUT: 
%    my3DMs:           a matrix of Pearson r values (correlation coefficients)
% OUTPUT: 
%    my3DCMs_fisherZ:  a matrix of fisher z values, approximately normally distributed and with a standard error of 1/sqrt(N-3)

% Remove all negative values due to their controversial interpretation [Qian J, et al. Positive connectivity predicts the dynamic intrinsic % topology of the human brain network. Front Syst Neurosci. 2018;12:38.]
% As said, this is controversial. I started removing negative values
% (uncomment below), but then I had reviewers saying I shouldn't do that...
% You can see and try what fits you best.. 


% % for s = 1: nSubj
% %    for r = 1:nRois 
% %        for c = 1:nRois
% %            if my3DCMs(r,c,s) < 0
% %                my3DCMs(r,c,s) = 0;
% %            end
% %        end
% %    end
% % end

% Fisher Z transform
for s=1:nSubj
    thisCM=my3DCMs(:,:,s);
    myNewCM=jFisherRtoZtransform(thisCM);
    my3DCMs_fisherZ(:,:,s)=myNewCM; 
end


% Make diagonale 0 (for better plotting, alternatively make it 1, for
% better plotting)
for s = 1: nSubj
   for r = 1:nRois 
     my3DCMs_fisherZ(r,r,s) = 0;
   end
end   

% Here we save the correlation matrices all in one .mat file. Which we can
% also just laod again, if you uncomment line 358. 
save(fullfile(CMpath,[atlas '_3DCMs_fisherZ.mat']),'my3DCMs_fisherZ');
load(fullfile(CMpath,[atlas '_3DCMs_fisherZ.mat']),'my3DCMs_fisherZ');
disp('****Build Functional Connectivity 3D Matrix: DONE****');

%% 6. Statistics
% Compare mean connectivity matrices between groups using t-tests or nonparametric Mann-U Whitney tests 

CM_all = mean(my3DCMs,3);
%Names for Yeo, we don't do that e.g. for AAL, because that would create a
%mess. 
switch atlas
    case 'Yeo'
        names = {'VisCen','VisPer','SomMotA','SomMotB','DorsAttnA','DorsAttnB','Sal/VenAttnA','Sal/VenAttnB',...
        'LimbicA','LimbicB','ContA','ContB - ECN','ContC','TempPar','DefaultA - DMN','DefaultB','DefaultC'};
    case 'AAL'
        names = [];
        for i = 1:90
            names{i} = num2str(i);
        end
end


[cmap] = cbrewer2('seq','Spectral',100); %https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
%figure; imagesc(CM_all); colorbar; colormap(cmap) ;
CMs1 = my3DCMs_fisherZ(:,:,1:N_HC);
figure; imagesc(mean(CMs1,3)); colorbar; colormap(flipud(cmap));
switch atlas 
    case 'Yeo'
        set(gca,'xtick',[1:17],'xticklabel',names,'FontSize',12);
        set(gca,'ytick',[1:17],'yticklabel',names,'FontSize',12);
        xtickangle(90);
end
axis('square','on');
title('HC Average Functional Connectivity');
CMs2 = my3DCMs_fisherZ(:,:,N_HC+1:end);
figure; imagesc(mean(CMs2,3)); colorbar; colormap(flipud(cmap)) ;
switch atlas 
    case 'Yeo'
        set(gca,'xtick',[1:17],'xticklabel',names,'FontSize',12);
        set(gca,'ytick',[1:17],'yticklabel',names,'FontSize',12);
        xtickangle(90);
end
axis('square','on');
title('FND Average Functional Connectivity');

% Seperate into FND and HCs - you can save it again, but not needed forlater
CMs1 = my3DCMs_fisherZ(:,:,1:N_HC);
CMs1(isnan(CMs1))=0;
save(fullfile(CMpath,[atlas '_3DCMs_fisherZ_HC.mat']),'CMs1');

CMs2 = my3DCMs_fisherZ(:,:,N_HC+1:end); 
CMs2(isnan(CMs2))=0;
save(fullfile(CMpath,[atlas '_3DCMs_fisherZ_FND.mat']),'CMs2');

%filename = fullfile(CMpath,'CorrelationMatrix_M8.xlsx');
%CMs2_full =[headers; num2cell(mean(CMs2,3))];
%xlswrite(filename,CMs2_full)

%% Testing

qVal=0.05; % alpha Level. If you use FDR correction (which you should), 0.05 is good. 

% Choose univariate test & correction for multiple comparisons
% This part of the code gives you a significance matrix (sig), a pVals
% matrix, a p_FDR matrix, those are the one we use for visualization. 

test='parametric'; % 'parametric or nonparametric
FDRcorr='yes'; % yes or no
switch test
    case 'parametric'
        switch FDRcorr
            case 'yes'
                %[pVals,sig,ci,stats]=jCorrmatSignificanceMulti2(CMs1,CMs2,qVal);
                [pVals,sig,ci,stats,p_FDR]=jCorrmatSignificanceMulti2_SW(CMs1,CMs2,qVal,nRois); %better FDR function, implemented by Samantha Feb 2022
            case 'no'
                [pVals,sig]=jCorrmatSignificanceMulti2_uncorr(CMs1,CMs2,qVal);
        end        
    case 'nonparametric'
        switch FDRcorr
            case 'yes'
                [pVals,h,stats]=ranksum2(CMs1,CMs2,qVal);
            case 'no'
                [pVals,h]=ranksum2_uncorr(CMs1,CMs2,qVal);
        end
end


%% Check outputs

% 1. Check if something is significant
max(max(sig)) % if 0 then nothing is significant

switch FDRcorr
    case 'yes'
        sig_mat = jVecToSymmetricMat(sig, nRois); % binary matrix of where it is significant
    case 'no'
        sig_mat = sig;
end

% 2. Visualise connections that are significantly different between your
% groups
newPvals=pVals;
newPvals(isnan(newPvals))=0;
%index=find(newPvals~=0);
index = nRois*nRois;
%3. Plot it --> this doesn't really help for understanding the
%significance!
figure; imagesc(newPvals); colorbar; colormap(flipud(hot)) ;

%4. Remove lower triangular part (matrix is symmetrical, better visible as such)
U = triu(newPvals,1);
% For better plotting
for r=1:index
    if U(r)== 0
        U(r)=1;
    end
end

[cmap2] = cbrewer2('seq','BuPu',50); %https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
% Plot if (different colours)
figure; imagesc(U); colorbar; colormap(flipud(cmap2)) ; %halve 
switch atlas 
    case 'Yeo'
        set(gca,'xtick',[1:17],'xticklabel',names);
        set(gca,'ytick',[1:17],'yticklabel',names);
        xtickangle(90);
end
axis('square','on');
title('HC - FND p-values');
%figure; imagesc(newPvals); colorbar; colormap(hot) ;


% Change all p-values higher than qVal to 1 (you can also just plot sig_mat)
for r=1:index%length(index)
    if newPvals(r)== 0.05
        newPvals(r)=1;
    else
    end
end


V = triu(newPvals,1);

% Change all values of 0 to 1 (for better plotting)
for r=1:index%length(index)
    if V(r) > qVal
        V(r)=1;
    else
    end
end

% plots plots plots
figure; imagesc(V); colorbar; colormap(flipud(cmap2)); %halve & shows you where it's significant (without FDR)

switch FDRcorr
    case 'yes'
        % FDR-thresholded p-value (binary)
        figure; imagesc(p_FDR); colorbar; %what's in yellow are those connections that survived FDR correction

        % Remove lower triangular part
        W = triu(p_FDR,1);
        %Plot it
        figure; imagesc(W); colorbar;colormap(flipud(cmap2)); %what's in white are those connections that survived FDR correction
end

% Close all figures
close all
clear sig

%% Visualization using BrainNetViewer

% Load codebook to know names regions that you'll be looking at
% This works for AAL only, can be adapted for other atlases, just change
% codebook, also adapt below the lobes etc. 
load(fullfile(atlasPath,atlas,'annotatedStructTemplate116codeBook.mat'));
myCB=codeBookStructTemplate.full;

% Find non-zero elements
[row, col] = find(W ~= 0); %For FDR corrected
[row, col] = find(sig_mat ~= 0); %For non-FDR corrected

for r=1:length(row)
        disp([myCB.name{row(r)} ' - ' myCB.name{col(r)} ' with mean FC in HC =' num2str(mean(CMs1(row(r),col(r)))) ' with mean FC in FND =' num2str(mean(CMs2(row(r),col(r)))) 'with p-value = ' num2str(pVals(row(r),col(r))) ]); % display
end


%% Make edge and node files for BrainNetViewer
FinalPath = '/Volumes/T7_Samantha/Analyses/Oxytocin/02_ConnectivityAnalyses/02_Plotting';
if ~exist(FinalPath),mkdir(FinalPath); end

cd(FinalPath);

% We make p_FDR diagonale 0, too

for r = 1:nRois
    p_FDR(r,r) = 0;
end

% Write edge file
pVals_new = 1 - pVals; %so that in brainnet viewer, you can only show singificant regions (you put threshold > 0.95, which corresponds then to 0.05)
for r = 1:nRois
    for c = 1:nRois
        if p_FDR(r,c) == 0
            pVals_new(r,c) = 0;
        end
    end
end
dlmwrite('P_FDR_HCvsFND.edge',pVals_new,'\t');


[row,col] = find(pVals_new(:,:) ~= 0);

for r = 1:nRois
    for c = 1:nRois
        if pVals_new(r,c) ~= 0
            pVals_new(r,c) = 1;
        end
    end
end
pVals_new(88,28)=2;
pVals_new(28,88)=2;
dlmwrite('P_FDR_HCvsFND.edge',pVals_new,'\t');

%Example creating Node (from BrainNetViewer)
% tmp = char(L);
% node = strcat(num2str(R),tmp);
% node(:,end-4:end-3) = char(' ');
% dlmwrite('test.node',node,'delimiter','');

% Making Nodes
node = zeros(90,6);
%Defining Coordingates
for r = 1:nRois
    node(r,1) = myCB.center{1,r}(1,1);
    node(r,2) = myCB.center{1,r}(2,1);
    node(r,3) = myCB.center{1,r}(3,1);
end
%Color brain according to lobe
for l = 1:16
    node(l,4) = 1; %frontal
end
for l=17:18
    node(l,4) = 2; %central
end
for l = 19: 28
    node(l,4) = 1; %frontal
end
for l = 29:42
    node(l,4) = 3; % limbic
end
for l =43:56
    node(l,4) = 4; % occipital
end
for l=57:70
    node(l,4) = 5; %parietal
end
for l = 71:78
    node(l,4) = 6; % subcortical
end
for l = 79:90
    node(l,4) = 7; %Temporal
end

% Add nodal degree (= number of edges connected to the node)
freq = zeros(90,1); % This is 90 because of AAL, if you use other atlas, adapt!

for k=1:nRois
  freq(k)=sum(pVals_new(:,k));
end

node(:,5) = freq(:);
cd(FinalPath);
% Add labels and convert to node file
tmp = char(myCB.sname);
node = strcat(num2str(node),tmp);
node(:,end-13:end-12) = char(' '); %longest character from the end is 13 long
dlmwrite('P_FDR_HCvsFND.node',node,'delimiter','');
 

%% Correlation with a scalar. 
% This has been used to correlate the Oxytocin Data connection-wise with
% the FC data. In this case, you have to load an excel file that has the
% same names of the p-codes (I had to change them below, in your case that
% won't be necessary anymore). If you understood the part above, this
% should be straight forward. 

T = readtable('/Volumes/T7_Samantha/Analyses/Oxytocin/00_Data/Liste_Oxytocin_Analysis_SamanthaWeber_T1.xlsx',"ReadRowNames", false);


for i = 1:length(oldValue)
    % Find the indices where the specified value occurs
    [row, col] = find(strcmp(T.ID, oldValue(i)));

    % Replace the value in the specified location
    T.ID{row, col} = newValue{i};

end

% Uncomment this if you want to use RowNames as true
% for i = 1:length(oldValue)
%     % Find the index of the row with the specified name
%     rowIndex = find(strcmp(T.Properties.RowNames, oldValue{i}));
% 
%     % Replace the row name
%     T.Properties.RowNames{rowIndex} = newValue{i};
% end

% -- Now we reorder our Oxytocin data so that it is in the same order like the fMRI data

% Find the indices of the rows in the desired order
[~, orderIndices] = ismember(mySubjects', T.ID);

% Remove rows with zero indices (not found in the original table)
orderIndices(orderIndices == 0) = [];

% Reorder the table based on the desired order
orderedTable = T(orderIndices, :); %Okei, we only have 117 left of the 124, we need to find out why. 
T = orderedTable;

% Find row names in the table that do not exist in the vector
missingInVector = setdiff(T.ID, mySubjects');
disp(missingInVector); %Those subjects make totally sense because they were all excluded for fMRI analyses. So all fine. 

% We have to remove these subjects now from the oxytocin data (because no fMRI exists
[~, rowsToRemove] = ismember(missingInVector, T.ID);
% Remove the specified rows
T(rowsToRemove, :) = [];


%% Calculate correlation between whole-brain FC and OT

% Initialize a matrix to store correlation results
corrResults_HC = zeros(size(CMs1,1), size(CMs1,1));
corrPVals_HC = zeros(size(CMs1,1), size(CMs1,1));

% FOR HC
% Loop through each connection of each subject and calculate correlation
x = T.OxytocinPg_ml((1:N_HC),:); % for computational speed we keep x out of the loop as it doesn't change, y will always be updated

for n_row = 1:size(CMs1,1)
    for n_col = 1:size(CMs1,1)
        %Extract 
        y = squeeze(CMs1(n_row, n_col,:));
        
        % Calculate the correlation coefficient and its p-value
        [correlationCoeff, pValue] = corrcoef(x, y);
        
        % the output has the following structure: 
        % [ correlation(x, x)  correlation(x, y)
        % correlation(y, x)  correlation(y, y) ]
        % This is why we need to access either [1,2] or [2,1]
        
        % Calculate the correlation between the scalar value and each cell of the matrix
        corrResults_HC(n_row, n_col) = correlationCoeff(1,2);
        corrPVals_HC(n_row, n_col) = pValue(1,2);
    end
end
% Replace NaN values with 0
corrPVals_HC(isnan(corrPVals_HC)) = 0;
corrResults_HC(isnan(corrResults_HC)) = 0;

% We make now a new matrix, that only shows those correlation coefficients
% which are statistically significant. We replace all > 0.05 with 1

% Threshold value
threshold = 0.05;
newValue = 1;

% Update cells in results where corresponding cells in pvals are smaller than the threshold
corrResults_HC(corrPVals_HC > threshold) = newValue; % Replace 'newValue' with the value you want to assign

% Let's visualize it
% Remove lower triangular part of correlation matrix (for better visability) and then plot
corrResults_HC_up = triu(corrResults_HC,1);

% Change all values of 0 to 1 (for better plotting)
for r=1:index%length(index)
    if corrResults_HC_up(r) == 1
        corrResults_HC_up(r)=0;
    else
    end
end

% plots plots plots
[cmap2] = cbrewer2('seq','BuPu',50); %https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
figure; imagesc(corrResults_HC_up); colorbar; colormap((cmap2));
axis('square','on');
title('HC Significant Correlations between Oxytocin and Functional Connectivity');

% FOR FND
% Loop through each connection of each subject and calculate correlation
corrResults_FND = zeros(size(CMs2,1), size(CMs2,1));
corrPVals_FND = zeros(size(CMs2,1), size(CMs2,1));

x = T.OxytocinPg_ml((N_HC+1:end),:); % for computational speed we keep x out of the loop as it doesn't change, y will always be updated

for n_row = 1:size(CMs2,1)
    for n_col = 1:size(CMs2,1)
        %Extract 
        y = squeeze(CMs2(n_row, n_col,:));
        
        % Calculate the correlation coefficient and its p-value
        [correlationCoeff, pValue] = corrcoef(x, y);
        
        % Calculate the correlation between the scalar value and each cell of the matrix
        corrResults_FND(n_row, n_col) = correlationCoeff(1,2);
        corrPVals_FND(n_row, n_col) = pValue(1,2);
    end
end
% Replace NaN values with 0
corrPVals_FND(isnan(corrPVals_FND)) = 0;
corrResults_FND(isnan(corrResults_FND)) = 0;

% Update cells in results where corresponding cells in pvals are smaller than the threshold
corrResults_FND(corrPVals_FND > threshold) = newValue; % Replace 'newValue' with the value you want to assign

% Let's visualize it
% Remove lower triangular part of correlation matrix (for better visability) and then plot
corrResults_FND_up = triu(corrResults_FND,1);

% Change all values of 0 to 1 (for better plotting)
for r=1:index%length(index)
    if corrResults_FND_up(r) == 1
        corrResults_FND_up(r)=0;
    else
    end
end

% plots plots plots
[cmap2] = cbrewer2('seq','BuPu',50); %https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
figure; imagesc(corrResults_FND_up); colorbar; colormap((cmap2));
axis('square','on');
title('FND Significant Correlations between Oxytocin and Functional Connectivity');


% Load codebook to know names regions that you'll be looking at
% This works for AAL only, can be adapted for other atlases, just change
% codebook, also adapt below the lobes etc. 
load(fullfile(atlasPath,atlas,'annotatedStructTemplate116codeBook.mat'));
myCB=codeBookStructTemplate.full;

% Find non-zero elements
[row, col] = find(corrResults_FND_up ~= 0);

for r=1:length(row)
        disp([myCB.name{row(r)} ' - ' myCB.name{col(r)} ' with Correlation Coefficient =' num2str(corrResults_FND_up(row(r),col(r))) ' with p =' num2str(corrPVals_FND(row(r),col(r))) ]); % display
end

%Now we only look at amygdala
corrResults_FND([1:40, 43:90],[1:40, 43:90]) = 0;
corrResults_FND_up([1:40, 43:90],[1:40, 43:90]) = 0;

% Find non-zero elements
[row, col] = find(corrResults_FND_up ~= 0);

for r=1:length(row)
        disp([myCB.name{row(r)} ' - ' myCB.name{col(r)} ' with Correlation Coefficient =' num2str(corrResults_FND_up(row(r),col(r))) ' with p =' num2str(corrPVals_FND(row(r),col(r))) ]); % display
end


%Same for HC
corrResults_HC([1:40, 43:90],[1:40, 43:90]) = 0;
corrResults_HC_up([1:40, 43:90],[1:40, 43:90]) = 0;

% Find non-zero elements
[row, col] = find(corrResults_HC_up ~= 0);

for r=1:length(row)
        disp([myCB.name{row(r)} ' - ' myCB.name{col(r)} ' with Correlation Coefficient =' num2str(corrResults_HC_up(row(r),col(r))) ' with p =' num2str(corrPVals_HC(row(r),col(r))) ]); % display
end


%% Make edge and node files for BrainNetViewer
FinalPath = '/Volumes/T7_Samantha/Analyses/Oxytocin/02_ConnectivityAnalyses/02_Plotting';
if ~exist(FinalPath),mkdir(FinalPath); end

cd(FinalPath);

% Write edge file
corrResults_FND_new = 1 - corrResults_FND; %so that in brainnet viewer, you can only show singificant regions
for r = 1:nRois
    for c = 1:nRois
        if corrResults_FND_new(r,c) ==1
            corrResults_FND_new(r,c) = 0;
        end
    end
end
dlmwrite('FND_Correlation_OT.edge',corrResults_FND_new,'\t');

corrResults_HC_new = 1 - corrResults_HC; %so that in brainnet viewer, you can only show singificant regions
for r = 1:nRois
    for c = 1:nRois
        if corrResults_HC_new(r,c) ==1
            corrResults_HC_new(r,c) = 0;
        end
    end
end
dlmwrite('HC_Correlation_OT.edge',corrResults_HC_new,'\t');



%Example creating Node (from BrainNetViewer)
% tmp = char(L);
% node = strcat(num2str(R),tmp);
% node(:,end-4:end-3) = char(' ');
% dlmwrite('test.node',node,'delimiter','');

% Making Nodes
node = zeros(90,6);
%Defining Coordingates
for r = 1:nRois
    node(r,1) = myCB.center{1,r}(1,1);
    node(r,2) = myCB.center{1,r}(2,1);
    node(r,3) = myCB.center{1,r}(3,1);
end
%Color brain according to lobe
for l = 1:16
    node(l,4) = 1; %frontal
end
for l=17:18
    node(l,4) = 2; %central
end
for l = 19: 28
    node(l,4) = 1; %frontal
end
for l = 29:42
    node(l,4) = 3; % limbic
end
for l =43:56
    node(l,4) = 4; % occipital
end
for l=57:70
    node(l,4) = 5; %parietal
end
for l = 71:78
    node(l,4) = 6; % subcortical
end
for l = 79:90
    node(l,4) = 7; %Temporal
end

%% Add nodal degree (= number of edges connected to the node)
freq = zeros(90,1);

for k=1:nRois
  freq(k)=sum(corrResults_FND_new(:,k));
end

node(:,5) = freq(:);
cd(FinalPath);
% Add labels and convert to node file
tmp = char(myCB.sname);
node = strcat(num2str(node),tmp);
node(:,end-13:end-12) = char(' '); %longest character from the end is 13 long
dlmwrite('FND_Correlation_OT.node',node,'delimiter','');

%HC
clear freq
freq = zeros(90,1);

for k=1:nRois
  freq(k)=sum(corrResults_HC_new(:,k));
end

node(:,5) = freq(:);
cd(FinalPath);
% Add labels and convert to node file
tmp = char(myCB.sname);
node = strcat(num2str(node),tmp);
node(:,end-13:end-12) = char(' '); %longest character from the end is 13 long
dlmwrite('HC_Correlation_OT.node',node,'delimiter','');


%% Show only highest 20 p-values
    myCorrs(:,1)=jUpperTriMatToVec(corrPVals_FND(:,:,1),1);
    [sortedImp,sortedNums]=sort(myCorrs(:,1),'ascend'); % sort mean importances in ascending order
    CorrMatrix=corrPVals_FND;
    numHighervalToShow=20; % decide how many connections with higher importance you want to explore

    % Make an importance mask only showing 20 most significant regions
    for r=numHighervalToShow+1 : length(sortedImp)
        mean_nums=find(corrPVals_FND==sortedImp(r)); % get number corresponding to highest importance in the mean importance matrix
        [I,J]=ind2sub(nRois,mean_nums(1)); % get row and column where the highest importance is
        corrPVals_FND(I,J) = 0;
        corrPVals_FND(J,I) = 0;
    end  

%corrPVals_FND = 1 - corrPVals_FND;
dlmwrite('FND_Correlation_OT_thresholded.edge',corrPVals_FND,'\t');


%% Alternatively we can remove some regions we are not interested in. 

% Let's say we are only interested in 
% amygdala = 41,42
% insula = 29,30
% ACC = 31,32
% PCC = 35, 36
% Hippocampus = 37, 38
% Orbitofrontal cortex = 5,6,9,10, 15,16, 25,26

% Indices of the rows and columns to keep
selectedRows = [41,42,29,30,31,32,35,36,37,38,5,6,9,10, 15,16, 25,26];
selectedCols = [41,42,29,30,31,32,35,36,37,38,5,6,9,10, 15,16, 25,26];

% Extract the selected rows and columns
newAdjMatrix_HC = CMs1(selectedRows, selectedCols,:);
newAdjMatrix_FND = CMs2(selectedRows, selectedCols,:);

%% Calculate correlation between whole-brain FC and OT

% Initialize a matrix to store correlation results
corrResults_HC = zeros(size(newAdjMatrix_HC,1), size(newAdjMatrix_HC,1));
corrPVals_HC = zeros(size(newAdjMatrix_HC,1), size(newAdjMatrix_HC,1));

names = {'Amygdala_L','Amygdala_R','Insula_L','Insula_R', 'ACC_L','ACC_R','PCC_L','PCC_R', 'Hippocampus_L','Hippocampus_R','OFC_L','OFC_R','OFC_L','OFC_R','OFC_L','OFC_R','OFC_L','OFC_R'};
% FOR HC
% Loop through each connection of each subject and calculate correlation
x = T.OxytocinPg_ml((1:N_HC),:); % for computational speed we keep x out of the loop as it doesn't change, y will always be updated

for n_row = 1:size(newAdjMatrix_HC,1)
    for n_col = 1:size(newAdjMatrix_HC,1)
        %Extract 
        y = squeeze(newAdjMatrix_HC(n_row, n_col,:));
        
        % Calculate the correlation coefficient and its p-value
        [correlationCoeff, pValue] = corrcoef(x, y);
        
        % the output has the following structure: 
        % [ correlation(x, x)  correlation(x, y)
        % correlation(y, x)  correlation(y, y) ]
        % This is why we need to access either [1,2] or [2,1]
        
        % Calculate the correlation between the scalar value and each cell of the matrix
        corrResults_HC(n_row, n_col) = correlationCoeff(1,2);
        corrPVals_HC(n_row, n_col) = pValue(1,2);
    end
end
% Replace NaN values with 0
corrPVals_HC(isnan(corrPVals_HC)) = 0;
corrResults_HC(isnan(corrResults_HC)) = 0;

% We make now a new matrix, that only shows those correlation coefficients
% which are statistically significant. We replace all > 0.05 with 1

% Threshold value
threshold = 0.05;
newValue = 1;

% Update cells in results where corresponding cells in pvals are smaller than the threshold
corrResults_HC(corrPVals_HC > threshold) = newValue; % Replace 'newValue' with the value you want to assign

% Let's visualize it
% Remove lower triangular part of correlation matrix (for better visability) and then plot
corrResults_HC_up = triu(corrResults_HC,1);

index = 18 * 18;
% Change all values of 0 to 1 (for better plotting)
for r=1:index%length(index)
    if corrResults_HC_up(r) == 1
        corrResults_HC_up(r)=0;
    else
    end
end

% plots plots plots
[cmap2] = cbrewer2('seq','BuPu',50); %https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
figure; imagesc(corrResults_HC_up); colorbar; colormap((cmap2));
axis('square','on');
set(gca,'xtick',[1:18],'xticklabel',names,'FontSize',12);
set(gca,'ytick',[1:18],'yticklabel',names,'FontSize',12);
xtickangle(90);
title('HC Significant Correlations between Oxytocin and Functional Connectivity');

% FOR FND
% Loop through each connection of each subject and calculate correlation
corrResults_FND = zeros(size(newAdjMatrix_FND,1), size(newAdjMatrix_FND,1));
corrPVals_FND = zeros(size(newAdjMatrix_FND,1), size(newAdjMatrix_FND,1));

x = T.OxytocinPg_ml((N_HC+1:end),:); % for computational speed we keep x out of the loop as it doesn't change, y will always be updated

for n_row = 1:size(newAdjMatrix_FND,1)
    for n_col = 1:size(newAdjMatrix_FND,1)
        %Extract 
        y = squeeze(newAdjMatrix_FND(n_row, n_col,:));
        
        % Calculate the correlation coefficient and its p-value
        [correlationCoeff, pValue] = corrcoef(x, y);
        
        % Calculate the correlation between the scalar value and each cell of the matrix
        corrResults_FND(n_row, n_col) = correlationCoeff(1,2);
        corrPVals_FND(n_row, n_col) = pValue(1,2);
    end
end
% Replace NaN values with 0
corrPVals_FND(isnan(corrPVals_FND)) = 0;
corrResults_FND(isnan(corrResults_FND)) = 0;

% Update cells in results where corresponding cells in pvals are smaller than the threshold
corrResults_FND(corrPVals_FND > threshold) = newValue; % Replace 'newValue' with the value you want to assign

% Let's visualize it
% Remove lower triangular part of correlation matrix (for better visability) and then plot
corrResults_FND_up = triu(corrResults_FND,1);

% Change all values of 0 to 1 (for better plotting)
for r=1:index%length(index)
    if corrResults_FND_up(r) == 1
        corrResults_FND_up(r)=0;
    else
    end
end

% plots plots plots
[cmap2] = cbrewer2('seq','BuPu',50); %https://bacteriophysics.web.illinois.edu/wp/wp-content/uploads/2019/02/figure-guide-ls-2016.pdf
figure; imagesc(corrResults_FND_up); colorbar; colormap((cmap2));
axis('square','on');
set(gca,'xtick',[1:18],'xticklabel',names,'FontSize',12);
set(gca,'ytick',[1:18],'yticklabel',names,'FontSize',12);
xtickangle(90);
title('FND Significant Correlations between Oxytocin and Functional Connectivity');

%% Visualization

% edges and nodes files can be used in BrainNetViewer for visualization. 
