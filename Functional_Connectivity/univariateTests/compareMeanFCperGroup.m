%clear

featPath = 'C:\Samantha\RestingState_Study\FNDstudy_Rest\04_Feature_Vectors';  % where your AAL_3DCMs.mat are saved during vCMextract_CHpreproc
% %%%featPath changed!%%%
%% set parameters
which_group='CD'; % OR_CD or patients
filtering ='bandpass';
switch which_group
    case 'OR_CD'
        N_OR=35; N_CD=23;
        cLabs=[zeros(1,N_OR) ones(1,N_CD)];
    case 'patients'
        cLabs=[0 0 1 0 1 0 1 0 0 0 1 1 0 1 0 0 1 1 1 0 0 1]; % 0=MOVEMENT 1=PARESIS (class labels)
end

atlas='AAL';
subband=4;
%features='wholeBrain'; % wholeBrain or networks
roi1 = 72; %% choose number corresponding to first region
roi2 = 41; %% choose number corresponding to second region

%% Load codebook to know number of regions that you'll be looking at
atlasPath='C:\Samantha\RestingState_Study\01_Scripts_Samantha\atlases';

switch atlas
    case 'AAL'
        load(fullfile(atlasPath,atlas,'annotatedStructTemplate116codeBook.mat'));
        myCB=codeBookStructTemplate.full;
    case 'Shirer'
        load(fullfile(atlasPath,atlas,'Shirer_codebook.mat'));
    case 'Hammers'
        load(fullfile(atlasPath,atlas,'Hammers_codebook.mat'));
end

load(fullfile(featPath,[atlas '_' filtering '_3DCMs_fisherZ.mat']));

%switch features
%    case 'wholeBrain'
        features_1=my3DCMs_fisherZ(:,:,cLabs==0);
        features_2=my3DCMs_fisherZ(:,:,cLabs==1);
%    case 'networks'
%        myRegions=double(45:49); %%%%%%% regions within network to list
%        features_1=my3DCMs_fisherZ(myRegions,myRegions,cLabs==0);
%        features_2=my3DCMs_fisherZ(myRegions,myRegions,cLabs==1);
%end

features_1(isnan(features_1))=0; features_2(isnan(features_2))=0;

% Get average connectivity for each group
mean_features_1=mean(features_1,3);
mean_features_2=mean(features_2,3);

mean_roi_connectivity_1=mean_features_1(roi1,roi2);
mean_roi_connectivity_2=mean_features_2(roi1,roi2);

disp(['mean connectivity between ' myCB.name{roi1} ' and ' myCB.name{roi2} ' ' num2str(mean_roi_connectivity_1,'%1.2f') ' in group 1 and ' num2str(mean_roi_connectivity_2,'%1.2f') ' in group 2']);
