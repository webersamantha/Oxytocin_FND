%% Compare mean connectivity matrices between groups using t-tests or nonparametric Mann-U Whitney tests %%

%clear
RootPath = 'C:\Samantha\RestingState_Study\FNDstudy_Rest';
scriptsPath = 'C:\Samantha\RestingState_Study\01_Scripts_Samantha';
atlasPath= fullfile(scriptsPath,'atlases');
featPath = fullfile(RootPath,'04_Feature_Vectors');
which_group='HR_CD'; % OR_CD or patients (patients vs controls or patients subgroups
atlas='AAL';
    prompt = {'On which Subband will you run the classification?'};
    title = 'Input Subband';
    dims = [1 35];
    definput = {'4'};
    answer = inputdlg(prompt,title,dims,definput);
    this_subband = str2num(answer{1});
    subband = num2str(this_subband);
N_OR=35;
filtering ='bandpass';
features='wholeBrain'; % wholeBrain or networks

load('C:\Samantha\RestingState_Study\FNDstudy_Rest\04_Feature_Vectors\AAL_bandpass_vectors_fisherZ.mat');
%load(fullfile(featPath,which_group,[atlas '_3DCM4s_fisherZ.mat']));

CMs1=my3DCMs_fisherZ(:,:,1:N_OR);
CMs2=my3DCMs_fisherZ(:,:,N_OR+1:end);
myRegions=[double(33:38)]; %%%%%%% regions within network to list
CMs1=my3DCMs_fisherZ(myRegions,myRegions,1:N_OR);
CMs2=my3DCMs_fisherZ(myRegions,myRegions,N_OR+1:end);

CMs1(isnan(CMs1))=0; CMs2(isnan(CMs2))=0;
qVal=0.05; % alpha

% Choose univariate test & correction for multiple comparisons
test='parametric'; % 'parametric or nonparametric
FDRcorr='yes'; % yes or no

%% Testing

switch test
    case 'parametric'
        switch FDRcorr
            case 'yes'
                [pVals,sig,ci,stats]=jCorrmatSignificanceMulti2(CMs1,CMs2,qVal);
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

% 2. Visualise connections that are significantly different between your
% groups
newPvals=pVals;
newPvals(isnan(newPvals))=0;
index=find(newPvals~=0);

% Change all p-values higher than qVal to 1
for r=1:length(index)
    if newPvals(r)>qVal
        newPvals(r)=1;
    end
end

figure; imagesc(newPvals); colorbar; % what's left in blue are your significantly different connections (with p-values lower than 0.05)


%% explore group means for significantly different correlations

meanCM_HC=mean(CMs1,3);
meanCM_PT=mean(CMs2,3);

% load results file containing allImportances variable (only with RF!)
%load('/Users/jweg/Desktop/Jen/All_rest/Resultats_All_rest_withoutold/classifResults_Shirer_wholeBrain_CM4_RF.mat')
load('D:\RestingState_Study\FNDstudy_Rest\05_Results_Classification\classifResults_perms_AAL_wholeBrain_bandpass_SVM_fisherZ.mat')
% clear workspace for unnecessary variables
clear subband X X_TR X_TE votes predLabels prediction_per_tree nFeats nFolds N_CD N_OR myAcc TEsIdx TRsIdx feats featureMatrix mvtParams myAccs myAUC myClassAccs myCls MYCLASSIFIER MYCPARAMS myStats allPredictions allRFfuncVals allTruths allVotes cLabs cLabs_TR cLabs_TE D DISPLAYPLOTS f 

% load codebook and get nRois
switch atlas
    case 'AAL'
        nRois=90;
        load(fullfile(atlasPath,atlas,'annotatedStructTemplate116codeBook.mat'));
    case 'Shirer'
        nRois=90;
        load(fullfile(atlasPath,atlas,'Shirer_codebook.mat'));
    case 'Hammers'
        nRois=83;
        load(fullfile(atlasPath,atlas,'Hammers_codebook.mat'));
end

meanImportances=mean(allImportances,2); % get mean importances across subjects
[sortedImp,sortedNums]=sort(meanImportances,'descend'); % sort mean importances in decreasing order
meanImpCM=jVecToSymmetricMat(meanImportances,nRois); % transform mean importance vector into mean importance matrix (nRegions x nRegions)

numHigherImpToShow=10; % decide how many connections with higher importance you want to explore

for r=1:numHigherImpToShow
    CM_nums=find(meanImpCM==sortedImp(r)); % get number corresponding to highest importance in the mean importance matrix
    [I,J]=ind2sub(nRois,CM_nums(1)); % get row and column where the highest importance is
    disp([myCB.name{I} ' - ' myCB.name{J} ' in controls mean=' num2str(meanCM_HC(I,J)) ' in patients mean=' num2str(meanCM_PT(I,J))]); % display
end
