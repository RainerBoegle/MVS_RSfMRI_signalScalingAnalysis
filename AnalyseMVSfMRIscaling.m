function MVSscaling = AnalyseMVSfMRIscaling()
%This function allows to analyse the distribution of voxel-wise scaling factors,
%i.e. the fraction of amplitudes, derived from dualreg analysis of groupICA results,
%between MRI2 = 3T MRI and MRI1 = 1.5T MRI.
%
%NB: The fraction that we look for should be MRI2/MRI1 = 2*sqrt(2)! 
%   (3T/1.5T = 2 & scaling of fMRI signal is sqrt(3T)/sqrt(1.5T) = sqrt(2))
%
% The script will ask the user to select the SPM.mat-file of a second-level
% ANOVA (for 2groups2MRIs&SubjectFactors) based on dualreg-outputs and will
% identify all important files and relationships from the SPM.mat.
%
%IT IS ASSUMED that the mask-files (as 0/1 values) for analysis are in the same folder as the SPM.mat-file!
%These include:
%       mask.img/hdr  <-- whole brain mask created by SPM
%       MaskA_IC.nii  <-- User-made: significant voxels of IC
%       MaskB_IC.nii  <-- User-made: significant voxels of IC(dualreg ANOVA AllEffectsContrast)
%       Mask_EffectOfMRI.nii  <-- User-made: significant voxels of dualreg ANOVA "Effect of MRI" or "MRI2>MRI1", depending on taste.
%
%ALL OTHER MASKs will be created by the script/functions!
%Analysis will include:
%       a) whole brain
%       b) significant voxels of IC (user selects which definition (MELODIC/dualregANOVA) is used)
%       c) significant voxels of ANOVA "Effects of MRI" (or "MRI2>MRI1")
%       d) whole brain WITHOUT significant voxels "Effects of MRI"
%       e) significant voxels of IC WITHOUT significant voxels "Effects of MRI"
%       f) whole brain WITHOUT the COMBINATION of significant voxels of IC AND significant voxels "Effects of MRI"
%
%WE EXPECT FOR:
%       a) some intersting mix of results per subject and per group
%       b) lots of voxels with scaling of 2*sqrt(2) but also others (if MRI effect coincides with IC map)
%       c) hopefully mostly voxels with scaling of 2*sqrt(2) (Hopefully very clearly different from the other distribution, -hopefully.)
%       d) whatever, but basically nothing at scaling of 2*sqrt(2), except some not
%          significant voxels that we missed somehow ie not in mask of "Effects of MRI".
%       e) hopefully not scaling of 2*sqrt(2), showing that the "Effects of MRI" has reasonable scaling!
%       f) a random rest without traces of scaling of 2, -hopefully.
%
%ROI-Analysis:
% Additionally the analysis will include a split of the "Effects of MRI" (or "MRI2>MRI1") 
% AND significant voxels of IC WITHOUT significant voxels "Effects of MRI" into AAL-Atlas regions 
% & also a ROI-ANalysis is possible given a ROI-NIFTI indicating the ROIs as integers 1:NROIs in the map.
%
%
%WHAT WE ANALYSE AND DISPLAY:
% Calculations in struct:
%              MVSscaling.
%                        .ScalingPerSubject(1:NVoxels_wholebrain,1:NSubj)  <-- scaling as is, ie can be negative. a
%                        NB: here we calculate the median scaling=MRI2/MRI1 ie for all combinations of runs per subject.
%
%
%Date/Version: 07.01.2015/4.1   (06.08.2014/1.0)
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)

% To Do:
% DONE: 1. check any writing of NIFTI should have a dt field of at least [16 0]!!! --> retry SPM Analysis
% 2. allow to import ROIs with description similar to AAL but selfmade.
%    --> this needs an extension of the output struct "MVSscaling" in addition to the field "MVSscaling.Masks.AALsplit."
%        we should also have a field "MVSscaling.Masks.ROIsplit." which uses the user input ROI mask from 1:max(ROI_Mask(:))
%        and the labels. Just like MVSscaling.Masks.AALsplit.EffectOfMRI = load([FunctionDir,filesep,'AAL_Atlas',filesep,'AAL_Labels.mat']);
% I.  make AALsplit generation into a function such that it can be used for each mask.
% II. add a function for importing ROI_Mask & Labels.
% III.change the plotting functions.

disp(' ');
disp('Starting analysis of MVSscaling of fMRI data...');

%% get inputs & main masks
MVSscaling = GetInputsMasks();

%% make AAL-Atlas splits
MVSscaling = MakeAALsplits(MVSscaling);

%% make ROI splits
MVSscaling = MakeROIsplits(MVSscaling);

%% For each Subject make all possible scaling calculations for all run combinations and calculate the median
MVSscaling.AmplitudesPerSubjectMRI = zeros(size(MVSscaling.Masks.AllMasks,1),MVSscaling.Design.NSubjs,MVSscaling.Design.N_MRIs);
MVSscaling.ScalingPerSubject       = zeros(size(MVSscaling.Masks.AllMasks,1),MVSscaling.Design.NSubjs);
MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum = zeros(MVSscaling.Design.NSubjs,1);
disp(' ');
disp('Creating fMRI scaling data per subject ...');
for IndSubj = 1:MVSscaling.Design.NSubjs
    disp(['Subj ',num2str(IndSubj)]);
    %average over runs per MRI
    for IndMRI = 1:MVSscaling.Design.N_MRIs
        MRIindex = find(MVSscaling.Design.SubjNrPerInput==IndSubj & MVSscaling.Design.MRI_NrPerInput==IndMRI);
        MVSscaling.AmplitudesPerSubjectMRI(:,IndSubj,IndMRI) = mean(MVSscaling.InputData.DataForAnalysis2D(:,MRIindex),2);
        clear MRIindex
    end
    %create scaling
    MVSscaling.ScalingPerSubject(:,IndSubj) = MVSscaling.AmplitudesPerSubjectMRI(:,IndSubj,2)./MVSscaling.AmplitudesPerSubjectMRI(:,IndSubj,1);
    %subject number according to ICA inputs & in accordance with eye movement analysis
    ANOVAInputsForCurrSubj = find(MVSscaling.Design.SubjNrPerInput==IndSubj);
    MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum(IndSubj) = MVSscaling.InputData.ICA.SubjNumPerANOVAInput(ANOVAInputsForCurrSubj(1));
    clear ANOVAInputsForCurrSubj
end

%% extra output of scaling data per subject
ScalingPerSubject.Design.SubjNumPerInput = MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum;
ScalingPerSubject.Design.GroupNrPerInput = MVSscaling.Design.GroupNrPerSubj;
ScalingPerSubject.Data2DMasked    = MVSscaling.ScalingPerSubject;
ScalingPerSubject.AbsData2DMasked = abs(MVSscaling.ScalingPerSubject);
ScalingPerSubject.MaskPath        = MVSscaling.Masks.MPaths{1};

%% save
disp(' ');
disp(['Saving results to "',MVSscaling.SPMdir,filesep,'ScalingAnalysis','".'])
if(~isdir([MVSscaling.SPMdir,filesep,'ScalingAnalysis']))
    mkdir([MVSscaling.SPMdir,filesep,'ScalingAnalysis']);
end
save([MVSscaling.SPMdir,filesep,'ScalingAnalysis',filesep,'MVSscaling_',MVSscaling.Masks.AllMasksName{2},'.mat'],'MVSscaling');
save([MVSscaling.SPMdir,filesep,'ScalingAnalysis',filesep,'MVS_ScalingPerSubject.mat'],'ScalingPerSubject');

%% Done.
disp(' ');
disp('DONE. (Check results using "DisplayMVSfMRIscalingResults()".)');

end

%% subfunctions
%% GetInputsMasks
function MVSscaling = GetInputsMasks()
% This function selects all input data and masks

%% get SPM.mat to extract information for analysis
MVSscaling.SPMmatPath = spm_select(1,'mat','Select SPM.mat for extraction of Information for analysis...',[],pwd,'^SPM',1);
load(MVSscaling.SPMmatPath);
MVSscaling.SPMdir = fileparts(MVSscaling.SPMmatPath);

%% get .filelist to extract subject information
MVSscaling.ICAfilelistPath = spm_select(1,'any','Select ".filelist" file specifying inputs to initial ICA...',[],pwd,'^.filelist',1);
MVSscaling.InputData.ICA.filelist = importdata(MVSscaling.ICAfilelistPath);

%% extract information from SPM struct & .filelist
%%Notes
% SPM.xX.X;  %DesignMatrix 
% SPM.xX.iH; %Effects columns
% SPM.xX.iC; %Covariates columns 
% SPM.xX.iB; %Subject columns 
% % SPM.xX.iG, %? unknown
% SPM.xX.name; %Name of columns/effects 
% SPM.xX.I; %grouping: 1stColumn: "Run" 2ndColumn: "Group" 3rdColumn: "MRI" 4thColumn: "SubjNr 1to27"
% SPM.xY.P; %file paths

MVSscaling.InputData.FilePaths = SPM.xY.P;
FileIndex    = zeros(length(MVSscaling.InputData.FilePaths),1);
FileIndexStr = cell( length(MVSscaling.InputData.FilePaths),1);
for Ind = 1:length(MVSscaling.InputData.FilePaths)
    [Base,File,Ext] = fileparts(MVSscaling.InputData.FilePaths{Ind});
    
    FileIndex(Ind)    = str2num(Ext(6:end));
    FileIndexStr{Ind} = Ext(6:end); 
end
StrPattern = 'Subj_p';
StartInds = strfind(MVSscaling.InputData.ICA.filelist,StrPattern);
MVSscaling.InputData.ICA.SubjNumPerICAInput = zeros(length(MVSscaling.InputData.ICA.filelist),1);
for IndInput = 1:length(MVSscaling.InputData.ICA.SubjNumPerICAInput)
    MVSscaling.InputData.ICA.SubjNumPerICAInput(IndInput) = str2num(MVSscaling.InputData.ICA.filelist{IndInput}(StartInds{IndInput}+length(StrPattern):StartInds{IndInput}+length(StrPattern)+2));
end
MVSscaling.InputData.ICA.SubjNumPerANOVAInput = MVSscaling.InputData.ICA.SubjNumPerICAInput(FileIndex); %Transform to design of ANOVA

[MVSscaling.InputData.BasePathInputs,MVSscaling.InputData.FNameICdualregInputs] = fileparts(SPM.xY.P{1});
MVSscaling.InputData.InputsNII_Path = [MVSscaling.InputData.BasePathInputs,filesep,MVSscaling.InputData.FNameICdualregInputs,'.nii']; %the NIFTI with all weigths per Subject & Run for the IC used in dualreg & subsequent ANOVA
if(~exist(MVSscaling.InputData.InputsNII_Path))
    MVSscaling.InputData.InputsNII_Path=spm_select(1,'image','Select DualReg 2ndStage 4D NIFTI manually...');
end

% get dualreg weights
MVSscaling.InputData.NII       = nifti(MVSscaling.InputData.InputsNII_Path);
MVSscaling.InputData.AllData4D = MVSscaling.InputData.NII.dat(:,:,:,:);
MVSscaling.InputData.AllData2D = reshape(MVSscaling.InputData.AllData4D,[],size(MVSscaling.InputData.AllData4D,4));
MVSscaling.InputData.mat       = MVSscaling.InputData.NII.mat;


%grouping vector for each level of each factor in the interaction Group*MRI
MVSscaling.Design.NSubjs              = max(SPM.xX.I(:,4));
MVSscaling.Design.SubjNrPerInput(:,1) = SPM.xX.I(:,4);
MVSscaling.Design.NGroups             = max(SPM.xX.I(:,2));
MVSscaling.Design.GroupNrPerInput(:,1)= SPM.xX.I(:,2);
if(max(SPM.xX.I(:,3))==2)
    MVSscaling.Design.N_MRIs          = max(SPM.xX.I(:,3));
else
    error('codes is not set up for more then two MRIs! Change for this application.');
end
MVSscaling.Design.MRI_NrPerInput(:,1) = SPM.xX.I(:,3);

MVSscaling.Design.GroupNrPerSubj = zeros(MVSscaling.Design.NSubjs,1);
for IndSubj = 1:MVSscaling.Design.NSubjs
    IndexSubj = find(MVSscaling.Design.SubjNrPerInput==IndSubj);
    MVSscaling.Design.GroupNrPerSubj(IndSubj) = MVSscaling.Design.GroupNrPerInput(IndexSubj(1));
    clear IndexSubj
end
if(any(MVSscaling.Design.GroupNrPerSubj==0))
    error('Some Subjects are still not assigned to a group.');
end

%% Check for the masks
%whole brain from SPM
if(exist([MVSscaling.SPMdir,filesep,'mask.img']))
    if(strcmp('mask.img',questdlg({'Do you want to use the Mask created by SPM ("SPMdir/mask.img") or select your own mask?'; ' '; 'NB: your mask has to fit the data dimensions and orientation, otherwise this will lead to an error!'},'Which Mask?','mask.img','MyMask','MyMask')))
        MVSscaling.Masks.MPaths{1} = [MVSscaling.SPMdir,filesep,'mask.img'];
    else
        MVSscaling.Masks.MPaths{1} = spm_select(1,'image','Select a mask image that is appropriate for the input scaling data (NO 4D-files only 3D!!!)...');
        [tmpbase,tmpfname,ext] = fileparts(MVSscaling.Masks.MPaths{1});
        if(~isempty(strfind(ext,',')))
            ext = ext(1:(strfind(ext,',')-1));
            MVSscaling.Masks.MPaths{1}=[tmpbase,filesep,tmpfname,ext]; %NB: this excludes 4D-files!!!
        end
    end
else
    error(['Could not find WHOLE BRAIN mask made by SPM "mask.img" (',[MVSscaling.SPMdir,filesep,'mask.hdr'],'"']);
end
NII_tmp = nifti(MVSscaling.Masks.MPaths{1});
if(sum((MVSscaling.InputData.mat-NII_tmp.mat).^2)~=0)
    error('Whole Brain NIFTI-MATRIX & DualRegInputData NIFTI-MATIX are NOT EQUAL! Orientation wrong?');
end
MVSscaling.Masks.WholeBrainRaw = NII_tmp.dat(:);
clear NII_tmp

%MaskA_IC.nii  <-- User-made: significant voxels of IC
if(exist([MVSscaling.SPMdir,filesep,'MaskA_IC.nii']))
    MVSscaling.Masks.MPaths{2} = [MVSscaling.SPMdir,filesep,'MaskA_IC.nii'];
else
    if(exist([MVSscaling.SPMdir,filesep,'MaskA_IC.img']))
        MVSscaling.Masks.MPaths{2} = [MVSscaling.SPMdir,filesep,'MaskA_IC.img'];
    else
        error(['Could not find MASK("significant voxels of IC") made by USER (',[MVSscaling.SPMdir,filesep,'MaskA_IC.nii/img'],'"']);
    end
end

%MaskB_IC.nii  <-- User-made: significant voxels of IC(dualreg ANOVA AllEffectsContrast)
if(exist([MVSscaling.SPMdir,filesep,'MaskB_IC.nii']))
    MVSscaling.Masks.MPaths{3} = [MVSscaling.SPMdir,filesep,'MaskB_IC.nii'];
else
    if(exist([MVSscaling.SPMdir,filesep,'MaskB_IC.img']))
        MVSscaling.Masks.MPaths{3} = [MVSscaling.SPMdir,filesep,'MaskB_IC.img'];
    else
        error(['Could not find MASK("significant voxels of IC(dualreg ANOVA AllEffectsContrast)") made by USER (',[MVSscaling.SPMdir,filesep,'MaskB_IC.nii/img'],'"']);
    end
end

%Mask_EffectOfMRI.nii  <-- User-made: significant voxels of dualreg ANOVA "Effect of MRI" or "MRI2>MRI1", depending on taste.
if(exist([MVSscaling.SPMdir,filesep,'Mask_EffectOfMRI.nii']))
    MVSscaling.Masks.MPaths{4} = [MVSscaling.SPMdir,filesep,'Mask_EffectOfMRI.nii'];
else
    if(exist([MVSscaling.SPMdir,filesep,'Mask_EffectOfMRI.img']))
        MVSscaling.Masks.MPaths{4} = [MVSscaling.SPMdir,filesep,'Mask_EffectOfMRI.img'];
    else
        error(['Could not find MASK("significant voxels of dualreg ANOVA "Effect of MRI"") made by USER (',[MVSscaling.SPMdir,filesep,'Mask_EffectOfMRI.nii/img'],'"']);
    end
end

%% FUTURE EXTENSION: get thresholded brain maps to make a plot of Scaling(lambda) over significance values!
%                  Are the most significant also those that scale with lambda==2?
%                  Is this only true for the Effect of MRI? If it is true at all?

%ThresMapA_IC.nii <-- User-made: zvals or tvals of significant voxels of IC
if(exist([MVSscaling.SPMdir,filesep,'ThresMapA_IC.nii']))
    MVSscaling.ThresMaps.MPaths{1} = [MVSscaling.SPMdir,filesep,'ThresMapA_IC.nii'];
else
    if(exist([MVSscaling.SPMdir,filesep,'ThresMapA_IC.img']))
        MVSscaling.ThresMaps.MPaths{1} = [MVSscaling.SPMdir,filesep,'ThresMapA_IC.img'];
    else
        error(['Could not find ThresMap("significant voxels of IC") made by USER (',[MVSscaling.SPMdir,filesep,'ThresMapA_IC.nii/img'],'"']);
    end
end

%ThresMapB_IC.nii <-- User-made: zvals or tvals of significant voxels of IC(dualreg ANOVA AllEffectsContrast)
if(exist([MVSscaling.SPMdir,filesep,'ThresMapB_IC.nii']))
    MVSscaling.ThresMaps.MPaths{2} = [MVSscaling.SPMdir,filesep,'ThresMapB_IC.nii'];
else
    if(exist([MVSscaling.SPMdir,filesep,'ThresMapB_IC.img']))
        MVSscaling.ThresMaps.MPaths{2} = [MVSscaling.SPMdir,filesep,'ThresMapB_IC.img'];
    else
        error(['Could not find ThresMap("significant voxels of IC(dualreg ANOVA AllEffectsContrast)") made by USER (',[MVSscaling.SPMdir,filesep,'ThresMapB_IC.nii/img'],'"']);
    end
end

%ThresMap_EffectOfMRI.nii  <-- User-made: significant voxels of dualreg ANOVA "Effect of MRI" or "MRI2>MRI1", depending on taste.
if(exist([MVSscaling.SPMdir,filesep,'ThresMap_EffectOfMRI.nii']))
    MVSscaling.ThresMaps.MPaths{3} = [MVSscaling.SPMdir,filesep,'ThresMap_EffectOfMRI.nii'];
else
    if(exist([MVSscaling.SPMdir,filesep,'ThresMap_EffectOfMRI.img']))
        MVSscaling.ThresMaps.MPaths{3} = [MVSscaling.SPMdir,filesep,'ThresMap_EffectOfMRI.img'];
    else
        error(['Could not find ThresMap("significant voxels of dualreg ANOVA "Effect of MRI"") made by USER (',[MVSscaling.SPMdir,filesep,'ThresMap_EffectOfMRI.nii/img'],'"']);
    end
end

%% create all necessary masks & get data for analysis & get ThresMaps
%data for analysis (whole brain) match to design in SPM! 
%NB: indices need to be rearranged, because NII goes from 1:108, but the design may be 1 2 5 6 7 8 ... 53 54 3 4 ...
MVSscaling.InputData.DataForAnalysis2D = MVSscaling.InputData.AllData2D(MVSscaling.Masks.WholeBrainRaw==1,FileIndex);

%masks
MVSscaling.Masks.AllMasks = ones(length(find(MVSscaling.Masks.WholeBrainRaw==1)),6); %all masks just like whole brain, then change them
MVSscaling.Masks.AllMasksName = cell(6,1);
%a) whole brain
MVSscaling.Masks.AllMasksName{1} = 'whole brain';
%b) significant voxels of IC (user selects which definition (MELODIC/dualregANOVA) is used)
choice = questdlg({'Which mask should be used to mark voxels that all belong to significant voxels of IC'; 'For those coming from MELODIC analysis including mixture model alternative hypothesis select "MaskA-IC"'; 'For those coming from the DualReg result for ANOVA "AllEffects of Interest" select "MaskB-IC"'},'Select IC mask','MaskA-IC','MaskB-IC','MaskA-IC');
switch(choice)
    case {'MaskA-IC'}
        NII_tmp = nifti(MVSscaling.Masks.MPaths{2});
        Mask_tmp= NII_tmp.dat(:).*MVSscaling.Masks.WholeBrainRaw;
        MVSscaling.Masks.AllMasks(:,2) = Mask_tmp(MVSscaling.Masks.WholeBrainRaw==1);
        MVSscaling.Masks.AllMasks(:,2) = (MVSscaling.Masks.AllMasks(:,2)>0); 
        clear NII_tmp Mask_tmp
        MVSscaling.Masks.AllMasksName{2} = 'IC(MELODIC)';
    case {'MaskB-IC'}
        NII_tmp = nifti(MVSscaling.Masks.MPaths{3});
        Mask_tmp= NII_tmp.dat(:).*MVSscaling.Masks.WholeBrainRaw;
        MVSscaling.Masks.AllMasks(:,2) = Mask_tmp(MVSscaling.Masks.WholeBrainRaw==1);
        MVSscaling.Masks.AllMasks(:,2) = (MVSscaling.Masks.AllMasks(:,2)>0);
        clear NII_tmp Mask_tmp
        MVSscaling.Masks.AllMasksName{2} = 'IC(DualRegANOVA)';
end
%c) significant voxels of ANOVA "Effects of MRI" (or "MRI2>MRI1")
NII_tmp = nifti(MVSscaling.Masks.MPaths{4});
Mask_tmp= NII_tmp.dat(:).*MVSscaling.Masks.WholeBrainRaw;
MVSscaling.Masks.AllMasks(:,3) = Mask_tmp(MVSscaling.Masks.WholeBrainRaw==1);
MVSscaling.Masks.AllMasks(:,3) = (MVSscaling.Masks.AllMasks(:,3)>0); 
clear NII_tmp Mask_tmp
MVSscaling.Masks.AllMasksName{3} = 'EffectOfMRI';
%d) whole brain WITHOUT significant voxels "Effects of MRI"
MVSscaling.Masks.AllMasks(:,4) = (~(MVSscaling.Masks.AllMasks(:,1).*MVSscaling.Masks.AllMasks(:,3)));
MVSscaling.Masks.AllMasksName{4} = 'whole brain WO EffectOfMRI';
%e) significant voxels of IC WITHOUT significant voxels "Effects of MRI"
MVSscaling.Masks.AllMasks(:,5) = (MVSscaling.Masks.AllMasks(:,2).*(~MVSscaling.Masks.AllMasks(:,3)));
MVSscaling.Masks.AllMasksName{5} = 'SignifVoxelsIC WO EffectOfMRI';
%f) whole brain WITHOUT the COMBINATION of significant voxels of IC AND significant voxels "Effects of MRI"
MVSscaling.Masks.AllMasks(:,6) = (~(MVSscaling.Masks.AllMasks(:,1).*(MVSscaling.Masks.AllMasks(:,2)+MVSscaling.Masks.AllMasks(:,3))));
MVSscaling.Masks.AllMasksName{6} = 'whole brain WO (IC&EffectOfMRI)';

%check
disp(' ');
for IndMask = 1:size(MVSscaling.Masks.AllMasks,2)
    disp([MVSscaling.Masks.AllMasksName{IndMask},': ',num2str(length(find(MVSscaling.Masks.AllMasks(:,IndMask)==1))),'Voxels']);
end

%% Thresholded statistics maps
%ThresMaps
MVSscaling.ThresMaps.AllMaps = zeros(length(find(MVSscaling.Masks.WholeBrainRaw==1)),3); %all ThresMaps just like whole brain, then change them
MVSscaling.Masks.AllMapsName = cell(3,1);
%i)  ZVals significant voxels of IC(MELODIC)
NII_tmp = nifti(MVSscaling.ThresMaps.MPaths{1});
Map_tmp = NII_tmp.dat(:).*MVSscaling.Masks.WholeBrainRaw;
MVSscaling.ThresMaps.AllMaps(:,1) = Map_tmp(MVSscaling.Masks.WholeBrainRaw==1);
clear NII_tmp Mask_tmp
MVSscaling.ThresMaps.AllMapsName{1} = 'ZVals-SignifVoxelsIC(MELODIC)';
%ii) ZVals significnt voxels of IC(dualregANOVA-AllEffects)
NII_tmp = nifti(MVSscaling.ThresMaps.MPaths{2});
Map_tmp = NII_tmp.dat(:).*MVSscaling.Masks.WholeBrainRaw;
MVSscaling.ThresMaps.AllMaps(:,2) = Map_tmp(MVSscaling.Masks.WholeBrainRaw==1);
clear NII_tmp Mask_tmp
MVSscaling.ThresMaps.AllMapsName{2} = 'ZVals-SignifVoxelsIC(DualRegANOVA)';
%iii)ZVals significant voxels of dualregANOVA-EffectOfMRI
NII_tmp = nifti(MVSscaling.ThresMaps.MPaths{3});
Map_tmp = NII_tmp.dat(:).*MVSscaling.Masks.WholeBrainRaw;
MVSscaling.ThresMaps.AllMaps(:,3) = Map_tmp(MVSscaling.Masks.WholeBrainRaw==1);
clear NII_tmp Mask_tmp
MVSscaling.ThresMaps.AllMapsName{3} = 'ZVals-SignifVoxels(EffectMRI)';

%check
disp(' ');
for IndThresMap = 1:size(MVSscaling.ThresMaps.AllMaps,2)
    disp([MVSscaling.ThresMaps.AllMapsName{IndThresMap},': ',num2str(length(find(MVSscaling.ThresMaps.AllMaps(:,IndThresMap)>0))),'Voxels']);
end

end

%% MakeAALsplits
function MVSscaling = MakeAALsplits(MVSscaling)
%% FUTURE EXTENSION(AAL): mark parts of Mask_EffectOfMRI according to AAL and make masks for these subparts to be analysed too.
% get AAL labels
[FunctionDir] = fileparts(mfilename('fullpath'));
MVSscaling.Masks.AALsplit.EffectOfMRI = load([FunctionDir,filesep,'AAL_Atlas',filesep,'AAL_Labels.mat']);

% make overlap with AAL & expand
AAL_splits_Dir = [MVSscaling.SPMdir,filesep,'ScalingAnalysis',filesep,'AAL_split',filesep,'EffectsOfMRI'];
[ExpandAtlasOverlap_OutputStruct] = MakeAALOverlap(MVSscaling.Masks.MPaths{4},AAL_splits_Dir,'EffectsOfMRI');

% % get mask files
% AAL_splits_Dir = spm_select(1,'dir','Select directory containing the AAL-splits (& recombinations) of "Effect of MRI"...');

%Expanded-Overlap with AAL
fName = spm_select('List',AAL_splits_Dir,'^Expanded.*.nii');
if(~isempty(fName))
    fName = cellstr(fName);
else
    error('fName is empty! Check AAL-split directory!');
end
if(length(fName)>1)
    [IndexSel,ok] = listdlg('ListString',fName,'Name','Which?','SelectionMode','single');
    if(ok)
        MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{1}    = [AAL_splits_Dir,filesep,fName{IndexSel}];
    else
        MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{1}    = [AAL_splits_Dir,filesep,fName{1}];
    end
else
    MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{1}    = [AAL_splits_Dir,filesep,fName{1}];
end
MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{1}= 'ExpandedOverlapAAL'; %assign names to this variable such that listdlg can be used to select which one to use in Display-function

%CombineBasic-Overlap with AAL
fName = spm_select('List',AAL_splits_Dir,'^CombineBasic.*.nii');
if(~isempty(fName))
    fName = cellstr(fName);
else
    error('fName is empty! Check AAL-split directory!');
end
if(length(fName)>1)
    [IndexSel,ok] = listdlg('ListString',fName,'Name','Which?','SelectionMode','single');
    if(ok)
        MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{2}= [AAL_splits_Dir,filesep,fName{IndexSel}];
    else
        MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{2}= [AAL_splits_Dir,filesep,fName{1}];
    end
else
    MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{2}    = [AAL_splits_Dir,filesep,fName{1}];
end
MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{2}= 'CombineBasicAAL';

%CombineGeneral-Overlap with AAL
fName = spm_select('List',AAL_splits_Dir,'^CombineGeneral.*.nii');
if(~isempty(fName))
    fName = cellstr(fName);
else
    error('fName is empty! Check AAL-split directory!');
end
if(length(fName)>1)
    [IndexSel,ok] = listdlg('ListString',fName,'Name','Which?','SelectionMode','single');
    if(ok)
        MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{3}= [AAL_splits_Dir,filesep,fName{IndexSel}];
    else
        MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{3}= [AAL_splits_Dir,filesep,fName{1}];
    end
else
    MVSscaling.Masks.AALsplit.EffectOfMRI.MPaths{3}    = [AAL_splits_Dir,filesep,fName{1}];
end
MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{3}= 'CombineGeneralAAL';


%% Get "EffectsOfMRI"-Mask AAL-splits and make masks for all such that user
%% can choose later in the display which split to look at.
MVSscaling.Masks.AALsplit.EffectOfMRI = MakeAllAALSplitMasks(MVSscaling.Masks.AALsplit.EffectOfMRI,MVSscaling.Masks.WholeBrainRaw);

%% FUTURE EXTENSION: 1.) output mask 5 as NIFTI, i.e. significant voxels of IC WITHOUT significant voxels "Effects of MRI" 
%%                   2.) Make AAL Overlap & Extension for this mask as well
% get AAL labels
[FunctionDir] = fileparts(mfilename('fullpath'));
MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI = load([FunctionDir,filesep,'AAL_Atlas',filesep,'AAL_Labels.mat']);

% output mask no 5, i.e. significant voxels of IC WITHOUT significant voxels "Effects of MRI" 
V_mask      = spm_vol(MVSscaling.Masks.MPaths{1});     %get whole brain mask as "blueprint"
NII_mask    = nifti(MVSscaling.Masks.MPaths{1});       %get volume
IndicesMask = find(MVSscaling.Masks.WholeBrainRaw~=0); %get indices in the volume belonging to whole brain mask
Y           = zeros(V_mask.dim);                       %setup empty volume to be filled with values of mask no 5, i.e. significant voxels of IC WITHOUT significant voxels "Effects of MRI"
for idx = 1:length(IndicesMask)
    Y(IndicesMask(idx)) = MVSscaling.Masks.AllMasks(idx,5); %put mask 5 into empty volume
end
V_file = V_mask; %prep vol struct for output of newly filled volume
V_file.fname = [MVSscaling.SPMdir,filesep,regexprep(regexprep(MVSscaling.Masks.AllMasksName{5},'SignifVoxelsIC',MVSscaling.Masks.AllMasksName{2}),' ','_'),'.nii'];
if(V_file.dt(1)<16) %unneccessary here but who cares
    V_file.dt(1) = 16; %not necessary but save
end
V_file = spm_write_vol(V_file, Y);

MVSscaling.Masks.MPaths{5} = V_file.fname;
clear V_mask NII_mask IndicesMask Y V_file

% make overlap with AAL & expand
AAL_splits_Dir = [MVSscaling.SPMdir,filesep,'ScalingAnalysis',filesep,'AAL_split',filesep,'IC_WITHOUT_EffectOfMRI'];
[ExpandAtlasOverlap_OutputStruct] = MakeAALOverlap(MVSscaling.Masks.MPaths{5},AAL_splits_Dir,'IC_WITHOUT_EffectOfMRI');


%Expanded-Overlap with AAL
fName = spm_select('List',AAL_splits_Dir,'^Expanded.*.nii');
if(~isempty(fName))
    fName = cellstr(fName);
else
    error('fName is empty! Check AAL-split directory!');
end
if(length(fName)>1)
    [IndexSel,ok] = listdlg('ListString',fName,'Name','Which?','SelectionMode','single');
    if(ok)
        MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{1}    = [AAL_splits_Dir,filesep,fName{IndexSel}];
    else
        MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{1}    = [AAL_splits_Dir,filesep,fName{1}];
    end
else
    MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{1}    = [AAL_splits_Dir,filesep,fName{1}];
end
MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames{1}= 'ExpandedOverlapAAL'; %assign names to this variable such that listdlg can be used to select which one to use in Display-function

%CombineBasic-Overlap with AAL
fName = spm_select('List',AAL_splits_Dir,'^CombineBasic.*.nii');
if(~isempty(fName))
    fName = cellstr(fName);
else
    error('fName is empty! Check AAL-split directory!');
end
if(length(fName)>1)
    [IndexSel,ok] = listdlg('ListString',fName,'Name','Which?','SelectionMode','single');
    if(ok)
        MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{2}= [AAL_splits_Dir,filesep,fName{IndexSel}];
    else
        MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{2}= [AAL_splits_Dir,filesep,fName{1}];
    end
else
    MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{2}    = [AAL_splits_Dir,filesep,fName{1}];
end
MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames{2}= 'CombineBasicAAL';

%CombineGeneral-Overlap with AAL
fName = spm_select('List',AAL_splits_Dir,'^CombineGeneral.*.nii');
if(~isempty(fName))
    fName = cellstr(fName);
else
    error('fName is empty! Check AAL-split directory!');
end
if(length(fName)>1)
    [IndexSel,ok] = listdlg('ListString',fName,'Name','Which?','SelectionMode','single');
    if(ok)
        MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{3}= [AAL_splits_Dir,filesep,fName{IndexSel}];
    else
        MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{3}= [AAL_splits_Dir,filesep,fName{1}];
    end
else
    MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.MPaths{3}    = [AAL_splits_Dir,filesep,fName{1}];
end
MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames{3}= 'CombineGeneralAAL';

%% make all masks for each split of IC_WITHOUT_EffectOfMRI
MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI = MakeAllAALSplitMasks(MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI,MVSscaling.Masks.WholeBrainRaw);

end

%% MakeAllAALSplitMasks
function SplitInput = MakeAllAALSplitMasks(SplitInput,WholeBrainRaw)
% This function makes all the masks as a logical matrix.

for IndSplit = 1:length(SplitInput.SplitNames)
    %get AAL split of EffectOfMRI
    NII_tmp                 = nifti(SplitInput.MPaths{IndSplit});
    tmp_AALsplit_Input      = round(NII_tmp.dat(:));
    tmp_AALsplit_Input(WholeBrainRaw~=1) = 0; %just to be sure.
    clear NII_tmp
    
    %make masks for each split
    UniqueAALsplitLabelNums  = unique(tmp_AALsplit_Input);
    if(any(UniqueAALsplitLabelNums==0))
        UniqueAALsplitLabelNums(UniqueAALsplitLabelNums==0) = [];
    end
    N_AALsplitMasks          = length(UniqueAALsplitLabelNums);
    SplitInput.Split(IndSplit).AllMasks       = zeros(length(find(WholeBrainRaw==1)),N_AALsplitMasks); %all masks that can be made from the split according to AAL.
    SplitInput.Split(IndSplit).AllMasksName   = cell(N_AALsplitMasks,1);
    SplitInput.Split(IndSplit).NVoxelAllMasks = zeros(N_AALsplitMasks,1);
    SplitInput.Split(IndSplit).LabelInds      = zeros(N_AALsplitMasks,1); %indices to keep track of what is in here.
    for IndSplitMasks = 1:length(UniqueAALsplitLabelNums)
        SplitInput.Split(IndSplit).LabelInds(IndSplitMasks) = UniqueAALsplitLabelNums(IndSplitMasks);
        SplitInput.Split(IndSplit).AllMasksName{IndSplitMasks} = SplitInput.Labels(IndSplit).Labels{UniqueAALsplitLabelNums(IndSplitMasks)};
        if(~isempty(tmp_AALsplit_Input==UniqueAALsplitLabelNums(IndSplitMasks)))
            SplitInput.Split(IndSplit).AllMasks(find(round(tmp_AALsplit_Input(WholeBrainRaw==1))==UniqueAALsplitLabelNums(IndSplitMasks)),IndSplitMasks) = 1;
        else
            error([SplitInput.Split(IndSplit).AllMasksName{IndSplitMasks},' is empty for label ',num2str(UniqueAALsplitLabelNums(IndSplitMasks)),' "',SplitInput.Labels(IndSplit).Labels{UniqueAALsplitLabelNums(IndSplitMasks)},'".']);
        end
        SplitInput.Split(IndSplit).NVoxelAllMasks(IndSplitMasks) = length(find(SplitInput.Split(IndSplit).AllMasks(:,IndSplitMasks)==1));
    end
end

end

%% MakeROIsplits
function MVSscaling = MakeROIsplits(MVSscaling)
% This function makes the matrix for the selection of voxels belonging to
% the ROIs 1:NROIs found in the NIFTI-mask.

%% select ROI-Mask NIFTIs
MVSscaling.Masks.ROIsplit.ROI_MaskPath = cell(2,1);
MVSscaling.Masks.ROIsplit.ROI_MaskPath{1} = spm_select(1,'image','Select ROI-masks for "Effect of MRI"...');
MVSscaling.Masks.ROIsplit.ROI_MaskPath{2} = spm_select(1,'image','Select ROI-masks for "IC-WITHOUT-Effect of MRI"...');

%% name of ROI-NIFTIs?
[tmp, fname1] = fileparts(MVSscaling.Masks.ROIsplit.ROI_MaskPath{1});
[tmp, fname2] = fileparts(MVSscaling.Masks.ROIsplit.ROI_MaskPath{2}); clear tmp
MVSscaling.Masks.ROIsplit.ROI_Name = inputdlg({'Name 1st ROIs: '; 'Name 2nd ROIs: '},'Name for ROIs?',1,{['"Effect of MRI" (',fname1,')'];['"IC-WITHOUT-Effect of MRI" (',fname2,')']});

%% select Labels
MVSscaling.Masks.ROIsplit.ROI_LabelsPath = cell(2,1);
MVSscaling.Masks.ROIsplit.ROI_LabelsPath{1} = spm_select(1,'mat','Select ROI-LABELS for "Effect of MRI"...');
MVSscaling.Masks.ROIsplit.ROI_LabelsPath{2} = spm_select(1,'mat','Select ROI-LABELS for "IC-WITHOUT-Effect of MRI"...');    

%% treat mask-files
%% Effects of MRI
MVSscaling.Masks.ROIsplit.EffectOfMRI = load(MVSscaling.Masks.ROIsplit.ROI_LabelsPath{1}); %must contain a variable "Labels"
if(~isfield(MVSscaling.Masks.ROIsplit.EffectOfMRI,'Labels'))
    disp('fieldnames(MVSscaling.Masks.ROIsplit.EffectOfMRI): ');
    fieldnames(MVSscaling.Masks.ROIsplit.EffectOfMRI)
    error('Could not find the required field "Labels"! Check MVSscaling structure ROIsplits.');
end
[MVSscaling.Masks.ROIsplit.EffectOfMRI.AllMasks,MVSscaling.Masks.ROIsplit.EffectOfMRI.NVoxelAllMasks] = MakeAllROIMasks(MVSscaling.Masks.ROIsplit.ROI_MaskPath{1},MVSscaling.Masks.WholeBrainRaw);
if(length(MVSscaling.Masks.ROIsplit.EffectOfMRI.NVoxelAllMasks)~=length(MVSscaling.Masks.ROIsplit.EffectOfMRI.Labels))
    error(['EffectOfMRI: Import of "',MVSscaling.Masks.ROIsplit.ROI_Name{1},'" resulted in ',num2str(length(MVSscaling.Masks.ROIsplit.EffectOfMRI.NVoxelAllMasks)),' masks in the ROI-NIFTI, but the number of Labels is different (==',num2str(length(MVSscaling.Masks.ROIsplit.EffectOfMRI.Labels)),').']);
end

%% IC WITHOUT EffectOfMRI
MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI = load(MVSscaling.Masks.ROIsplit.ROI_LabelsPath{2}); %must contain a variable "Labels"
if(~isfield(MVSscaling.Masks.ROIsplit.EffectOfMRI,'Labels'))
    disp('fieldnames(MVSscaling.Masks.ROIsplit.EffectOfMRI): ');
    fieldnames(MVSscaling.Masks.ROIsplit.EffectOfMRI)
    error('Could not find the required field "Labels"! Check MVSscaling structure ROIsplits.');
end
[MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.AllMasks,MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.NVoxelAllMasks] = MakeAllROIMasks(MVSscaling.Masks.ROIsplit.ROI_MaskPath{2},MVSscaling.Masks.WholeBrainRaw);
if(length(MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.NVoxelAllMasks)~=length(MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.Labels))
    error(['IC_WITHOUT_EffectOfMRI: Import of "',MVSscaling.Masks.ROIsplit.ROI_Name{2},'" resulted in ',num2str(length(MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.NVoxelAllMasks)),' masks in the ROI-NIFTI, but the number of Labels is different (==',num2str(length(MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.Labels)),').']);
end

end

%% MakeAllROIMasks
function [AllMasks,NVoxelAllMasks] = MakeAllROIMasks(ROI_MaskPath,WholeBrainRaw)
% This function makes all the masks as a logical matrix.

%get ROI-NIFTI
NII_tmp        = nifti(ROI_MaskPath);
tmp_Input      = round(NII_tmp.dat(:));
tmp_Input(WholeBrainRaw~=1) = 0; %just to be sure.
clear NII_tmp

%make masks for each split
UniqueROILabelNums  = unique(tmp_Input);
if(any(UniqueROILabelNums==0))
    UniqueROILabelNums(UniqueROILabelNums==0) = [];
end
N_ROIMasks = length(UniqueROILabelNums);

%% init masks
AllMasks       = zeros(length(find(WholeBrainRaw==1)),N_ROIMasks);
NVoxelAllMasks = zeros(N_ROIMasks,1);

%% create masks
for IndROIMasks = 1:length(UniqueROILabelNums)
    if(~isempty(tmp_Input==UniqueROILabelNums(IndROIMasks)))
        AllMasks(find(round(tmp_Input(WholeBrainRaw==1))==UniqueROILabelNums(IndROIMasks)),IndROIMasks) = 1;
    else
        error(['ROI ',num2str(IndROIMasks),' is empty!']);
    end
    NVoxelAllMasks(IndROIMasks) = length(find(AllMasks(:,IndROIMasks)==1));
end

end