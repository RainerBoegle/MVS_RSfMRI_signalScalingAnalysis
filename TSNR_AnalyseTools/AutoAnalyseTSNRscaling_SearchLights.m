%% Ask User if fMRI 4D-Data should be read and temporal SNR calculated 
%% OR if temporal SNR files should be read
ChoiceDataType = questdlg('Do you want to start analysis on the basis of fMRI time series (4D-NIFTIs) or already calculated TSNR images (3D-NIFTI)?','fMRI 4D-Data or TSNR 3D-Data?','4D fMRI timeseries','3D TSNR images','Load TSNR_Design.mat','3D TSNR images');

%% if 4D fMRI data --> make TSNR files
%% Ask user to select 4D fMRI data files (use spm_select(Inf,'any') to select the 4D files not each frame of each single 4D file.
%% Ask user to select output directory for TSNR files
%% Calculate TSNR & generate TSNR_Design.filelist (save in dir as TSNR_Design.mat)
if(strcmp(ChoiceDataType,'4D fMRI timeseries'))
    fMRI4Dfiles = cellstr(spm_select(Inf,'any','Select all 4D fMRI-files of all subjects...'));
    TSNRoutdir  = spm_select(1,'dir','Select output directory for TSNR-files...');
    
    %% select Mask (whole brain)
    MaskPath = spm_select(1,'image','Select the mask to be used for whole brain data extraction...');
    
    %% precess fMRI data & collect TSNR paths
    TSNR_Design.filelist = cell(length(fMRI4Dfiles),1);
    for IndFile = 1:length(fMRI4Dfiles)
        [TSNR,Mu,Stdev,Mask2D,V_out,OutputNIIpath] = CalcTempStats(fMRI4Dfiles{IndFile},MaskPath,TSNRoutdir); clear TSNR Mu Stdev Mask2D V_out
        TSNR_Design.filelist{IndFile} = OutputNIIpath{1}; %the TSNR file
    end
    
    %% save TSNR_Design in TSNRoutdir
    save([TSNRoutdir,filesep,'TSNR_Design.mat'],'TSNR_Design');
    
elseif(strcmp(ChoiceDataType,'3D TSNR images'))
    %% if 3D TSNR data --> select files or TSNR_Design.mat
    TSNR_Design.filelist = cellstr(spm_select(Inf,'any','Select all 3D TSNR-files of all subjects...'));
    TSNRoutdir  = spm_select(1,'dir','Select directory for saving TSNR_Design.mat...');
    
    %% select Mask (whole brain)
    MaskPath = spm_select(1,'image','Select the mask to be used for whole brain data extraction...');
    
    %% save TSNR_Design in TSNRoutdir
    save([TSNRoutdir,filesep,'TSNR_Design.mat'],'TSNR_Design');
    
else
    TSNR_DesignPath = spm_select(1,'mat','Select TSNR_Design.mat...');
    load(TSNR_DesignPath);
    
    %% select Mask (whole brain)
    MaskPath = spm_select(1,'image','Select the mask to be used for whole brain data extraction...');
    
end

%% select directory for saving results
ResultsBaseDir = spm_select(1,'dir','Select BASE output directory for results (additional directories will be created)...');

%% Output ParameterEstimate NIFTI-files and show overlay?
ChoiceParamEstNIFTIfiles = questdlg('Create all ParameterEstimate NIFTI-files for Overlay after performing TSNRscalingTest?','ParamEst NIFTI-files?','Yes','No','Yes, but do not display now.','Yes');

%% if TSNR_Design.SubjNrMRInrRunNr does not exist (or files selected) get subject nr, MRI nr and Run nr from filenames
%% suggest search strings (regexp), -always two: e.g. for subject: Subj_p0\d\d  &   \d\d\d  --> the first one finds the entry in the file names the second one extracts the subject number from the results.
%% do this for MRI nr and Run nr as well --> then save in field TSNR_Design.SubjNrMRInrRunNr matrix (nFiles-x-3) --> output to file TSNR_Design.mat
if(~isfield(TSNR_Design,'SubjNrMRInrRunNr'))
    H = helpdlg({'Need to select subject, MRI and Run numbers to form design specifications.'; 'This will be done by specifying the strings that extract the subject related string, the MRI related string and the Run number related string from the filename of the TSNR files.';'In the next step the numbers will be extracted from these strings.';'Therefore for each such string there are two regexp suggestions necessary.'},'How to select subjects, MRIs & Runs?');
    uiwait(H);
    answer_subj = inputdlg({'Subject regexp string: ';'Subject number regexp string: '},'Subj Number selection',1,{'Subj_p\d\d\d';'\d\d\d'});
    answer_MRI  = inputdlg({'MRI regexp string: ';    'MRI number regexp string: '},    'MRI Number selection', 1,{'MRI\d';       '\d'});
    answer_Run  = inputdlg({'Run regexp string: ';    'Run number regexp string: '},    'Run Number selection', 1,{'rsfMRI_\d';   '\d'});
    
    TSNR_Design.SubjNrMRInrRunNr = zeros(length(TSNR_Design.filelist),3);
    TSNR_Design.SubjNrMRInrRunNr(:,1) = cellfun(@str2num,cellfun((@(x)x),regexp(cellfun((@(x)x),regexp(TSNR_Design.filelist,answer_subj{1},'match')),answer_subj{2},'match')));
    TSNR_Design.SubjNrMRInrRunNr(:,2) = cellfun(@str2num,cellfun((@(x)x),regexp(cellfun((@(x)x),regexp(TSNR_Design.filelist, answer_MRI{1},'match')), answer_MRI{2},'match')));
    TSNR_Design.SubjNrMRInrRunNr(:,3) = cellfun(@str2num,cellfun((@(x)x),regexp(cellfun((@(x)x),regexp(TSNR_Design.filelist, answer_Run{1},'match')), answer_Run{2},'match')));
    
    TSNR_Design.UniqueSubjNrs = unique(TSNR_Design.SubjNrMRInrRunNr(:,1));
    TSNR_Design.UniqueMRInrs  = unique(TSNR_Design.SubjNrMRInrRunNr(:,2));
    
    TSNR_Design.NSubjs = length(unique(TSNR_Design.SubjNrMRInrRunNr(:,1)));
    TSNR_Design.N_MRIs = length(unique(TSNR_Design.SubjNrMRInrRunNr(:,2)));
end

%% define SearchLight or load?
if(strcmp('Define NEW',questdlg('Define a new searchlight or select an existing definition?','SL definition?','Define NEW','Load definition','Load definition')))
    SLmaskPath = spm_select(1,'image','Select SEARCHLIGHT brain mask, i.e. all voxels that are allowed searchlight centers...');
    
    answer_SL_NHood = inputdlg({'NHood= '},'SLdef',1,{'1'});
    NHood = eval(answer_SL_NHood{1});
    
    [tmp,WBmaskName] = fileparts(MaskPath);
    [tmp,SLmaskName] = fileparts(SLmaskPath); clear tmp %backwards compatible
    answer_SLname    = inputdlg({'SLdef FileName= '},'SLdef',1,{['SLdef_NHood',num2str(NHood),'_crop_',SLmaskName,'_in_',WBmaskName]});
    if(strcmp('Yes',questdlg(['Save SearchLight definition "',answer_SLname{1},'.mat" in ResultsBaseDir directory or select a directory?'],'Where to put SLdef.mat?','Yes','Select','Yes')))
        SLdefOutPath = [ResultsBaseDir,answer_SLname{1},'.mat'];
    else
        SLdefOutPath = [spm_select(1,'dir',['Select output directory for SLdefinition "',answer_SLname{1},'.mat"...']),answer_SLname{1},'.mat'];
    end
    
    RestrictAlsoToStartInds = 1; %restrict to start inds coming from SLmask
    CreateNewSLdef = 1; %this switch will lead to the searchlight being created later on, right now we want to define everything and don't have the nerve to wait for this...
else
    SLdefPath = spm_select(1,'mat','Select searchlight definition *.mat file...');
    load(SLdefPath);
    answer_SLname    = cell(1,1);
    [tmp,answer_SLname{1}] = fileparts(SLdefPath); clear tmp
    if(exist('SLight','var'))
        CreateNewSLdef = 0; %this switch will lead to the searchlight being created later on, right now we want to define everything and don't have the nerve to wait for this...
    else
        error('Searchlight definition not loaded OR loaded definition is not saved as structure with name "SLight".');
    end
end
TSNRscalingTest_OutDir = [ResultsBaseDir,filesep,'Results_',answer_SLname{1}];

%% FUTURE: %% use the same searchlight in TSNRscaling.ScalingPerSubject{2} calculation or another one?
%% FUTURE: %% select or define other if necessary

%% setup permutation test (FUTURE!)

%% create or load SLight definition
if(CreateNewSLdef)
    SLight = GenerateSLight(MaskPath,SLmaskPath,NHood,RestrictAlsoToStartInds,SLdefOutPath);
end

%% calculate TSNRscaling.ScalingPerSubject{1} as the simple fraction just like MVSscaling from amplitudes (average over Runs first)
%% calculate TSNRscaling.ScalingPerSubject{2} using searchlight to collect data in space from both runs (of a subject & MRI) and then use nanmedian to get TSNR for searchlight; --> project back and then calculate lambda{2} as fraction
% TSNRscaling.TSNR{1&2}  TSNRscaling.ScalingPerSubject & TSNRscaling.Design = TSNR_Design (TSNR_Design.filelist & TSNR_Design.SubjNrMRInrRunNr)
% TSNRscaling.Design   = TSNR_Design;
% TSNRscaling.Mask.MaskNII_FilePath = MaskPath;
% TSNRscaling.SLight   = SLight;
% TSNRscaling.TSNR              = cell(2,1); %TSNR{1or2}   will be (NVoxelWholeBrain,NSubj,NMRI); -i.e. average over Runs OR in case of TSNR{2} nanmedian over space (searchlight) and Runs
% TSNRscaling.ScalingPerSubject = cell(2,1); %lambda{1or2} will be (NVoxelWholeBrain,NSubj);      -i.e. fraction of MRI==2 divided by MRI==1. 

TSNRscaling = CalcTSNRscaling(TSNR_Design,MaskPath,SLight);
disp(['Saving TSNRscaling.mat in ResultsBaseDir "',ResultsBaseDir,'".']);
save([ResultsBaseDir,filesep,'TSNRscaling.mat'],'TSNRscaling');

%% set up scaling test variables for parameter estimate test
conf = 0.95; %conficence interval for median (as in boxplot
ParamEstLambdaTest.mH0               = 2*sqrt(2); %the expected median, ie. null hypothesis is H0(m(Lambda)=2*sqrt(2))
ParamEstLambdaTest.mConfInt          = conf; %the confidence interval (0.95)
ParamEstLambdaTest.mConfIntIQRfactor = (1/(norminv(.75)-norminv(.25))) * sqrt(pi/2) * (norminv(1-(1-conf)/2)+(norminv(1-(1-conf)/2)/sqrt(2)))/2; %factor before the IQR (inter quartile range, i.e. 0.25 to 0.75 range of data) that gives the confidence interval (symmetric assumption)
ParamEstLambdaTest.SignRankTest.DimInfoStr = {'All(i.e. in each cell) have dimensions (NVoxel,3).'; 'The second dimension indicates:'; '1==median subject; 2==aggregate data; 3==median searchlight.'};
ParamEstLambdaTest.SignRankTest.zVals      = []; %(NVoxel,3);   %init later %z-score values from signed rank test for median subject, aggregate data and median searchlight
ParamEstLambdaTest.SignRankTest.pVals      = []; %(NVoxel,3);   %init later %   p    values from signed rank test for median subject, aggregate data and median searchlight
ParamEstLambdaTest.SignRankTest.signedrank = []; %(NVoxel,3);   %init later %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
ParamEstLambdaTest.MedianQrtCIwidth_DimInfoStr = {'MedianQrtCIwidth(i.e. in each cell) has dimensions (NVoxel,4,3).'; ' '; 'The second dimension(1:4) indicates: '; '1==median; 2==1stQuartile (25%point); 3==3rdQuartile (75%point); 4==Confidence-Interval width (CIwidth=CIfactor*(Qrt3-Qrt1)/sqrt(NInputData))'; ' '; 'The third dimension indicates:'; '1==median subject; 2==aggregate data; 3==median searchlight.'};
ParamEstLambdaTest.MedianQrtCIwidth            = []; %(NVoxel,4,3); %init later %parameter estimates for median subject, aggregate data and median searchlight (that is the last dim); MiddleDim: 1==Median,2==1stQuartile,3==3rdQuartile,4==CIwidth 

%% RunTSNRscalingTest.m to do MedS, Agg & MedSL estimation and testing as in RunMVSscalingTest.m
%% output (&input) TSNRscalingTest in TSNRscalingTest.mat
TSNRscalingTest                    = TSNRscaling;
TSNRscalingTest.ParamEstLambdaTest = ParamEstLambdaTest; %parameters and fields for performing test and storing results

TSNRscalingTest = RunTSNRscalingTest(TSNRscalingTest); %run test
if(~exist(TSNRscalingTest_OutDir,'dir'))
    mkdir(TSNRscalingTest_OutDir);
end
disp(['Saving TSNRscalingTest.mat in TSNRscalingTest_OutDir "',TSNRscalingTest_OutDir,'".']);
save([TSNRscalingTest_OutDir,filesep,'TSNRscalingTest.mat'],'TSNRscalingTest');

%% RunTSNRscalingCorrPermTest like RunMVScorrPermTest.m (FUTURE!)

%% DONE
disp('DONE.');
hdone = helpdlg('Done with analysis.','Done.');
uiwait(hdone,3);

%% Create ParameterEstimate NIFTI files for overlay? 
%% output the results of the scaling test just like the other results as Raw and Binned files for overlay. Then use the common overlay function
switch(ChoiceParamEstNIFTIfiles)
    case {'Yes','Yes, but do not display now.'}
        Limits = [0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)];
        Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1];
        if(strcmp('Yes',ChoiceParamEstNIFTIfiles))
            SlicesForDisplay = [-48; -42; -36; -30; -26; -22; -18; -12; -4; 0; +2; +6; +10; 14; +16; +24; +30; +36; +42; +48];
        else
            SlicesForDisplay = NaN(1);
        end
        
        MakeOverlayParamEstLambda(TSNRscalingTest,Limits,Colors,TSNRscalingTest_OutDir,SlicesForDisplay);
        if(~strcmp('Yes',ChoiceParamEstNIFTIfiles))
            disp('Use "DisplayOverlayParamEstLambda.m" to display the created ParameterEstimate NIFTI-files.');
        end
    otherwise
        disp('Use "MakeOverlayParamEstLambda.m" to make ParameterEstimate NIFTI-files and display them as colored overlays.');
end
disp('DONE.');
