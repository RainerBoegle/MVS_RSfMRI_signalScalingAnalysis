% This script allows to automatic analysis of dual regression stage 2 data
% from two MRIs for ScalingValues Lambda with median 2*sqrt(2) (MVS present) 
% and correlations with covariates using a permutation test.
% This is done using searchlights making each voxel and it's neighbors into
% a ROI and prejecting the results back.
% In this way we can look at the median subject, the aggregate and the median searchlight data.
% All as described in (Boegle et al., 2015) [hopefully submitted, reviewed and published soon...].
%
%
%V1.5
%Date: V1.5(20.09.2015): include adjusted version of amplitudes (e.g. divided by stdev of residuals or such). V1.1(03.09.2015): adjustments for allowing searchlight averaging of amplitudes before scaling analysis in addition to original method. V1.0(22.7.2015): initial implementation based on test script for analysis of scaling data.
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% settings
if(strcmp('Yes',questdlg('Display basic overview figures?','DispOverview?','Yes','No','No')))
    DispOverview = 1; %display overview figures
else
    DispOverview = 0; %Do NOT display overview figures
end

%% select folder with dual regression stage 2 amplitudes/weights
DualRegDir = spm_select(1,'dir','Select dual-regression directory containing stage2 amplitudes...');

%% are there adjusted amplitudes?
TmpAdjFiles = spm_select('List',DualRegDir,'^Adj.*.nii');
if(~isempty(TmpAdjFiles))
    if(strcmp('Use adjusted amplitudes',questdlg({['This folder "',DualRegDir,'"']; 'seems to contain ADJUSTED DUAL REGRESSION AMPLITUDES as well as the regular dual regression amplitudes.'; ' '; 'Do you want to use the regular amplitudes or the adjusted amplitudes?'},'Which amplitudes?','Use regular amplitudes','Use adjusted amplitudes','Use adjusted amplitudes')))
        ChoiceAdjAmps = questdlg('Which kind of adjusted amplitudes do you want to use? Adjusted by standard deviation of the residuals after dual regression or the standard deviation of the full timeseries?','Which adjusted amplitudes?','AdjStdevRes','AdjStdevInputs','Define','AdjStdevRes');
        switch(ChoiceAdjAmps)
            case {'AdjStdevRes','AdjStdevInputs'}
                PrefixStr = [ChoiceAdjAmps,'_'];
                disp(['Using ADJUSTED dual regression amplitudes "',PrefixStr,'"...']);
            case 'Define'
                answer = inputdlg({'Define prefix for "dr_stage2_ic????.nii"-file (i.e. "PREFIX"_dr_stage2_ic????): '},'Prefix?',1,{'Adj'});
                if(~isempty(answer{1}))
                    ChoiceAdjAmps = answer{1};
                    PrefixStr     = [answer{1},'_'];
                    disp(['Using ADJUSTED dual regression amplitudes "',PrefixStr,'"...']);
                else
                    disp('Using regular dual regression amplitudes, NOT the adjusted...');
                    PrefixStr = '';
                end
        end
    else
        disp('Using regular dual regression amplitudes...');
        ChoiceAdjAmps = '';
        PrefixStr     = '';
    end
else
    disp('Only regular dual regression amplitudes present. Will use these...');
    ChoiceAdjAmps = '';
    PrefixStr     = '';
end

%% search for files either (MVSscaling_IC#.mat) or (dualreg stage2; .filelist; mask.nii) in folder AND ask for covariates struct 
%MVSscaling_IC#.mat files
MVSscaling_ICFiles = spm_select('List',DualRegDir,['^',ChoiceAdjAmps,'MVSscaling_IC.*.mat']); %The files themselves
if(~isempty(MVSscaling_ICFiles))
    if(strcmp('Yes',questdlg(['Use ',ChoiceAdjAmps,'MVSscaling_IC#.mat files (',num2str(length(MVSscaling_ICFiles)),') that are already in the dual-regression folder, or overwrite?'],'Overwrite?','Yes','Overwrite','Overwrite')))
        UseMVSscalingFiles = 1;
        MVSscaling_ICFiles = cellstr(MVSscaling_ICFiles); %as cellstring
        
        MVSscaling_ICFilePaths = cell(size(MVSscaling_ICFiles)); %full path
        for IndIC = 1:length(MVSscaling_ICFiles)
            MVSscaling_ICFilePaths{IndIC} = [DualRegDir,filesep,MVSscaling_ICFiles{IndIC}];
        end
        NICs = length(MVSscaling_ICFiles);
    else
        disp(['Will OVERWRITE ',ChoiceAdjAmps,'MVSscaling_IC#.mat files found in "',DualRegDir,'".']);
        UseMVSscalingFiles = 0;
    end
else
    UseMVSscalingFiles = 0;
    if(isempty(PrefixStr))
        disp(['No ',ChoiceAdjAmps,'MVSscaling_IC#.mat files found in "',DualRegDir,'". --> Will have to create them from "',PrefixStr,'"dual-regression stage 2 files...']);
    else
        disp(['No ',ChoiceAdjAmps,'MVSscaling_IC#.mat files found in "',DualRegDir,'". --> Will have to create them from "',PrefixStr,'"dual-regression stage 2 files...']);
    end
end

%dualreg stage2; .filelist; mask.nii
if(~UseMVSscalingFiles)
    %dualreg stage2
    DualRegStage2Files = spm_select('List',DualRegDir,['^',PrefixStr,'dr_stage2_ic.*.nii']); %The files themselves
    if(~isempty(DualRegStage2Files))
        DualRegStage2Files = cellstr(DualRegStage2Files); %as cellstring
        
        DualRegStage2FilePaths = cell(size(DualRegStage2Files)); %full path
        for IndIC = 1:length(DualRegStage2Files)
            DualRegStage2FilePaths{IndIC} = [DualRegDir,filesep,DualRegStage2Files{IndIC}];
        end
    else
        error(['No dual-regression stage2 files found in "',DualRegDir,'".']);
    end
    NICs = length(DualRegStage2Files);
    
    %.filelist
    FileListPath = spm_select('List',DualRegDir,'^.filelist');
    if(isempty(FileListPath))
        FListAutoSelect = 0;
        FileListPath = spm_select(1,'any','Select ".filelist" from MELODIC analysis...');
    else
        FListAutoSelect = 1;
        FileListPath = [DualRegDir,filesep,FileListPath]; %make full path
    end
    
    %mask.nii (MELODIC)
    MELODICmaskPath = spm_select('List',DualRegDir,'^mask.nii'); %The mask nifti file
    if(isempty(MELODICmaskPath))
        MELODICmaskAutoSelect = 0;
        MELODICmaskPath = spm_select(1,'image','Select whole brain mask from MELODIC analysis...');
        MELODICmaskPath = MELODICmaskPath(1:end-2);
    else
        MELODICmaskAutoSelect = 1;
        MELODICmaskPath = [DualRegDir,filesep,MELODICmaskPath];
    end
else
    FListAutoSelect = 0;
    MELODICmaskAutoSelect = 0;
end

%% list options, e.g. select ICs of interest with option to process all and so on
% 1.select ICs or take all
% 2.select .filelist or take the one found?
% 3.select mask.nii or take the one found?
% 4.define or load searchlight definition?
%   --> if define then ask for SLmask that should be in same resolution and orientation as mask.nii, ie. the whole brain mask. (SLmask specified which voxels can be center voxels)
%       Also ask for NHood of course. NB: here we always restrict to SLmask! (which should be a tighter brain/gray matter/whatever mask)
% 5.get covariates for eye movement based prediction of scaling
% 6.select covariates *.mat file
%   --> suggest permutation test 
%       NB: test only for corr_perm.*sign(corr_org) >= abs(corr_org)
%           such that if corr_perm has other sign, it is not a failure of the permtest
%           AND this will work for any sign that corr_org has.
%           --> i.e. if corr_org is negative then we only want "more negative" correlations to fail the test
%           --> i.e. if corr_org is positive then we only want "more positive" correlations to fail the test
%     --> suggest parameters (NPermTestFirstRun = 100; NPermTestFinal = 10000; NFailSkip = 5;)
%         NB: Make a first run over all voxels doing 1000 permutation tests at maximum.
%             Those that pass, i.e. NFail<=1 (p<=0.001) get a second (final) run (should be much less than
%             whole brain) with 10000 permutation tests.
%             TO SPEED THINGS UP EVEN MORE, stop doing permutation tests for each voxel that has more then 5 failed tests,
%             because these are not significant anyways and can not get significant.
%             I.e. if we have done 40 permutation tests and have reached 5 failed ones, we don't need to do the remaining 60 (or 9960) other tests,
%             because it is clearly not gonna significant anyways.
%


%1.which ICs?
if(UseMVSscalingFiles)
    IC_SelStr = cell(length(MVSscaling_ICFiles),1);
    for Ind = 1:length(IC_SelStr)
        IC_SelStr{Ind} = MVSscaling_ICFiles{Ind}((length([ChoiceAdjAmps,'MVSscaling_IC'])+1):(length(MVSscaling_ICFiles{Ind})-4));
    end
else
    IC_SelStr = cellstr(num2str((1:NICs)'));
end
[ICSelection,UseSomeICsOnly] = listdlg('ListString',cellstr(num2str((1:NICs)')),'SelectionMode','multiple','InitialValue',1:NICs,'Name','ICs?','PromptString','Select ICs of interest','CancelString','Use All');
if(~UseSomeICsOnly) %i.e. use all
    ICSelection = 1:NICs;
else
    if(isempty(ICSelection))
        error('No ICs of interest selected!');
    end
end

if(~UseMVSscalingFiles)
    %2.use this filelist? (if autoselected)
    if(FListAutoSelect)
        if(strcmp('No, select another one now.',questdlg('Do you want to use the ".filelist" that was found in the dualreg directory or select another one?','???".filelist"???','Yes, use this one.','No, select another one now.','Yes, use this one.')))
            FileListPath = spm_select(1,'any','Select ".filelist" from MELODIC analysis...');
        end
    end
    %get filelist INFORMATION
    FileListInfo = Determine_SubjNrScannerRun_from_FileList(FileListPath);
    
    %3.select mask.nii or take the one found? (if autoselected)
    if(MELODICmaskAutoSelect)
        if(strcmp('No, select another one now.',questdlg('Do you want to use the "mask.nii" nifti (WHOLE BRAIN MASK) that was found in the dualreg directory or select another one?','"mask.nii"?','Yes, use this one.','No, select another one now.','Yes, use this one.')))
            MELODICmaskPath = spm_select(1,'image','Select whole brain mask from MELODIC analysis...');
            MELODICmaskPath = MELODICmaskPath(1:end-2);
        end
    end
else
    load(MVSscaling_ICFilePaths{1}); %use the first as template
    FileListInfo    = MVSscaling.Design;
    MELODICmaskPath = MVSscaling.Mask.MaskNII_FilePath;
end

%4.define or load searchlight definition?
if(strcmp('Define NEW',questdlg('Define a new searchlight or select an existing definition?','SL definition?','Define NEW','Load definition','Load definition')))
    SLmaskPath = spm_select(1,'image','Select SEARCHLIGHT brain mask, i.e. all voxels that are allowed searchlight centers...');
    
    answer_SL_NHood = inputdlg({'NHood= '},'SLdef',1,{'1'});
    NHood = eval(answer_SL_NHood{1});
    
    [tmp,WBmaskName] = fileparts(MELODICmaskPath);
    [tmp,SLmaskName] = fileparts(SLmaskPath); clear tmp %backwards compatible
    answer_SLname   = inputdlg({'SLdef FileName= '},'SLdef',1,{['SLdef_NHood',num2str(NHood),'_crop_',SLmaskName,'_in_',WBmaskName]});
    if(strcmp('Yes',questdlg(['Save SearchLight definition "',answer_SLname{1},'.mat" in dualreg directory or select a directory?'],'Where to put SLdef.mat?','Yes','Select','Yes')))
        SLdefOutPath = [DualRegDir,filesep,answer_SLname{1},'.mat'];
    else
        SLdefOutPath = [spm_select(1,'dir',['Select output directory for SLdefinition "',answer_SLname{1},'.mat"...']),filesep,answer_SLname{1},'.mat'];
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
%4b.Output ParameterEstimate NIFTI-files and show overlay?
ChoiceParamEstNIFTIfiles = questdlg('Create all ParameterEstimate NIFTI-files for Overlay after performing MVSscalingTest?','ParamEst NIFTI-files?','Yes','No','Yes, but do not display now.','Yes');

%5.Make eye movement based prediction for scaling values?
if(strcmp('Yes',questdlg('Try to predict scaling values from eye movements? (sqrt(2)*lambda_EyeMovementsScaling)','Lambda=?=sqrt(2)*lambda','Yes','No','Yes')))
    UsePredictLambdaFromlambda = 1;
    lambda_CovarsFile = spm_select(1,'mat','Select Covars.mat file...');
    TmpCovarsStruct   = GetCovarsViaSubjNrs(lambda_CovarsFile,FileListInfo);
    if(~isempty(find(~cellfun(@isempty,strfind(TmpCovarsStruct.CovarName,'lambda')))))
        lambda_EyeMovements = TmpCovarsStruct.Regressors(:,find(~cellfun(@isempty,strfind(TmpCovarsStruct.CovarName,'lambda'))));
    else
        error(['lambda (scaling of eye movements) not found in "',lambda_CovarsFile,'".']);
    end
else
    UsePredictLambdaFromlambda = 0;
end 

%6.Use covariates with permtest?
choiceCovarPermtest = questdlg({'Use covariates for correlation analysis (on scaling values & amplitudes)?'; 'Covariates and permutation test can be defined here but applied later, i.e. fist do normal statistics (fast) then later run permutation test (slow) in another session or even another machine (need to transfer data in that case).'},'Covars & Permtest-Corr?','Yes','No','Define, but skip for now.','Yes');
switch(choiceCovarPermtest)
    case 'Yes'
        DoPermTest           = 1; %set it up
        ContinueWithPermTest = 1; %run it in this session
    case 'Define, but skip for now.'
        DoPermTest           = 1; %set it up
        ContinueWithPermTest = 0; %run it in another session...
    case 'No'
        DoPermTest           = 0; %Don't do a permtest --> no setup
        ContinueWithPermTest = 0; %don't do it...
end
if(DoPermTest) %need to define it
    if(UsePredictLambdaFromlambda)
        if(strcmp('Yes',questdlg({'Use the same covariates *.mat file as for the prediction of Lambda from eye movement scaling value lambda?'; ['i.e. "',lambda_CovarsFile,'".']},'Reuse Covar-File?','Yes','No, select another one.','Yes')))
            PermTest.CovarsMatFilePath = lambda_CovarsFile;
        else
            PermTest.CovarsMatFilePath = spm_select(1,'mat','Select Covars.mat file...');
        end
    else
        PermTest.CovarsMatFilePath = spm_select(1,'mat','Select Covars.mat file...');
    end
    PermTest.CovarsStruct          = GetCovarsViaSubjNrs(PermTest.CovarsMatFilePath,FileListInfo);
    hinfo = helpdlg({'NPermTest is the number of tests to perform in total (if NFail<NFailSkip)';'NFailSkip is the maximum number of tests that may fail before we skip the current searchlight altogether.'; ' '; 'I.e. if NPermTest= 10000 and NFailSkip= 5, then it will be tried to do 10000 tests for each searchlight, but each one with 5 or more failed test will just be skipped, even if 10000 test have not been reached, e.g. if only 500 tests have been done by then we will still skip to the next.'; 'This saves time, alot of time actually!'},'Info PermTest');
    uiwait(hinfo);
    answer_perm = inputdlg({'NPermTests= ';'NFailSkip= '},'PermTest',1,{'10000';'5'});
    PermTest.NPermTest = eval(answer_perm{1});
    PermTest.NFailSkip = eval(answer_perm{2});
    
    disp(' ');
    PermTest.PermMatrix = zeros(FileListInfo.NSubjs,PermTest.NPermTest); %the permutations for the covariates
    for IndTest = 1:PermTest.NPermTest
        PermTest.PermMatrix(:,IndTest) = randperm(FileListInfo.NSubjs);
    end
    PermTest.PermMatrix_Rho_pVals_InfoStr = {'Dimension is (NPermTest, MCovar, Rho==1/P==2, Pearson==1/Spearman==2).'; 'Use "InspectPermMatrix_Rho_pVals.m" to see the results of a correlation of the original covariates with their PERMUTED counterparts.'};
    PermTest.PermMatrix_Rho_pVals         = []; %init later (NPermTest, MCovar, Rho==1/P==2, Pearson==1/Spearman==2). Use "InspectPermMatrix_Rho_pVals.m" to see the results.
    %Results init
    PermTest.Results.DimInfoStr = {'All(i.e. in each cell) have dimensions (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2).'; '[Except ".CurrVoxInd" which just keeps track if the current voxel in the analysis, such that continuation is possible.]'};
    PermTest.Results.Rho        = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %correlation at voxel for UNPERMUTED covariates
    PermTest.Results.pCorr      = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %p value from corr-function of correlation at voxel for UNPERMUTED covariates
    PermTest.Results.NFail      = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %the number of failed tests
    PermTest.Results.NTests     = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %the number of test that were actually done.
    PermTest.Results.pPermTest  = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %NFail./NTests %assume equal number of failed tests until full number test, i.e rate stays the same. This should not be a big problem if NFail is >=5, because then (roughly) we would have all rejected tests anyways.
    PermTest.Results.CurrVoxInd = []; %(NVoxel,1) current voxel that we are working on, if empty then we start at the first, if any number then we continue at that one.
else
    PermTest = [];
end

%% check if output directory for MVSscalingTest exists already (based on searchlight that is used)
if(~exist([DualRegDir,filesep,'Results',ChoiceAdjAmps,'_',answer_SLname{1}],'dir'))
    MVSscalingTest_OutDir = [DualRegDir,filesep,'Results',ChoiceAdjAmps,'_',answer_SLname{1}];
    mkdir(MVSscalingTest_OutDir);
else
    if(strcmp('Yes',questdlg(['Overwrite results in Folder "',[DualRegDir,filesep,'Results',ChoiceAdjAmps,'_',answer_SLname{1}],'"?'],'Overwrite?','Yes','No','Yes')))
        MVSscalingTest_OutDir = [DualRegDir,filesep,'Results',ChoiceAdjAmps,'_',answer_SLname{1}];
    else
        MVSscalingTest_OutDir = spm_select(1,'dir','Select output directory for MVSscalingTest-analysis...');
    end
end


%% Do the processing given all the above settings %%
%% Do the processing given all the above settings %%
%% Do the processing given all the above settings %%

%% create/define or load searchlight definition
if(CreateNewSLdef)
    SLight = GenerateSLight(MELODICmaskPath,SLmaskPath,NHood,RestrictAlsoToStartInds,SLdefOutPath);
end

%% create MVSscaling (do all other things that are necessary for this e.g. using data from steps 2. & 3.)
disp(' ');
if(UseMVSscalingFiles)
    disp([ChoiceAdjAmps,'MVSscaling struct will be loaded from previous results stored in dual-regression directory "',DualRegDir,'"...']);
else
    for IndIC = 1:length(ICSelection)
        disp(['Creating "MVSscaling"-struct for IC ',num2str(ICSelection(IndIC),'%02g')]);
        MVSscaling = CreateAmplitudeScalingData(DualRegStage2FilePaths{ICSelection(IndIC)},MELODICmaskPath,FileListInfo,SLight);
        
        %% save MVSscaling.mat file
        disp('saving...');
        save([DualRegDir,filesep,ChoiceAdjAmps,'MVSscaling_IC',num2str(ICSelection(IndIC),'%02g'),'.mat'],'MVSscaling');
        disp(' ');
    end
    disp(['Done creating all "',ChoiceAdjAmps,'MVSscaling_IC#.mat" files.']);
end

%% setup ParamEstLambdaTest-struct for MVSscalingTest-struct (statistics for parameter estimate and H0)
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
if(UsePredictLambdaFromlambda)
    ParamEstLambdaTest.LambdaFromEyeMovements.lambda_EyeMovements = lambda_EyeMovements;
    ParamEstLambdaTest.LambdaFromEyeMovements.mH0                 = lambda_EyeMovements.*sqrt(2); %the expected median, ie. null hypothesis is H0(m(Lambda)=lambda_EyeMovements*sqrt(2))
    % continue here...
else
    ParamEstLambdaTest.LambdaFromEyeMovements = [];
end

%% testing of Lambda values
clear MVSscaling 
disp(' ');
for IndIC = 1:length(ICSelection)
    %% assign
    disp(['Creating "',ChoiceAdjAmps,'MVSscalingTest"-struct for IC ',num2str(ICSelection(IndIC),'%02g'),' and testing Lambda values.']);
    load([DualRegDir,filesep,ChoiceAdjAmps,'MVSscaling_IC',num2str(ICSelection(IndIC),'%02g'),'.mat']);
    MVSscalingTest                   = MVSscaling;
    %old(V1.0) now SLight is already included in MVSscaling struct: MVSscalingTest.SLight            = SLight;
    MVSscalingTest.PermTest          = PermTest;
    MVSscalingTest.ParamEstLambdaTest= ParamEstLambdaTest;
    
    %% do the basic testing, i.e. PARAMETER ESTIMATION & SIGNED-RANK TEST OF LAMBDA, WITHOUT PERMUTATION TEST and write to MVSscalingTest
    MVSscalingTest = RunMVSscalingTest(MVSscalingTest);
    
    %% save MVSscalingTest to MVSscalingTest.mat
    disp('saving...');
    save([MVSscalingTest_OutDir,filesep,ChoiceAdjAmps,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat'],'MVSscalingTest');
    disp(' ');
    
    %% Display Overview
    if(DispOverview)
        [H,ax] = DispOverview_ParamEstMVSscalingTest(MVSscalingTest);
    end
    clear MVSscaling MVSscalingTest
end
disp(['Done with testing median(Lambda) and creating all "',ChoiceAdjAmps,'MVSscalingTest_IC#.mat" files.']);
    
%% continue WITH PERMUTATION TEST? (i.e. can be restarted on another machine if wanted or if failed.)
%input need to be covariates mat file and all parameters
%function will inspect MVSscalingTest struct for determining at which voxel the permutation test is was before.
%if at a certain voxel then we continue from there.
%NB: permutations for each test are produced first and saved in MVSscalingTest.PermTest fields
%    this way continuing will be simple

if(DoPermTest)
    if(~ContinueWithPermTest)
        h=helpdlg({'Permutation test will not be done now, but it is always possible to do it later based on saved setup and data.'; ' '; '(Never put off till tomorrow what you can do on the day after tomorrow just as well. -Mark Twain)'},'Do PermTest later...');
        uiwait(h,5);
    else
        %do permtest
        try
            clear MVSscaling MVSscalingTest
        end
        disp(' ');
        disp('Performing permutation test for correlation of covariates with Lambda values.');
        for IndIC = 1:length(ICSelection)
            %% load
            disp(['Loading "',ChoiceAdjAmps,'MVSscalingTest"-struct for IC ',num2str(ICSelection(IndIC),'%02g'),' for PermTest.']);
            load([MVSscalingTest_OutDir,filesep,ChoiceAdjAmps,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat']);
            
            %% do permtest
            MVSscalingTest = RunMVScorrPermTest(MVSscalingTest,[MVSscalingTest_OutDir,filesep,ChoiceAdjAmps,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat']);
            
            %% save MVSscalingTest to MVSscalingTest.mat
            disp('saving...');
            save([MVSscalingTest_OutDir,filesep,ChoiceAdjAmps,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat'],'MVSscalingTest');
            disp(' ');
            clear MVSscaling MVSscalingTest
        end
    end
end



%% DONE
disp('DONE.');
hdone = helpdlg('Done with analysis.','Done.');
uiwait(hdone,3);

%% Create ParameterEstimate NIFTI files for overlay?
switch(ChoiceParamEstNIFTIfiles)
    case {'Yes','Yes, but do not display now.'}
        Limits = [     0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2);    10/4*sqrt(2); 3*sqrt(2)];
        Colors = [0 0 .5;       0 0 1;       0 1 1;      0 .5 0;       0 1 0;  .5 1 0;       1 1 0;      1 .5 0;       1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;     1 1 1];
        if(strcmp('Yes',ChoiceParamEstNIFTIfiles))
            SlicesForDisplay = [-48; -42; -36; -30; -26; -22; -18; -12; -4; 0; +2; +6; +10; 14; +16; +24; +30; +36; +42; +48];
        else
            SlicesForDisplay = NaN(1);
        end
        clear MVSscalingTest
        for IndIC = 1:length(ICSelection)
            %% load
            disp(['Loading "',ChoiceAdjAmps,'MVSscalingTest"-struct for IC ',num2str(ICSelection(IndIC),'%02g'),' for ParameterEstimate NIFTI-file generation.']);
            load([MVSscalingTest_OutDir,filesep,ChoiceAdjAmps,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat']);
            
            MakeOverlayParamEstLambda(MVSscalingTest,Limits,Colors,MVSscalingTest_OutDir,SlicesForDisplay);
            clear MVSscalingTest
        end
        if(~strcmp('Yes',ChoiceParamEstNIFTIfiles))
            disp('Use "DisplayOverlayParamEstLambda.m" to display the created ParameterEstimate NIFTI-files.');
        end
    otherwise
        disp('Use "MakeOverlayParamEstLambda.m" to make ParameterEstimate NIFTI-files and display them as colored overlays.');
end
disp('DONE.');
