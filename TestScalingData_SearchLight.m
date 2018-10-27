%% Get MVSscaling struct
ScalingAnaMatPath = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
load(ScalingAnaMatPath);

[BasePath, FName, ext] = fileparts(ScalingAnaMatPath);
MVSscalingTest.InputPath = ScalingAnaMatPath;

%% check for eye movement data
if(exist([BasePath,filesep,'EyeMovementData',filesep,'EyeMovementDataForPlots.mat'],'file'))
    EyeMovementData = load([BasePath,filesep,'EyeMovementData',filesep,'EyeMovementDataForPlots.mat']);
    
    if(strcmp('Yes',questdlg('Use permtest on correlation?','use permtest?','Yes','No','Yes')))
        UsePermTest = 1;
        NPermAnswer = inputdlg({'N Permtest= '},'N Permtest',1,{'1000'});
        NPermTest   = eval(NPermAnswer{1});
    else
        UsePermTest = 0;
        NPermTest   = 0;
    end
end
MVSscalingTest.StatisticsSettings.Corr.UsePermTest = UsePermTest;
MVSscalingTest.StatisticsSettings.Corr.NPermTest   = NPermTest;

% assign subjects' data according to coding of IDs and ANOVA design
IDsubjEyeMovementData = zeros(size(EyeMovementData.InfoStruct.ID_Subj)); %number assign to subjects from eye movement data --> assign from ID_Subj
for IndSubj = 1:length(EyeMovementData.InfoStruct.ID_Subj)
    IDsubjEyeMovementData(IndSubj) = str2num(EyeMovementData.InfoStruct.ID_Subj{IndSubj}(2:end));
end
ScalingFactorsLambda = zeros(size(MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum)); %the scaling factors lambda that we will use for correlation --> need to assign right
SPVcorrData          = zeros(length(MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum),2); %the SPV data that we will use for correlation --> need to assign right
for Ind = 1:length(MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum)
    CorrespondingSubjInd = find(IDsubjEyeMovementData==MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum(Ind)); %the index that fits the subject number relative to eye movement data for the fMRI data that we use here.
    ScalingFactorsLambda(Ind) = EyeMovementData.ScalingFactorsLambda(CorrespondingSubjInd); %this should be the correct scaling data for the correlation analysis
    SPVcorrData(Ind,:) = EyeMovementData.SPVdata(CorrespondingSubjInd,[2,6]); %this should be the correct SPV data for the correlation analysis
end

figure(42); 
boxplot(ScalingFactorsLambda,'notch','on','labels',{'Lambda'}); title('Scaling data for correlation'); hold on
plot(ones(size(ScalingFactorsLambda))-(rand(size(ScalingFactorsLambda))-0.5)./10,ScalingFactorsLambda,'ko'); ylim([0 6.125]); 

%% define Searchlight (NHood and such)
if(strcmp('Generate',questdlg('Generate or load searchlight definition?','SL?','Generate','Load','Generate')))
    answer_SLight = inputdlg({'SearchLight extend NHood= '},'SLight def',1,{'1'});
    SLight = GenerateSLight([],[],eval(answer_SLight{1}),[],[BasePath,filesep,'TestScalingDataSEARCHLIGHT_NHood',num2str(eval(answer_SLight{1}))]);
else
    load(spm_select(1,'mat','Select appropriate searchlight definition *.mat'));
end
MVSscalingTest.SLight = SLight; %assign SLight

%% create directory for output
OutDir = [BasePath,filesep,'TestScalingDataSEARCHLIGHT_NHood',num2str(MVSscalingTest.SLight.NHood)];
if(~exist(OutDir))
    mkdir(OutDir);
end

%% smooth data by applying median of searchlight on each subject and write median to center voxel?
if(strcmp('Yes',questdlg('Apply SearchLight smoothing per subject?','SLightSmooth?','Yes','No','Yes')))    
    UseSLsmooth = 1;
    MVSscaling.ScalingPerSubject_org = MVSscaling.ScalingPerSubject; %for safety
    MVSscaling.ScalingPerSubject     = abs(MVSscaling.ScalingPerSubject); %prep before SLight smooth
    if(strcmp('Always',questdlg('ALWAYS change center value to median of the searchlight or only if center value is an outlier(like boxplot)?','When change?','Always','Outlier','Outlier')))
        [MVSscaling.ScalingPerSubject,MVSscaling.SLightSmooth]  = SLsmooth(MVSscaling.ScalingPerSubject,[],'AssignMedianSL');
    else
        if(strcmp('IfOutlierAssignMedianSL',questdlg('OUTLIERS(like boxplot): Remove outliers (assign it as NaN) or assign nanmedian of the searchlight to center value?','When change?','IfOutlierAssignNaN','IfOutlierAssignMedianSL','IfOutlierAssignMedianSL')))
            [MVSscaling.ScalingPerSubject,MVSscaling.SLightSmooth]  = SLsmooth(MVSscaling.ScalingPerSubject,[],'IfOutlierAssignMedianSL');
        else
            [MVSscaling.ScalingPerSubject,MVSscaling.SLightSmooth]  = SLsmooth(MVSscaling.ScalingPerSubject,[],'IfOutlierAssignNaN');
        end
    end
    save([OutDir,filesep,'MVSscaling_SLsmooth.mat'],'MVSscaling');
else
    UseSLsmooth = 0;
end
MVSscalingTest.UseSLsmooth = UseSLsmooth;

%% define test statistics (need to be robust so I will test the median somehow???)
answer_SetupStats = inputdlg({'ExpectedMedian= '},'Setup stats',1,{'2*sqrt(2)'});
ExpectedMedian = eval(answer_SetupStats{1});
ChoiceTest = questdlg({'Do you want to test if the data''s median is:'; ' '; ['A. "different" from expected median(',answer_SetupStats{1},')']; ['B. "higher" than the expected median(',answer_SetupStats{1},')']; 'OR ';  ['C. "lower" than the expected median(',answer_SetupStats{1},')']},'Which test?','Different','Higher','Lower','Higher');
switch(ChoiceTest)
    case 'Different'
        TailType = 'both';
    case 'Higher'
        TailType = 'right';
    case 'Lower'
        TailType = 'left';
end
MVSscalingTest.StatisticsSettings.answer_SetupStats = answer_SetupStats;
MVSscalingTest.StatisticsSettings.ExpectedMedian    = ExpectedMedian;
MVSscalingTest.StatisticsSettings.ChoiceTest        = ChoiceTest;
MVSscalingTest.StatisticsSettings.TailType          = TailType;

%% create directory for output
OutDir = [BasePath,filesep,'TestScalingDataSEARCHLIGHT_NHood',num2str(MVSscalingTest.SLight.NHood),filesep,'Med',ChoiceTest,regexprep(answer_SetupStats{1},'*|\(|\)','')];
if(~exist(OutDir))
    mkdir(OutDir);
end

%% remove outliers
if(strcmp('Yes',questdlg('Do you want to remove strong outliers from the data and track this in masks?','Track STRONG outliers?','Yes','No','No')))
    MVSscalingTest.StatisticsSettings.RemOutliers.Use = 1;
    MVSscalingTest.StatisticsSettings.RemOutliers.method = 'Threshold';
    ansOutlierThres = [];
    while(isempty(ansOutlierThres))
        ansOutlierThres = inputdlg({'Threshold for STRONG Outliers (High): ';'Threshold for STRONG Outliers (Low): '},'ThresOutliers?',1,{'10';'0.1'});
    end
    MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi = eval(ansOutlierThres{1});
    MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo = eval(ansOutlierThres{2});
    
    if(~isempty(eval(ansOutlierThres{1}))&&~isempty(eval(ansOutlierThres{2})))
        OulierStr = ['RemOutliersThres_Hi',regexprep(ansOutlierThres{1},'*|\(|\)',''),'Lo',regexprep(ansOutlierThres{2},'*|\(|\)','')];
    else
        if(~isempty(eval(ansOutlierThres{1})))
            OulierStr = ['RemOutliersThres_Hi',regexprep(ansOutlierThres{1},'*|\(|\)','')];
        else
            OulierStr = ['RemOutliersThres_Lo',regexprep(ansOutlierThres{2},'*|\(|\)','')];
        end
    end
else
    MVSscalingTest.StatisticsSettings.RemOutliers.Use = 0;
    OulierStr = '';
end
    
%% create directory for output
if(~isempty(OulierStr))
    OutDir = [OutDir,filesep,OulierStr];
    if(~exist(OutDir))
        mkdir(OutDir);
    end
end

%% also try testing the aggregated data in each searchlight? I.e. totally different approach where all data from all voxels in searchlight and from all subjects is aggregated and then tested for effect.
MVSscalingTest.Results.DataQuality.UsableData     = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as all data usable
MVSscalingTest.Results.DataQuality.SomeMissingData= zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as no missing data
MVSscalingTest.Results.DataQuality.OutliersGeneral= zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as no outliers generally
MVSscalingTest.Results.DataQuality.OutliersCenter = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as no outliers at center

MVSscalingTest.Results.Aggregate.MedianQrtCI = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),5); %init not significant
MVSscalingTest.Results.Aggregate.SkewKurt    = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init 
MVSscalingTest.Results.Aggregate.p           = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.h           = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.signedrank  = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.zval        = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.p_sign      = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.h_sign      = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.sign        = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.zval_sign   = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.Aggregate.CIoverlap   = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant

MVSscalingTest.Results.MedianSubj.MedianQrtCI = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),5); %init not significant
MVSscalingTest.Results.MedianSubj.SkewKurt    = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init 
MVSscalingTest.Results.MedianSubj.p           = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.h           = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.signedrank  = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.zval        = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.p_sign      = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.h_sign      = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.sign        = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.zval_sign   = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSubj.CIoverlap   = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant

MVSscalingTest.Results.MedianSLight.MedianQrtCI = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),5); %init not significant
MVSscalingTest.Results.MedianSLight.SkewKurt    = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init 
MVSscalingTest.Results.MedianSLight.p           = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.h           = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.signedrank  = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.zval        = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.p_sign      = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.h_sign      = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.sign        = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.zval_sign   = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.CIoverlap   = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init not significant
MVSscalingTest.Results.MedianSLight.Corr.Rho    = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.Corr.P      = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.Corr.Pperm  = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.Corr.Pperm2 = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AmpCorr.Rho    = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AmpCorr.P      = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm  = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm2 = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2 
MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Rho    = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %abs AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AbsAmpCorr.P      = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %abs AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm  = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %abs AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2
MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm2 = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),2); %abs AMPLITUDE correlation init not significant %normal corr==1 %rank corr==2 

%% run SLight analysis
if(size(MVSscaling.ScalingPerSubject,1)<length(MVSscalingTest.SLight.SLightIndsInMaskCell))
    error('more searchlight centers than voxels in scaling analysis data!');
end
if(~UsePermTest)
    H_waitbar = waitbar(0,['SearchLight-test centered at every voxel starting...']);
end
%use autosave
if(UsePermTest&&exist([OutDir,filesep,'autosave.mat'],'file'))
    if(strcmp('Yes',questdlg('Continue analysis from "autosave" position?','Continue from autosave?','Yes','No','Yes')))
        load([OutDir,filesep,'autosave.mat']);
        disp(['Restarting from autosave.mat at IndSLCenterVox= ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell))]);
        SLCenterVoxInds = IndSLCenterVox:length(MVSscalingTest.SLight.SLightIndsInMaskCell);
    else
        SLCenterVoxInds = 1:length(MVSscalingTest.SLight.SLightIndsInMaskCell);
    end
else
    SLCenterVoxInds = 1:length(MVSscalingTest.SLight.SLightIndsInMaskCell);
end
for IndSLCenterVox = SLCenterVoxInds
    CurrInds = MVSscalingTest.SLight.SLightIndsInMaskCell{IndSLCenterVox};
    DataOrg   = abs(MVSscaling.ScalingPerSubject(CurrInds,:));
    AmpDataOrg= MVSscaling.AmplitudesPerSubjectMRI(CurrInds,:,:);
    if(MVSscalingTest.StatisticsSettings.RemOutliers.Use==1) %mask outliers and check if center voxel is a outlier --> mask
        Data = DataOrg;
        if(any(isnan(Data))) %originally some NaNs contained
            MVSscalingTest.Results.DataQuality.SomeMissingData(IndSLCenterVox) = 1;
            MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 4; %mark original NaNs with 4
            if(isnan(Data(1))) %THE CENTER VOXEL!
                MVSscalingTest.Results.DataQuality.OutliersCenter(IndSLCenterVox)  = 5; %mark original NaNs of center voxel with 5
            end
            %test high threshold
            if(~isempty(MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi))
                if(any(DataOrg>=MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi))
                    MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 2; %mark STRONG Outliers with 2 if above high threshold
                    if(Data(1)>=MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi)
                        MVSscalingTest.Results.DataQuality.OutliersCenter(IndSLCenterVox)  = 2; %mark STRONG Outliers with 2 if above high threshold
                    end
                    Data(DataOrg>=MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi) = NaN;
                end
            end
            %test low threshold
            if(~isempty(MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo))
                if(any(DataOrg<=MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo))
                    if(MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox)==2) %also high outliers!
                        MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 3; %mark STRONG Outliers with 3 if under low and above high threshold are present
                    else
                        MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 1; %mark STRONG Outliers with 1 if under low threshold
                    end
                    if(Data(1)<=MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo)
                        MVSscalingTest.Results.DataQuality.OutliersCenter(IndSLCenterVox)  = 1; %mark STRONG Outliers with 1 if under low threshold
                    end
                    Data(DataOrg<=MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo) = NaN;
                end
            end
        else
            %no original NaNs but maybe strong outliers
            %test high threshold
            if(~isempty(MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi))
                if(any(DataOrg>=MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi))
                    MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 2; %mark STRONG Outliers with 2 if above high threshold
                    if(Data(1)>=MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi)
                        MVSscalingTest.Results.DataQuality.OutliersCenter(IndSLCenterVox)  = 2; %mark STRONG Outliers with 2 if above high threshold
                    end
                    Data(DataOrg>=MVSscalingTest.StatisticsSettings.RemOutliers.ThresHi) = NaN;
                end
            end
            %test low threshold
            if(~isempty(MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo))
                if(any(DataOrg<=MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo))
                    if(MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox)==2) %also high outliers!
                        MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 3; %mark STRONG Outliers with 3 if under low and above high threshold are present
                    else
                        MVSscalingTest.Results.DataQuality.OutliersGeneral(IndSLCenterVox) = 1; %mark STRONG Outliers with 1 if under low threshold
                    end
                    if(Data(1)<=MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo)
                        MVSscalingTest.Results.DataQuality.OutliersCenter(IndSLCenterVox)  = 1; %mark STRONG Outliers with 1 if under low threshold
                    end
                    Data(DataOrg<=MVSscalingTest.StatisticsSettings.RemOutliers.ThresLo) = NaN;
                end
            end
        end
    else
        Data = DataOrg; %just use data as it is
        if(any(isnan(Data))) %originally some NaNs contained
            MVSscalingTest.Results.DataQuality.SomeMissingData(IndSLCenterVox) = 1;
        end
    end
    MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,1:3)= quantile(                  Data(:), [.5 .25 .75]);
    MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,1:3)= quantile(squeeze(nanmedian(Data,2)),[.5 .25 .75]);
    MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,1:3)= quantile(squeeze(nanmedian(Data,1)),[.5 .25 .75]);
    
    %correlation with lambda for medianSLight, i.e. subject data in Searchlight
    %pearson corr==1
    [RHO1,PVal1] = corr(squeeze(nanmedian(Data,1))',ScalingFactorsLambda); %normal corr==1
    MVSscalingTest.Results.MedianSLight.Corr.Rho(IndSLCenterVox,1) = RHO1(1);
    MVSscalingTest.Results.MedianSLight.Corr.P(IndSLCenterVox,1)   = PVal1(1);
    
    %spearman rank corr==2
    [RHO2,PVal2] = corr(squeeze(nanmedian(Data,1))',ScalingFactorsLambda,'type','Spearman'); %rank corr==2
    MVSscalingTest.Results.MedianSLight.Corr.Rho(IndSLCenterVox,2) = RHO2(1);
    MVSscalingTest.Results.MedianSLight.Corr.P(IndSLCenterVox,2)   = PVal2(1);
    
    %% amplitude data corr with SPVs
    AmpData = squeeze(nanmedian(AmpDataOrg,1));
    %pearson corr==1
    [RHO1_AmpData,PVal1_AmpData] = corr(AmpData(:),SPVcorrData(:)); %normal corr==1
    MVSscalingTest.Results.MedianSLight.AmpCorr.Rho(IndSLCenterVox,1) = RHO1_AmpData(1);
    MVSscalingTest.Results.MedianSLight.AmpCorr.P(IndSLCenterVox,1)   = PVal1_AmpData(1);
    
    %spearman rank corr==2
    [RHO2_AmpData,PVal2_AmpData] = corr(AmpData(:),SPVcorrData(:),'type','Spearman'); %rank corr==2
    MVSscalingTest.Results.MedianSLight.AmpCorr.Rho(IndSLCenterVox,2) = RHO2_AmpData(1);
    MVSscalingTest.Results.MedianSLight.AmpCorr.P(IndSLCenterVox,2)   = PVal2_AmpData(1);
    
    %% abs amplitude data corr with SPVs
    AbsAmpData = squeeze(nanmedian(abs(AmpDataOrg),1));
    %pearson corr==1
    [RHO1_AbsAmpData,PVal1_AbsAmpData] = corr(AbsAmpData(:),SPVcorrData(:)); %normal corr==1
    MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Rho(IndSLCenterVox,1) = RHO1_AbsAmpData(1);
    MVSscalingTest.Results.MedianSLight.AbsAmpCorr.P(IndSLCenterVox,1)   = PVal1_AbsAmpData(1);
    
    %spearman rank corr==2
    [RHO2_AbsAmpData,PVal2_AbsAmpData] = corr(AbsAmpData(:),SPVcorrData(:),'type','Spearman'); %rank corr==2
    MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Rho(IndSLCenterVox,2) = RHO2_AbsAmpData(1);
    MVSscalingTest.Results.MedianSLight.AbsAmpCorr.P(IndSLCenterVox,2)   = PVal2_AbsAmpData(1);
    
    
    %% PERM-Test for correlations
    if(UsePermTest)
%         H_waitbarPTest = waitbar(0,['Perm-Test(',num2str(IndSLCenterVox),') ',num2str(0),'% done...']);
        NLargerPermTest1_1= 0;
        NLargerPermTest1_2= 0;
        NLargerPermTest2_1= 0;
        NLargerPermTest2_2= 0;
        MaxCovarPermAbsCorr = 0;
        
        NLargerSPVPermTest1_1= 0;
        NLargerSPVPermTest1_2= 0;
        NLargerSPVPermTest2_1= 0;
        NLargerSPVPermTest2_2= 0;
        NLargerSPVabsAmpDataPermTest1_1= 0;
        NLargerSPVabsAmpDataPermTest1_2= 0;
        NLargerSPVabsAmpDataPermTest2_1= 0;
        NLargerSPVabsAmpDataPermTest2_2= 0;
        MaxCovarPermAbsCorr_SPV = 0;
        for IndPermTest = 1:NPermTest
            %scaling data lambda
            IndicesPermutation = randperm(length(ScalingFactorsLambda)); %for the permutation of subjects
            R = corr(ScalingFactorsLambda,ScalingFactorsLambda(IndicesPermutation));
            MaxCovarPermAbsCorr = max([MaxCovarPermAbsCorr; abs(R(1))]); %maximum correlation between covariates and permutation of them --> test later
            RHOpermtest1 = corr(squeeze(nanmedian(Data,1))',ScalingFactorsLambda(IndicesPermutation)); %normal corr==1
            RHOpermtest2 = corr(squeeze(nanmedian(Data,1))',ScalingFactorsLambda(IndicesPermutation),'type','Spearman'); %rank corr==2
            %pearson
            if(abs(RHO1(1))<=abs(RHOpermtest1(1))) %ABS
                NLargerPermTest1_2 = NLargerPermTest1_2+1;
            end
            if(RHO1(1)<=0)
                if(RHO1(1)>RHOpermtest1(1))
                    NLargerPermTest1_1 = NLargerPermTest1_1+1;
                end
            else
                if(RHO1(1)<RHOpermtest1(1))
                    NLargerPermTest1_1 = NLargerPermTest1_1+1;
                end
            end
            %spearman
            if(abs(RHO2(1))<=abs(RHOpermtest2(1))) %ABS
                NLargerPermTest2_2 = NLargerPermTest2_2+1;
            end
            if(RHO2(1)<=0)
                if(RHO2(1)>RHOpermtest2(1))
                    NLargerPermTest2_1 = NLargerPermTest2_1+1;
                end
            else
                if(RHO2(1)<RHOpermtest2(1))
                    NLargerPermTest2_1 = NLargerPermTest2_1+1;
                end
            end
%             try
%                 H_waitbarPTest = waitbar(IndPermTest/NPermTest,H_waitbarPTest,['Perm-Test(',num2str(IndSLCenterVox),') ',num2str(IndPermTest*100/NPermTest),'% done...']);
%             catch
%                 H_waitbarPTest = waitbar(IndPermTest/NPermTest,['Perm-Test(',num2str(IndSLCenterVox),') ',num2str(IndPermTest*100/NPermTest),'% done...']);
%             end

            %SPV corr permtest
            IndicesPermutation_SPVs = randperm(length(SPVcorrData(:))); %for the permutation of subjects & mris --> make difference more pronounced because SPVs are quite similar
            R_SPVs = corr(SPVcorrData(:),SPVcorrData(IndicesPermutation_SPVs(:)));
            MaxCovarPermAbsCorr_SPV = max([MaxCovarPermAbsCorr_SPV; abs(R_SPVs(1))]); %maximum correlation between covariates and permutation of them --> test later
            
            %amplitude data
            RHO_AmpData_permtest1 = corr(AmpData(:),SPVcorrData(IndicesPermutation_SPVs(:))); %normal corr==1
            RHO_AmpData_permtest2 = corr(AmpData(:),SPVcorrData(IndicesPermutation_SPVs(:)),'type','Spearman'); %rank corr==2
            %pearson
            if(abs(RHO1_AmpData(1))<=abs(RHO_AmpData_permtest1(1))) %ABS
                NLargerSPVPermTest1_2 = NLargerSPVPermTest1_2+1;
            end
            if(RHO1_AmpData(1)<=0)
                if(RHO1_AmpData(1)>RHO_AmpData_permtest1(1))
                    NLargerSPVPermTest1_1 = NLargerSPVPermTest1_1+1;
                end
            else
                if(RHO1_AmpData(1)<RHO_AmpData_permtest1(1))
                    NLargerSPVPermTest1_1 = NLargerSPVPermTest1_1+1;
                end
            end
            %spearman
            if(abs(RHO2_AmpData(1))<=abs(RHO_AmpData_permtest2(1))) %ABS
                NLargerSPVPermTest2_2 = NLargerSPVPermTest2_2+1;
            end
            if(RHO2_AmpData(1)<=0)
                if(RHO2_AmpData(1)>RHO_AmpData_permtest2(1))
                    NLargerSPVPermTest2_1 = NLargerSPVPermTest2_1+1;
                end
            else
                if(RHO2_AmpData(1)<RHO_AmpData_permtest2(1))
                    NLargerSPVPermTest2_1 = NLargerSPVPermTest2_1+1;
                end
            end

            %ABS amplitude data
            RHO_AbsAmpData_permtest1 = corr(AbsAmpData(:),SPVcorrData(IndicesPermutation_SPVs(:))); %normal corr==1
            RHO_AbsAmpData_permtest2 = corr(AbsAmpData(:),SPVcorrData(IndicesPermutation_SPVs(:)),'type','Spearman'); %rank corr==2
            %pearson
            if(abs(RHO1_AbsAmpData(1))<=abs(RHO_AbsAmpData_permtest1(1))) %ABS
                NLargerSPVabsAmpDataPermTest1_2 = NLargerSPVabsAmpDataPermTest1_2+1;
            end
            if(RHO1_AbsAmpData(1)<=0)
                if(RHO1_AbsAmpData(1)>RHO_AbsAmpData_permtest1(1))
                    NLargerSPVabsAmpDataPermTest1_1 = NLargerSPVabsAmpDataPermTest1_1+1;
                end
            else
                if(RHO1_AbsAmpData(1)<RHO_AbsAmpData_permtest1(1))
                    NLargerSPVabsAmpDataPermTest1_1 = NLargerSPVabsAmpDataPermTest1_1+1;
                end
            end
            %spearman
            if(abs(RHO2_AbsAmpData(1))<=abs(RHO_AbsAmpData_permtest2(1))) %ABS
                NLargerSPVabsAmpDataPermTest2_2 = NLargerSPVabsAmpDataPermTest2_2+1;
            end
            if(RHO2_AbsAmpData(1)<=0)
                if(RHO2_AbsAmpData(1)>RHO_AbsAmpData_permtest2(1))
                    NLargerSPVabsAmpDataPermTest2_1 = NLargerSPVabsAmpDataPermTest2_1+1;
                end
            else
                if(RHO2_AbsAmpData(1)<RHO_AbsAmpData_permtest2(1))
                    NLargerSPVabsAmpDataPermTest2_1 = NLargerSPVabsAmpDataPermTest2_1+1;
                end
            end
        end
        %assign results
        %scaling data lambda
        MVSscalingTest.Results.MedianSLight.Corr.Pperm( IndSLCenterVox,1) = NLargerPermTest1_1/NPermTest;
        MVSscalingTest.Results.MedianSLight.Corr.Pperm2(IndSLCenterVox,1) = NLargerPermTest1_2/NPermTest;
        MVSscalingTest.Results.MedianSLight.Corr.Pperm( IndSLCenterVox,2) = NLargerPermTest2_1/NPermTest;
        MVSscalingTest.Results.MedianSLight.Corr.Pperm2(IndSLCenterVox,2) = NLargerPermTest2_2/NPermTest;
        %display and tests on the go
        if(NLargerPermTest1_1==0||NLargerPermTest1_2==0||NLargerPermTest2_1==0||NLargerPermTest2_2==0)
            disp(['SL ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),': LAMBDA-Correlation is highly significant! (',num2str(NLargerPermTest1_1),'&',num2str(NLargerPermTest1_2),' or ',num2str(NLargerPermTest2_1),'&',num2str(NLargerPermTest2_2),') [MaxCovarPermAbsCorr=',num2str(MaxCovarPermAbsCorr),']']);
        else
            if(NLargerPermTest1_1==1||NLargerPermTest1_2==1||NLargerPermTest2_1==1||NLargerPermTest2_2==1)
                disp(['SL ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),': LAMBDA-Correlation is significant! (',num2str(NLargerPermTest1_1),'&',num2str(NLargerPermTest1_2),' or ',num2str(NLargerPermTest2_1),'&',num2str(NLargerPermTest2_2),') [MaxCovarPermAbsCorr=',num2str(MaxCovarPermAbsCorr),']']);
            end
        end
        
        
        %SPV data corr with AmpData
        MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm( IndSLCenterVox,1) = NLargerSPVPermTest1_1/NPermTest;
        MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm2(IndSLCenterVox,1) = NLargerSPVPermTest1_2/NPermTest;
        MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm( IndSLCenterVox,2) = NLargerSPVPermTest2_1/NPermTest;
        MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm2(IndSLCenterVox,2) = NLargerSPVPermTest2_2/NPermTest;
        %display and tests on the go
        if(NLargerSPVPermTest1_1==0||NLargerSPVPermTest1_2==0||NLargerSPVPermTest2_1==0||NLargerSPVPermTest2_2==0)
            disp(['SL ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),': SPV-Correlation-AmpData is highly significant! (',num2str(NLargerSPVPermTest1_1),'&',num2str(NLargerSPVPermTest1_2),' or ',num2str(NLargerSPVPermTest2_1),'&',num2str(NLargerSPVPermTest2_2),') [MaxCovarPermAbsCorr_SPV=',num2str(MaxCovarPermAbsCorr_SPV),']']);
        else
            if(NLargerSPVPermTest1_1==1||NLargerSPVPermTest1_2==1||NLargerSPVPermTest2_1==1||NLargerSPVPermTest2_2==1)
                disp(['SL ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),': SPV-Correlation-AmpData is significant! (',num2str(NLargerSPVPermTest1_1),'&',num2str(NLargerSPVPermTest1_2),' or ',num2str(NLargerSPVPermTest2_1),'&',num2str(NLargerSPVPermTest2_2),') [MaxCovarPermAbsCorr_SPV=',num2str(MaxCovarPermAbsCorr_SPV),']']);
            end
        end
        %SPV data corr with AbsAmpData
        MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm( IndSLCenterVox,1) = NLargerSPVabsAmpDataPermTest1_1/NPermTest;
        MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm2(IndSLCenterVox,1) = NLargerSPVabsAmpDataPermTest1_2/NPermTest;
        MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm( IndSLCenterVox,2) = NLargerSPVabsAmpDataPermTest2_1/NPermTest;
        MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm2(IndSLCenterVox,2) = NLargerSPVabsAmpDataPermTest2_2/NPermTest;
        %display and tests on the go
        if(NLargerSPVabsAmpDataPermTest1_1==0||NLargerSPVabsAmpDataPermTest1_2==0||NLargerSPVabsAmpDataPermTest2_1==0||NLargerSPVabsAmpDataPermTest2_2==0)
            disp(['SL ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),': SPV-Correlation-AmpData is highly significant! (',num2str(NLargerSPVabsAmpDataPermTest1_1),'&',num2str(NLargerSPVabsAmpDataPermTest1_2),' or ',num2str(NLargerSPVabsAmpDataPermTest2_1),'&',num2str(NLargerSPVabsAmpDataPermTest2_2),') [MaxCovarPermAbsCorr_SPV=',num2str(MaxCovarPermAbsCorr_SPV),']']);
        else
            if(NLargerSPVabsAmpDataPermTest1_1==1||NLargerSPVabsAmpDataPermTest1_2==1||NLargerSPVabsAmpDataPermTest2_1==1||NLargerSPVabsAmpDataPermTest2_2==1)
                disp(['SL ',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),': SPV-Correlation-AmpData is significant! (',num2str(NLargerSPVabsAmpDataPermTest1_1),'&',num2str(NLargerSPVabsAmpDataPermTest1_2),' or ',num2str(NLargerSPVabsAmpDataPermTest2_1),'&',num2str(NLargerSPVabsAmpDataPermTest2_2),') [MaxCovarPermAbsCorr_SPV=',num2str(MaxCovarPermAbsCorr_SPV),']']);
            end
        end
    end
    
    
    %skewness
    MVSscalingTest.Results.Aggregate.SkewKurt(   IndSLCenterVox,1)    = skewness(                  Data(:), 0); %skewness bias corrected (needs more then 3 data values, but that should be alright.)
    MVSscalingTest.Results.MedianSubj.SkewKurt(  IndSLCenterVox,1)    = skewness(squeeze(nanmedian(Data,2)),0); %skewness bias corrected (needs more then 3 data values, but that should be alright.)
    MVSscalingTest.Results.MedianSLight.SkewKurt(IndSLCenterVox,1)    = skewness(squeeze(nanmedian(Data,1)),0); %skewness bias corrected (needs more then 3 data values, but that should be alright.)
    %kurtosis
    MVSscalingTest.Results.Aggregate.SkewKurt(   IndSLCenterVox,2)    = kurtosis(                  Data(:), 0); %kurtosis bias corrected (needs more then 3 data values, but that should be alright.)
    MVSscalingTest.Results.MedianSubj.SkewKurt(  IndSLCenterVox,2)    = kurtosis(squeeze(nanmedian(Data,2)),0); %kurtosis bias corrected (needs more then 3 data values, but that should be alright.)
    MVSscalingTest.Results.MedianSLight.SkewKurt(IndSLCenterVox,2)    = kurtosis(squeeze(nanmedian(Data,1)),0); %kurtosis bias corrected (needs more then 3 data values, but that should be alright.)
    %CI(1/2) = median -/+ 1.57*(Qrt3-Qrt1)/sqrt(length(InputData)); %--> Qrt3-Qrt1 == iQR
    MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,4)  = MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,1)-1.57*(MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,3)-MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,2))/sqrt(length(Data(:)));
    MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,5)  = MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,1)+1.57*(MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,3)-MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,2))/sqrt(length(Data(:)));
    MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,4)  = MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,1)-1.57*(MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,3)-MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,2))/sqrt(length(squeeze(median(Data,2))));
    MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,5)  = MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,1)+1.57*(MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,3)-MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,2))/sqrt(length(squeeze(median(Data,2))));
    MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)  = MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,1)-1.57*(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,3)-MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,2))/sqrt(length(squeeze(median(Data,1))));
    MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)  = MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,1)+1.57*(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,3)-MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,2))/sqrt(length(squeeze(median(Data,1))));
    try
        MVSscalingTest.Results.Type = 'NewSignRankFun';
        % check that some data is left
        if(all(isnan(Data(:)))||all(isnan(nanmedian(Data,2)))||all(isnan(nanmedian(Data,1))))
            MVSscalingTest.Results.DataQuality.UsableData(IndSLCenterVox) = 0; %note this in the mask for all possibilities
            if(all(isnan(Data(:)))) %the complete searchlight
                %% signrank-test
                MVSscalingTest.Results.Aggregate.p(IndSLCenterVox) = NaN;
                MVSscalingTest.Results.Aggregate.h(IndSLCenterVox) = NaN;
                statsAggregate.signedrank = NaN;
                statsAggregate.zval       = NaN;
                
                MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox) = NaN;
                MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox) = NaN;
                statsMedianSubj.signedrank = NaN;
                statsMedianSubj.zval       = NaN;
                
                MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox) = NaN;
                MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox) = NaN;
                statsMedianSLight.signedrank = NaN;
                statsMedianSLight.zval       = NaN;
                
                %% sign-test
                MVSscalingTest.Results.Aggregate.p_sign(IndSLCenterVox) = NaN;
                MVSscalingTest.Results.Aggregate.h_sign(IndSLCenterVox) = NaN;
                statsAggregate_sign.sign = NaN;
                statsAggregate_sign.zval = NaN;
                
                MVSscalingTest.Results.MedianSubj.p_sign(IndSLCenterVox) = NaN;
                MVSscalingTest.Results.MedianSubj.h_sign(IndSLCenterVox) = NaN;
                statsMedianSubj_sign.sign = NaN;
                statsMedianSubj_sign.zval = NaN;
                
                MVSscalingTest.Results.MedianSLight.p_sign(IndSLCenterVox) = NaN;
                MVSscalingTest.Results.MedianSLight.h_sign(IndSLCenterVox) = NaN;
                statsMedianSLight_sign.sign = NaN;
                statsMedianSLight_sign.zval = NaN;
            else %maybe only median subject or median searchlight
                %% signrank-test & %% sign-test
                [MVSscalingTest.Results.Aggregate.p(IndSLCenterVox),        MVSscalingTest.Results.Aggregate.h(IndSLCenterVox),        statsAggregate]     = signrank(Data(:), MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
                [MVSscalingTest.Results.Aggregate.p_sign(IndSLCenterVox),   MVSscalingTest.Results.Aggregate.h_sign(IndSLCenterVox),   statsAggregate_sign]= signtest(Data(:), MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
                if(all(isnan(nanmedian(Data,2))))
                    %% signrank-test
                    MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox) = NaN;
                    MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox) = NaN;
                    statsMedianSubj.signedrank = NaN;
                    statsMedianSubj.zval       = NaN;
                    %% sign-test
                    MVSscalingTest.Results.MedianSubj.p_sign(IndSLCenterVox) = NaN;
                    MVSscalingTest.Results.MedianSubj.h_sign(IndSLCenterVox) = NaN;
                    statsMedianSubj_sign.sign = NaN;
                    statsMedianSubj_sign.zval = NaN;
                else
                    %% signrank-test & %% sign-test
                    [MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox),       MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox),       statsMedianSubj]     = signrank(squeeze(nanmedian(Data,2)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
                    [MVSscalingTest.Results.MedianSubj.p_sign(IndSLCenterVox),  MVSscalingTest.Results.MedianSubj.h_sign(IndSLCenterVox),  statsMedianSubj_sign]= signtest(squeeze(nanmedian(Data,2)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
                    if(all(isnan(nanmedian(Data,1))))
                        %% signrank-test
                        MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox) = NaN;
                        MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox) = NaN;
                        statsMedianSLight.signedrank = NaN;
                        statsMedianSLight.zval       = NaN;
                        %% sign-test
                        MVSscalingTest.Results.MedianSLight.p_sign(IndSLCenterVox) = NaN;
                        MVSscalingTest.Results.MedianSLight.h_sign(IndSLCenterVox) = NaN;
                        statsMedianSLight_sign.sign = NaN;
                        statsMedianSLight_sign.zval = NaN;
                    else
                        %% signrank-test & %% sign-test
                        [MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox),     MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox),     statsMedianSLight]     = signrank(squeeze(nanmedian(Data,1)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
                        [MVSscalingTest.Results.MedianSLight.p_sign(IndSLCenterVox),MVSscalingTest.Results.MedianSLight.h_sign(IndSLCenterVox),statsMedianSLight_sign]= signtest(squeeze(nanmedian(Data,1)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
                    end
                end
            end
        else
            %% signrank-test
            [MVSscalingTest.Results.Aggregate.p(IndSLCenterVox),   MVSscalingTest.Results.Aggregate.h(IndSLCenterVox),   statsAggregate]    = signrank(                  Data(:), MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
            [MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox),  MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox),  statsMedianSubj]   = signrank(squeeze(nanmedian(Data,2)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
            [MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox),MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox),statsMedianSLight] = signrank(squeeze(nanmedian(Data,1)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
            
            %% sign-test
            [MVSscalingTest.Results.Aggregate.p_sign(IndSLCenterVox),   MVSscalingTest.Results.Aggregate.h_sign(IndSLCenterVox),   statsAggregate_sign]    = signtest(                  Data(:), MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
            [MVSscalingTest.Results.MedianSubj.p_sign(IndSLCenterVox),  MVSscalingTest.Results.MedianSubj.h_sign(IndSLCenterVox),  statsMedianSubj_sign]   = signtest(squeeze(nanmedian(Data,2)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
            [MVSscalingTest.Results.MedianSLight.p_sign(IndSLCenterVox),MVSscalingTest.Results.MedianSLight.h_sign(IndSLCenterVox),statsMedianSLight_sign] = signtest(squeeze(nanmedian(Data,1)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
        end
        
        %% signrank-test
        MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox) = statsAggregate.signedrank;
        if(isfield(statsAggregate,'zval'))
            MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox)   = statsAggregate.zval;
        end
        MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox) = statsMedianSubj.signedrank;
        if(isfield(statsMedianSubj,'zval'))
            MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox)   = statsMedianSubj.zval;
        end
        MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox) = statsMedianSLight.signedrank;
        if(isfield(statsMedianSLight,'zval'))
            MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox)   = statsMedianSLight.zval;
        end
        
        %% sign-test (NOT RANKSIGN, but SIGN!)
        MVSscalingTest.Results.Aggregate.sign(IndSLCenterVox) = statsAggregate_sign.sign;
        if(isfield(statsAggregate_sign,'zval'))
            if(~isnan(statsAggregate_sign.zval))
                MVSscalingTest.Results.Aggregate.zval_sign(IndSLCenterVox)   = statsAggregate_sign.zval;
            else
                MVSscalingTest.Results.Aggregate.zval_sign(IndSLCenterVox)   = Inf;
            end
        end
        MVSscalingTest.Results.MedianSubj.sign(IndSLCenterVox) = statsMedianSubj_sign.sign;
        if(isfield(statsMedianSubj_sign,'zval'))
            if(~isnan(statsMedianSubj_sign.zval))
                MVSscalingTest.Results.MedianSubj.zval_sign(IndSLCenterVox)   = statsMedianSubj_sign.zval;
            else
                MVSscalingTest.Results.MedianSubj.zval_sign(IndSLCenterVox)   = Inf;
            end
        end
        MVSscalingTest.Results.MedianSLight.sign(IndSLCenterVox) = statsMedianSLight_sign.sign;
        if(isfield(statsMedianSLight_sign,'zval'))
            if(~isnan(statsMedianSLight_sign.zval))
                MVSscalingTest.Results.MedianSLight.zval_sign(IndSLCenterVox)   = statsMedianSLight_sign.zval;
            else
                MVSscalingTest.Results.MedianSLight.zval_sign(IndSLCenterVox)   = Inf;
            end
        end
        switch(ChoiceTest)
            case 'Higher'
                %statsAggregate
                if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H0
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                else
                    if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H1
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when below because of the "Higher" test-case
                    else %accept H0
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                    end
                end
                %statsMedianSubj
                if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                else
                    if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when below because of the "Higher" test-case
                    else %accept H0
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                    end
                end
                %statsMedianSLight
                if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                else
                    if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when below because of the "Higher" test-case
                    else %accept H0
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                    end
                end         
            case 'Lower'
                %statsAggregate
                if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H0
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                else
                    if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H1
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when above because of the "Lower" test-case
                    else
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                    end
                end
                %statsMedianSubj
                if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                else
                    if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 1; %in this case we consider the overlap also when above because of the "Lower" test-case
                    else %accept H0
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                    end
                end
                %statsMedianSLight
                if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                else
                    if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when above because of the "Lower" test-case
                    else %accept H0
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                    end
                end
            otherwise
                %statsAggregate
                if((MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian)||(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian))
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %Not overlapping!
                else
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 1; %Overlapping!
                end
                %statsMedianSubj
                if((MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian)||(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian))
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %Not overlapping!
                else
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 1; %Overlapping!
                end
                %statsMedianSLight
                if((MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian)||(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian))
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %Not overlapping!
                else
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 1; %Overlapping!
                end
        end
    catch CATCH_stats
        MVSscalingTest.Results.Type = 'OrgSignRankFun+MyChanges';
        [MVSscalingTest.Results.Aggregate.p(IndSLCenterVox),   MVSscalingTest.Results.Aggregate.h(IndSLCenterVox),   statsAggregate]    = signrank(                  Data(:), ExpectedMedian);
        [MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox),  MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox),  statsMedianSubj]   = signrank(squeeze(nanmedian(Data,2)),ExpectedMedian);
        [MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox),MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox),statsMedianSLight] = signrank(squeeze(nanmedian(Data,1)),ExpectedMedian);
        switch(ChoiceTest)
            case 'Higher'
                %statsAggregate
                if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H0
                    MVSscalingTest.Results.Aggregate.p(IndSLCenterVox)         = 1;
                    MVSscalingTest.Results.Aggregate.h(IndSLCenterVox)         = 0;
                    MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox)= Inf;
                    MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox)      = 0;
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                else
                    if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H1
                        MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox) = statsAggregate.signedrank;
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when below because of the "Higher" test-case
                        if(isfield(statsAggregate,'zval'))
                            MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox) = MVSscalingTest.Results.Aggregate.h(IndSLCenterVox).*abs(statsAggregate.zval);
                        end
                    else %accept H0
                        MVSscalingTest.Results.Aggregate.p(IndSLCenterVox)         = 1;
                        MVSscalingTest.Results.Aggregate.h(IndSLCenterVox)         = 0;
                        MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox)= Inf;
                        MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox)      = 0;
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                    end
                end
                %statsMedianSubj
                if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox)         = 1;
                    MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox)         = 0;
                    MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox)= Inf;
                    MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox)      = 0;
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                else
                    if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox) = statsMedianSubj.signedrank;
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when below because of the "Higher" test-case
                        if(isfield(statsMedianSubj,'zval'))
                            MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox) = MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox).*abs(statsMedianSubj.zval);
                        end
                    else %accept H0
                        MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox)         = 1;
                        MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox)         = 0;
                        MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox)= Inf;
                        MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox)      = 0;
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                    end
                end
                %statsMedianSLight
                if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox)         = 1;
                    MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox)         = 0;
                    MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox)= Inf;
                    MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox)      = 0;
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                else
                    if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox) = statsMedianSLight.signedrank;
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when below because of the "Higher" test-case
                        if(isfield(statsMedianSLight,'zval'))
                            MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox) = MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox).*abs(statsMedianSLight.zval);
                        end
                    else %accept H0
                        MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox)         = 1;
                        MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox)         = 0;
                        MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox)= Inf;
                        MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox)      = 0;
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when below because of the "Higher" test-case
                    end
                end
            case 'Lower'
                %statsAggregate
                if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H0
                    MVSscalingTest.Results.Aggregate.p(IndSLCenterVox)         = 1;
                    MVSscalingTest.Results.Aggregate.h(IndSLCenterVox)         = 0;
                    MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox)= Inf;
                    MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox)      = 0;
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                else
                    if(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H1
                        MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox) = statsAggregate.signedrank;
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when above because of the "Lower" test-case
                        if(isfield(statsAggregate,'zval'))
                            MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox) = MVSscalingTest.Aggregate.h(IndSLCenterVox).*abs(statsAggregate.zval);
                        end
                    else
                        MVSscalingTest.Results.Aggregate.p(IndSLCenterVox)         = 1;
                        MVSscalingTest.Results.Aggregate.h(IndSLCenterVox)         = 0;
                        MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox)= Inf;
                        MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox)      = 0;
                        MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                    end
                end
                %statsMedianSubj
                if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox)         = 1;
                    MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox)         = 0;
                    MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox)= Inf;
                    MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox)      = 0;
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                else
                    if(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox) = statsMedianSubj.signedrank;
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 1; %in this case we consider the overlap also when above because of the "Lower" test-case
                        if(isfield(statsMedianSubj,'zval'))
                            MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox) = MVSscalingTest.MedianSubj.h(IndSLCenterVox).*abs(statsMedianSubj.zval);
                        end
                    else %accept H0
                        MVSscalingTest.Results.MedianSubj.p(IndSLCenterVox)         = 1;
                        MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox)         = 0;
                        MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox)= Inf;
                        MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox)      = 0;
                        MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                    end
                end
                %statsMedianSLight
                if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian) %accept H0
                    MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox)         = 1;
                    MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox)         = 0;
                    MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox)= Inf;
                    MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox)      = 0;
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                else
                    if(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian) %accept H1
                        MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox) = statsMedianSLight.signedrank;
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox)= 1; %in this case we consider the overlap also when above because of the "Lower" test-case
                        if(isfield(statsMedianSLight,'zval'))
                            MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox) = MVSscalingTest.MedianSLight.h(IndSLCenterVox).*abs(statsMedianSLight.zval);
                        end
                    else %accept H0
                        MVSscalingTest.Results.MedianSLight.p(IndSLCenterVox)         = 1;
                        MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox)         = 0;
                        MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox)= Inf;
                        MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox)      = 0;
                        MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %in this case we consider the overlap also when above because of the "Lower" test-case
                    end
                end
            otherwise
                %statsAggregate
                MVSscalingTest.Results.Aggregate.signedrank(IndSLCenterVox)= statsAggregate.signedrank;
                if(isfield(statsAggregate,'zval'))
                    MVSscalingTest.Results.Aggregate.zval(IndSLCenterVox)  = MVSscalingTest.Results.Aggregate.h(IndSLCenterVox).*abs(statsAggregate.zval);
                end
                if((MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian)||(MVSscalingTest.Results.Aggregate.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian))
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 0; %Not overlapping!
                else
                    MVSscalingTest.Results.Aggregate.CIoverlap(IndSLCenterVox) = 1; %Overlapping!
                end
                %statsMedianSubj
                MVSscalingTest.Results.MedianSubj.signedrank(IndSLCenterVox)= statsMedianSubj.signedrank;
                if(isfield(statsMedianSubj,'zval'))
                    MVSscalingTest.Results.MedianSubj.zval(IndSLCenterVox)  = MVSscalingTest.Results.MedianSubj.h(IndSLCenterVox).*abs(statsMedianSubj.zval);
                end
                if((MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian)||(MVSscalingTest.Results.MedianSubj.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian))
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 0; %Not overlapping!
                else
                    MVSscalingTest.Results.MedianSubj.CIoverlap(IndSLCenterVox) = 1; %Overlapping!
                end
                %statsMedianSLight
                MVSscalingTest.Results.MedianSLight.signedrank(IndSLCenterVox)= statsMedianSLight.signedrank;
                if(isfield(statsMedianSLight,'zval'))
                    MVSscalingTest.Results.MedianSLight.zval(IndSLCenterVox)  = MVSscalingTest.Results.MedianSLight.h(IndSLCenterVox).*abs(statsMedianSLight.zval);
                end
                if((MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)>ExpectedMedian)||(MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)<ExpectedMedian))
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 0; %Not overlapping!
                else
                    MVSscalingTest.Results.MedianSLight.CIoverlap(IndSLCenterVox) = 1; %Overlapping!
                end
        end
    end
    if(~UsePermTest)
        try
            H_waitbar = waitbar(IndSLCenterVox/length(MVSscalingTest.SLight.SLightIndsInMaskCell),H_waitbar,['SearchLight-test centered at every voxel ',num2str(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell)),'% done...']);
        catch
            H_waitbar = waitbar(IndSLCenterVox/length(MVSscalingTest.SLight.SLightIndsInMaskCell),['SearchLight-test centered at every voxel ',num2str(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell)),'% done...']);
        end
    end
    if(mod(IndSLCenterVox,1000)==0) %save a copy
        disp(['autosave: "',OutDir,'".']);
        save([OutDir,filesep,'autosave.mat'],'*');
    end
end
if(~UsePermTest)
    try
        H_waitbar = waitbar(IndSLCenterVox/length(MVSscalingTest.SLight.SLightIndsInMaskCell),H_waitbar,'SearchLight-test centered at every voxel is finished!');
    catch
        H_waitbar = waitbar(IndSLCenterVox/length(MVSscalingTest.SLight.SLightIndsInMaskCell),'SearchLight-test centered at every voxel is finished!');
    end
    uiwait(H_waitbar,2);
    close(H_waitbar);
end

%% save the results
[basedir,fnameMask]   = fileparts(MVSscalingTest.SLight.V_SLmask.fname);
[basedir,fnameWBMask] = fileparts(MVSscalingTest.SLight.V_WholeBrainMask.fname);
if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
    try
        if(UseSLsmooth)
            save([OutDir,filesep,'sResultsSL',num2str(MVSscalingTest.SLight.NHood),'_Crop',regexprep(fnameMask,' ','_'),'_W_',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        else
            save([OutDir,filesep,'ResultsSL',num2str(MVSscalingTest.SLight.NHood),'_Crop', regexprep(fnameMask,' ','_'),'_W_',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        end
    catch
        CurrDir = pwd;
        cd(OutDir); 
        if(UseSLsmooth)
            save(['sResultsSL',num2str(MVSscalingTest.SLight.NHood),'_Crop',               regexprep(fnameMask,' ','_'),'_W_',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        else
            save(['ResultsSL',num2str(MVSscalingTest.SLight.NHood),'_Crop',                regexprep(fnameMask,' ','_'),'_W_',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        end
        cd(CurrDir);
    end
else
    try
        if(UseSLsmooth)
            save([OutDir,filesep,'sResultsSL',num2str(MVSscalingTest.SLight.NHood),'_All', regexprep(fnameMask,' ','_'),'_In',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        else
            save([OutDir,filesep,'ResultsSL',num2str(MVSscalingTest.SLight.NHood),'_All',  regexprep(fnameMask,' ','_'),'_In',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        end
    catch
        CurrDir = pwd;
        cd(OutDir);
        if(UseSLsmooth)
            save(['sResultsSL',num2str(MVSscalingTest.SLight.NHood),'_All',                regexprep(fnameMask,' ','_'),'_In',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        else
            save(['ResultsSL',num2str(MVSscalingTest.SLight.NHood),'_All',                 regexprep(fnameMask,' ','_'),'_In',regexprep(fnameWBMask,' ','_'),'.mat'],'MVSscalingTest');
        end
        cd(CurrDir);
    end
end

%% clean up autosave
try
    delete([OutDir,filesep,'autosave.mat']);
end

%% inform user 
h = helpdlg({['Done with gathering stats (using ',MVSscalingTest.Results.Type,').']; 'Next will be writing out to NIFTI.'},'Done.');
uiwait(h,1);
try
    close(h);
end

%% write out the statistics in the way that is appropriate.
[V_out]=WriteNII_SLResults(MVSscalingTest,OutDir);
% %% write out a mask
% V_out = spm_vol(MVSscaling.Masks.MPaths{1});
% if(V_out.dt(1)<16)
%     V_out.dt(1) = 16; %not necessary but save
% end
% V_out.fname = [OutDir,filesep,'ResultsScalingTestSEARCHLIGHT_MaskSignif_Median',ChoiceTest,regexprep(answer_SetupStats{1},'*',''),'.nii'];
% Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
% Y(MVSscaling.Masks.WholeBrainRaw~=0) = MVSscalingTest.h;
% Y = reshape(Y,V_out.dim);
%
% V_out = spm_write_vol(V_out,Y);
%
% %% save inverse of mask?
% if(strcmp('Yes',questdlg('Save INVERSE of mask?','Save InverseMask?','Yes','No','No')))
%     V_out.fname = [OutDir,filesep,'ResultsScalingTestSEARCHLIGHT_INVERSEMaskSignif_Median',ChoiceTest,regexprep(answer_SetupStats{1},'*',''),'.nii'];
%     Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
%     Y(MVSscaling.Masks.WholeBrainRaw~=0) = ~MVSscalingTest.h;
%     Y = reshape(Y,V_out.dim);
%
%     V_out = spm_write_vol(V_out,Y);
% end
%
% %% write out zvals
% V_out.fname = [OutDir,filesep,'ResultsScalingTestSEARCHLIGHT_zVals_Median',ChoiceTest,regexprep(answer_SetupStats{1},'*',''),'.nii'];
% Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
% Y(MVSscaling.Masks.WholeBrainRaw~=0) = MVSscalingTest.zval;
% Y = reshape(Y,V_out.dim);
%
% V_out = spm_write_vol(V_out,Y);


%% Done.
disp(' ');
disp('Done.');