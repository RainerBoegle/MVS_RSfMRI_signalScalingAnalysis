function MVSscalingTest = RunMVScorrPermTest(MVSscalingTest,BackupSavePath)
% This function does the correlation with the covariates and permutation test for significance.
% The input data for this test is the median searchlight of the scaling data or the amplitude data,
% as correlations need to be done over subjects.
%
%Usage:
%       MVSscalingTest = RunMVScorrPermTest(MVSscalingTest,BackupSavePath);
%      
%
%V1.1
%Date: V1.1(04.09.2015): update to process data prepared in two ways (Voxel-Wise scaling values & Searchlight-averaged amplitudes turned into scaling data), both analysed with searchlight here as before. V1.0(26.7.2015) (initial implementation based on test script for analysis of scaling data.)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% settings
DisplayPermTestCorrsOfCovars = 0; %diplay boxplots of covariates correlations with permuted covariates --> test if covars are sufficiently random to be viable for permtest and correlation test.
PercentSaveStep              = 5; %save every 5% of searchlight progress

%% init IC message string
if(isfield(MVSscalingTest,'ICnum'))
    if(isnumeric(MVSscalingTest.ICnum))
        ICinfoStr = ['IC ',num2str(MVSscalingTest.ICnum),': '];
    else
        ICinfoStr =  'IC <unknown>: ';
    end
else
    warning('MATLAB:struct:missing','Field "ICnum" missing in structure "MVSscalingTest"');
    ICinfoStr = 'IC <???>: ';
end

%% check also NHood of searchlight for degeneracy ie NHood == 0
if(MVSscalingTest.SLight.NHood==0)
    ICinfoStr = [ICinfoStr,'[NHood==0] '];
end

%% reminder
% PermTest.CovarsMatFilePath = spm_select(1,'mat','Select Covars.mat file...');
% PermTest.CovarsStruct      = GetCovarsViaSubjNrs(PermTest.CovarsMatFilePath,FileListInfo);
% PermTest.CovarsStruct.Regressors
% PermTest.CovarsStruct.CovarName
%
% PermTest.NPermTest = eval(answer_perm{1});
% PermTest.NFailSkip = eval(answer_perm{2});
% PermTest.PermMatrix = zeros(FileListInfo.NSubjs,PermTest.NPermTest); %the permutations for the covariates
% PermTest.PermMatrix_Rho_pVals_InfoStr = {'Dimension is (NPermTest, MCovar, Rho==1/P==2, Pearson==1/Spearman==2).'; 'Use "InspectPermMatrix_Rho_pVals.m" to see the results of a correlation of the original covariates with their PERMUTED counterparts.'};
% PermTest.PermMatrix_Rho_pVals         = []; %init later (NPermTest, MCovar, Rho==1/P==2, Pearson==1/Spearman==2). Use "InspectPermMatrix_Rho_pVals.m" to see the results.
% 
% PermTest.Results.DimInfoStr = {'All have dimensions (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2).'; '[Except ".CurrVoxInd" which just keeps track if the current voxel in the analysis, such that continuation is possible.]'};
% PermTest.Results.Rho        = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %correlation at voxel for UNPERMUTED covariates
% PermTest.Results.pCorr      = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %p value from corr-function of correlation at voxel for UNPERMUTED covariates
% PermTest.Results.NFail      = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %the number of failed tests
% PermTest.Results.NTests     = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %the number of test that were actually done.
% PermTest.Results.pPermTest  = []; %init later (NVoxel,MCovar,2==Pearson&SpearmanCorr,3==Lambda&Amplitude1&Amplitude2) %NFail./NTests %assume equal number of failed tests until full number test, i.e rate stays the same. This should not be a big problem if NFail is >=5, because then (roughly) we would have all rejected tests anyways.
% PermTest.Results.CurrVoxInd = []; %(NVoxel,1) current voxel that we are working on, if empty then we start at the first, if any number then we continue at that one.

%% permtest even set up?
if(~isfield(MVSscalingTest,'PermTest')) %no permtest set up?!!! --> quit!
    warning('MATLAB:missingfield','No permutation test was set up?! (field is missing) Bailing out without any processing...');
    return;
else
    if(isempty(MVSscalingTest.PermTest)) %no permtest set up?!!! --> quit!
        warning('MATLAB:emptyfield','No permutation test was set up! Bailing out without any processing...');
        return;
    else
        PermTest = MVSscalingTest.PermTest;
    end
end

%% estimate permutation matrix correlations with original covariates and plot them, such that user can get an idea if covariates are appropriate for permutation test
if(isempty(PermTest.PermMatrix_Rho_pVals))
    disp('Before permutation Test: ');
    disp('Checking sufficient difference in covariates under permutation of subject...');
    reverseStr = '';
    PermTest.PermMatrix_Rho_pVals = zeros(PermTest.NPermTest,size(PermTest.CovarsStruct.Regressors,2),2,2); %the correlations & corresponding p-values (Pearson/Spearman) for the original UNPERMUTED covariates with their permuted counterparts --> use this to inspect if permutations are too similar to original data, i.e. data is not useful for permutation test. Use the function "InspectPermMatrix_Rho_pVals.m".
    for IndPermTest = 1:PermTest.NPermTest
        %disp(['PermMatrix_Rho_pVals: ',num2str(floor(IndPermTest*1000/PermTest.NPermTest)/10),'% done.']);
        msg = sprintf('PermMatrix_Rho_pVals percent done: %3.1f ...', floor(IndPermTest*1000/PermTest.NPermTest)/10); %Don't forget this semicolon
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        [Rho,P] = corr(PermTest.CovarsStruct.Regressors,PermTest.CovarsStruct.Regressors(PermTest.PermMatrix(:,IndPermTest),:));
        PermTest.PermMatrix_Rho_pVals(IndPermTest,:,1,1) = diag(Rho);
        PermTest.PermMatrix_Rho_pVals(IndPermTest,:,2,1) = diag(P);
        [Rho,P] = corr(PermTest.CovarsStruct.Regressors,PermTest.CovarsStruct.Regressors(PermTest.PermMatrix(:,IndPermTest),:),'type','Spearman');
        PermTest.PermMatrix_Rho_pVals(IndPermTest,:,1,2) = diag(Rho);
        PermTest.PermMatrix_Rho_pVals(IndPermTest,:,2,2) = diag(P);
        clear Rho P
    end
    MVSscalingTest.PermTest = PermTest; %assign back
    fprintf('Done.\n');
end
if(DisplayPermTestCorrsOfCovars)
    [H,ax1,ax2,ax3,ax4] = InspectPermMatrix_Rho_pVals(MVSscalingTest);
end

%% init to prepare for processing 
ApproachesInfo = MVSscalingTest.ApproachesInfo;
if(isempty(PermTest.Results.CurrVoxInd)) %initial without having been run EVER before.
    PermTest.Results.CurrVoxInd     = zeros(size(ApproachesInfo,1),1); %(NVoxel,1) current voxel that we are working on, if empty then we start at the first, if any number then we continue at that one.
    PermTest.DataQuality.UsableData =  cell(size(ApproachesInfo,1),1); %if median searchlight ONLY CONTAINS NaNs then this in not usable data. --> keep track of this. (but assume in the beginning that everything is usable)
    PermTest.Results.Rho            =  cell(size(ApproachesInfo,1),1); %correlation at voxel for UNPERMUTED covariates
    PermTest.Results.pCorr          =  cell(size(ApproachesInfo,1),1); %p value from corr-function of correlation at voxel for UNPERMUTED covariates
    PermTest.Results.NFail          =  cell(size(ApproachesInfo,1),1); %the number of failed tests
    PermTest.Results.NTests         =  cell(size(ApproachesInfo,1),1); %the number of test that were actually done.
    PermTest.Results.pPermTest      =  cell(size(ApproachesInfo,1),1); %NFail./NTests %assume equal number of failed tests until full number test, i.e rate stays the same. This should not be a big problem if NFail is >=5, because then (roughly) we would have all rejected tests anyways.
end %no else part needed at this point

%% do the processing & loop over approaches
disp('------------------------------------------------------------------');
disp(ICinfoStr);
for IndApproach = 1:size(ApproachesInfo,1)
    disp(['Data is from approach: ',ApproachesInfo{IndApproach,1},': ',ApproachesInfo{IndApproach,2}]);
    %% init fields
    if(PermTest.Results.CurrVoxInd(IndApproach)==0) %for this APPROACH we did not do the tests yet.
        PermTest.DataQuality.UsableData{IndApproach} = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1);  %if median searchlight ONLY CONTAINS NaNs then this in not usable data. --> keep track of this. (but assume in the beginning that everything is usable)
        PermTest.Results.Rho{IndApproach}            = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),length(PermTest.CovarsStruct.CovarName),2,3); %correlation at voxel for UNPERMUTED covariates
        PermTest.Results.pCorr{IndApproach}          = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),length(PermTest.CovarsStruct.CovarName),2,3); %p value from corr-function of correlation at voxel for UNPERMUTED covariates
        PermTest.Results.NFail{IndApproach}          = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),length(PermTest.CovarsStruct.CovarName),2,3); %the number of failed tests
        PermTest.Results.NTests{IndApproach}         = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),length(PermTest.CovarsStruct.CovarName),2,3); %the number of test that were actually done.
        PermTest.Results.pPermTest{IndApproach}      = zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),length(PermTest.CovarsStruct.CovarName),2,3); %NFail./NTests %assume equal number of failed tests until full number test, i.e rate stays the same. This should not be a big problem if NFail is >=5, because then (roughly) we would have all rejected tests anyways.
        PermTest.Results.CurrVoxInd(IndApproach)     = 1; %(NVoxel,1) current voxel that we are working on, if empty then we start at the first, if any number then we continue at that one.
    else
        %continuing from PermTest.Results.CurrVoxInd(IndApproach) %NB: assuming that everything is initialized already
        disp(['Continuing correlation & permutation test from voxel #',num2str(PermTest.Results.CurrVoxInd(IndApproach)),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),'...']);
        %NB: assuming that everything is initialized already
    end
    
    %% do basic correlation for all voxels & permutation test for all searchlights
    Regressors = PermTest.CovarsStruct.Regressors; %Covariates for correlation.
    
    NextPercentDone = floor(PermTest.Results.CurrVoxInd(IndApproach)*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell)); %init %for indicating the progress 1%-stepwise
    PrevPercentSave = NextPercentDone-mod(NextPercentDone,PercentSaveStep);
    reverseStr = '';
    for IndSLCenterVox = PermTest.Results.CurrVoxInd(IndApproach):length(MVSscalingTest.SLight.SLightIndsInMaskCell)
        if(NextPercentDone==0)%start
            msg = sprintf('PermTests of Lambda & Amplitude distribution done for %03.0f percent of searchlights.',floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            NextPercentDone = 1;
            PrevPercentSave = 0;
        else
            if(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))==NextPercentDone)
                msg = sprintf('PermTests of Lambda & Amplitude distribution done for %03.0f percent of searchlights.',floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                NextPercentDone = NextPercentDone+1;
            end
        end
        
        %% get Data, i.e. the median searchlight (from aggregate data using nanmedian) and combine for correlation
        CurrInds   = MVSscalingTest.SLight.SLightIndsInMaskCell{IndSLCenterVox};
        if(MVSscalingTest.SLight.NHood>0)
            DataLambda = nanmedian(          abs(MVSscalingTest.ScalingPerSubject{IndApproach}(CurrInds,:))  ,1); %median searchlight
            DataAmp1   = nanmedian(squeeze(MVSscalingTest.AmplitudesPerSubjectMRI{IndApproach}(CurrInds,:,1)),1); %median searchlight
            DataAmp2   = nanmedian(squeeze(MVSscalingTest.AmplitudesPerSubjectMRI{IndApproach}(CurrInds,:,2)),1); %median searchlight
        else
            DataLambda = squeeze(  abs(MVSscalingTest.ScalingPerSubject{IndApproach}(CurrInds,:)));  %for single voxel, i.e. NHood==0 it is the "aggregate data" that we have to use
            DataAmp1   = squeeze(MVSscalingTest.AmplitudesPerSubjectMRI{IndApproach}(CurrInds,:,1)); %for single voxel, i.e. NHood==0 it is the "aggregate data" that we have to use
            DataAmp2   = squeeze(MVSscalingTest.AmplitudesPerSubjectMRI{IndApproach}(CurrInds,:,2)); %for single voxel, i.e. NHood==0 it is the "aggregate data" that we have to use
        end
        if(all(isnan(DataLambda))||all(isnan(DataAmp1))||all(isnan(DataAmp2))) %if DataLambda is okay then we should not have a problem with the amplitudes, but safety is best.
            PermTest.DataQuality.UsableData{IndApproach}(IndSLCenterVox) = 0; %note this in the mask for all possibilities
        end
        
        Data = [DataLambda(:),DataAmp1(:),DataAmp2(:)];
        %Pearson correlation
        [RhoOrg1,pOrg1] = corr(Data,Regressors);
        PermTest.Results.Rho{IndApproach}(IndSLCenterVox,:,1,1)   = RhoOrg1(1,:); %RHO: Lambda         correalations with all covars
        PermTest.Results.Rho{IndApproach}(IndSLCenterVox,:,1,2)   = RhoOrg1(2,:); %RHO: Amplitude-MRI1 correalations with all covars
        PermTest.Results.Rho{IndApproach}(IndSLCenterVox,:,1,3)   = RhoOrg1(3,:); %RHO: Amplitude-MRI2 correalations with all covars
        
        PermTest.Results.pCorr{IndApproach}(IndSLCenterVox,:,1,1) = pOrg1(1,:); %pCorr: Lambda         correalations with all covars
        PermTest.Results.pCorr{IndApproach}(IndSLCenterVox,:,1,2) = pOrg1(2,:); %pCorr: Amplitude-MRI1 correalations with all covars
        PermTest.Results.pCorr{IndApproach}(IndSLCenterVox,:,1,3) = pOrg1(3,:); %pCorr: Amplitude-MRI2 correalations with all covars
        
        %Spearman correlation
        [RhoOrg2,pOrg2] = corr(Data,Regressors,'type','Spearman');
        PermTest.Results.Rho{IndApproach}(IndSLCenterVox,:,2,1)   = RhoOrg2(1,:); %RHO: Lambda         correalations with all covars
        PermTest.Results.Rho{IndApproach}(IndSLCenterVox,:,2,2)   = RhoOrg2(2,:); %RHO: Amplitude-MRI1 correalations with all covars
        PermTest.Results.Rho{IndApproach}(IndSLCenterVox,:,2,3)   = RhoOrg2(3,:); %RHO: Amplitude-MRI2 correalations with all covars
        
        PermTest.Results.pCorr{IndApproach}(IndSLCenterVox,:,2,1) = pOrg2(1,:); %pCorr: Lambda         correalations with all covars
        PermTest.Results.pCorr{IndApproach}(IndSLCenterVox,:,2,2) = pOrg2(2,:); %pCorr: Amplitude-MRI1 correalations with all covars
        PermTest.Results.pCorr{IndApproach}(IndSLCenterVox,:,2,3) = pOrg2(3,:); %pCorr: Amplitude-MRI2 correalations with all covars
        
        %% permtest (check if all NFail>NFailSkip --> skip,with message)
        NextPercentPermTestDone = 0;
        fprintf('\n');
        reverseStr = '';
        disp(['Starting permutation test... [SL#',num2str(IndSLCenterVox),'of',num2str(length(MVSscalingTest.SLight.SLightIndsInMaskCell)),']']);
        for IndPerm = 1:PermTest.NPermTest
            if(NextPercentPermTestDone==0)%start
                msg = sprintf('PermTests percent done: %03.0f...', floor(IndPerm*100/PermTest.NPermTest)); %Don't forget this semicolon
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                NextPercentPermTestDone = 1;
            else
                if(floor(IndPerm*100/PermTest.NPermTest)==NextPercentPermTestDone)
                    msg = sprintf('PermTests percent done: %03.0f...', floor(IndPerm*100/PermTest.NPermTest)); %Don't forget this semicolon
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    NextPercentPermTestDone = NextPercentPermTestDone+1;
                end
            end
            
            %Pearson correlation
            RhoPerm1 = corr(Data,Regressors(PermTest.PermMatrix(:,IndPerm),:));
            PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,1,1) = PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,1,1)+((sign(RhoOrg1(1,:)).*RhoPerm1(1,:))>=abs(RhoOrg1(1,:))); %Lambda:         fail only if equal or "higher" correlation IN THE SAME DIRECTION!!!
            PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,1,2) = PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,1,2)+((sign(RhoOrg1(2,:)).*RhoPerm1(2,:))>=abs(RhoOrg1(2,:))); %Amplitude-MRI1: fail only if equal or "higher" correlation IN THE SAME DIRECTION!!!
            PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,1,3) = PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,1,3)+((sign(RhoOrg1(3,:)).*RhoPerm1(3,:))>=abs(RhoOrg1(3,:))); %Amplitude-MRI2: fail only if equal or "higher" correlation IN THE SAME DIRECTION!!!
            
            PermTest.Results.NTests{IndApproach}(IndSLCenterVox,:,1,:) = PermTest.Results.NTests{IndApproach}(IndSLCenterVox,:,1,:)+1; %updata NTests
            
            %Spearman correlation
            RhoPerm2 = corr(Data,Regressors(PermTest.PermMatrix(:,IndPerm),:),'type','Spearman');
            PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,2,1) = PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,2,1)+((sign(RhoOrg2(1,:)).*RhoPerm2(1,:))>=abs(RhoOrg2(1,:))); %Lambda:         fail only if equal or "higher" correlation IN THE SAME DIRECTION!!!
            PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,2,2) = PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,2,2)+((sign(RhoOrg2(2,:)).*RhoPerm2(2,:))>=abs(RhoOrg2(2,:))); %Amplitude-MRI1: fail only if equal or "higher" correlation IN THE SAME DIRECTION!!!
            PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,2,3) = PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,2,3)+((sign(RhoOrg2(3,:)).*RhoPerm2(3,:))>=abs(RhoOrg2(3,:))); %Amplitude-MRI2: fail only if equal or "higher" correlation IN THE SAME DIRECTION!!!
            
            PermTest.Results.NTests{IndApproach}(IndSLCenterVox,:,2,:) = PermTest.Results.NTests{IndApproach}(IndSLCenterVox,:,2,:)+1; %updata NTests
            
            %% check if NFailSkip is reached for all MCovars,2CorrelationTests&3DataTypes --> we won't do 10000 tests for most voxels even though some could just be significant for a few or one covariate in one type of test for one certain data-variable
            CurrNFail = squeeze(PermTest.Results.NFail{IndApproach}(IndSLCenterVox,:,:,:));
            if(all(CurrNFail(:)>=PermTest.NFailSkip))
                disp(['NFailSkip reached after ',num2str(IndPerm),'of',num2str(PermTest.NPermTest),' possible permutations! Skipping to next searchlight.']);
                break;
            end
        end
        disp(' ');
        PermTest.Results.CurrVoxInd(IndApproach) = IndSLCenterVox;
        
        %% backup results
        MVSscalingTest.PermTest = PermTest;
        if(PrevPercentSave+PercentSaveStep==NextPercentDone)
            disp('Backing up results so far...');
            save(BackupSavePath,'MVSscalingTest');
            PrevPercentSave = PrevPercentSave+PercentSaveStep; %next interval at which backup is done.
        end
    end
    %% calculate pVals for permtest and assign back
    PermTest.Results.pPermTest{IndApproach} = PermTest.Results.NFail{IndApproach}./PermTest.Results.NTests{IndApproach}; %pVals (approximate)
    fprintf('DONE.\n');
end
MVSscalingTest.PermTest    = PermTest; %assign back
fprintf('ALL DONE.\n');

end