function TSNRscalingTest = RunTSNRscalingTest(TSNRscalingTest)
% This function tests the scaling values of the TSNRs between field strength (i.e. between MRIs),
% similar as in the MVSscalingTest function for the DualReg Amplitudes.
%
%Usage:
%       TSNRscalingTest = RunTSNRscalingTest(TSNRscalingTest);
%
%
%V1.0
%Date: V1.0(18.09.2015): initial implementation based on MVSscalingTest.m.
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% Remember
% TSNRscalingTest.Design  = TSNR_Design;
% TSNRscalingTest.Design.filelist
% TSNRscalingTest.Design.SubjNrMRInrRunNr
% TSNRscalingTest.Design.UniqueSubjNrs
% TSNRscalingTest.Design.UniqueMRInrs
% TSNRscalingTest.Design.NSubjs
% TSNRscalingTest.Design.N_MRIs
%
% TSNRscalingTest.Mask.MaskNII_FilePath= MaskPath;
% ...
% TSNRscalingTest.SLight  = SLight;
% TSNRscalingTest.AmplitudesPerSubjectMRI = cell(2,1); %AmplitudesPerSubjectMRI{1or2} will be (NVoxelWholeBrain,NSubj,NMRI); -i.e. average over Runs OR in case of TSNR{2} nanmedian over space (searchlight) and Runs
% TSNRscalingTest.ScalingPerSubject       = cell(2,1); %ScalingPerSubject{1or2}       will be (NVoxelWholeBrain,NSubj);      -i.e. fraction of MRI==2 divided by MRI==1. 


% TSNRscalingTest.ParamEstLambdaTest.mH0               = 2*sqrt(2); %the expected median, ie. null hypothesis is H0(m(Lambda)=2*sqrt(2))
% TSNRscalingTest.ParamEstLambdaTest.mConfInt          = conf; %the confidence interval (0.95)
% TSNRscalingTest.ParamEstLambdaTest.mConfIntIQRfactor = (1/(norminv(.75)-norminv(.25))) * sqrt(pi/2) * (norminv(1-(1-conf)/2)+(norminv(1-(1-conf)/2)/sqrt(2)))/2; %factor before the IQR (inter quartile range, i.e. 0.25 to 0.75 range of data) that gives the confidence interval (symmetric assumption)
% TSNRscalingTest.ParamEstLambdaTest.SignRankTest.zVals      = []; %(NVoxel,3);   %init later %z-score values from signed rank test for median subject, aggregate data and median searchlight
% TSNRscalingTest.ParamEstLambdaTest.SignRankTest.pVals      = []; %(NVoxel,3);   %init later %   p    values from signed rank test for median subject, aggregate data and median searchlight
% TSNRscalingTest.ParamEstLambdaTest.SignRankTest.signedrank = []; %(NVoxel,3);   %init later %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
% TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth        = []; %(NVoxel,4,3); %init later %parameter estimates for median subject, aggregate data and median searchlight (that is the last dim); MiddleDim: 1==Median,2==1stQuartile,3==3rdQuartile,4==CIwidth 

%% init to prepare for processing 
ApproachesInfo = TSNRscalingTest.ApproachesInfo;
disp('------------------------------------------------------------------');
disp('Parameter estimate & testing of Lambda distribution from TSNR values...');
TSNRscalingTest.DataQuality.UsableData                     = cell(size(ApproachesInfo,1),1);  %if any of the aggregate or median subject or median searchlight ONLY CONTAIN NaNs then this in not usable data. --> keep track of this. (but assume in the beginning that everything is usable)
TSNRscalingTest.ParamEstLambdaTest.SignRankTest.zVals      = cell(size(ApproachesInfo,1),1);   %(NVoxel,3);   %z-score values from signed rank test for median subject, aggregate data and median searchlight
TSNRscalingTest.ParamEstLambdaTest.SignRankTest.pVals      = cell(size(ApproachesInfo,1),1);   %(NVoxel,3);   %   p    values from signed rank test for median subject, aggregate data and median searchlight
TSNRscalingTest.ParamEstLambdaTest.SignRankTest.signedrank = cell(size(ApproachesInfo,1),1);   %(NVoxel,3);   %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth        = cell(size(ApproachesInfo,1),1);

%% do the processing & loop over approaches
for IndApproach = 1:size(ApproachesInfo,1)
    disp(['Data is from approach: ',ApproachesInfo{IndApproach,1},': ',ApproachesInfo{IndApproach,2}]);
    %% init
    NextPercentDone = 0; %init %for indicating the progress 1%-stepwise
    TSNRscalingTest.DataQuality.UsableData{IndApproach}                     = ones(length(TSNRscalingTest.SLight.SLightIndsInMaskCell),1);  %if any of the aggregate or median subject or median searchlight ONLY CONTAIN NaNs then this in not usable data. --> keep track of this. (but assume in the beginning that everything is usable)
    TSNRscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}      =  NaN(length(TSNRscalingTest.SLight.SLightIndsInMaskCell),3);   %(NVoxel,3);   %z-score values from signed rank test for median subject, aggregate data and median searchlight
    TSNRscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}      =  NaN(length(TSNRscalingTest.SLight.SLightIndsInMaskCell),3);   %(NVoxel,3);   %   p    values from signed rank test for median subject, aggregate data and median searchlight
    TSNRscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach} =  NaN(length(TSNRscalingTest.SLight.SLightIndsInMaskCell),3);   %(NVoxel,3);   %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
    TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}        =  NaN(length(TSNRscalingTest.SLight.SLightIndsInMaskCell),4,3); %(NVoxel,4,3); %parameter estimates for median subject, aggregate data and median searchlight (that is the last dim); MiddleDim: 1==Median,2==1stQuartile,3==3rdQuartile,4==CIwidth
    
    reverseStr = '';
    for IndSLCenterVox = 1:length(TSNRscalingTest.SLight.SLightIndsInMaskCell)
        if(NextPercentDone==0)%start
            %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(TSNRscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
            msg = sprintf('      %03.0f percent done...', floor(IndSLCenterVox*100/length(TSNRscalingTest.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            NextPercentDone = 1;
        else
            if(floor(IndSLCenterVox*100/length(TSNRscalingTest.SLight.SLightIndsInMaskCell))==NextPercentDone)
                %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(TSNRscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
                msg = sprintf('      %03.0f percent done...', floor(IndSLCenterVox*100/length(TSNRscalingTest.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                NextPercentDone = NextPercentDone+1;
            end
        end
        
        %% get Data, i.e. the aggregate data for the searchlight
        CurrInds = TSNRscalingTest.SLight.SLightIndsInMaskCell{IndSLCenterVox};
        Data     = abs(TSNRscalingTest.ScalingPerSubject{IndApproach}(CurrInds,:));
        if(TSNRscalingTest.SLight.NHood==0)
            if(all(isnan(Data(:))))
                TSNRscalingTest.DataQuality.UsableData{IndApproach}(IndSLCenterVox) = 0; %note this in the mask for all possibilities
                msg = 'skipping one...';
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                continue;
            end
        else
            if(all(isnan(Data(:)))||all(isnan(nanmedian(Data,2)))||all(isnan(nanmedian(Data,1))))
                TSNRscalingTest.DataQuality.UsableData{IndApproach}(IndSLCenterVox) = 0; %note this in the mask for all possibilities
            end
        end
        
        %% get parameter estimates
        TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,1:3,2)     = quantile(                  Data(:), [.5 .25 .75]);
        if(TSNRscalingTest.SLight.NHood>0)
            TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,1:3,1) = quantile(squeeze(nanmedian(Data,2)),[.5 .25 .75]);
            TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,1:3,3) = quantile(squeeze(nanmedian(Data,1)),[.5 .25 .75]);
        end
        
        %CIwidth = mConfIntIQRfactor*(Qrt3-Qrt1)/sqrt(length(InputData)); %--> Qrt3-Qrt1 == iQR   %OLD: CI(1/2) = median -/+ 1.57*(Qrt3-Qrt1)/sqrt(length(InputData)); %--> Qrt3-Qrt1 == iQR
        TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,4,2)     = TSNRscalingTest.ParamEstLambdaTest.mConfIntIQRfactor*(TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,3,2)-TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,2,2))/sqrt(length(Data(:)));
        if(TSNRscalingTest.SLight.NHood>0)
            TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,4,1) = TSNRscalingTest.ParamEstLambdaTest.mConfIntIQRfactor*(TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,3,1)-TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,2,1))/sqrt(length(squeeze(nanmedian(Data,2))));
            TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,4,3) = TSNRscalingTest.ParamEstLambdaTest.mConfIntIQRfactor*(TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,3,3)-TSNRscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,2,3))/sqrt(length(squeeze(nanmedian(Data,1))));
        end
        
        %% signrank-test
        [p_Agg,   h_Agg,   stats_Agg]       = signrank(                  Data(:), TSNRscalingTest.ParamEstLambdaTest.mH0);
        if(TSNRscalingTest.SLight.NHood>0)
            [p_MedS,  h_MedS,  stats_MedS]  = signrank(squeeze(nanmedian(Data,2)),TSNRscalingTest.ParamEstLambdaTest.mH0);
            [p_MedSL, h_MedSL, stats_MedSL] = signrank(squeeze(nanmedian(Data,1)),TSNRscalingTest.ParamEstLambdaTest.mH0);
        end
        
        %% assign p vals
        TSNRscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}(IndSLCenterVox,2)     = p_Agg;
        if(TSNRscalingTest.SLight.NHood>0)
            TSNRscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}(IndSLCenterVox,1) = p_MedS;
            TSNRscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}(IndSLCenterVox,3) = p_MedSL;
        end
        
        %% assign z vals
        if(isfield(stats_Agg,'zval'))
            TSNRscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}(IndSLCenterVox,2)     = stats_Agg.zval;
        end
        if(TSNRscalingTest.SLight.NHood>0)
            if(isfield(stats_MedS,'zval'))
                TSNRscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}(IndSLCenterVox,1) = stats_MedS.zval;
            end
            if(isfield(stats_MedSL,'zval'))
                TSNRscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}(IndSLCenterVox,3) = stats_MedSL.zval;
            end
        end
        
        %% assign signed rank test statistic value
        TSNRscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach}(IndSLCenterVox,2)     = stats_Agg.signedrank;
        if(TSNRscalingTest.SLight.NHood>0)
            TSNRscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach}(IndSLCenterVox,1) = stats_MedS.signedrank;
            TSNRscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach}(IndSLCenterVox,3) = stats_MedSL.signedrank;
        end
        
        %% sign-test NOT USED RIGHT NOW
        % [MVSscalingTest.Results.Aggregate.p_sign(IndSLCenterVox),   MVSscalingTest.Results.Aggregate.h_sign(IndSLCenterVox),   statsAggregate_sign]    = signtest(                  Data(:), MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
        % [MVSscalingTest.Results.MedianSubj.p_sign(IndSLCenterVox),  MVSscalingTest.Results.MedianSubj.h_sign(IndSLCenterVox),  statsMedianSubj_sign]   = signtest(squeeze(nanmedian(Data,2)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
        % [MVSscalingTest.Results.MedianSLight.p_sign(IndSLCenterVox),MVSscalingTest.Results.MedianSLight.h_sign(IndSLCenterVox),statsMedianSLight_sign] = signtest(squeeze(nanmedian(Data,1)),MVSscalingTest.StatisticsSettings.ExpectedMedian,'tail',MVSscalingTest.StatisticsSettings.TailType);
        
    end
    fprintf('DONE.\n');
end
fprintf('ALL DONE.\n');

end