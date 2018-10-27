function MVSscalingTest = RunMVSscalingTest(MVSscalingTest)
% This function performes the scaling values parameter estimate, i.e. the
% median scaling value for the median subject, the aggregate and the median
% SLight data, and also determines the quartiles and Confidence Interval width.
%
% Furthermore, this function performs a Wilcoxon signed rank test on the
% same data as above (the median subject, the aggregate and the median
% SLight data) with the null hypothesis H0(\Lambda=2*sqrt(2)).
% This can later be used to determine where this hypothesis is accepted and
% where it is rejected. 
% (I.e. make a d=log2(d*)-values plot for d*(zVal,Cutoff)=Cutoff/abs(zVal).)
% See (Boegle et al., 2015) [hopefully submitted, reviewed and published soon...].
%
%Usage:
%       MVSscalingTest = RunMVSscalingTest(MVSscalingTest);
%
%V1.1
%Date: V1.1(04.09.2015): update to process data prepared in two ways (Voxel-Wise scaling values & Searchlight-averaged amplitudes turned into scaling data), both analysed with searchlight here as before. V1.0(24.7.2015) (initial implementation based on test script for analysis of scaling data.)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

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
% MVSscalingTest.ParamEstLambdaTest.mH0               = 2*sqrt(2); %the expected median, ie. null hypothesis is H0(m(Lambda)=2*sqrt(2))
% MVSscalingTest.ParamEstLambdaTest.mConfInt          = conf; %the confidence interval (0.95)
% MVSscalingTest.ParamEstLambdaTest.mConfIntIQRfactor = (1/(norminv(.75)-norminv(.25))) * sqrt(pi/2) * (norminv(1-(1-conf)/2)+(norminv(1-(1-conf)/2)/sqrt(2)))/2; %factor before the IQR (inter quartile range, i.e. 0.25 to 0.75 range of data) that gives the confidence interval (symmetric assumption)
% MVSscalingTest.ParamEstLambdaTest.SignRankTest.zVals      = []; %(NVoxel,3);   %init later %z-score values from signed rank test for median subject, aggregate data and median searchlight
% MVSscalingTest.ParamEstLambdaTest.SignRankTest.pVals      = []; %(NVoxel,3);   %init later %   p    values from signed rank test for median subject, aggregate data and median searchlight
% MVSscalingTest.ParamEstLambdaTest.SignRankTest.signedrank = []; %(NVoxel,3);   %init later %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
% MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth        = []; %(NVoxel,4,3); %init later %parameter estimates for median subject, aggregate data and median searchlight (that is the last dim); MiddleDim: 1==Median,2==1stQuartile,3==3rdQuartile,4==CIwidth 

%% init to prepare for processing 
ApproachesInfo = MVSscalingTest.ApproachesInfo;
disp('------------------------------------------------------------------');
disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution...']);
MVSscalingTest.DataQuality.UsableData                     = cell(size(ApproachesInfo,1),1);  %if any of the aggregate or median subject or median searchlight ONLY CONTAIN NaNs then this in not usable data. --> keep track of this. (but assume in the beginning that everything is usable)
MVSscalingTest.ParamEstLambdaTest.SignRankTest.zVals      = cell(size(ApproachesInfo,1),1);   %(NVoxel,3);   %z-score values from signed rank test for median subject, aggregate data and median searchlight
MVSscalingTest.ParamEstLambdaTest.SignRankTest.pVals      = cell(size(ApproachesInfo,1),1);   %(NVoxel,3);   %   p    values from signed rank test for median subject, aggregate data and median searchlight
MVSscalingTest.ParamEstLambdaTest.SignRankTest.signedrank = cell(size(ApproachesInfo,1),1);   %(NVoxel,3);   %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth        = cell(size(ApproachesInfo,1),1);

%% do the processing & loop over approaches
for IndApproach = 1:size(ApproachesInfo,1)
    disp(['Data is from approach: ',ApproachesInfo{IndApproach,1},': ',ApproachesInfo{IndApproach,2}]);
    %% init
    NextPercentDone = 0; %init %for indicating the progress 1%-stepwise
    MVSscalingTest.DataQuality.UsableData{IndApproach}                     = ones(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1);  %if any of the aggregate or median subject or median searchlight ONLY CONTAIN NaNs then this in not usable data. --> keep track of this. (but assume in the beginning that everything is usable)
    MVSscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}      =  NaN(length(MVSscalingTest.SLight.SLightIndsInMaskCell),3);   %(NVoxel,3);   %z-score values from signed rank test for median subject, aggregate data and median searchlight
    MVSscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}      =  NaN(length(MVSscalingTest.SLight.SLightIndsInMaskCell),3);   %(NVoxel,3);   %   p    values from signed rank test for median subject, aggregate data and median searchlight
    MVSscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach} =  NaN(length(MVSscalingTest.SLight.SLightIndsInMaskCell),3);   %(NVoxel,3);   %signed rank test statistic values from signed rank test for median subject, aggregate data and median searchlight
    MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}        =  NaN(length(MVSscalingTest.SLight.SLightIndsInMaskCell),4,3); %(NVoxel,4,3); %parameter estimates for median subject, aggregate data and median searchlight (that is the last dim); MiddleDim: 1==Median,2==1stQuartile,3==3rdQuartile,4==CIwidth
    
    reverseStr = '';
    for IndSLCenterVox = 1:length(MVSscalingTest.SLight.SLightIndsInMaskCell)
        if(NextPercentDone==0)%start
            %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
            msg = sprintf('      %03.0f percent done...', floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            NextPercentDone = 1;
        else
            if(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))==NextPercentDone)
                %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
                msg = sprintf('      %03.0f percent done...', floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                NextPercentDone = NextPercentDone+1;
            end
        end
        
        %% get Data, i.e. the aggregate data for the searchlight
        CurrInds = MVSscalingTest.SLight.SLightIndsInMaskCell{IndSLCenterVox};
        Data     = abs(MVSscalingTest.ScalingPerSubject{IndApproach}(CurrInds,:));
        if(MVSscalingTest.SLight.NHood==0)
            if(all(isnan(Data(:))))
                MVSscalingTest.DataQuality.UsableData{IndApproach}(IndSLCenterVox) = 0; %note this in the mask for all possibilities
                msg = 'skipping one...';
                fprintf([reverseStr, msg]);
                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                continue;
            end
        else
            if(all(isnan(Data(:)))||all(isnan(nanmedian(Data,2)))||all(isnan(nanmedian(Data,1))))
                MVSscalingTest.DataQuality.UsableData{IndApproach}(IndSLCenterVox) = 0; %note this in the mask for all possibilities
            end
        end
        
        %% get parameter estimates
        MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,1:3,2)     = quantile(                  Data(:), [.5 .25 .75]);
        if(MVSscalingTest.SLight.NHood>0)
            MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,1:3,1) = quantile(squeeze(nanmedian(Data,2)),[.5 .25 .75]);
            MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,1:3,3) = quantile(squeeze(nanmedian(Data,1)),[.5 .25 .75]);
        end
        
        %CIwidth = mConfIntIQRfactor*(Qrt3-Qrt1)/sqrt(length(InputData)); %--> Qrt3-Qrt1 == iQR   %OLD: CI(1/2) = median -/+ 1.57*(Qrt3-Qrt1)/sqrt(length(InputData)); %--> Qrt3-Qrt1 == iQR
        MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,4,2)     = MVSscalingTest.ParamEstLambdaTest.mConfIntIQRfactor*(MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,3,2)-MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,2,2))/sqrt(length(Data(:)));
        if(MVSscalingTest.SLight.NHood>0)
            MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,4,1) = MVSscalingTest.ParamEstLambdaTest.mConfIntIQRfactor*(MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,3,1)-MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,2,1))/sqrt(length(squeeze(nanmedian(Data,2))));
            MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,4,3) = MVSscalingTest.ParamEstLambdaTest.mConfIntIQRfactor*(MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,3,3)-MVSscalingTest.ParamEstLambdaTest.MedianQrtCIwidth{IndApproach}(IndSLCenterVox,2,3))/sqrt(length(squeeze(nanmedian(Data,1))));
        end
        
        %% signrank-test
        [p_Agg,   h_Agg,   stats_Agg]       = signrank(                  Data(:), MVSscalingTest.ParamEstLambdaTest.mH0);
        if(MVSscalingTest.SLight.NHood>0)
            [p_MedS,  h_MedS,  stats_MedS]  = signrank(squeeze(nanmedian(Data,2)),MVSscalingTest.ParamEstLambdaTest.mH0);
            [p_MedSL, h_MedSL, stats_MedSL] = signrank(squeeze(nanmedian(Data,1)),MVSscalingTest.ParamEstLambdaTest.mH0);
        end
        
        %% assign p vals
        MVSscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}(IndSLCenterVox,2)     = p_Agg;
        if(MVSscalingTest.SLight.NHood>0)
            MVSscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}(IndSLCenterVox,1) = p_MedS;
            MVSscalingTest.ParamEstLambdaTest.SignRankTest.pVals{IndApproach}(IndSLCenterVox,3) = p_MedSL;
        end
        
        %% assign z vals
        if(isfield(stats_Agg,'zval'))
            MVSscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}(IndSLCenterVox,2)     = stats_Agg.zval;
        end
        if(MVSscalingTest.SLight.NHood>0)
            if(isfield(stats_MedS,'zval'))
                MVSscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}(IndSLCenterVox,1) = stats_MedS.zval;
            end
            if(isfield(stats_MedSL,'zval'))
                MVSscalingTest.ParamEstLambdaTest.SignRankTest.zVals{IndApproach}(IndSLCenterVox,3) = stats_MedSL.zval;
            end
        end
        
        %% assign signed rank test statistic value
        MVSscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach}(IndSLCenterVox,2)     = stats_Agg.signedrank;
        if(MVSscalingTest.SLight.NHood>0)
            MVSscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach}(IndSLCenterVox,1) = stats_MedS.signedrank;
            MVSscalingTest.ParamEstLambdaTest.SignRankTest.signedrank{IndApproach}(IndSLCenterVox,3) = stats_MedSL.signedrank;
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