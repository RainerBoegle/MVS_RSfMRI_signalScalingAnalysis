%% get data
load('MVSscaling_IC(MELODIC).mat');

%% find outliers ORIGINAL DATA
[DataNew,IsOutlier] = FindOutliersLikeBB(MVSscaling.ScalingPerSubject);
[DataNewAll,IsOutlierAll] = FindOutliersLikeBB(MVSscaling.ScalingPerSubject(:));
IsOutlierCombinedFindFun = abs(IsOutlier)+abs(reshape(IsOutlierAll,size(MVSscaling.ScalingPerSubject)));
IsOriginalNaNs = double(isnan(MVSscaling.ScalingPerSubject));
IsOutlierCombined = abs(IsOutlier)+abs(reshape(IsOutlierAll,size(MVSscaling.ScalingPerSubject)))+IsOriginalNaNs;

%% make new data without outliers
DataOrg = MVSscaling.ScalingPerSubject;
%DataOrg(IsOutlierCombined~=0) = NaN;

%% test for normality
jbtestPerSubj.h          = zeros(size(DataOrg,2),1);
jbtestPerSubj.p          = ones(size(DataOrg,2),1);
jbtestPerSubj.h_outliers = zeros(size(DataOrg,2),1);
jbtestPerSubj.p_outliers = ones(size(DataOrg,2),1);
for IndSubj = 1:size(DataOrg,2)
    [jbtestPerSubj.h_outliers(IndSubj),jbtestPerSubj.p_outliers(IndSubj)] = jbtest(DataOrg(IsOutlierCombined(:,IndSubj)~=0,IndSubj));
    [jbtestPerSubj.h(IndSubj),         jbtestPerSubj.p(IndSubj)]          = jbtest(DataOrg(IsOutlierCombined(:,IndSubj)==0,IndSubj));
end

%% do histogram per subject
Hist_xvals = [0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [0, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist(abs(DataOrg(IsOutlierCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['Outliers IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']);
    subplot(2,2,2); hist(abs(DataOrg(IsOutlierCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['Outliers IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist(abs(DataOrg(IsOutlierCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']);
    subplot(2,2,4); hist(abs(DataOrg(IsOutlierCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% try hist without abs
Hist_xvals = [-1000:2:-100,-100:1:-50,-50:.5:-25,-25:.25:-10,-10:.1:0,0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [-6, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist((DataOrg(IsOutlierCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['Outliers IndSubj= ',num2str(IndSubj),' (DataOrg)']);
    subplot(2,2,2); hist((DataOrg(IsOutlierCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['Outliers IndSubj= ',num2str(IndSubj),' (DataOrg)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist((DataOrg(IsOutlierCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataOrg)']);
    subplot(2,2,4); hist((DataOrg(IsOutlierCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataOrg)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% as above (starting after loading the data, but do outlier detection on abs of data.
%% make new data without outliers
DataAbs = abs(MVSscaling.ScalingPerSubject);
%DataOrg(IsOutlierCombined~=0) = NaN;
%% find outliers ORIGINAL DATA
[DataNew,IsOutlier_DataAbs] = FindOutliersLikeBB(DataAbs);
[DataNewAll,IsOutlierAll_DataAbs] = FindOutliersLikeBB(DataAbs(:));
IsOutlierCombinedFindFun_DataAbs = abs(IsOutlier_DataAbs)+abs(reshape(IsOutlierAll_DataAbs,size(DataAbs)));
IsOriginalNaNs_DataAbs = double(isnan(DataAbs));
IsOutlierCombined_DataAbs = abs(IsOutlier_DataAbs)+abs(reshape(IsOutlierAll_DataAbs,size(DataAbs)))+IsOriginalNaNs_DataAbs;


%% test for normality
jbtestPerSubj.DataAbs.h          = zeros(size(DataAbs,2),1);
jbtestPerSubj.DataAbs.p          = ones(size(DataAbs,2),1);
jbtestPerSubj.DataAbs.h_outliers = zeros(size(DataAbs,2),1);
jbtestPerSubj.DataAbs.p_outliers = ones(size(DataAbs,2),1);
for IndSubj = 1:size(DataOrg,2)
    [jbtestPerSubj.DataAbs.h_outliers(IndSubj),jbtestPerSubj.DataAbs.p_outliers(IndSubj)] = jbtest(DataAbs(IsOutlierCombined_DataAbs(:,IndSubj)~=0,IndSubj));
    [jbtestPerSubj.DataAbs.h(IndSubj),         jbtestPerSubj.DataAbs.p(IndSubj)]          = jbtest(DataAbs(IsOutlierCombined_DataAbs(:,IndSubj)==0,IndSubj));
end

%% do histogram per subject
Hist_xvals = [0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [0, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist(abs(DataAbs(IsOutlierCombined_DataAbs(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersDataAbs IndSubj= ',num2str(IndSubj),' (DataAbs)']);
    subplot(2,2,2); hist(abs(DataAbs(IsOutlierCombined_DataAbs(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersDataAbs IndSubj= ',num2str(IndSubj),' (DataAbs)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist(abs(DataAbs(IsOutlierCombined_DataAbs(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataAbs)']);
    subplot(2,2,4); hist(abs(DataAbs(IsOutlierCombined_DataAbs(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataAbs)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% combine judgements about outliers of both approaches
IsOutlierBothCombined = IsOutlierCombined_DataAbs + IsOutlierCombined;

%% do histogram per subject
Hist_xvals = [0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [0, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist(abs(DataAbs(IsOutlierBothCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersBothCombined IndSubj= ',num2str(IndSubj),' (DataAbs)']);
    subplot(2,2,2); hist(abs(DataAbs(IsOutlierBothCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersBothCombined IndSubj= ',num2str(IndSubj),' (DataAbs)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist(abs(DataAbs(IsOutlierBothCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataAbs)']);
    subplot(2,2,4); hist(abs(DataAbs(IsOutlierBothCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataAbs)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% also do avove plots again
%% do histogram per subject
Hist_xvals = [0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [0, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist(abs(DataOrg(IsOutlierBothCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersBothCombined IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']);
    subplot(2,2,2); hist(abs(DataOrg(IsOutlierBothCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersBothCombined IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist(abs(DataOrg(IsOutlierBothCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']);
    subplot(2,2,4); hist(abs(DataOrg(IsOutlierBothCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataOrg)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% try hist without abs
Hist_xvals = [-1000:2:-100,-100:1:-50,-50:.5:-25,-25:.25:-10,-10:.1:0,0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [-6, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist((DataOrg(IsOutlierBothCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersBothCombined IndSubj= ',num2str(IndSubj),' (DataOrg)']);
    subplot(2,2,2); hist((DataOrg(IsOutlierBothCombined(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['OutliersBothCombined IndSubj= ',num2str(IndSubj),' (DataOrg)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist((DataOrg(IsOutlierBothCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataOrg)']);
    subplot(2,2,4); hist((DataOrg(IsOutlierBothCombined(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataOrg)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end


%% New try alltogether
%% What if I throw out everything over a threshold for the absolute of the data, i.e. symmetrically
%% and then apply the above outlier methods and look at the data again?
%% The initial threshold could be quite high, e.g. 8-12...
Threshold = 10;
DataThres = DataOrg;
DataThres(DataAbs>Threshold) = NaN;

%% find outliers
[DataNew,IsOutlier_DataThres] = FindOutliersLikeBB(DataThres);
[DataNewAll,IsOutlierAll_DataThres] = FindOutliersLikeBB(DataThres(:));
IsOutlierCombinedFindFun_DataThres = abs(IsOutlier_DataThres)+abs(reshape(IsOutlierAll_DataThres,size(DataThres)));
IsOriginalNaNs_DataThres = double(isnan(DataThres));
IsOutlierCombined_DataThres = abs(IsOutlier_DataThres)+abs(reshape(IsOutlierAll_DataThres,size(DataThres)))+IsOriginalNaNs_DataThres;

%% do histogram plots
Hist_xvals = [0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [0, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist(abs(DataThres(IsOutlierCombined_DataThres(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']);
    subplot(2,2,2); hist(abs(DataThres(IsOutlierCombined_DataThres(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist(abs(DataThres(IsOutlierCombined_DataThres(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']);
    subplot(2,2,4); hist(abs(DataThres(IsOutlierCombined_DataThres(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% try hist without abs
Hist_xvals = [-1000:2:-100,-100:1:-50,-50:.5:-25,-25:.25:-10,-10:.1:0,0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [-6, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist((DataThres(IsOutlierCombined_DataThres(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres IndSubj= ',num2str(IndSubj),' (DataThres)']);
    subplot(2,2,2); hist((DataThres(IsOutlierCombined_DataThres(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres IndSubj= ',num2str(IndSubj),' (DataThres)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist((DataThres(IsOutlierCombined_DataThres(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataThres)']);
    subplot(2,2,4); hist((DataThres(IsOutlierCombined_DataThres(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataThres)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% find outliers abs DataThres as well and combine with others
[DataNew,IsOutlier_DataThresAbs] = FindOutliersLikeBB(abs(DataThres));
[DataNewAll,IsOutlierAll_DataThresAbs] = FindOutliersLikeBB(abs(DataThres(:)));
IsOutlierCombinedFindFun_DataThresAbs = abs(IsOutlier_DataThresAbs)+abs(reshape(IsOutlierAll_DataThresAbs,size(DataThres)));
IsOutlierCombined_DataThres_n_DataThresAbs = abs(IsOutlier_DataThresAbs)+abs(reshape(IsOutlierAll_DataThresAbs,size(DataThres)))+abs(IsOutlier_DataThres)+abs(reshape(IsOutlierAll_DataThres,size(DataThres)))+IsOriginalNaNs_DataThres;

%% do histogram plots using new combination of outliers
Hist_xvals = [0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [0, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist(abs(DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres_n_DataThresAbs IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']);
    subplot(2,2,2); hist(abs(DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres_n_DataThresAbs IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist(abs(DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']);
    subplot(2,2,4); hist(abs(DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (useAbsDataThres)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% try hist without abs
Hist_xvals = [-1000:2:-100,-100:1:-50,-50:.5:-25,-25:.25:-10,-10:.1:0,0:0.1:10,10:0.25:25,25:.5:50,50:1:100,100:2:1000];
XLimitSubPlot2 = [-6, 6];

for IndSubj = 1:size(DataOrg,2)
    figure(42); clf;
    subplot(2,2,1); hist((DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres_n_DataThresAbs IndSubj= ',num2str(IndSubj),' (DataThres)']);
    subplot(2,2,2); hist((DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)~=0,IndSubj)),Hist_xvals); title(['IsOutlierCombined_DataThres_n_DataThresAbs IndSubj= ',num2str(IndSubj),' (DataThres)']); xlim(XLimitSubPlot2);
    subplot(2,2,3); hist((DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataThres)']);
    subplot(2,2,4); hist((DataThres(IsOutlierCombined_DataThres_n_DataThresAbs(:,IndSubj)==0,IndSubj)),Hist_xvals); title(['IndSubj= ',num2str(IndSubj),' (DataThres)']); xlim(XLimitSubPlot2);
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% New try HISTOGRAM based
%% Run Histogram symmectrically around zero from -max(abs(Data)) to max(abs(Data))
%% select a specific step, e.g. 0.1 or 0.2 --> look at distribution of counts (just to get a feeling for it)
%% throw out those bins that are below a threshold or derived value from distribution
%% for those bins that are still there create the regions of values that are still allowed and only keep them from the original data.
%% i.e. create the IsOutlier mask from this.
%% ==> check out the data afterwards. orgininal-outliers and abs(original-outliers)
%% try the above outlier algoriths on this data as well to see what the split is then.
Hist_xStepDeNoise = [0.1,.25,  .5,   1,   10,   100, 1000]; %step
Hist_LimitsDeNoise= [ 10, 50, 100, 500, 1000, 10000,  Inf];
Hist_xvalsDeNoise = cell(size(DataOrg,2),1); %for each subject
Hist_CountsDeNoise= cell(size(DataOrg,2),1); %for each subject
for IndSubj = 1:size(DataOrg,2)
    %create centers within limits
    for IndStep = 1:length(Hist_LimitsDeNoise)
        if(isfinite(Hist_LimitsDeNoise(IndStep)))
            if(Hist_LimitsDeNoise(IndStep)<max(abs(DataOrg(:,IndSubj))))
                if(IndStep==1)
                    Hist_xvalsDeNoise{IndSubj}  = [Hist_xvalsDeNoise{IndSubj}, -Hist_LimitsDeNoise(IndStep):Hist_xStepDeNoise(IndStep):Hist_LimitsDeNoise(IndStep)];
                else
                    Hist_xvalsDeNoise{IndSubj}  = [Hist_xvalsDeNoise{IndSubj}, -Hist_LimitsDeNoise(IndStep):Hist_xStepDeNoise(IndStep):-Hist_LimitsDeNoise(IndStep-1),Hist_LimitsDeNoise(IndStep-1):Hist_xStepDeNoise(IndStep):Hist_LimitsDeNoise(IndStep)];
                end
            else
                if(IndStep==1)
                    Hist_xvalsDeNoise{IndSubj}  = [Hist_xvalsDeNoise{IndSubj}, -max(abs(DataOrg(:,IndSubj))):Hist_xStepDeNoise(IndStep):max(abs(DataOrg(:,IndSubj)))];
                else
                    Hist_xvalsDeNoise{IndSubj}  = [Hist_xvalsDeNoise{IndSubj}, -max(abs(DataOrg(:,IndSubj))):Hist_xStepDeNoise(IndStep):-Hist_LimitsDeNoise(IndStep-1),Hist_LimitsDeNoise(IndStep-1):Hist_xStepDeNoise(IndStep):max(abs(DataOrg(:,IndSubj)))];
                end
                disp(['IndSub= ',num2str(IndSubj),',IndStep= ',num2str(IndStep),' stopping bincenter creation at maximum. (current stepwidth= ',num2str(Hist_xStepDeNoise(IndStep)),')']);
                break; %we are done!
            end
        else
            Hist_xvalsDeNoise{IndSubj}  = [Hist_xvalsDeNoise{IndSubj}, -max(abs(DataOrg(:,IndSubj))):Hist_xStepDeNoise(IndStep):-Hist_LimitsDeNoise(IndStep-1),Hist_LimitsDeNoise(IndStep-1):Hist_xStepDeNoise(IndStep):max(abs(DataOrg(:,IndSubj)))];
            if(IndStep~=length(Hist_LimitsDeNoise)) %this would be an error, but we can avoid further danger by breaking out of the loop.
                disp('WTF! breaking out of steps-loop.');
                break;
            else
                disp(['IndSub= ',num2str(IndSubj),' stopping bincenter creation at maximum by force.']);
            end
        end
    end
    %order values and remove the redundant ones
    Hist_xvalsDeNoise{IndSubj} = unique(Hist_xvalsDeNoise{IndSubj})';
    Hist_CountsDeNoise{IndSubj} = hist(DataOrg(:,IndSubj),Hist_xvalsDeNoise{IndSubj});
    
    figure(81); clf;
    subplot(1,2,1); bar(Hist_xvalsDeNoise{IndSubj},Hist_CountsDeNoise{IndSubj}); title(['IndSubj= ',num2str(IndSubj),' (Supp==[',num2str(min(Hist_xvalsDeNoise{IndSubj})),':',num2str(max(Hist_xvalsDeNoise{IndSubj})),'])']); set(gca,'YScale','log');
    subplot(1,2,2); boxplot(Hist_CountsDeNoise{IndSubj}(:)); title(['IndSubj= ',num2str(IndSubj),' (Supp==[',num2str(min(Hist_xvalsDeNoise{IndSubj})),':',num2str(max(Hist_xvalsDeNoise{IndSubj})),'])']); set(gca,'YScale','log');
    
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end

%% Hist for abs data
%% try the above outlier algoriths on this data as well to see what the split is then.
Hist_xStepDeNoise_DataAbs = [0.05, 0.1, .25,  .5,    1,    10,    100, 1000]; %step
Hist_LimitsDeNoise_DataAbs= [  10,  50, 100, 500, 1000, 10000, 100000,  Inf];
Hist_xvalsDeNoise_DataAbs = cell(size(DataOrg,2),1); %for each subject
Hist_xvalGrDeNoise_DataAbs= cell(size(DataOrg,2),1); %for each subject
Hist_CountsDeNoise_DataAbs= cell(size(DataOrg,2),1); %for each subject
for IndSubj = 1:size(DataAbs,2)
    %create centers within limits
    for IndStep = 1:length(Hist_LimitsDeNoise_DataAbs)
        if(isfinite(Hist_LimitsDeNoise_DataAbs(IndStep)))
            if(Hist_LimitsDeNoise_DataAbs(IndStep)<max(abs(DataOrg(:,IndSubj))))
                if(IndStep==1)
                    Hist_xvalsDeNoise_DataAbs{IndSubj}  = [Hist_xvalsDeNoise_DataAbs{IndSubj}, 0:Hist_xStepDeNoise_DataAbs(IndStep):Hist_LimitsDeNoise_DataAbs(IndStep)];
                    Hist_xvalGrDeNoise_DataAbs{IndSubj} = [Hist_xvalGrDeNoise_DataAbs{IndSubj},IndStep.*ones(size(0:Hist_xStepDeNoise_DataAbs(IndStep):Hist_LimitsDeNoise_DataAbs(IndStep)))];
                else
                    Hist_xvalsDeNoise_DataAbs{IndSubj}  = [Hist_xvalsDeNoise_DataAbs{IndSubj}, Hist_LimitsDeNoise_DataAbs(IndStep-1):Hist_xStepDeNoise_DataAbs(IndStep):Hist_LimitsDeNoise_DataAbs(IndStep)];
                    Hist_xvalGrDeNoise_DataAbs{IndSubj} = [Hist_xvalGrDeNoise_DataAbs{IndSubj},IndStep.*ones(size(Hist_LimitsDeNoise_DataAbs(IndStep-1):Hist_xStepDeNoise_DataAbs(IndStep):Hist_LimitsDeNoise_DataAbs(IndStep)))];
                end
            else
                if(IndStep==1)
                    Hist_xvalsDeNoise_DataAbs{IndSubj}  = [Hist_xvalsDeNoise_DataAbs{IndSubj}, 0:Hist_xStepDeNoise_DataAbs(IndStep):max(abs(DataOrg(:,IndSubj)))];
                    Hist_xvalGrDeNoise_DataAbs{IndSubj} = [Hist_xvalGrDeNoise_DataAbs{IndSubj},IndStep.*ones(size(0:Hist_xStepDeNoise_DataAbs(IndStep):max(abs(DataOrg(:,IndSubj)))))];
                else
                    Hist_xvalsDeNoise_DataAbs{IndSubj}  = [Hist_xvalsDeNoise_DataAbs{IndSubj}, Hist_LimitsDeNoise_DataAbs(IndStep-1):Hist_xStepDeNoise_DataAbs(IndStep):max(abs(DataOrg(:,IndSubj)))];
                    Hist_xvalGrDeNoise_DataAbs{IndSubj} = [Hist_xvalGrDeNoise_DataAbs{IndSubj},IndStep.*ones(size(Hist_LimitsDeNoise_DataAbs(IndStep-1):Hist_xStepDeNoise_DataAbs(IndStep):max(abs(DataOrg(:,IndSubj)))))];
                end
                disp(['IndSub= ',num2str(IndSubj),',IndStep= ',num2str(IndStep),' stopping bincenter creation at maximum. (current stepwidth= ',num2str(Hist_xStepDeNoise_DataAbs(IndStep)),')']);
                break; %we are done!
            end
        else
            Hist_xvalsDeNoise_DataAbs{IndSubj}  = [Hist_xvalsDeNoise_DataAbs{IndSubj}, Hist_LimitsDeNoise_DataAbs(IndStep-1):Hist_xStepDeNoise_DataAbs(IndStep):max(abs(DataOrg(:,IndSubj)))];
            Hist_xvalGrDeNoise_DataAbs{IndSubj} = [Hist_xvalGrDeNoise_DataAbs{IndSubj},IndStep.*ones(size(Hist_LimitsDeNoise_DataAbs(IndStep-1):Hist_xStepDeNoise_DataAbs(IndStep):max(abs(DataOrg(:,IndSubj)))))];
            if(IndStep~=length(Hist_LimitsDeNoise_DataAbs)) %this would be an error, but we can avoid further danger by breaking out of the loop.
                disp('WTF! breaking out of steps-loop.');
                break;
            else
                disp(['IndSub= ',num2str(IndSubj),' stopping bincenter creation at maximum by force.']);
            end
        end
    end
    %order values and remove the redundant ones
    [Hist_xvalsDeNoise_DataAbs{IndSubj},TrafoInds] = unique(Hist_xvalsDeNoise_DataAbs{IndSubj}(:)');
    Hist_xvalGrDeNoise_DataAbs{IndSubj} = Hist_xvalGrDeNoise_DataAbs{IndSubj}(TrafoInds);
    Hist_CountsDeNoise_DataAbs{IndSubj} = hist(DataAbs(:,IndSubj),Hist_xvalsDeNoise_DataAbs{IndSubj});
    
    figure(81); clf;
    subplot(1,2,1); bar(Hist_xvalsDeNoise_DataAbs{IndSubj},Hist_CountsDeNoise_DataAbs{IndSubj}); title(['DataAbs: IndSubj= ',num2str(IndSubj),' (Supp==[',num2str(min(Hist_xvalsDeNoise_DataAbs{IndSubj})),':',num2str(max(Hist_xvalsDeNoise_DataAbs{IndSubj})),'])']); set(gca,'YScale','log');
    subplot(1,2,2); boxplot(Hist_CountsDeNoise_DataAbs{IndSubj}(:),Hist_xvalGrDeNoise_DataAbs{IndSubj}(:)); title(['DataAbs: IndSubj= ',num2str(IndSubj),' (Supp==[',num2str(min(Hist_xvalsDeNoise_DataAbs{IndSubj})),':',num2str(max(Hist_xvalsDeNoise_DataAbs{IndSubj})),'])']); set(gca,'YScale','log');
    
    h = helpdlg(['Subj ',num2str(IndSubj)],'CurrSubj');
    uiwait(h);
end