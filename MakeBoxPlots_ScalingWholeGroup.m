function [DataForPlotting,Handles] = MakeBoxPlots_ScalingWholeGroup(varargin)
% This function makes boxplots for the whole group regarding the rsfMRI scaling-values.
%
% This includes:
%               Boxplot 1: The main masks "EffectOfMRI", IC-Map, WholeBrain and all exclusions
%               Boxplot 2: The split of mask "EffectOfMRI"
%               Boxplot 3: The split of mask IC-Map
% and in each boxplot the subplots(3,1,1:3) OR subplots(1,3,1:3) depending on orientation that is choosen:
%               Subplot(3,1,1): "VOXELS-MEDIAN SUBJECT" - Distribution of Voxel-Values(per Masks) for the MEDIAN SUBJECT. 
%               Subplot(3,1,2): "AGGREGATE DATA"        - Distribution of ALL Voxel-Values(per Masks) of ALL SUBJECTs COLLECTED TOGETHER. 
%               Subplot(3,1,3): "SUBJECTS-MEDIAN VOXEL" - Distribution of Subject-Values(per Masks) for the MEDIAN VOXEL.
%
% The boxplots will be repeated differently in Boxplots 11, 12, and 13 which correspond to Boxplots 1, 2, and 3 respectively.
% The difference in these will be that the user can choose
%       if the median values (2nd quartile==50%) should be replaced by 1st(25%) or 3rd(75%) quartile. 
%       Or if median, 1st quartile and 2nd quartile or any other mix should be collected(aggregated) and plotted together.
%       Or if the mean of that collected mix should be plotted.
%       The default will be to use 1st(25%) and 3rd(75%)-quartile together without averaging(mean) them.
% Note: This means that the AGGREGATE DATA will be the same as in the other box plots but the other subplots change accordingly.
%       This allows a closer inspection of the behavior of scaling values for all voxels and subjects.
%
%Usage:
%      [DataForPlotting,Handles] = MakeBoxPlots_ScalingWholeGroup(varargin);
%      [DataForPlotting,Handles] = MakeBoxPlots_ScalingWholeGroup(DataForPlotting); %automatically plot data from a previous plot configuration, i.e. for a given split, sorting of areas according to reference value of median and so on
%      [DataForPlotting,Handles] = MakeBoxPlots_ScalingWholeGroup();                %ask user to select *.mat-file to load MVSscaling-results struct and then ask for parameters for plotting. All settings and data for plotting will be output in "DataForPlotting"-struct. 
%
%
% NB: scaling-values == magnitude of the fraction of amplitudes for higher field strength(MRI 2) divided by amplitudes for lower field strength(MRI 1).
%     amplitudes     == voxel-wise results from the 2nd-stage of dual-regression analysis of group ICA for the preprocessed data.
%
% NB2: The structure needed for setting up the plots can be made using the function "AnalyseMVSfMRIscaling.m".
%      The function will need the SPM-ANOVA(using the dual-regression data as input) SPM-structure (SPM.mat)
%      and ".filelist" from MELODIC groupICA to create the scaling values in the first place.
%      See the help there.
%
%Date: 16.10.2014 (V2.5(Include log2-scaling; including more sorting options+some plot edits) adapted from V5 of DisplayMVSfMRIscalingResults.m)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% Fix Values
FSizeSplits = 12; %this only applies to vertical plot of the splits
UseFSizeSplitsHor = 1; %if(1) use FSizeSplits also for horizontal boxplot 

%% get MVSscaling struct
if(nargin==0) %in case user has not input the plotting struct
    if(evalin('base','exist(''MVSscalingMAT_Path'',''var'')'))
        MVSscalingMAT_Path = evalin('base','MVSscalingMAT_Path');
        choice_load = questdlg({'Load data OR Rerun previous display?'; ['"',MVSscalingMAT_Path,'"']},'Load OR Rerun','Load','Rerun','Quit','Rerun');
        switch(choice_load)
            case 'Quit'
                return;
            case 'Load'
                MVSscalingMAT_Path = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
            case 'Rerun'
                disp(['Will display data from "',MVSscalingMAT_Path,'"']);
        end
    else
        MVSscalingMAT_Path = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
    end
    assignin('base','MVSscalingMAT_Path',MVSscalingMAT_Path);
    
    %% get/create masks and get data for plotting
    DataForPlotting = PrepDataForPlotting(MVSscalingMAT_Path);
    if(strcmp(DataForPlotting.Choice_split,'Quit'))
        Handles = [];
        return; %Quit
    end
else
    %we got an input, assume it is DataForPlotting and just plot
    DataForPlotting = varargin{1};
end

%% use log2 scaling on Data-axis?
if(strcmp('Yes',questdlg('Scale Y-Axis (i.e. data-axis) with log2?','Y-Scale?','Yes','No','No')))
    UseLog2 = 1;
else
    UseLog2 = 0;
end

%% ylim?
switch(questdlg('Add YLIM?','YLIM?','Yes','No','Yes')) %YLimit in the generalized sense, in case plot-arrangement is changed, we also need to use xlim instead.
    case 'Yes'
        UseYLim = 1;
        answer_YLim = inputdlg({'YLim = '},'Ylim?',1,{'[0 6]'});
        YLimits = sort(eval(answer_YLim{1}));
    otherwise
        UseYLim = 0;
        YLimits = [];
end

%% expected results lines?
switch(questdlg('Use DEFAULT theory lines?','DEFAULT theory?','Yes','No','Yes'))
    case 'Yes'
        TheoryVals = [2*sqrt(2),2,sqrt(2),1];
    otherwise
        answer_Theory = inputdlg({'\lambda MVS-EFFECT = ','\lambda MVS-EFFECT(red) = ','\lambda Constant(fMRI)EFFECT = ','\lambda NoEffect = '},'Theory?',1,{'2*sqrt(2)','2','sqrt(2)','1'});
        TheoryVals(1) = eval(answer_Theory{1}); %MVS-Effect
        TheoryVals(2) = eval(answer_Theory{2}); %MVS-Effect(reduced)
        TheoryVals(3) = eval(answer_Theory{3}); %fMRI-scaling
        TheoryVals(4) = eval(answer_Theory{4}); %No-Effect~random scaling
end
if(UseLog2)
    PlotTheoryVals = log2(TheoryVals);
else
    PlotTheoryVals = TheoryVals;
end

answer_cols = inputdlg({'Color 1: ','Color 2: ','Color 3: ','Color 4: '},'Colors TheoryVals?',1,{'g:','r-','r--','b-'});
if(~isempty(answer_cols))
    ColsLines = answer_cols;
else
    ColsLines = {'g:','r-','r--','b-'};
end

%% boxplot 1: Scaling Values MainMasks subplots: median-subject; aggregate of voxels&subjects; median-voxel
%% Resort data?
choice_resort = questdlg({'Do you want to resort the displayed subregion data according to'; 'the absolute difference of the median from a reference "\lambda"?'; ' '; 'Otherwise the areas are listed according to size(descending)'},'Resort data by abs(median(data)-RefLambda)?','Yes','No','Yes');
switch(choice_resort)
    case 'Yes'
        CheckResort = 1;
    otherwise
        CheckResort = 0;
end

if(CheckResort)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer      = inputdlg({'\lambda-reference= '},'Main Masks',1,{num2str(TheoryVals(1))},options);
    LambdaRef   = eval(answer{1});
    if(UseLog2)
        LambdaRef = log2(LambdaRef);
    end
    
    choice_resort_all = questdlg('Resort EACH or resort all according to MEDIAN-SUBJECT?','Resort EACH or according to MEDIAN-SUBJECT?','EACH','MEDIAN-SUBJECT','EACH');
    switch(choice_resort_all)
        case 'EACH'
            ResortInds = 1:3;
        case 'MEDIAN-SUBJECT'
            ResortInds = 1;
    end
else
    ResortInds= [];
    LambdaRef = [];
end

%% prep masks
MasksForDisp        = DataForPlotting.AllMasks;
MasksNamesForDisp   = DataForPlotting.AllMasksNames;

IndFix = []; %don't keep any fixed if sorting

%% get data & plot
%median-subject
[BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
if(UseLog2)
    BPData = log2(BPData);
end

Handles{1} = figure(1); clf;
subplot(3,1,1); plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
subplot(3,1,1); plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
subplot(3,1,1); plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
subplot(3,1,1); plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['MAIN-MASKS: Scaling-Values for Voxels of MEDIAN-SUBJECT [',InfoStr,']']);
subplot(3,1,1); boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
if(UseYLim)
    ylim(YLimits);
end
if(UseLog2)
    ylabel('log2-scaled');
end

%aggregate-data
if(any(ResortInds==2)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
else
    %according to median-subject if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

subplot(3,1,2); plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
subplot(3,1,2); plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
subplot(3,1,2); plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
subplot(3,1,2); plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['MAIN-MASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
subplot(3,1,2); boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
if(UseYLim)
    ylim(YLimits);
end
if(UseLog2)
    ylabel('log2-scaled');
end


%median-voxel
if(any(ResortInds==3)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianVoxel',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
else
    %according to median-subject
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianVoxel',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

subplot(3,1,3); plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
subplot(3,1,3); plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),'r-'); hold('on');
subplot(3,1,3); plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
subplot(3,1,3); plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['MAIN-MASKS: Scaling-Values of Subjects for MEDIAN-VOXEL [',InfoStr,']']);
subplot(3,1,3); boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
if(UseYLim)
    ylim(YLimits);
end
if(UseLog2)
    ylabel('log2-scaled');
end


clear BPData BPGrouping BPLabels InfoStr IndicesSort %for cleanup

%% ask if a sorting (currently based on median-subject) of either "EffectsOfMRI" or "IC\EffectsOfMRI" should be used to resort the other,
%% as best as possible and resorting the remaining accoring to setting of reference sorting.
Choice_Global_Sort = questdlg('Should a global resorting be used? Should the basis for the global resorting be the MEDIAN-SUBJECT results for "EffectsOfMRI" or "IC\EffectsOfMRI"?','GlobalSorting?','NO global sorting','Use EffectsOfMRI','Use IC\EffectsOfMRI','Use EffectsOfMRI');
switch(Choice_Global_Sort)
    case 'NO global sorting'
        UseGlobalSortSplits = 0;
    case 'Use EffectsOfMRI'
        UseGlobalSortSplits = 1;
        ResortInds = 1;
        %% prep masks                                                       DataForPlotting.EffectOfMRI_split.MasksToBeUsed
        MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)+1);
        MasksNamesForDisp   = cell(length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)+1,1);
        
        MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,3);
        MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{3}; %EffectsOfMRI WHOLE-MASK
        for Ind=1:length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)
            MasksForDisp(:,1+Ind)   = DataForPlotting.EffectOfMRI_split.AllMasks(:,DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp(Ind));
            MasksNamesForDisp{1+Ind}= DataForPlotting.EffectOfMRI_split.AllMasksNames{DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp(Ind)};
        end
        
        IndFix = 1; %fix first when sorting
        
        %% get data 
        if(CheckResort)
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer      = inputdlg({'\lambda-reference= '},'"EffectsOfMRI"',1,{num2str(TheoryVals(1))},options);
            LambdaRef   = eval(answer{1});
            if(UseLog2)
                LambdaRef = log2(LambdaRef);
            end
        else
            %ResortInds= [];
            LambdaRef = [];
        end
        
        %median-subject
        [BPData,BPGrouping,BPLabels] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
        if(UseGlobalSortSplits)
            GlobalSort_Names = BPLabels;
        end
        clear MasksForDisp MasksNamesForDisp IndFix LambdaRef BPData BPGrouping BPLabels %clean up to avoid pain later

    case 'Use IC\EffectsOfMRI'
        UseGlobalSortSplits = 1;
        ResortInds = 1;
        %% prep masks                                                       DataForPlotting.ICmap_split.MasksToBeUsed
        MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)+1);
        MasksNamesForDisp   = cell(length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)+1,1);
        
        MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,5);
        MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{5}; %IC WITHOUT EffectsOfMRI (WHOLE-MASK)
        for Ind=1:length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)
            MasksForDisp(:,1+Ind)   = DataForPlotting.ICmap_split.AllMasks(:,DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp(Ind));
            MasksNamesForDisp{1+Ind}= DataForPlotting.ICmap_split.AllMasksNames{DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp(Ind)};
        end
        
        IndFix = 1; %fix first when sorting
        
        %% get data 
        if(CheckResort)
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer      = inputdlg({'\lambda-reference= '},'"IC-EffectsOfMRI"',1,{num2str(TheoryVals(3))},options);
            LambdaRef   = eval(answer{1});
            if(UseLog2)
                LambdaRef = log2(LambdaRef);
            end
        else
            %ResortInds= [];
            LambdaRef = [];
        end
        [BPData,BPGrouping,BPLabels] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
        if(UseGlobalSortSplits)
            GlobalSort_Names = BPLabels;
        end
        clear MasksForDisp MasksNamesForDisp IndFix LambdaRef BPData BPGrouping BPLabels %clean up to avoid pain later
end

%% for splits: Boxplots Vertical or Horizontal? --> changes also subplots
choice_plotarr = questdlg({'For the splits it could be advantageous to use a different plot style than "Horizontal"(default-boxplot).'; 'Should the boxplots&subplots be done "Horizontal"(Default) OR "Vertical"?'; 'I.e. where should the labels be?'},'Plot-arrangement?','Horizontal','Vertical','Horizontal');
switch(choice_plotarr) %NB: boxplot function notes arrangement from the point of the data, i.e. opposite to our scheme of where to put the labels, therefore this switch.
    case 'Horizontal'
        choice_plotarr = 'vertical'; %default
    case 'Vertical'
        choice_plotarr = 'horizontal';
    otherwise
        choice_plotarr = 'vertical'; %default
end

%% boxplot 2: Scaling Values split "EffectsOfMRI" subplots: median-subject; aggregate of voxels&subjects; median-voxel
%% Resort data?
if(CheckResort)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer      = inputdlg({'\lambda-reference= '},'"EffectsOfMRI"',1,{num2str(TheoryVals(1))},options);
    LambdaRef   = eval(answer{1});
    if(UseLog2)
        LambdaRef = log2(LambdaRef);
    end

    
    if(~UseGlobalSortSplits)
        choice_resort_all = questdlg('Resort EACH or resort all according to MEDIAN-SUBJECT?','Resort EACH or according to MEDIAN-SUBJECT?','EACH','MEDIAN-SUBJECT','EACH');
        switch(choice_resort_all)
            case 'EACH'
                ResortInds = 1:3;
            case 'MEDIAN-SUBJECT'
                ResortInds = 1;
        end
    end
else
    ResortInds= [];
    LambdaRef = [];
end

%% prep masks                                                       DataForPlotting.EffectOfMRI_split.MasksToBeUsed
MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)+1);
MasksNamesForDisp   = cell(length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)+1,1);

MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,3);
MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{3}; %EffectsOfMRI WHOLE-MASK
for Ind=1:length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)
    MasksForDisp(:,1+Ind)   = DataForPlotting.EffectOfMRI_split.AllMasks(:,DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp(Ind));
    MasksNamesForDisp{1+Ind}= DataForPlotting.EffectOfMRI_split.AllMasksNames{DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp(Ind)};
end

IndFix = 1; %fix first when sorting

%% get data & plot
%median-subject
if(UseGlobalSortSplits)
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,GlobalSort_Names);
else
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
end
if(UseLog2)
    BPData = log2(BPData);
end

Handles{2} = figure(2); clf;
if(strcmpi(choice_plotarr,'vertical')) %default
    subplot(3,1,1);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of MEDIAN-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseFSizeSplitsHor)
        set(gca,'XTick',1:length(BPLabels),'XTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    end
    if(UseYLim)
        ylim(YLimits);
    end
else %'horizontal'
    subplot(1,3,1);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of MEDIAN-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
end
if(UseLog2)
    ylabel('log2-scaled');
end


%aggregate-data
if(any(ResortInds==2)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
else
    %according to median-subject if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,2); 
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseFSizeSplitsHor)
        set(gca,'XTick',1:length(BPLabels),'XTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    end
    if(UseYLim)
        ylim(YLimits);
    end
else %'horizontal'
    subplot(1,3,2);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
%     set(gca,'YTickLabel',{' '});
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    set(findobj(gca,'Type','text'),'FontSize',FSizeSplits)
end
if(UseLog2)
    ylabel('log2-scaled');
end



%median-voxel
if(any(ResortInds==3)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianVoxel',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
else
    %according to median-subject
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianVoxel',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,3);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for MEDIAN-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr,'labels',BPLabels); hold('on');
    if(UseFSizeSplitsHor)
        set(gca,'XTick',1:length(BPLabels),'XTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    end
    if(UseYLim)
        ylim(YLimits);
    end
else %'horizontal'
    subplot(1,3,3);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for MEDIAN-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
%     set(gca,'YTickLabel',{' '});
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    set(findobj(gca,'Type','text'),'FontSize',FSizeSplits)
end
if(UseLog2)
    ylabel('log2-scaled');
end

clear BPData BPGrouping BPLabels InfoStr IndicesSort %for cleanup


%% boxplot 3: Scaling Values split "IC-EffectsOfMRI" subplots: median-subject; aggregate of voxels&subjects; median-voxel
%% Resort data?
if(CheckResort)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer      = inputdlg({'\lambda-reference= '},'"IC-EffectsOfMRI"',1,{num2str(TheoryVals(3))},options);
    LambdaRef   = eval(answer{1});
    if(UseLog2)
        LambdaRef = log2(LambdaRef);
    end
    
    if(UseGlobalSortSplits)
        ResortInds = 1;
    else
        choice_resort_all = questdlg('Resort EACH or resort all according to MEDIAN-SUBJECT?','Resort EACH or according to MEDIAN-SUBJECT?','EACH','MEDIAN-SUBJECT','EACH');
        switch(choice_resort_all)
            case 'EACH'
                ResortInds = 1:3;
            case 'MEDIAN-SUBJECT'
                ResortInds = 1;
        end
    end
else
    ResortInds= [];
    LambdaRef = [];
end

%% prep masks                                                       DataForPlotting.ICmap_split.MasksToBeUsed
MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)+1);
MasksNamesForDisp   = cell(length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)+1,1);

MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,5);
MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{5}; %IC WITHOUT EffectsOfMRI (WHOLE-MASK)
for Ind=1:length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)
    MasksForDisp(:,1+Ind)   = DataForPlotting.ICmap_split.AllMasks(:,DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp(Ind));
    MasksNamesForDisp{1+Ind}= DataForPlotting.ICmap_split.AllMasksNames{DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp(Ind)};
end

IndFix = 1; %fix first when sorting

%% get data & plot
%median-subject
if(UseGlobalSortSplits)
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,GlobalSort_Names);
else
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianSubj',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
end
if(UseLog2)
    BPData = log2(BPData);
end

Handles{3} = figure(3); clf;
if(strcmpi(choice_plotarr,'vertical'))  %default (i.e. what should be called HORIZONTAL!!!)
    subplot(3,1,1);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of MEDIAN-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseFSizeSplitsHor)
        set(gca,'XTick',1:length(BPLabels),'XTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    end
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,1);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of MEDIAN-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
    set(findobj(gca,'Type','text'),'FontSize',FSizeSplits)
end



%aggregate-data
if(any(ResortInds==2)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
else
    %according to median-subject if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,2);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseFSizeSplitsHor)
        set(gca,'XTick',1:length(BPLabels),'XTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    end
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,2);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
    set(findobj(gca,'Type','text'),'FontSize',FSizeSplits)
end


%median-voxel
if(any(ResortInds==3)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianVoxel',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix);
else
    %according to median-subject
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'MedianVoxel',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,3);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for MEDIAN-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseFSizeSplitsHor)
        set(gca,'XTick',1:length(BPLabels),'XTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    end
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,3);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for MEDIAN-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
    set(findobj(gca,'Type','text'),'FontSize',FSizeSplits)
end


clear BPData BPGrouping BPLabels InfoStr IndicesSort %for cleanup


%% boxplot 11, 12, and 13: Ask user which values to mix or average
%which quartiles should be used for plotting?
ListChoice = {'1stQuartile(25%)';'2ndQuartile(50%,i.e. median)';'3rdQuartile(75%)'};
[MixIndices,ok] = listdlg('ListString',ListChoice,'InitialValue',[1 3],'PromptString','Select Quartiles for plotting.');
if(~ok)
    disp('Done (plotting stopped by user)');
    disp(' ');
    return;
end
        
%should these be averaged(mean)? if(length(QuartInds)>1)
if(length(MixIndices)>1)
    ChoiceAvMix = questdlg('Average(mean) the quartiles?','Average?','Yes','No','No');
    switch(ChoiceAvMix)
        case 'Yes'
            AverageMix = 1;
        otherwise
            AverageMix = 0;
    end
else
    AverageMix = 0;
end
InfoStrMix = [];
for Ind = 1:length(MixIndices)
    if(Ind~=length(MixIndices))
        InfoStrMix = [InfoStrMix,ListChoice{MixIndices(Ind)},'&'];
    else
        InfoStrMix = [InfoStrMix,ListChoice{MixIndices(Ind)},'"'];
    end
end
if(AverageMix)
    InfoStrMix = ['"Average(mean) of ',InfoStrMix];
else
    InfoStrMix = ['"Aggregate of ',InfoStrMix];
end


%% boxplot 11: Scaling Values MainMasks subplots: median-subject; aggregate of voxels&subjects; median-voxel
%% Resort data?
choice_resort = questdlg({'Do you want to resort the displayed subregion data according to'; 'the absolute difference of the median from a reference "\lambda"?'; ' '; 'Otherwise the areas are listed according to size(descending)'},'Resort data by abs(median(data)-RefLambda)?','Yes','No','Yes');
switch(choice_resort)
    case 'Yes'
        CheckResort = 1;
    otherwise
        CheckResort = 0;
end

if(CheckResort)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer      = inputdlg({'\lambda-reference= '},'Main Masks',1,{num2str(TheoryVals(1))},options);
    LambdaRef   = eval(answer{1});
    if(UseLog2)
        LambdaRef = log2(LambdaRef);
    end
    
    choice_resort_all = questdlg('Resort EACH or resort all according to MEDIAN-SUBJECT?','Resort EACH or according to MEDIAN-SUBJECT?','EACH','MEDIAN-SUBJECT','EACH');
    switch(choice_resort_all)
        case 'EACH'
            ResortInds = 1:3;
        case 'MEDIAN-SUBJECT'
            ResortInds = 1;
    end
else
    ResortInds= [];
    LambdaRef = [];
end

%% prep masks
MasksForDisp        = DataForPlotting.AllMasks;
MasksNamesForDisp   = DataForPlotting.AllMasksNames;

IndFix = []; %don't keep any fixed if sorting

%% get data & plot
%median-subject
[BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesSubjMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
if(UseLog2)
    BPData = log2(BPData);
end

Handles{4} = figure(11); clf;
subplot(3,1,1); 
plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['MAIN-MASKS: Scaling-Values for Voxels of ',InfoStrMix,'-SUBJECT [',InfoStr,']']);
boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
if(UseYLim)
    ylim(YLimits);
end
if(UseLog2)
    ylabel('log2-scaled');
end

%aggregate-data
if(any(ResortInds==2)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
else
    %according to median-subject if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort,MixIndices,AverageMix);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

subplot(3,1,2); 
plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['MAIN-MASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
if(UseYLim)
    ylim(YLimits);
end
if(UseLog2)
    ylabel('log2-scaled');
end

%median-voxel
if(any(ResortInds==3)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesVoxelMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
else
    %according to median-subject
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesVoxelMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort,MixIndices,AverageMix);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

subplot(3,1,3); 
plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['MAIN-MASKS: Scaling-Values of Subjects for ',InfoStrMix,'-VOXEL [',InfoStr,']']);
boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
if(UseYLim)
    ylim(YLimits);
end
if(UseLog2)
    ylabel('log2-scaled');
end

clear BPData BPGrouping BPLabels InfoStr IndicesSort %for cleanup

%% boxplot 12: Scaling Values split "EffectsOfMRI" subplots: median-subject; aggregate of voxels&subjects; median-voxel
%% Resort data?
if(CheckResort)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer      = inputdlg({'\lambda-reference= '},'"EffectsOfMRI"',1,{num2str(TheoryVals(1))},options);
    LambdaRef   = eval(answer{1});
    if(UseLog2)
        LambdaRef = log2(LambdaRef);
    end
    
    if(UseGlobalSortSplits)
        ResortInds = 1;
    else
        choice_resort_all = questdlg('Resort EACH or resort all according to MEDIAN-SUBJECT?','Resort EACH or according to MEDIAN-SUBJECT?','EACH','MEDIAN-SUBJECT','EACH');
        switch(choice_resort_all)
            case 'EACH'
                ResortInds = 1:3;
            case 'MEDIAN-SUBJECT'
                ResortInds = 1;
        end
    end
else
    ResortInds= [];
    LambdaRef = [];
end

%% prep masks                                                       DataForPlotting.EffectOfMRI_split.MasksToBeUsed
MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)+1);
MasksNamesForDisp   = cell(length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)+1,1);

MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,3);
MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{3}; %EffectsOfMRI WHOLE-MASK
for Ind=1:length(DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp)
    MasksForDisp(:,1+Ind)   = DataForPlotting.EffectOfMRI_split.AllMasks(:,DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp(Ind));
    MasksNamesForDisp{1+Ind}= DataForPlotting.EffectOfMRI_split.AllMasksNames{DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp(Ind)};
end

IndFix = 1; %fix first when sorting

%% get data & plot
%median-subject
if(UseGlobalSortSplits)
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesSubjMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,GlobalSort_Names,MixIndices,AverageMix);
else
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesSubjMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
end
if(UseLog2)
    BPData = log2(BPData);
end

Handles{5} = figure(12); clf;
if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,1);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of ',InfoStrMix,'-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,1);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of ',InfoStrMix,'-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
end


%aggregate-data
if(any(ResortInds==2)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
else
    %according to median-subject if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort,MixIndices,AverageMix);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,2);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,2);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
end


%median-voxel
if(any(ResortInds==3)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesVoxelMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
else
    %according to median-subject
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesVoxelMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort,MixIndices,AverageMix);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,3);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for ',InfoStrMix,'-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,3);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['"EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for ',InfoStrMix,'-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
end


clear BPData BPGrouping BPLabels InfoStr IndicesSort %for cleanup


%% boxplot 13: Scaling Values split "IC-EffectsOfMRI" subplots: median-subject; aggregate of voxels&subjects; median-voxel
%% Resort data?
if(CheckResort)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer      = inputdlg({'\lambda-reference= '},'"IC-EffectsOfMRI"',1,{num2str(TheoryVals(3))},options);
    LambdaRef   = eval(answer{1});
    if(UseLog2)
        LambdaRef = log2(LambdaRef);
    end
    
    if(UseGlobalSortSplits)
        ResortInds = 1;
    else
        choice_resort_all = questdlg('Resort EACH or resort all according to MEDIAN-SUBJECT?','Resort EACH or according to MEDIAN-SUBJECT?','EACH','MEDIAN-SUBJECT','EACH');
        switch(choice_resort_all)
            case 'EACH'
                ResortInds = 1:3;
            case 'MEDIAN-SUBJECT'
                ResortInds = 1;
        end
    end
else
    ResortInds= [];
    LambdaRef = [];
end

%% prep masks                                                       DataForPlotting.ICmap_split.MasksToBeUsed
MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)+1);
MasksNamesForDisp   = cell(length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)+1,1);

MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,5);
MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{5}; %IC WITHOUT EffectsOfMRI (WHOLE-MASK)
for Ind=1:length(DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp)
    MasksForDisp(:,1+Ind)   = DataForPlotting.ICmap_split.AllMasks(:,DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp(Ind));
    MasksNamesForDisp{1+Ind}= DataForPlotting.ICmap_split.AllMasksNames{DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp(Ind)};
end

IndFix = 1; %fix first when sorting

%% get data & plot
%median-subject
if(UseGlobalSortSplits)
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesSubjMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,GlobalSort_Names,MixIndices,AverageMix);
else
    [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesSubjMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
end
if(UseLog2)
    BPData = log2(BPData);
end

Handles{6} = figure(13); clf;
if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,1);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of ',InfoStrMix,'-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,1);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for Voxels of ',InfoStrMix,'-SUBJECT [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
end


%aggregate-data
if(any(ResortInds==2)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
else
    %according to median-subject if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'Aggregate',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort,MixIndices,AverageMix);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,2);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');    
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,2);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values for AGGREGATE-DATA(Voxels&Subjects) [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
end

%median-voxel
if(any(ResortInds==3)) %resort each if CheckResort==1
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesVoxelMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,[],MixIndices,AverageMix);
else
    %according to median-subject
    [BPData,BPGrouping,BPLabels,InfoStr] = AssembleBPinputs_n_Resort(DataForPlotting.ScalingData,'QuartilesVoxelMix',MasksNamesForDisp,MasksForDisp,CheckResort,LambdaRef,IndFix,IndicesSort,MixIndices,AverageMix);
    if(CheckResort)
        InfoStr = ['FromMedianSubject-',InfoStr];
    end
end
if(UseLog2)
    BPData = log2(BPData);
end

if(strcmpi(choice_plotarr,'vertical'))  %default
    subplot(3,1,3);
    plot(1:length(BPLabels),PlotTheoryVals(1).*ones(length(BPLabels),1),ColsLines{1},'LineWidth',2); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(2).*ones(length(BPLabels),1),ColsLines{2}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(3).*ones(length(BPLabels),1),ColsLines{3}); hold('on');
    plot(1:length(BPLabels),PlotTheoryVals(4).*ones(length(BPLabels),1),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for ',InfoStrMix,'-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    if(UseYLim)
        ylim(YLimits);
    end
    if(UseLog2)
        ylabel('log2-scaled');
    end
else %'horizontal'
    subplot(1,3,3);
    plot(PlotTheoryVals(1).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{1},'LineWidth',2); hold('on');
    plot(PlotTheoryVals(2).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{2}); hold('on');
    plot(PlotTheoryVals(3).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{3}); hold('on');
    plot(PlotTheoryVals(4).*ones(length(BPLabels),1),1:length(BPLabels),ColsLines{4}); hold('on'); title(['ICmap WO "EffectsOfMRI"-SplitMASKS: Scaling-Values of Subjects for ',InfoStrMix,'-VOXEL [',InfoStr,']']);
    boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on','orientation',choice_plotarr); hold('on');
    set(gca,'YTick',1:length(BPLabels),'YTickLabel',BPLabels,'FontSize',FSizeSplits); %this fixes a bug in matlab boxplot function
    if(UseYLim)
        xlim(YLimits);
    end
    if(UseLog2)
        xlabel('log2-scaled');
    end
end


clear BPData BPGrouping BPLabels InfoStr IndicesSort %for cleanup


%% Done
disp('Done');
disp(' ');

end

%% Subfunctions
function DataForPlotting = PrepDataForPlotting(MVSscalingMAT_Path)
% This function assembles data for plotting 
% (Sorting according to median value of distribution is done by assemble function.)

%% load MVSscaling-struct
load(MVSscalingMAT_Path);

%% magnitude of scaling data input --> scaling data for plotting
if(1) %use magnitude of the fraction of amplitudes for MRI 2 to MRI 1
DataForPlotting.ScalingData = abs(MVSscaling.ScalingPerSubject);
else
    DataForPlotting.ScalingData = (MVSscaling.ScalingPerSubject);
end

%% get main masks
DataForPlotting.AllMasks      = MVSscaling.Masks.AllMasks;
DataForPlotting.AllMasksNames = MVSscaling.Masks.AllMasksName;

%% Do you want to look at the AAL-splits or ROI-splits
Choice_split = questdlg({'Which split do you want to look at?'; ' '; 'ROI-splits or AAL-splits?'},'ROI-splits or AAL-splits?','ROI-splits','AAL-splits','Quit','ROI-splits');
DataForPlotting.Choice_split = Choice_split;
switch(Choice_split)
    case 'ROI-splits'
        %% assign split masks
        DataForPlotting.EffectOfMRI_split.AllMasks      = MVSscaling.Masks.ROIsplit.EffectOfMRI.AllMasks;
        DataForPlotting.EffectOfMRI_split.AllMasksNames = MVSscaling.Masks.ROIsplit.EffectOfMRI.Labels;
        DataForPlotting.EffectOfMRI_split.NVoxelAllMasks= MVSscaling.Masks.ROIsplit.EffectOfMRI.NVoxelAllMasks;
        DataForPlotting.EffectOfMRI_split.MasksToBeUsed.IndsMasksDisp = 1:size(MVSscaling.Masks.ROIsplit.EffectOfMRI.AllMasks,2); %NThres; PThres; IndsMasksDisp; InfoStr;
        
        %% assign split masks
        DataForPlotting.ICmap_split.AllMasks      = MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.AllMasks;
        DataForPlotting.ICmap_split.AllMasksNames = MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.Labels;
        DataForPlotting.ICmap_split.NVoxelAllMasks= MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.NVoxelAllMasks;
        DataForPlotting.ICmap_split.MasksToBeUsed.IndsMasksDisp = 1:size(MVSscaling.Masks.ROIsplit.IC_WITHOUT_EffectOfMRI.AllMasks,2); %NThres; PThres; IndsMasksDisp; InfoStr;
        
    case 'AAL-splits'
        %% decide on the split and exclusion due to voxel number after merging (if applicable)
        %% pick a split of mask "EffectsOfMRI" to use for plotting
        H = helpdlg('Select split: EffectsOfMRI','Select split: EffectsOfMRI');
        uiwait(H);
        [IndSplitEffectsOfMRI,ok] = listdlg('ListString',MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames,'SelectionMode','single','Name','Select split: EffectsOfMRI');
        if(~ok)
            IndSplitEffectsOfMRI = 1; %default to one if none is picked.
        end
        
        %% get masks for split and ask for merging if necessary
        [EffectOfMRI_split] = CombineSubparts(MVSscaling.Masks.AALsplit.EffectOfMRI,IndSplitEffectsOfMRI,'EffectsOfMRI');
        
        %% assign split masks
        DataForPlotting.EffectOfMRI_split.IndSplit      = IndSplitEffectsOfMRI;
        DataForPlotting.EffectOfMRI_split.SplitName     = MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{IndSplitEffectsOfMRI};
        DataForPlotting.EffectOfMRI_split.AllMasks      = EffectOfMRI_split.AllMasks;
        DataForPlotting.EffectOfMRI_split.AllMasksNames = EffectOfMRI_split.AllMasksNames;
        DataForPlotting.EffectOfMRI_split.NVoxelAllMasks= EffectOfMRI_split.NVoxelAllMasks;
        DataForPlotting.EffectOfMRI_split.MasksToBeUsed = EffectOfMRI_split.MasksToBeUsed; %NThres; PThres; IndsMasksDisp; InfoStr;
        
        
        %% pick a split of mask ICmap WITHOUT "EffectsOfMRI" to use for plotting
        H = helpdlg('Select split: IC-Map WO EffectsOfMRI','Select split: IC-Map WO EffectsOfMRI');
        uiwait(H);
        [IndSplitIC,ok] = listdlg('ListString',MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames,'SelectionMode','single','Name','Select split: IC-Map');
        if(~ok)
            IndSplitIC = 1; %default to one if none is picked.
        end
        
        %% get masks for split and ask for merging if necessary
        [ICmap_split] = CombineSubparts(MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI,IndSplitIC,'ICMap-EffectsOfMRI');
        
        %% assign split masks
        DataForPlotting.ICmap_split.IndSplit      = IndSplitEffectsOfMRI;
        DataForPlotting.ICmap_split.SplitName     = MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames{IndSplitIC};
        DataForPlotting.ICmap_split.AllMasks      = ICmap_split.AllMasks;
        DataForPlotting.ICmap_split.AllMasksNames = ICmap_split.AllMasksNames;
        DataForPlotting.ICmap_split.NVoxelAllMasks= ICmap_split.NVoxelAllMasks;
        DataForPlotting.ICmap_split.MasksToBeUsed = ICmap_split.MasksToBeUsed; %NThres; PThres; IndsMasksDisp; InfoStr;
    case 'Quit'
        return;
end

end


function [OutputStruct] = CombineSubparts(MaskSubpartSplitStruct,IndSplitToBeUsed,MaskName)
% This function checks if subparts need to be mergered and allows user to
% do this manually and suggest mergers and even exclude subparts due to
% their size.

%% ask user if subparts should be combined.
% DataForPlotting.EffectOfMRI_split.AllMasks      = MVSscaling.Masks.AALsplit.EffectOfMRI.AllMasks;
% DataForPlotting.EffectOfMRI_split.AllMasksNames = MVSscaling.Masks.AALsplit.EffectOfMRI.AllMasksName;
% DataForPlotting.EffectOfMRI_split.NVoxelAllMasks= MVSscaling.Masks.AALsplit.EffectOfMRI.NVoxelAllMasks;
NMaskMax = 5; %let's just say 5 is a lot because then the figure will have six boxes, so let's at least ask if user wants to combine masks
if(length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)>=NMaskMax) 
    InfoDlgTxt = cell(length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)+6,1);
    InfoDlgTxt{1} = ['"',MaskName,'" has been split(',MaskSubpartSplitStruct.SplitNames{IndSplitToBeUsed},') into ',num2str(length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)),' subparts,'];
    InfoDlgTxt{2} = ['This might be a bit too much to look at at once. '];
    InfoDlgTxt{3} = ['If you want to combine some of them, click "combine", otherwise "skip".'];
    InfoDlgTxt{4} = ['[NB: You can still exclude small subparts from display in a later step instead of combining them with others.]'];
    InfoDlgTxt{5} =  ' ';
    InfoDlgTxt{6} = ['These are all the included subparts and their sizes: '];
    for Ind = 1:length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)
        InfoDlgTxt{6+Ind} = [num2str(Ind),'.',MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{Ind},': ',num2str(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)./sum(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks).*100),'%(',num2str(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)),' Voxels) of the whole Mask(',num2str(sum(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks)),' Voxels)'];
    end
    
    choice = questdlg(InfoDlgTxt,'Combine subparts "EffectOfMRI"?','combine','skip','skip');
    switch(choice)
        case 'combine'
            %% Inform user that mergers can be choosen as far as possible
            %% and then an automatic algorithm will try to combine the rest.
            InfoManualSuggestionTxt{ 1} = 'You should try to make some manual selection';
            InfoManualSuggestionTxt{ 2} = 'of mergers or selections of single areas.';
            InfoManualSuggestionTxt{ 3} = ' ';
            InfoManualSuggestionTxt{ 4} = 'This will help the algorithm that comes next,';
            InfoManualSuggestionTxt{ 5} = 'which can find possible merging of the left-over areas.';
            InfoManualSuggestionTxt{ 6} = ' ';
            InfoManualSuggestionTxt{ 7} = 'Usage: ';
            InfoManualSuggestionTxt{ 8} = '      Just choose one or multiple areas to be merged.';
            InfoManualSuggestionTxt{ 9} = '      If you are sufficiently satisfied you can quit';
            InfoManualSuggestionTxt{10} = '      and let the automatic algorithm try to find';
            InfoManualSuggestionTxt{11} = '      some more mergers for the remaining areas.';
                        
            H_InfoManualSuggestion = helpdlg(InfoManualSuggestionTxt,'Help for ManualSuggestMergers');
            uiwait(H_InfoManualSuggestion,2);
            
            %% do manual merging
            SuggestionsInit = ManualMergers(MaskSubpartSplitStruct,IndSplitToBeUsed);
            
            %% Inform user about inputs for suggestion
            InfoSuggestionTxt{ 1} = 'You will be asked next to define the Limits for the merging of Subparts/Areas of the AAL-split.';
            InfoSuggestionTxt{ 2} = ' ';
            InfoSuggestionTxt{ 3} = 'Note that too excessive merging may compromise judgement of results!';
            InfoSuggestionTxt{ 4} = ' ';
            InfoSuggestionTxt{ 5} = 'You can make the following Inputs: ';
            InfoSuggestionTxt{ 6} = '             1.) The search range in number of neighbors per voxel.';
            InfoSuggestionTxt{ 7} = '                 The best choice is usually nearest neightbors, i.e. N=1.';
            InfoSuggestionTxt{ 8} = '          2.&3.) The other two parameters prevent runaway effects.';
            InfoSuggestionTxt{ 9} = '                 If they are reached a stop of merging is performed.';
            InfoSuggestionTxt{10} = '                 Best choice are N=8 Areas & N=100 Voxels maximum.';
            InfoSuggestionTxt{11} = '                 [NB: if NAreas<=2 ignore NVoxelMax, to allow one merge.]';
            InfoSuggestionTxt{12} = ' ';
            InfoSuggestionTxt{13} = '   HINT   HINT   HINT:     KEEP THIS WINDOW OPEN     :HINT   HINT   HINT   ';
            
            H_InfoSuggestion = helpdlg(InfoSuggestionTxt,'Help for SuggestMergers');
            uiwait(H_InfoSuggestion,2);
            
            %% SearchRange: ask user for input
            prompt   ={'SearchRange3D[VoxelNeighbors]: ','NAreaMaxPerMerger: ','NVoxelMaxPerMerger: '};
            dlg_title='Limits for combining subparts';
            num_lines=1;
            defAns   ={'1','8','100'};
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer   = inputdlg(prompt,dlg_title,num_lines,defAns,options);
            clear prompt dlg_title num_lines defAns options %clean up or consequences are bad!
            SearchRange = eval(answer{1}); %how many neighborhoods to search
            NAreasMax   = eval(answer{2}); %how many areas a maximally allowed per merger
            NVoxelsMax  = eval(answer{3}); %maximum number of voxels in a merger
            %% suggest subparts to be merged
            Suggestions = SuggestMergers(MaskSubpartSplitStruct,IndSplitToBeUsed,SearchRange,NAreasMax,NVoxelsMax,SuggestionsInit);
                        
            %% show user which subparts are available
            InfoDlgTxt = cell(length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)+3,1);
            InfoDlgTxt{1} = ['The following subparts are available: '];
            for Ind = 1:length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)
                InfoDlgTxt{1+Ind} = [num2str(Ind),'.',MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{Ind},': ',num2str(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)./sum(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks).*100),'%(',num2str(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)),' Voxels) of the whole Mask(',num2str(sum(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks)),' Voxels)'];
            end
            InfoDlgTxt{1+length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)+1} = ' ';
            InfoDlgTxt{1+length(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName)+2} = '       HINT:   KEEP THIS WINDOW OPEN   :HINT';
            h_subparts=helpdlg(InfoDlgTxt,'Available subparts');
            uiwait(h_subparts,1.5); %wait a second and a half
            
            %% ask user for input using suggestions
            prompt   ={['Subparts Combinations: ']};
            dlg_title='List of subparts to combine';
            num_lines=1;
            if(Suggestions.Info.NAreasBeforeMerge>=Suggestions.Info.NAreasAfterMerge) %suggest has worked out without fatal error put a list together
                ChoiceStr = []; %init
                for IndMerger = 1:length(Suggestions.IndsMergers)
                    for IndIndex = 1:length(Suggestions.IndsMergers{IndMerger})
                        if(IndIndex==1)
                            if(IndIndex~=length(Suggestions.IndsMergers{IndMerger}))
                                ChoiceStr = [ChoiceStr,'[',num2str(Suggestions.IndsMergers{IndMerger}(IndIndex)),', ']; %first subpart of this merger
                            else
                                ChoiceStr = [ChoiceStr,'[',num2str(Suggestions.IndsMergers{IndMerger}(IndIndex)),']']; %only one subpart, i.e. finish
                            end
                        else
                            if(IndIndex~=length(Suggestions.IndsMergers{IndMerger}))
                                ChoiceStr = [ChoiceStr,num2str(Suggestions.IndsMergers{IndMerger}(IndIndex)),', ']; %more subparts to follow for this merger
                            else
                                ChoiceStr = [ChoiceStr,num2str(Suggestions.IndsMergers{IndMerger}(IndIndex)),']']; %last subpart of this merger, i.e. finish
                            end
                        end
                    end
                    if(IndMerger~=length(Suggestions.IndsMergers))
                        ChoiceStr = [ChoiceStr,'; ']; %there are further mergers to come so let's put a line break/next line in cell
                    end
                end %done with adding mergers
                defAns ={['{',ChoiceStr,'}']}; %put curly brakets around it for eval to make it into a cell.
            else %leave choice to user
                defAns ={'{[Choice1, ..., ChoiceN]; [Choice2, ..., ChoiceM]; Choice42}'};
            end
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer   = inputdlg(prompt,dlg_title,num_lines,defAns,options);
            clear prompt dlg_title num_lines defAns options %clean up or consequences are bad!
            SubpartsCombine = eval(answer{1});
            
            AllMasks      = zeros(size(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasks,1),length(SubpartsCombine));
            AllMasksNames = cell(length(SubpartsCombine),1);
            NVoxelAllMasks= zeros(length(SubpartsCombine),1);
            for Ind = 1:length(SubpartsCombine)
                tmpMask     = zeros(size(AllMasks,1),1);
                tmpMaskName = [];
                for IndComb = 1:length(SubpartsCombine{Ind})
                    IndMask = find(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).LabelInds==SubpartsCombine{Ind}(IndComb));
                    tmpMask = tmpMask+MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasks(:,IndMask);
                    if(IndComb==1)
                        tmpMaskName = [tmpMaskName,MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                    else
                        tmpMaskName = [tmpMaskName,'+',MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                    end
                end
                AllMasks(:,Ind)    = (tmpMask>0); %just to be sure include logical for failsafe.
                AllMasksNames{Ind} = tmpMaskName;
                NVoxelAllMasks(Ind)= length(find(AllMasks(:,Ind)));
            end
            InfoCombDlgTxt = cell(length(SubpartsCombine)+1,1);
            InfoCombDlgTxt{1} = ['Will combine the following subparts: '];
            for Ind = 1:length(SubpartsCombine)
                CombTxt = [];
                for IndComb = 1:length(SubpartsCombine{Ind})
                    IndMask = find(MaskSubpartSplitStruct.Split(IndSplitToBeUsed).LabelInds==SubpartsCombine{Ind}(IndComb));
                    if(IndComb==1)
                        CombTxt = [CombTxt,num2str(Ind),'.(',num2str(SubpartsCombine{Ind}),'): "',MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                    else
                        if(IndComb~=length(SubpartsCombine{Ind}))
                            CombTxt = [CombTxt,'+',MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                        else
                            CombTxt = [CombTxt,'+',MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName{IndMask},'" (',num2str(NVoxelAllMasks(Ind)),'Voxels of ',num2str(sum(NVoxelAllMasks)),'Voxels in total.)'];
                        end
                    end
                end
                InfoCombDlgTxt{1+Ind} = CombTxt;
            end
            h_comb=helpdlg(InfoCombDlgTxt,'Combined subparts'); 
        otherwise
            AllMasks      = MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasks;
            AllMasksNames = MaskSubpartSplitStruct.Split(IndSplitToBeUsed).AllMasksName;
            NVoxelAllMasks= MaskSubpartSplitStruct.Split(IndSplitToBeUsed).NVoxelAllMasks;
    end
end

%% assemble outputs
OutputStruct.AllMasks      = AllMasks;
OutputStruct.AllMasksNames = AllMasksNames;
OutputStruct.NVoxelAllMasks= NVoxelAllMasks;

%% exclude subparts that are too small?
%check if there are many splits with only few voxels i.e. <5% or 1% of total number of voxels
%OR if the maximum subpart/submask is smaller than 10%
%OR if there are just more than 5 subparts.
if((length(find((OutputStruct.NVoxelAllMasks./sum(OutputStruct.NVoxelAllMasks))<0.05))>6) || (length(find((OutputStruct.NVoxelAllMasks./sum(OutputStruct.NVoxelAllMasks))<0.01))>3) || (max(OutputStruct.NVoxelAllMasks(:))./sum(OutputStruct.NVoxelAllMasks)<0.1) || length(OutputStruct.NVoxelAllMasks)>5)
    choice = questdlg({['There are many subparts of the Mask "',MaskName,'" that are small!']; ...
        [num2str(length(find((OutputStruct.NVoxelAllMasks./sum(OutputStruct.NVoxelAllMasks))<0.05))),'of',num2str(length(OutputStruct.NVoxelAllMasks)),' are smaller than 5%(',num2str(round(sum(OutputStruct.NVoxelAllMasks).*0.05)),'Voxels) of the whole Mask(',num2str(sum(OutputStruct.NVoxelAllMasks)),'Voxels).']; ...
        [num2str(length(find((OutputStruct.NVoxelAllMasks./sum(OutputStruct.NVoxelAllMasks))<0.01))),'of',num2str(length(OutputStruct.NVoxelAllMasks)),' are smaller than 1%(',num2str(round(sum(OutputStruct.NVoxelAllMasks).*0.01)),'Voxels) of the whole Mask(',num2str(sum(OutputStruct.NVoxelAllMasks)),'Voxels).']; ...
        ['The biggest subpart is ',num2str(max(OutputStruct.NVoxelAllMasks(:))./sum(OutputStruct.NVoxelAllMasks).*100),'%(',num2str(max(OutputStruct.NVoxelAllMasks(:))),'Voxels) of the whole Mask(',num2str(sum(OutputStruct.NVoxelAllMasks)),'Voxels).']; ...
        ' ';
        'Do you want to limit the minimum size of the subparts that are displayed?'; ...
        'Options are: '; ...
        'Yes --> "enter percentage cutoff".';...
        ['T<5% --> "Exclude subparts smaller 5%(',num2str(round(sum(OutputStruct.NVoxelAllMasks).*0.05)),'Voxels) of the whole Mask(',num2str(sum(OutputStruct.NVoxelAllMasks)),'Voxels).']; ...
        'No --> "do not exclude any subpart". (However, areas will be sorted according to size.)';
        }, 'Exclude subparts by size?','Yes','T<5%','No','T<5%');
    switch(choice)
        case 'Yes'
            %ask user for input
            prompt   ={['Cutoff Percentage: [1%==',num2str(round(sum(OutputStruct.NVoxelAllMasks).*0.01)),'Voxels]']};
            dlg_title=['Exclude Subparts of Mask "',MaskName,'"'];
            num_lines=1;
            defAns   ={'1'};
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer   = inputdlg(prompt,dlg_title,num_lines,defAns,options);
            
            OutputStruct.MasksToBeUsed.NThres = round(sum(OutputStruct.NVoxelAllMasks).*(eval(answer{1})/100));
            OutputStruct.MasksToBeUsed.PThres = [answer{1},'%(',num2str(OutputStruct.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(OutputStruct.NVoxelAllMasks)),'Total)'];
            IndsMasksDisp = find(OutputStruct.NVoxelAllMasks>OutputStruct.MasksToBeUsed.NThres);
            [tmp,SortInds] = sort(OutputStruct.NVoxelAllMasks(IndsMasksDisp),'descend');
            OutputStruct.MasksToBeUsed.IndsMasksDisp = IndsMasksDisp(SortInds);
            
            OutputStruct.MasksToBeUsed.InfoStr= ['Only subparts (sorted by size) of Mask "',MaskName,'" are included that are bigger than ',OutputStruct.MasksToBeUsed.PThres];
        case 'T<5%'
            OutputStruct.MasksToBeUsed.NThres = round(sum(OutputStruct.NVoxelAllMasks).*0.05);
            OutputStruct.MasksToBeUsed.PThres = ['5%(',num2str(OutputStruct.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(OutputStruct.NVoxelAllMasks)),'Total)'];
            IndsMasksDisp = find(OutputStruct.NVoxelAllMasks>OutputStruct.MasksToBeUsed.NThres);
            [tmp,SortInds] = sort(OutputStruct.NVoxelAllMasks(IndsMasksDisp),'descend');
            OutputStruct.MasksToBeUsed.IndsMasksDisp = IndsMasksDisp(SortInds);
            
            OutputStruct.MasksToBeUsed.InfoStr= ['Only subparts (sorted by size) of Mask "',MaskName,'" are included that are bigger than ',OutputStruct.MasksToBeUsed.PThres];
        case 'No'
            OutputStruct.MasksToBeUsed.NThres = 0;
            OutputStruct.MasksToBeUsed.PThres = ['0%(',num2str(OutputStruct.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(OutputStruct.NVoxelAllMasks)),'Total)'];
            [tmp,SortInds] = sort(OutputStruct.NVoxelAllMasks,'descend');
            OutputStruct.MasksToBeUsed.IndsMasksDisp = SortInds;
            
            OutputStruct.MasksToBeUsed.InfoStr= ['All subparts (sorted by size) of Mask "',MaskName,'" are included.'];
    end
else
    OutputStruct.MasksToBeUsed.NThres = 0;
    OutputStruct.MasksToBeUsed.PThres = ['0%(',num2str(OutputStruct.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(OutputStruct.NVoxelAllMasks)),'Total)'];
    [tmp,SortInds] = sort(OutputStruct.NVoxelAllMasks,'descend');
    OutputStruct.MasksToBeUsed.IndsMasksDisp = SortInds;
    
    OutputStruct.MasksToBeUsed.InfoStr= ['All subparts (sorted by size) of Mask "',MaskName,'" are included. (Sorted by size)'];
end

end

function [BPData,BPGrouping,BPLabels,InfoStr,IndicesSort] = AssembleBPinputs_n_Resort(ScalingData,PlotType,AllMasksNames,AllMasks,CheckResort,LambdaRef,IndFix,varargin)
% Assemble data for Boxplot.
% Resort the distributions in boxplot according to absolute deviation from reference median.
% IN CASE GLOBAL SORTING IS USED INDICESSORT WILL BE A CELLSTR WITH THE NAMES FOR MATCHING.

if(nargin>7)
    IndicesSort = varargin{1};
    if(isempty(IndicesSort))
        clear IndicesSort %none given so let's sort for this instance too.
    end
    if(nargin>8)
        MixIndices = varargin{2};
        AverageMix = varargin{3};
    end
end

%% do resorting or just leave data as is
SubregionMedians = zeros(length(AllMasksNames),1);
BPData     = [];
BPGrouping = [];
BPLabels   = cell(length(AllMasksNames),1);
for IndMask = 1:length(AllMasksNames)
    DataInputTmp = ScalingData(AllMasks(:,IndMask)>0,:); %pick the part of the data that is in the mask, i.e. mask~=0
    if(1) %use all data
    switch(PlotType) 
        case 'Aggregate'
            DataTmp = DataInputTmp;
        case 'MedianSubj'
            DataTmp = median(DataInputTmp,2);
        case 'MedianVoxel'
            DataTmp = squeeze(median(DataInputTmp,1));
        case 'QuartilesSubjMix'
            DataTmp = quantile(DataInputTmp,[.25 .50 .75],2);
            if(AverageMix)
                DataTmp = mean(DataTmp(:,MixIndices),2);
            else
                DataTmp =      DataTmp(:,MixIndices);
            end
        case 'QuartilesVoxelMix'
            DataTmp = permute(quantile(DataInputTmp,[.25 .50 .75],1),[2 1]);
            if(AverageMix)
                DataTmp = mean(DataTmp(:,MixIndices),2);
            else
                DataTmp =      DataTmp(:,MixIndices);
            end
    end
    else %remove negative as outliers
        switch(PlotType)
        case 'Aggregate'
            DataTmp = DataInputTmp;
            DataTmp(DataInputTmp<0) = [];
        case 'MedianSubj'
            DataTmp = median(DataInputTmp,2);
            DataTmp(DataTmp<0) = [];
        case 'MedianVoxel'
            DataTmp = squeeze(median(DataInputTmp,1));
        case 'QuartilesSubjMix'
            DataTmp = quantile(DataInputTmp,[.25 .50 .75],2);
            if(AverageMix)
                DataTmp = mean(DataTmp(:,MixIndices),2);
            else
                DataTmp =      DataTmp(:,MixIndices);
            end
            DataTmp(DataTmp<0) = [];
        case 'QuartilesVoxelMix'
            DataTmp = permute(quantile(DataInputTmp,[.25 .50 .75],1),[2 1]);
            if(AverageMix)
                DataTmp = mean(DataTmp(:,MixIndices),2);
            else
                DataTmp =      DataTmp(:,MixIndices);
            end
            DataTmp(DataTmp<0) = [];
        end
    end
    
    if(IndFix==IndMask)
        SubregionMedians(IndMask) = 0;
        %BugFix
        if(1)
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},'SignifVoxelsIC WO EffectOfMRI','DMNunmod');
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},'EffectOfMRI','DMNmod');
        else
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},'SignifVoxels','');
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},' WO ','\');
        end
    else
        SubregionMedians(IndMask) = abs(median(DataTmp(:))-LambdaRef); %for sorting even if not used.
        %BigFix2
        if(1)
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},'ParaHippocampal','ParaHipp.');
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},'Hippocampus','Hipp.');
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},'Cingulum ','Cing.');
            AllMasksNames{IndMask} = regexprep(AllMasksNames{IndMask},' ','');
        end
    end
    BPData     = [BPData;DataTmp(:)]; %just add the data that is in the mask==1 parts otherwise this gets too big.
    BPGrouping = [BPGrouping; IndMask.*ones(length(DataTmp(:)),1)]; %label this data with the index of the mask
    BPLabels{IndMask} = AllMasksNames{IndMask};
end

if(CheckResort)
    if(~exist('IndicesSort','var'))
        [SortedMedians,IndicesSort] = sort(SubregionMedians);
    else
        if(iscellstr(IndicesSort)) %this means we have to apply global sorting
            LabelsToMatch = IndicesSort;
            clear IndicesSort %remove to avoid pain later
            IndicesSortLeft(:,1) = 1:length(BPLabels); %start with original indices which will be moved around
            IndicesSort    = []; %fill this
            for IndLMatch = 1:length(LabelsToMatch)
                CurLabelToMatch = LabelsToMatch{IndLMatch};
                MatchStr{1} = regexprep(CurLabelToMatch,'_',' '); %first, use this to look for FULL STRING WITHOUT UNDERSCORE BUT WHITESPACE
                Indices     = strfind(MatchStr{1},' ');
                if(~isempty(Indices))
                    MatchStr{2} = CurLabelToMatch(1:(Indices(1)-1)); %second, use this to look for FIRST PART OF STRING
                end
                
                %go over all masks not assigned yet and check if they match
                if(~isempty(IndicesSortLeft)) %redundancy check
                    LastInitialEndIndex = IndicesSortLeft(end); %take this for checking if we didn't already go through the whole list.
                    NSearch = 0;
                    while(~isempty(IndicesSortLeft))
                        %                 for IndMask = 1:length(IndicesSortLeft)
                        [startIndex] = regexpi(BPLabels{IndicesSortLeft(1)}, MatchStr{1}); %full string match?
                        if(length(MatchStr)>1)
                            if(isempty(startIndex)) %full string not found so let's try the shorter version if this is the second round
                                if(NSearch~=0) %just to make sure we don't make some stupid mistake like assigning Temporal Inf the place of Temporal Mid because we weren't patient enough to first keep looking for it.
                                    [startIndex] = regexpi(BPLabels{IndicesSortLeft(1)}, MatchStr{2}); %short sting match?
                                end
                            end
                        end
                        %check if there was any match
                        if(~isempty(startIndex)) %match!
                            IndicesSort       = [IndicesSort; IndicesSortLeft(1)]; %add indice NB must be first one, otherwise it is moved to end.
                            IndicesSortLeft(1)= []; %don't look for this one any more.
                            break; %don't continue because we already found one
                        else %no match!
                            if(IndicesSortLeft(1)==LastInitialEndIndex) %we are through the list already once and have no match so let's break out of this and try next one.
                                %if we are already more than once throught --> break out of loop
                                if(NSearch==0)
                                    NSearch = NSearch + 1;
                                    %--> put at the end
                                    IndicesSortLeft(end+1) = IndicesSortLeft(1);
                                    IndicesSortLeft(1)     = [];
                                else
                                    break;
                                end
                            else %--> put at the end
                                IndicesSortLeft(end+1) = IndicesSortLeft(1);
                                IndicesSortLeft(1)     = [];
                            end
                        end
                    end
                else
                    break; %there are none left to match so let's stop
                end
            end
            if(~isempty(IndicesSortLeft)) %sort them according to median
                [SortedMedians,IndicesSortTmp] = sort(SubregionMedians(IndicesSortLeft));
                IndicesSortLeft = IndicesSortLeft(IndicesSortTmp);
            end
            IndicesSort = [IndicesSort; IndicesSortLeft];
            %% check that IndFix is at first place
            IndicesSort = [IndFix; IndicesSort(IndicesSort~=IndFix)];
        end
    end
    BPGroupingTmp = BPGrouping;
    BPLabelsTmp   = BPLabels;
    for IndMask = 1:length(AllMasksNames)
        BPGrouping(BPGroupingTmp==IndicesSort(IndMask)) = IndMask;
        BPLabels{IndMask} = BPLabelsTmp{IndicesSort(IndMask)};
    end
    InfoStr = ['Sorted[abs(median(Data(MASK))-',num2str(LambdaRef),'),ascending]'];
else
    IndicesSort = 1:length(AllMasksNames);
    InfoStr     = 'Sorted[size(Mask),ascending]';
end

end