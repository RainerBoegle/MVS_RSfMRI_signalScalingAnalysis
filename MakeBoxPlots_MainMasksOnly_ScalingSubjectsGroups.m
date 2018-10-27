function [DataForPlotting,Handles] = MakeBoxPlots_MainMasksOnly_ScalingSubjectsGroups(varargin)
% This function makes boxplots for each main mask containing the subjects of
% each group regarding the rsfMRI scaling-values. I.e. a check of the distributions
% that will show more the deviations over the group/subjects but not so easy to see the commonality.
%
% This includes:
%               Boxplot 1:
%                           subplot(1): The main mask "EffectOfMRI"
%                           subplot(2): The main mask IC-Map\"EffectOfMRI", i.e. the exclusion
%
%
%Usage:
%      [DataForPlotting,Handles] = MakeBoxPlots_MainMasksOnly_ScalingSubjectsGroups(varargin);
%      [DataForPlotting,Handles] = MakeBoxPlots_MainMasksOnly_ScalingSubjectsGroups(DataForPlotting); %automatically plot data from a previous plot configuration, i.e. for a given split, sorting of areas according to reference value of median and so on
%      [DataForPlotting,Handles] = MakeBoxPlots_MainMasksOnly_ScalingSubjectsGroups();                %ask user to select *.mat-file to load MVSscaling-results struct and then ask for parameters for plotting. All settings and data for plotting will be output in "DataForPlotting"-struct. 
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
%Date: 25.09.2014 (V1 adapted from V5 of DisplayMVSfMRIscalingResults.m)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)


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
else
    %we got an input, assume it is DataForPlotting and just plot
    DataForPlotting = varargin{1};
end


%% boxplot 1: Scaling Values MainMasks 
if(~isempty(DataForPlotting.Covars))
    [CovarsSelection,UseCovars] = listdlg('ListString',DataForPlotting.Covars.Name,'CancelString','No Covars','SelectionMode','single');
    if(UseCovars)
        CovarName = DataForPlotting.Covars.Name{CovarsSelection};
    end
else
    UseCovars = 0;
end

%% ylim?
switch(questdlg('Add YLIM?','YLIM?','Yes','No','Yes'))
    case 'Yes'
        UseYLim = 1;
        if(UseCovars)
            answer_YLim = inputdlg({'YLim = '},'Ylim?',1,{['[',num2str(floor(min(DataForPlotting.Covars.RegForPlot(:,CovarsSelection)))-1),' ',num2str(ceil(max(DataForPlotting.Covars.RegForPlot(:,CovarsSelection)))+1),']']});
        else
            answer_YLim = inputdlg({'YLim = '},'Ylim?',1,{'[0 6]'});
        end
        YLimits = sort(eval(answer_YLim{1}));
    otherwise
        UseYLim = 0;
        YLimits = [];
end


Handles{1} = figure(1); clf;
for IndSubplot = 1:length(DataForPlotting.Plot)
    try
        BPLabels   = DataForPlotting.Plot(IndSubplot).BPLabels;
        BPData     = DataForPlotting.Plot(IndSubplot).BPData;
        BPGrouping = DataForPlotting.Plot(IndSubplot).BPGrouping;
        Title      = DataForPlotting.Plot(IndSubplot).Title;
    catch CATCH_InterpretDataForPlotting
        assignin('base','CATCH_InterpretDataForPlotting',CATCH_InterpretDataForPlotting);
        disp('The structure "DataForPlotting" might not belong to this function/can not be interpreted.');
        disp('Try another function, maybe.');
        return;
    end
    
    
    subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(1:length(BPLabels),2.*ones(length(BPLabels),1),'r-'); hold('on');
    subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(1:length(BPLabels),sqrt(2).*ones(length(BPLabels),1),'r--'); hold('on');
    subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(1:length(BPLabels),1.*ones(length(BPLabels),1),'b-'); hold('on'); title(Title);
    if(UseCovars)
        subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(1:length(BPLabels),DataForPlotting.Covars.RegForPlot(:,CovarsSelection),'kx','MarkerSize',10); hold('on'); 
        subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(1:length(BPLabels),DataForPlotting.Covars.RegForPlot(:,CovarsSelection),'ko','MarkerSize',10); hold('on'); 
    end
    subplot(length(DataForPlotting.Plot),1,IndSubplot); boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); hold('on');
    if(UseYLim)
        ylim(YLimits);
    end
end

%% scatter plot Covar over rsfMRI data (if(UseCovars))
if(UseCovars)
    ChoiceScatterPlot = questdlg('Plot Covars over rsfMRI Quartiles or rsfMRI Quartiles over Covars?','Plot choice?','Covars-x-Quartiles','Quartiles-x-Covars','Quartiles-x-Covars');    
    
    if(strcmp(ChoiceScatterPlot,'Quartiles-x-Covars'))
        Handles{2} = figure(2); clf;
        x = DataForPlotting.Covars.RegForPlotSubj(:,CovarsSelection);
    else
        Handles{2} = figure(3); clf;
        y = DataForPlotting.Covars.RegForPlotSubj(:,CovarsSelection);
    end
    for IndSubplot = 1:length(DataForPlotting.Plot)
        LegendStr  = cell(size(DataForPlotting.Plot(IndSubplot).QrtPerSubj,2)+3,1);
        for IndQrt = 1:size(DataForPlotting.Plot(IndSubplot).QrtPerSubj,2)
            if(strcmp(ChoiceScatterPlot,'Quartiles-x-Covars'))
                y = DataForPlotting.Plot(IndSubplot).QrtPerSubj(:,IndQrt);
            else
                x = DataForPlotting.Plot(IndSubplot).QrtPerSubj(:,IndQrt);
                subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,      2.*ones(size(x)),'r-'); hold('on');
                subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,sqrt(2).*ones(size(x)),'r:'); hold('on');
                subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,      1.*ones(size(x)),'b-'); hold('on');
            end
            switch(IndQrt)
                case 1
                    Col = 'g';
                case 2 %median
                    Col = 'r';
                case 3
                    Col = 'b';
                otherwise
                    Col = 'y';
            end
            subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,y,[Col,'x'],'MarkerSize',10); hold('on');
            LegendStr{IndQrt} = [num2str(IndQrt),'.Qrt (',ChoiceScatterPlot,')'];
        end
        if(strcmp(ChoiceScatterPlot,'Quartiles-x-Covars'))
            subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,      2.*ones(size(x)),'r-'); hold('on');
            subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,sqrt(2).*ones(size(x)),'r:'); hold('on');
            subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(x,      1.*ones(size(x)),'b-'); hold('on');
        end
        LegendStr{size(DataForPlotting.Plot(IndSubplot).QrtPerSubj,2)+1} = '$$\lambda= 2$$';
        LegendStr{size(DataForPlotting.Plot(IndSubplot).QrtPerSubj,2)+2} = '$$\lambda= \sqrt{2}$$';
        LegendStr{size(DataForPlotting.Plot(IndSubplot).QrtPerSubj,2)+3} = '$$\lambda= 1$$';
            
        title([regexprep(ChoiceScatterPlot,'Covars',['Covar(',CovarName,')']),': ',Title]);
        h = legend(LegendStr);
        set(h,'Interpreter','latex') 
        if(UseYLim)
            ylim(YLimits);
        end
    end
end


%% Diff plot
if(UseCovars)
    IndQrt = 2; %median
    
    Handles{3} = figure(Handles{2}+1); clf;
    for IndSubplot = 1:length(DataForPlotting.Plot)
        Data = DataForPlotting.Plot(IndSubplot).QrtPerSubj(:,IndQrt);
        AbsDiff = abs(Data - DataForPlotting.Covars.RegForPlotSubj(:,CovarsSelection));
        subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(AbsDiff,'kx','MarkerSize',10); title(['AbsDiff to ',CovarName])
    end
    
    Handles{4} = figure(Handles{3}+1); clf;
    for IndSubplot = 1:length(DataForPlotting.Plot)
        Data = DataForPlotting.Plot(IndSubplot).QrtPerSubj(:,IndQrt);
        AbsDiff = abs(Data - DataForPlotting.Covars.RegForPlotSubj(:,CovarsSelection));
        subplot(length(DataForPlotting.Plot),1,IndSubplot); plot(Data,AbsDiff,'kx','MarkerSize',10); title(['AbsDiff to ',CovarName,' over MedianPerSubject'])
    end
end
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

%% define which subplots to make
[Selection,ok] = listdlg('ListString',DataForPlotting.AllMasksNames,'InitialValue',[3, 5],'CancelString','Defaults');
if(~ok)
    Selection = [3, 5];
end
NSubplots = length(Selection); %set number of subplots

%% collect boxplot data
DataForPlotting.SubjIDnum     = MVSscaling.InputData.ICA.EyeMovAnaSubjNumPerANOVASubjNum;
DataForPlotting.GroupNrPerSubj= MVSscaling.Design.GroupNrPerSubj;
DataForPlotting.NGroups       = MVSscaling.Design.NGroups;
DataForPlotting.NSubjs        = MVSscaling.Design.NSubjs;
DataForPlotting.NSubjPerGroup = zeros(DataForPlotting.NGroups,1);
for IndGroup = 1:DataForPlotting.NGroups
    DataForPlotting.NSubjPerGroup(IndGroup) = length(find(DataForPlotting.GroupNrPerSubj==IndGroup));
end

for IndSubplot = 1:NSubplots    
    IndMask = Selection(IndSubplot); %pick mask from selection
    Mask = DataForPlotting.AllMasks(:,IndMask); %assign mask
    Title = ['ScalingValues for Mask "',DataForPlotting.AllMasksNames{IndMask},'" for all individual subjects per group.']; %assign title of plot
    
    BPData = [];
    BPGrouping = [];
    BPLabels   = cell((1+DataForPlotting.NGroups+DataForPlotting.NSubjs),1);
    QrtPerSubj = zeros(DataForPlotting.NSubjs,3);
    for IndGroup = 1:DataForPlotting.NGroups %go over groups
        %% add median of all subjects
        Index = 1;
        Data = median(DataForPlotting.ScalingData(Mask~=0,:),2);
        BPData     = [BPData; Data(:)];
        BPGrouping = [BPGrouping; Index.*ones(length(Data(:)),1)];
        BPLabels{Index} = ['MEDIAN(AllSubj)'];
        
        %% add median subject for the group        
        SubjInds = find(DataForPlotting.GroupNrPerSubj==IndGroup); %pick subjects belonging to group
        
        Index = 1+IndGroup;
        Data = median(DataForPlotting.ScalingData(Mask~=0,SubjInds),2);
        BPData     = [BPData; Data(:)];
        BPGrouping = [BPGrouping; Index.*ones(length(Data(:)),1)];
        BPLabels{Index} = ['MEDIAN(G',num2str(IndGroup),')'];
        
        %% add subjects of the group
        for IndSubj = 1:length(SubjInds) %assign subjects to boxplot
            if(IndGroup==1)
                Index     = 1+DataForPlotting.NGroups+IndSubj;
                IndexSubj = IndSubj;
            else
                Index     = 1+DataForPlotting.NGroups+sum(DataForPlotting.NSubjPerGroup(1:(IndGroup-1)))+IndSubj;
                IndexSubj = sum(DataForPlotting.NSubjPerGroup(1:(IndGroup-1)))+IndSubj;
            end
            Data = DataForPlotting.ScalingData(Mask~=0,SubjInds(IndSubj));
            QrtPerSubj(IndexSubj,:) = quantile(Data,[.25 .5 .75]);
            BPData     = [BPData; Data(:)];
            BPGrouping = [BPGrouping; Index.*ones(length(Data(:)),1)];
            BPLabels{Index} = ['G',num2str(IndGroup),'-S',num2str(SubjInds(IndSubj))];
        end
    end
    DataForPlotting.Plot(IndSubplot).QrtPerSubj = QrtPerSubj;
    DataForPlotting.Plot(IndSubplot).BPLabels   = BPLabels;
    DataForPlotting.Plot(IndSubplot).BPData     = BPData;
    DataForPlotting.Plot(IndSubplot).BPGrouping = BPGrouping;
    DataForPlotting.Plot(IndSubplot).Title      = Title;
end

%% add covariate for display?
CovarsPath = spm_select(1,'mat','Select covariate for plot...');
if(~isempty(CovarsPath))
    load(CovarsPath);
    
    DataForPlotting.CovarsOrg = Covars;
    
    NSubj = length(DataForPlotting.SubjIDnum);
    for IndSubj = 1:NSubj
        IndexSubjMatch = find(Covars.SubjNum==DataForPlotting.SubjIDnum(IndSubj));
        DataForPlotting.Covars.Regressors(IndSubj,:) = DataForPlotting.CovarsOrg.Regressors(IndexSubjMatch,:);
    end
    DataForPlotting.Covars.Name = DataForPlotting.CovarsOrg.Name;
    
    DataForPlotting.Covars.RegForPlot     = zeros((1+DataForPlotting.NGroups+DataForPlotting.NSubjs),length(DataForPlotting.Covars.Name));
    DataForPlotting.Covars.RegForPlotSubj = zeros(                           DataForPlotting.NSubjs, length(DataForPlotting.Covars.Name));
    %% median of all subjects
    DataForPlotting.Covars.RegForPlot(1,:) = squeeze(median(DataForPlotting.Covars.Regressors,1));
    for IndGroup = 1:DataForPlotting.NGroups %go over groups
        %% add median subject for the group        
        SubjInds = find(DataForPlotting.GroupNrPerSubj==IndGroup); %pick subjects belonging to group
        
        DataForPlotting.Covars.RegForPlot(1+IndGroup,:) = squeeze(median(DataForPlotting.Covars.Regressors(SubjInds,:),1));
        %% add subjects of the group
        for IndSubj = 1:length(SubjInds) %assign subjects to boxplot
            if(IndGroup==1)
                Index     = 1+DataForPlotting.NGroups+IndSubj;
                IndexSubj = IndSubj; 
            else
                Index     = 1+DataForPlotting.NGroups+sum(DataForPlotting.NSubjPerGroup(1:(IndGroup-1)))+IndSubj;
                IndexSubj = sum(DataForPlotting.NSubjPerGroup(1:(IndGroup-1)))+IndSubj;
            end
            DataForPlotting.Covars.RegForPlot(Index,:)         = squeeze(DataForPlotting.Covars.Regressors(SubjInds(IndSubj),:));
            DataForPlotting.Covars.RegForPlotSubj(IndexSubj,:) = squeeze(DataForPlotting.Covars.Regressors(SubjInds(IndSubj),:));
        end
    end    
else
    DataForPlotting.Covars = [];
end

%% decide on the split and exclusion due to voxel number after merging (if applicable)
% %% pick a split of mask "EffectsOfMRI" to use for plotting
% H = helpdlg('Select split: EffectsOfMRI','Select split: EffectsOfMRI');
% uiwait(H);
% [IndSplitEffectsOfMRI,ok] = listdlg('ListString',MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames,'SelectionMode','single','Name','Select split: EffectsOfMRI');
% if(~ok)
%     IndSplitEffectsOfMRI = 1; %default to one if none is picked.
% end
% 
% %% get masks for split and ask for merging if necessary
% [EffectOfMRI_AALsplit] = CombineSubparts(MVSscaling.Masks.AALsplit.EffectOfMRI,IndSplitEffectsOfMRI,'EffectsOfMRI');
% 
% %% assign split masks
% DataForPlotting.EffectOfMRI_AALsplit.IndSplit      = IndSplitEffectsOfMRI;
% DataForPlotting.EffectOfMRI_AALsplit.SplitName     = MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{IndSplitEffectsOfMRI};
% DataForPlotting.EffectOfMRI_AALsplit.AllMasks      = EffectOfMRI_AALsplit.AllMasks;
% DataForPlotting.EffectOfMRI_AALsplit.AllMasksNames = EffectOfMRI_AALsplit.AllMasksNames;
% DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks= EffectOfMRI_AALsplit.NVoxelAllMasks;
% DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed = EffectOfMRI_AALsplit.MasksToBeUsed; %NThres; PThres; IndsMasksDisp; InfoStr;
% 
% 
% %% pick a split of mask ICmap WITHOUT "EffectsOfMRI" to use for plotting
% H = helpdlg('Select split: IC-Map WO EffectsOfMRI','Select split: IC-Map WO EffectsOfMRI');
% uiwait(H);
% [IndSplitIC,ok] = listdlg('ListString',MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames,'SelectionMode','single','Name','Select split: IC-Map');
% if(~ok)
%     IndSplitIC = 1; %default to one if none is picked.
% end
% 
% %% get masks for split and ask for merging if necessary
% [ICmap_AALsplit] = CombineSubparts(MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI,IndSplitIC,'ICMap-EffectsOfMRI');
% 
% %% assign split masks
% DataForPlotting.ICmap_AALsplit.IndSplit      = IndSplitEffectsOfMRI;
% DataForPlotting.ICmap_AALsplit.SplitName     = MVSscaling.Masks.AALsplit.IC_WITHOUT_EffectOfMRI.SplitNames{IndSplitIC};
% DataForPlotting.ICmap_AALsplit.AllMasks      = ICmap_AALsplit.AllMasks;
% DataForPlotting.ICmap_AALsplit.AllMasksNames = ICmap_AALsplit.AllMasksNames;
% DataForPlotting.ICmap_AALsplit.NVoxelAllMasks= ICmap_AALsplit.NVoxelAllMasks;
% DataForPlotting.ICmap_AALsplit.MasksToBeUsed = ICmap_AALsplit.MasksToBeUsed; %NThres; PThres; IndsMasksDisp; InfoStr;


end


