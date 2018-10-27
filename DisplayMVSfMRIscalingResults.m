function [DataForPlotting] = DisplayMVSfMRIscalingResults(varargin)
% Function for displaying results of scaling mask-based analysis.
%
%USAGE:
%       [DataForPlotting] = DisplayMVSfMRIscalingResults(varargin);
%        DataForPlotting  = DisplayMVSfMRIscalingResults();                % setup plots from file and output all settings and data needed for plots
%        DataForPlotting  = DisplayMVSfMRIscalingResults(DataForPlotting); % reuse data for making plots without having to input details of plots again.
%
%Date: 29.08.2014 (V5)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% to do
H = helpdlg('Split also IC and IC-EffectOfMRI into subparts and show that they also do not modulate like EffectOfMRI, or only because they belong to EffectOfMRI.','one more split');
uiwait(H);

%% are there figures open? If yes, do you want to close them?
CurrentFigHandles = findobj('Type','figure'); %all open figures that are not hidden. %figHandles = findall(0,'Type','figure'); %this even gets all hidden figure handles
if(~isempty(CurrentFigHandles))
    choice_CFigs = questdlg(['There are currently ',num2str(length(CurrentFigHandles)),' Figures opened. Do you want to close them?'],'Close Figures?','Yes','No','Quit','Yes');
    switch(choice_CFigs)
        case 'Quit'
            return;
        case 'Yes'
            close(CurrentFigHandles);
    end
end

if(nargin==0) %load data of scaling analysis & assemble data for plotting
    %% load results file
    % choice_load = questdlg('Load data OR Rerun display?','Load OR Rerun','Load','Rerun','Quit','Rerun');
    % switch(choice_load)
    %     case 'Quit'
    %         return;
    %     case 'Load'
    %         MVSscalingMAT_Path = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...');
    %     case 'Rerun'
    %         disp(['Will display data from "',MVSscalingMat_Path,'"']);
    % end
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
    
    %% Group results or Single Subjects?
    choice_plotting = questdlg('Plot Groups Average results, Groups Aggregate results OR each single Subject?','Groups OR Subjects?','GroupsAverage','GroupsAggregate','SingleSubjects','GroupsAverage');
    if(isempty(choice_plotting))
        return;
    end

    %% Get Data
    DataForPlotting = GetDataForPlotting(MVSscalingMAT_Path,choice_plotting);
    
    % DataForPlotting.PlotType      = 'SingleSubjects'; %or 'GroupsAverage';
    % DataForPlotting.ScalingData
    % DataForPlotting.AmplitudesData
    % DataForPlotting.DataTitle
    %
    %
    % DataForPlotting.AllMasks
    % DataForPlotting.AllMasksNames
    % DataForPlotting.ThresMaps
    % DataForPlotting.ThresMapsNames
else
    %we got an input, assume it is DataForPlotting and just plot
    DataForPlotting = varargin{1};
end

%% Plot Data
%% hist
XHist = [-0.01,0:0.125:4]; %[-1:0.125:4]; %[-10:0.1:10];
for IndMask = 1:length(DataForPlotting.AllMasksNames)
    if(strcmpi(DataForPlotting.PlotType,'GroupsAverage'))
        for IndData = 1:length(DataForPlotting.DataTitle)
            HistFigNum = 10+IndData-1;  %ONE GROUP PER FIGURE with all masks displayed
            switch(IndMask)
                case {1,2,3,4,5}
                    IndSubplot = IndMask;
                case 6
                    IndSubplot = 7;
            end
            figure(HistFigNum); subplot(3,3,IndSubplot); ...
                hist(DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,IndData),XHist); title(['histogram scaling[Average] (Mask: ',DataForPlotting.AllMasksNames{IndMask},')'])
            ylim([0, 500]);
            xlim([-0.01,3]);
            set(figure(HistFigNum),'Name',[DataForPlotting.DataTitle{IndData}]);
        end
    else
        if(strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
            for IndData = 1:length(DataForPlotting.DataTitle)
                HistFigNum = 10+IndData-1;  %ONE GROUP PER FIGURE with all masks displayed
                switch(IndMask)
                    case {1,2,3,4,5}
                        IndSubplot = IndMask;
                    case 6
                        IndSubplot = 7;
                end
                SubjInds = DataForPlotting.DataGrouping{IndData};
                DataTmp  = DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,SubjInds);
                figure(HistFigNum); subplot(3,3,IndSubplot); ...
                    hist(DataTmp(:),XHist); title(['histogram scaling[Aggregate] (Mask: ',DataForPlotting.AllMasksNames{IndMask},')'])
                ylim([0, 500]);
                xlim([-0.01,3]);
                set(figure(HistFigNum),'Name',[DataForPlotting.DataTitle{IndData}]);
            end
        else
            if(strcmpi(DataForPlotting.PlotType,'SingleSubjects'))
                for IndData = 1:length(DataForPlotting.DataTitle)
                    if(IndData<=DataForPlotting.NSubjPerGroup(1))
                        HistFigNum = 100+IndMask; %ONE MASK PER FIGURE with all subjects per group displayed (sorted by group)
                        CurrSubPlotLength = DataForPlotting.NSubjPerGroup(1);
                        IndSubplot = IndData;
                    else
                        HistFigNum = 110+IndMask; %ONE MASK PER FIGURE with all subjects per group displayed (sorted by group)
                        CurrSubPlotLength = DataForPlotting.NSubjPerGroup(2);
                        IndSubplot = IndData-DataForPlotting.NSubjPerGroup(1);
                    end
                    
                    figure(HistFigNum); subplot(CurrSubPlotLength,1,IndSubplot); ...
                        hist(DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,IndData),XHist); title(['histogram scaling ',DataForPlotting.DataTitle{IndData}])
                    
                    if(IndSubplot~=CurrSubPlotLength)
                        set(gca,'xtick',[])
                    end
                    ylim([0, 100]);
                    xlim([-0.01,3]);
                    set(figure(HistFigNum),'Name',[DataForPlotting.AllMasksNames{IndMask},'-MASK']);
                end
            end
        end
    end
end
hold off

%% boxplot 1: subplot per mask comparing groups
BPFigNum1 = 20;
for IndMask = 1:length(DataForPlotting.AllMasksNames)
    BPData     = [];
    BPGrouping = [];
    BPLabels   = cell(length(DataForPlotting.DataTitle),1);
    
    if(strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
        for IndData = 1:length(DataForPlotting.DataTitle)
            SubjInds = DataForPlotting.DataGrouping{IndData};
            DataTmp  = DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,SubjInds);
            BPData     = [BPData;DataTmp(:)];
            BPGrouping = [BPGrouping; IndData.*ones(length(DataTmp(:)),1)];
            BPLabels{IndData} = DataForPlotting.DataTitle{IndData};
        end
        InfoStr = '[Aggregate]';
    else
        if(strcmpi(DataForPlotting.PlotType,'GroupsAverage'))
            InfoStr = '[Average]';
        else
            InfoStr = []; %single subjects
        end
        for IndData = 1:length(DataForPlotting.DataTitle)
            BPData     = [BPData;DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,IndData)];
            BPGrouping = [BPGrouping; IndData.*ones(length(DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,IndData)),1)];
            BPLabels{IndData} = DataForPlotting.DataTitle{IndData};
        end
    end
    TitleString = ['boxplot scaling',InfoStr,' (Mask: ',DataForPlotting.AllMasksNames{IndMask},')'];
    figure(BPFigNum1); subplot(length(DataForPlotting.AllMasksNames),1,IndMask); ...
        plot(1:length(BPLabels),2.*ones(length(BPLabels),1),'r-'); hold on
        plot(1:length(BPLabels),sqrt(2).*ones(length(BPLabels),1),'r--'); hold on
        plot(1:length(BPLabels),1.*ones(length(BPLabels),1),'b-'); hold on
        if(strcmpi(DataForPlotting.PlotType,'SingleSubjects'))
            if(IndMask~=length(DataForPlotting.AllMasksNames))
                boxplot(BPData,BPGrouping,'notch','on'); title(TitleString); hold on
                set(gca,'xtick',[])
            else
                boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleString); hold off
                set(gca,'XTick',1:max(BPGrouping(:)))
                set(gca,'XTickLabel',BPLabels)
            end
        else
            if(IndMask~=length(DataForPlotting.AllMasksNames))
                boxplot(BPData,BPGrouping,'notch','on'); title(TitleString); hold off
                set(gca,'xtick',[])
            else
                boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleString); hold off
                set(gca,'XTick',1:max(BPGrouping(:)))
                set(gca,'XTickLabel',BPLabels)
            end
        end
    ylim([0, 3]);
    if(IndMask==1)
        legend('\lambda = 3T/1.5T = 2','\lambda = sqrt(3T/1.5T) = sqrt(2)','\lambda = 1','Data','Location','Best');
    end
    set(figure(BPFigNum1),'Name',['boxplot (groups) scaling ',DataForPlotting.PlotType]);
end
hold off

%% boxplot 2: subplot per group comparing masks

choice_resort = questdlg({'Do you want to resort the displayed subregion data, i.e. aside from the complete mask "EffectsOfMRI",'; 'according to the absolute difference of the median from "2==\lambda"?'},'Resort data by abs(median(data)-2)?','Yes','No','Yes');
switch(choice_resort)
    case 'Yes'
        CheckResort = 1;
    otherwise
        CheckResort = 0;
end

BPFigNum2 = 21;
for IndData = 1:length(DataForPlotting.DataTitle)
    BPData     = [];
    BPGrouping = [];
    BPLabels   = cell(length(DataForPlotting.AllMasksNames),1);
    
    if(strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
        SubjInds = DataForPlotting.DataGrouping{IndData};
        SubregionMedians = zeros(length(DataForPlotting.AllMasksNames),1);
        for IndMask = 1:length(DataForPlotting.AllMasksNames)
            DataTmp = DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,SubjInds);
            SubregionMedians(IndMask) = abs(median(DataTmp(:))-2); %for sorting
            BPData     = [BPData;DataTmp(:)];
            BPGrouping = [BPGrouping; IndMask.*ones(length(DataTmp(:)),1)];
            BPLabels{IndMask} = DataForPlotting.AllMasksNames{IndMask};
        end
        InfoStr = '[Aggregate]';
        if(CheckResort)
            [SortedMedians,IndicesSort] = sort(SubregionMedians);
            BPGroupingTmp = BPGrouping;
            BPLabelsTmp   = BPLabels;
            for IndMask = 1:length(DataForPlotting.AllMasksNames)
                BPGrouping(BPGroupingTmp==IndicesSort(IndMask)) = IndMask;
                BPLabels{IndMask} = BPLabelsTmp{IndicesSort(IndMask)};
            end 
        end
    else
        if(strcmpi(DataForPlotting.PlotType,'GroupsAverage'))
            InfoStr = '[Average]';
        else
            InfoStr = []; %single subjects
        end
        SubregionMedians = zeros(length(DataForPlotting.AllMasksNames),1);
        for IndMask = 1:length(DataForPlotting.AllMasksNames)
            DataTmp = DataForPlotting.ScalingData(DataForPlotting.AllMasks(:,IndMask)>0,IndData);
            SubregionMedians(IndMask) = abs(median(DataTmp(:))-2); %for sorting
            BPData     = [BPData; DataTmp(:)];
            BPGrouping = [BPGrouping; IndMask.*ones(length(DataTmp(:)),1)];
            BPLabels{IndMask} = DataForPlotting.AllMasksNames{IndMask};
        end
        if(CheckResort)
            [SortedMedians,IndicesSort] = sort(SubregionMedians);
            BPGroupingTmp = BPGrouping;
            BPLabelsTmp   = BPLabels;
            for IndMask = 1:length(DataForPlotting.AllMasksNames)
                BPGrouping(BPGroupingTmp==IndicesSort(IndMask)) = IndMask;
                BPLabels{IndMask} = BPLabelsTmp{IndicesSort(IndMask)};
            end 
        end
    end
    if(CheckResort)
        TitleString = ['boxplot scaling',InfoStr,'(sorted by ABS(median-2)) (',DataForPlotting.DataTitle{IndData},')'];
    else
        TitleString = ['boxplot scaling',InfoStr,' (',DataForPlotting.DataTitle{IndData},')'];
    end
    figure(BPFigNum2); subplot(length(DataForPlotting.DataTitle),1,IndData); ...
        plot(1:length(BPLabels),2.*ones(length(BPLabels),1),'r-'); hold on
        plot(1:length(BPLabels),sqrt(2).*ones(length(BPLabels),1),'r--'); hold on
        plot(1:length(BPLabels),1.*ones(length(BPLabels),1),'b-'); hold on
        if(strcmpi(DataForPlotting.PlotType,'SingleSubjects'))
            if(IndData~=length(DataForPlotting.DataTitle))
                boxplot(BPData,BPGrouping,'notch','on'); title(TitleString); hold on
                set(gca,'xtick',[])
            else
                boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleString); hold off
            end
        else
            boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleString); hold off
        end
    ylim([0, 4]);
    if(IndData==3)
        legend('\lambda = 3T/1.5T = 2','\lambda = sqrt(3T/1.5T) = sqrt(2)','\lambda = 1','Data','Location','Best');
    end
    set(figure(BPFigNum2),'Name',['boxplot (masks) scaling ',DataForPlotting.PlotType]);
end
hold off

%% boxplot (subregions of MaskEffectsOfMRI) - 1

%prep masks
MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp)+1);
MasksNamesForDisp   = cell(length(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp)+1,1);

MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,3);
MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{3}; %EffectsOfMRI WHOLE-MASK
for Ind=1:length(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp)
    MasksForDisp(:,1+Ind)   = DataForPlotting.EffectOfMRI_AALsplit.AllMasks(:,DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp(Ind));
    MasksNamesForDisp{1+Ind}= DataForPlotting.EffectOfMRI_AALsplit.AllMasksNames{DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp(Ind)};
end

BPsubFigNum = 22;      
if(strcmpi(DataForPlotting.PlotType,'GroupsAverage'))
    for IndGroup = 1:length(DataForPlotting.DataTitle)
        TitleStr = [DataForPlotting.DataTitle{IndGroup},'[Average]: ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr];
        
        BPData    = [];
        BPGrouping= [];
        BPLabels  = MasksNamesForDisp;
        SubregionMedians = zeros(length(MasksNamesForDisp),1);
        for IndMask = 1:length(MasksNamesForDisp)
            DataTmp    = DataForPlotting.ScalingData(MasksForDisp(:,IndMask)>0,IndGroup);
            BPData     = [BPData;DataForPlotting.ScalingData(MasksForDisp(:,IndMask)>0,IndGroup)];
            BPGrouping = [BPGrouping; IndMask.*ones(length(DataForPlotting.ScalingData(MasksForDisp(:,IndMask)>0,IndGroup)),1)];
            if(IndMask==1)
                SubregionMedians(IndMask) = 0; %trick to ensure that WholeMask "EffectsOfMRI" is still first. ;)
            else
                SubregionMedians(IndMask) = abs(median(DataTmp(:))-2);
            end
        end
        if(CheckResort)
            [SortedMedians,IndicesSort] = sort(SubregionMedians);
            BPGroupingTmp = BPGrouping;
            BPLabelsTmp   = BPLabels;
            for IndMask = 1:length(MasksNamesForDisp)
                BPGrouping(BPGroupingTmp==IndicesSort(IndMask)) = IndMask;
                BPLabels{IndMask} = BPLabelsTmp{IndicesSort(IndMask)};
            end 
            TitleStr =  regexprep(TitleStr, '(sorted by size)', '(sorted by ABS(median-2))');
        end
        figure(BPsubFigNum); subplot(length(DataForPlotting.DataTitle),1,IndGroup); ...
        plot(1:length(BPLabels),2.*ones(length(BPLabels),1),'r-'); hold on
        plot(1:length(BPLabels),sqrt(2).*ones(length(BPLabels),1),'r--'); hold on
        plot(1:length(BPLabels),1.*ones(length(BPLabels),1),'b-'); hold on
        boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleStr); hold off
        ylim([0, 6]);
    end
    set(figure(BPsubFigNum),'Name',['boxplot[Average] (subregions-EffectOfMRI) scaling ',DataForPlotting.PlotType]);
else
    if(strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
        for IndGroup = 1:length(DataForPlotting.DataTitle)
            TitleStr = [DataForPlotting.DataTitle{IndGroup},'[Aggregate]: ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr];
            SubjInds = DataForPlotting.DataGrouping{IndGroup};
            
            BPData    = [];
            BPGrouping= [];
            BPLabels  = MasksNamesForDisp;
            SubregionMedians = zeros(length(MasksNamesForDisp),1);
            for IndMask = 1:length(MasksNamesForDisp)
                DataTmp = DataForPlotting.ScalingData(MasksForDisp(:,IndMask)>0,SubjInds);
                BPData     = [BPData;DataTmp(:)];
                BPGrouping = [BPGrouping; IndMask.*ones(length(DataTmp(:)),1)];
                if(IndMask==1)
                    SubregionMedians(IndMask) = 0; %trick to ensure that WholeMask "EffectsOfMRI" is still first. ;)
                else
                    SubregionMedians(IndMask) = abs(median(DataTmp(:))-2);
                end
            end
            if(CheckResort)
                [SortedMedians,IndicesSort] = sort(SubregionMedians);
                BPGroupingTmp = BPGrouping;
                BPLabelsTmp   = BPLabels;
                for IndMask = 1:length(MasksNamesForDisp)
                    BPGrouping(BPGroupingTmp==IndicesSort(IndMask)) = IndMask;
                    BPLabels{IndMask} = BPLabelsTmp{IndicesSort(IndMask)};
                end
                TitleStr =  regexprep(TitleStr, '(sorted by size)', '(sorted by ABS(median-2))');
            end
            figure(BPsubFigNum); subplot(length(DataForPlotting.DataTitle),1,IndGroup); ...
            plot(1:length(BPLabels),2.*ones(length(BPLabels),1),'r-'); hold on
            plot(1:length(BPLabels),sqrt(2).*ones(length(BPLabels),1),'r--'); hold on
            plot(1:length(BPLabels),1.*ones(length(BPLabels),1),'b-'); hold on
            boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleStr); hold off
            ylim([0, 6]);
        end
        set(figure(BPsubFigNum),'Name',['boxplot[Aggregate] (subregions-EffectOfMRI) scaling ',DataForPlotting.PlotType]);
    else
        disp('Single subjects plots are not implemented yet, but if you do that');
        disp('please plot them together depending on the group that they belong to,');
        disp('i.e. one group per figure and one subplots for each subject of that group.');
    end
end
hold off

%% boxplot (subregions of MaskEffectsOfMRI) - 2: ONLY FOR AGGREGATE DATA 
%% For this plot the median per region is calculated and the distribution
%% of subjects is shown, as an alternative to showing the distribution of
%% voxels for the median subject (GroupAverage) or the collection of all
%% voxels and subjects.

UseQuart = 1;
IndQuart = 2; %1:3;%all %2;%median %1; %3; %
UseAvQuart = 0;%if(1) --> average the selected quarts

%prep masks
MasksForDisp        = zeros(size(DataForPlotting.AllMasks,1),length(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp)+1);
MasksNamesForDisp   = cell(length(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp)+1,1);

MasksForDisp(:,1)   = DataForPlotting.AllMasks(:,3);
MasksNamesForDisp{1}= DataForPlotting.AllMasksNames{3}; %EffectsOfMRI WHOLE-MASK
for Ind=1:length(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp)
    MasksForDisp(:,1+Ind)   = DataForPlotting.EffectOfMRI_AALsplit.AllMasks(:,DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp(Ind));
    MasksNamesForDisp{1+Ind}= DataForPlotting.EffectOfMRI_AALsplit.AllMasksNames{DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp(Ind)};
end


BPsubFigNum2 = 23;    
if(strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
    for IndGroup = 1:length(DataForPlotting.DataTitle)
        SubjInds = DataForPlotting.DataGrouping{IndGroup};
        
        BPData    = [];
        BPGrouping= [];
        BPLabels  = MasksNamesForDisp;
        SubregionMedians = zeros(length(MasksNamesForDisp),1);
        for IndMask = 1:length(MasksNamesForDisp)
            if(UseQuart)
                DataTmp = quantile(DataForPlotting.ScalingData(MasksForDisp(:,IndMask)>0,SubjInds),[.25 .50 .75],1);
                if(UseAvQuart)
                    DataTmp = mean(DataTmp(IndQuart,:),1); %average over selected quartiles
                    TitleStr = [DataForPlotting.DataTitle{IndGroup},'[Aggregate mean(quartiles(Area))]: ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr];
                else
                    DataTmp = DataTmp(IndQuart,:); %use all selected quartiles
                    TitleStr = [DataForPlotting.DataTitle{IndGroup},'[Aggregate quartiles(Area)]: ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr];
                end
                FigName  = ['boxplot[Aggregate] quartiles(subregions-EffectOfMRI scaling) per subject ',DataForPlotting.PlotType];
            else
                DataTmp = median(DataForPlotting.ScalingData(MasksForDisp(:,IndMask)>0,SubjInds),1);
                TitleStr = [DataForPlotting.DataTitle{IndGroup},'[Aggregate median(Area)]: ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr];
                FigName  = ['boxplot[Aggregate] median(subregions-EffectOfMRI scaling) per subject ',DataForPlotting.PlotType];
            end
            BPData     = [BPData;DataTmp(:)];
            BPGrouping = [BPGrouping; IndMask.*ones(length(DataTmp(:)),1)];
            if(IndMask==1)
                SubregionMedians(IndMask) = 0; %trick to ensure that WholeMask "EffectsOfMRI" is still first. ;)
            else
                SubregionMedians(IndMask) = abs(median(DataTmp(:))-2);
            end
        end
        if(CheckResort)
            [SortedMedians,IndicesSort] = sort(SubregionMedians);
            BPGroupingTmp = BPGrouping;
            BPLabelsTmp   = BPLabels;
            for IndMask = 1:length(MasksNamesForDisp)
                BPGrouping(BPGroupingTmp==IndicesSort(IndMask)) = IndMask;
                BPLabels{IndMask} = BPLabelsTmp{IndicesSort(IndMask)};
            end
            TitleStr =  regexprep(TitleStr, '(sorted by size)', '(sorted by ABS(median-2))');
        end
        
        figure(BPsubFigNum2); subplot(length(DataForPlotting.DataTitle),1,IndGroup); ...
        plot(1:length(BPLabels),2.*ones(length(BPLabels),1),'r-'); hold on
        plot(1:length(BPLabels),sqrt(2).*ones(length(BPLabels),1),'r--'); hold on
        plot(1:length(BPLabels),1.*ones(length(BPLabels),1),'b-'); hold on
        boxplot(BPData,BPGrouping,'labels',BPLabels,'notch','on'); title(TitleStr); hold off
        ylim([0, 6]);
    end
    set(figure(BPsubFigNum2),'Name',FigName);
end

%% stop?
choice_stop = questdlg('Stop now?','Stop?','Stop','Continue','Stop');
switch(choice_stop)
    case 'Stop'
        return;
end

%% scaling over significance values
for IndData = 1:length(DataForPlotting.DataTitle)
    if(~strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
        if(strcmpi(DataForPlotting.PlotType,'GroupsAverage'))
            SCFigNum = 30+IndData-1;
            InfoStr = '[Average]';
        else
            if(strcmpi(DataForPlotting.PlotType,'SingleSubjects'))
                SCFigNum = 300+IndData;
            end
            InfoStr = [];
        end
        for IndZValsMap = 1:length(DataForPlotting.ThresMapsNames)
            figure(SCFigNum); subplot(1,length(DataForPlotting.ThresMapsNames),IndZValsMap); ...
                plot(DataForPlotting.ThresMaps(:,IndZValsMap),2.*ones(size(DataForPlotting.ScalingData(:,IndData))),'r-'); hold on
            plot(DataForPlotting.ThresMaps(:,IndZValsMap),sqrt(2).*ones(size(DataForPlotting.ScalingData(:,IndData))),'r--'); hold on
            plot(DataForPlotting.ThresMaps(:,IndZValsMap),1.*ones(size(DataForPlotting.ScalingData(:,IndData))),'b-'); hold on
            plot(DataForPlotting.ThresMaps(:,IndZValsMap),DataForPlotting.ScalingData(:,IndData),'kx'); title(['Scatter plot Scaling',InfoStr,'-over-',DataForPlotting.ThresMapsNames{IndZValsMap},'-Map']); hold off
            ylim([0, 4]);
            legend('\lambda = 3T/1.5T = 2','\lambda = sqrt(3T/1.5T) = sqrt(2)','\lambda = 1','Data');
        end
    else
        InfoStr = '[Aggregate]';
        SCFigNum = 30+IndData-1;
        SubjInds = find(DataForPlotting.DataGrouping{IndGroup});
        for IndZValsMap = 1:length(DataForPlotting.ThresMapsNames)
            figure(SCFigNum); subplot(1,length(DataForPlotting.ThresMapsNames),IndZValsMap); ...
            plot(DataForPlotting.ThresMaps(:,IndZValsMap),2.*ones(size(DataForPlotting.ScalingData(:,SubjInds(1)))),'r-'); hold on
            plot(DataForPlotting.ThresMaps(:,IndZValsMap),sqrt(2).*ones(size(DataForPlotting.ScalingData(:,SubjInds(1)))),'r--'); hold on
            plot(DataForPlotting.ThresMaps(:,IndZValsMap),1.*ones(size(DataForPlotting.ScalingData(:,SubjInds(1)))),'b-'); hold on
            for IndSubj = 1:length(SubjInds)
                plot(DataForPlotting.ThresMaps(:,IndZValsMap),DataForPlotting.ScalingData(:,SubjInds(IndSubj)),'kx'); title(['Scatter plot Scaling',InfoStr,'-over-',DataForPlotting.ThresMapsNames{IndZValsMap},'-Map']); hold on
            end
            ylim([0, 4]);
            legend('\lambda = 3T/1.5T = 2','\lambda = sqrt(3T/1.5T) = sqrt(2)','\lambda = 1','Data');
        end
    end
    set(figure(SCFigNum),'Name',['SC',InfoStr,'-x-ZVals ',DataForPlotting.DataTitle{IndData}]);
end
hold off

%% scatterplot MRI2 over MRI1
% DataForPlotting.AmplitudesData(Voxels,MRIs,Data)
for IndData = 1:length(DataForPlotting.DataTitle)
    if(~strcmpi(DataForPlotting.PlotType,'GroupsAggregate'))
        if(strcmpi(DataForPlotting.PlotType,'GroupsAverage'))
            InfoStr = '[Average]';
            SCMRIsFigNum = 40+IndData-1;
        else
            if(strcmpi(DataForPlotting.PlotType,'SingleSubjects'))
                InfoStr = [];
                SCMRIsFigNum = 400+IndData;
            end
        end
        for IndMask = 1:length(DataForPlotting.AllMasksNames)
            switch(IndMask)
                case {1,2,3,4,5}
                    IndSubplot = IndMask;
                case 6
                    IndSubplot = 7;
            end
            figure(SCMRIsFigNum); subplot(3,3,IndSubplot); ...
            plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),   DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,2,IndData),'kx'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
            plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),2.*DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),'r-'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
            plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),sqrt(2).*DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),'r--'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
            plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),1.*DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,IndData),'b-'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold off
            if(IndSubplot==3)
                legend('Data','\lambda = 3T/1.5T = 2','\lambda = sqrt(3T/1.5T) = sqrt(2)','\lambda = 1','Location','Best');
            end
        end
    else
        InfoStr = '[Aggregate]';
        SCMRIsFigNum = 40+IndData-1;
        SubjInds = find(DataForPlotting.DataGrouping{IndGroup});
        for IndSubj = 1:length(SubjInds)
            for IndMask = 1:length(DataForPlotting.AllMasksNames)
                switch(IndMask)
                    case {1,2,3,4,5}
                        IndSubplot = IndMask;
                    case 6
                        IndSubplot = 7;
                end
                figure(SCMRIsFigNum); subplot(3,3,IndSubplot); ...
                plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),   DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,2,SubjInds(IndSubj)),'kx'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
                plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),2.*DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),'r-'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
                plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),sqrt(2).*DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),'r--'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
                plot(DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),1.*DataForPlotting.AmplitudesData(DataForPlotting.AllMasks(:,IndMask)>0,1,SubjInds(IndSubj)),'b-'); title([DataForPlotting.AllMasksNames{IndMask},'-Mask ',InfoStr]); hold on
                if(IndSubplot==3)
                    legend('Data','\lambda = 3T/1.5T = 2','\lambda = sqrt(3T/1.5T) = sqrt(2)','\lambda = 1','Location','Best');
                end
            end
        end
    end
    set(figure(SCMRIsFigNum),'Name',['MRI2-x-MRI1 Amplitude',InfoStr,' (',DataForPlotting.DataTitle{IndData},')']);
end
hold off

%% Done.
disp(' ');
disp('DONE.  [Hint: You can replot this using "[DataForPlotting] = DisplayMVSfMRIscalingResults(DataForPlotting)"]');
disp(' ');
end


%% SUBFUNCTIONS

%% get data
function DataForPlotting = GetDataForPlotting(MVSscalingMAT_Path,choice_plotting)
%% load data such that plotting will work
load(MVSscalingMAT_Path);

%% FUTURE EXTENSION: add choice of preproc e.g. median or mean & abs or not
choice_abs = questdlg('Take absolute value of data?','ABS(Data)?','Yes','No','Yes');
switch(choice_abs)
    case 'Yes'
        UseAbs = 1;
    otherwise
        UseAbs = 0;
end
DataForPlotting.UseAbs = UseAbs;

if(strcmp(choice_plotting,'GroupsAverage')) %how to average?
    choice_average = questdlg({'How to average the data? Mean or Median?'; 'NB: Median is more robust against outliers!'},'Mean OR Median?','Mean','Median','Median');
else
    choice_average = [];
end
DataForPlotting.choice_average = choice_average;

%% pick a split of mask "EffectsOfMRI" to use for plotting
[IndSplitToBeUsed,ok] = listdlg('ListString',MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames,'SelectionMode','single','Name','Select split');
if(~ok)
    IndSplitToBeUsed = 1; %default to one if none is picked.
end

%% FUTURE EXTENSION: ask user if subparts should be combined.
% DataForPlotting.EffectOfMRI_AALsplit.AllMasks      = MVSscaling.Masks.AALsplit.EffectOfMRI.AllMasks;
% DataForPlotting.EffectOfMRI_AALsplit.AllMasksNames = MVSscaling.Masks.AALsplit.EffectOfMRI.AllMasksName;
% DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks= MVSscaling.Masks.AALsplit.EffectOfMRI.NVoxelAllMasks;
NMaskMax = 5; %let's just say 5 is a lot because then the figure will have six boxes, so let's at least ask if user wants to combine masks
if(length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)>=NMaskMax) 
    InfoDlgTxt = cell(length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)+6,1);
    InfoDlgTxt{1} = ['"EffectsOfMRI" has been split(',MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{IndSplitToBeUsed},') into ',num2str(length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)),' subparts,'];
    InfoDlgTxt{2} = ['This might be a bit too much to look at at once. '];
    InfoDlgTxt{3} = ['If you want to combine some of them, click "combine", otherwise "skip".'];
    InfoDlgTxt{4} = ['[NB: You can still exclude small subparts from display in a later step instead of combining them with others.]'];
    InfoDlgTxt{5} =  ' ';
    InfoDlgTxt{6} = ['These are all the included subparts and their sizes: '];
    for Ind = 1:length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)
        InfoDlgTxt{6+Ind} = [num2str(Ind),'.',MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{Ind},': ',num2str(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)./sum(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks).*100),'%(',num2str(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)),' Voxels) of the whole Mask(',num2str(sum(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks)),' Voxels)'];
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
            SuggestionsInit = ManualMergers(MVSscaling.Masks.AALsplit.EffectOfMRI,IndSplitToBeUsed);
            
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
            Suggestions = SuggestMergers(MVSscaling.Masks.AALsplit.EffectOfMRI,IndSplitToBeUsed,SearchRange,NAreasMax,NVoxelsMax,SuggestionsInit);
                        
            %% show user which subparts are available
            InfoDlgTxt = cell(length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)+3,1);
            InfoDlgTxt{1} = ['The following subparts are available: '];
            for Ind = 1:length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)
                InfoDlgTxt{1+Ind} = [num2str(Ind),'.',MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{Ind},': ',num2str(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)./sum(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks).*100),'%(',num2str(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks(Ind)),' Voxels) of the whole Mask(',num2str(sum(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks)),' Voxels)'];
            end
            InfoDlgTxt{1+length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)+1} = ' ';
            InfoDlgTxt{1+length(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName)+2} = '       HINT:   KEEP THIS WINDOW OPEN   :HINT';
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
            
            AllMasks      = zeros(size(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasks,1),length(SubpartsCombine));
            AllMasksNames = cell(length(SubpartsCombine),1);
            NVoxelAllMasks= zeros(length(SubpartsCombine),1);
            for Ind = 1:length(SubpartsCombine)
                tmpMask     = zeros(size(AllMasks,1),1);
                tmpMaskName = [];
                for IndComb = 1:length(SubpartsCombine{Ind})
                    IndMask = find(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).LabelInds==SubpartsCombine{Ind}(IndComb));
                    tmpMask = tmpMask+MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasks(:,IndMask);
                    if(IndComb==1)
                        tmpMaskName = [tmpMaskName,MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                    else
                        tmpMaskName = [tmpMaskName,'+',MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
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
                    IndMask = find(MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).LabelInds==SubpartsCombine{Ind}(IndComb));
                    if(IndComb==1)
                        CombTxt = [CombTxt,num2str(Ind),'.(',num2str(SubpartsCombine{Ind}),'): "',MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                    else
                        if(IndComb~=length(SubpartsCombine{Ind}))
                            CombTxt = [CombTxt,'+',MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{IndMask}];
                        else
                            CombTxt = [CombTxt,'+',MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName{IndMask},'" (',num2str(NVoxelAllMasks(Ind)),'Voxels of ',num2str(sum(NVoxelAllMasks)),'Voxels in total.)'];
                        end
                    end
                end
                InfoCombDlgTxt{1+Ind} = CombTxt;
            end
            h_comb=helpdlg(InfoCombDlgTxt,'Combined subparts'); 
        otherwise
            AllMasks      = MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasks;
            AllMasksNames = MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).AllMasksName;
            NVoxelAllMasks= MVSscaling.Masks.AALsplit.EffectOfMRI.Split(IndSplitToBeUsed).NVoxelAllMasks;
    end
end

%% prep masks
DataForPlotting.AllMasks      = MVSscaling.Masks.AllMasks;
DataForPlotting.AllMasksNames = MVSscaling.Masks.AllMasksName;

DataForPlotting.ThresMaps     = MVSscaling.ThresMaps.AllMaps;
DataForPlotting.ThresMapsNames= MVSscaling.ThresMaps.AllMapsName;

DataForPlotting.EffectOfMRI_AALsplit.IndSplit      = IndSplitToBeUsed;
DataForPlotting.EffectOfMRI_AALsplit.SplitName     = MVSscaling.Masks.AALsplit.EffectOfMRI.SplitNames{IndSplitToBeUsed};
DataForPlotting.EffectOfMRI_AALsplit.AllMasks      = AllMasks;
DataForPlotting.EffectOfMRI_AALsplit.AllMasksNames = AllMasksNames;
DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks= NVoxelAllMasks;

%% exclude subparts that are too small?
%check if there are many splits with only few voxels i.e. <5% or 1% of total number of voxels
%OR if the maximum subpart/submask is smaller than 10%
%OR if there are just more than 5 subparts.
if((length(find((DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks./sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks))<0.05))>6) || (length(find((DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks./sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks))<0.01))>3) || (max(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks(:))./sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)<0.1) || length(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)>5)
    choice = questdlg({'There are many subparts of the Mask "EffectsOfMRI" that are small!'; ...
        [num2str(length(find((DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks./sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks))<0.05))),'of',num2str(length(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),' are smaller than 5%(',num2str(round(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*0.05)),'Voxels) of the whole Mask(',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Voxels).']; ...
        [num2str(length(find((DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks./sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks))<0.01))),'of',num2str(length(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),' are smaller than 1%(',num2str(round(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*0.01)),'Voxels) of the whole Mask(',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Voxels).']; ...
        ['The biggest subpart is ',num2str(max(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks(:))./sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*100),'%(',num2str(max(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks(:))),'Voxels) of the whole Mask(',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Voxels).']; ...
        ' ';
        'Do you want to limit the minimum size of the subparts that are displayed?'; ...
        'Options are: '; ...
        'Yes --> "enter percentage cutoff".';...
        ['T<5% --> "Exclude subparts smaller 5%(',num2str(round(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*0.05)),'Voxels) of the whole Mask(',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Voxels).']; ...
        'No --> "do not exclude any subpart". (However, areas will be sorted according to size.)';
        }, 'Exclude subparts by size?','Yes','T<5%','No','T<5%');
    switch(choice)
        case 'Yes'
            %ask user for input
            prompt   ={['Cutoff Percentage: [1%==',num2str(round(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*0.01)),'Voxels]']};
            dlg_title='Exclude Subparts of Mask "EffectsOfMRI"';
            num_lines=1;
            defAns   ={'1'};
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer   = inputdlg(prompt,dlg_title,num_lines,defAns,options);
            
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres = round(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*(eval(answer{1})/100));
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.PThres = [answer{1},'%(',num2str(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Total)'];
            IndsMasksDisp = find(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks>DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres);
            [tmp,SortInds] = sort(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks(IndsMasksDisp),'descend');
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp = IndsMasksDisp(SortInds);
            
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr= ['Only subparts (sorted by size) of Mask "EffectsOfMRI" are included that are bigger than ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.PThres];
        case 'T<5%'
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres = round(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks).*0.05);
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.PThres = ['5%(',num2str(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Total)'];
            IndsMasksDisp = find(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks>DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres);
            [tmp,SortInds] = sort(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks(IndsMasksDisp),'descend');
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp = IndsMasksDisp(SortInds);
            
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr= ['Only subparts (sorted by size) of Mask "EffectsOfMRI" are included that are bigger than ',DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.PThres];
        case 'No'
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres = 0;
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.PThres = ['0%(',num2str(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Total)'];
            [tmp,SortInds] = sort(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks,'descend');
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp = SortInds;
            
            DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr= 'All subparts (sorted by size) of Mask "EffectsOfMRI" are included.';
    end
else
    DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres = 0;
    DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.PThres = ['0%(',num2str(DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.NThres),'VoxelsOf',num2str(sum(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks)),'Total)'];
    [tmp,SortInds] = sort(DataForPlotting.EffectOfMRI_AALsplit.NVoxelAllMasks,'descend');
    DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.IndsMasksDisp = SortInds;
    
    DataForPlotting.EffectOfMRI_AALsplit.MasksToBeUsed.InfoStr= 'All subparts (sorted by size) of Mask "EffectsOfMRI" are included. (Sorted by size)';
end

%% assemble data according to choice
switch(choice_plotting)
    case 'GroupsAverage'
        DataForPlotting.PlotType      = 'GroupsAverage';
        for IndGroup = 1:MVSscaling.Design.NGroups
            DataForPlotting.NSubjPerGroup(IndGroup,1) = length(find(MVSscaling.Design.GroupNrPerSubj==IndGroup));
        end
        DataForPlotting.ScalingData   = zeros(size(MVSscaling.ScalingPerSubject,1),3);
        DataForPlotting.AmplitudesData= zeros(size(MVSscaling.ScalingPerSubject,1),MVSscaling.Design.N_MRIs,3);
        DataForPlotting.DataTitle     = cell(3,1);
        
        for IndGroup = 0:2
            if(IndGroup==0)
                SubjIndices = (MVSscaling.Design.GroupNrPerSubj>0);
                DataForPlotting.DataTitle{IndGroup+1} = 'AllSubjects';
            else
                SubjIndices = (MVSscaling.Design.GroupNrPerSubj==IndGroup);
                DataForPlotting.DataTitle{IndGroup+1} = ['Group ',num2str(IndGroup)];
            end
            
            if(UseAbs)
                InputData = abs(MVSscaling.ScalingPerSubject(:,SubjIndices));
            else
                InputData =     MVSscaling.ScalingPerSubject(:,SubjIndices);
            end
            if(strcmp(choice_average,'Median'))
                DataForPlotting.ScalingData(:,IndGroup+1)     = median(InputData,2);
                DataForPlotting.AmplitudesData(:,:,IndGroup+1)= squeeze(median(MVSscaling.AmplitudesPerSubjectMRI(:,SubjIndices,:),2));
            else
                DataForPlotting.ScalingData(:,IndGroup+1)     = mean(InputData,2);
                DataForPlotting.AmplitudesData(:,:,IndGroup+1)= squeeze(mean(MVSscaling.AmplitudesPerSubjectMRI(:,SubjIndices,:),2));
            end
        end
    case 'GroupsAggregate'
        DataForPlotting.PlotType      = 'GroupsAggregate';
        for IndGroup = 1:MVSscaling.Design.NGroups
            DataForPlotting.NSubjPerGroup(IndGroup,1) = length(find(MVSscaling.Design.GroupNrPerSubj==IndGroup));
        end
        DataForPlotting.DataTitle     = cell(3,1);
        DataForPlotting.DataGrouping  = cell(3,1);
        
        for IndGroup = 0:2
            if(IndGroup==0)
                SubjIndices = (MVSscaling.Design.GroupNrPerSubj>0);
                DataForPlotting.DataTitle{IndGroup+1} = 'AllSubjects';
            else
                SubjIndices = (MVSscaling.Design.GroupNrPerSubj==IndGroup);
                DataForPlotting.DataTitle{IndGroup+1} = ['Group ',num2str(IndGroup)];
            end
            
            DataForPlotting.DataGrouping{IndGroup+1}= SubjIndices;
        end
        
        if(UseAbs)
            DataForPlotting.ScalingData             = abs(MVSscaling.ScalingPerSubject);
        else
            DataForPlotting.ScalingData             =     MVSscaling.ScalingPerSubject;
        end
        DataForPlotting.AmplitudesData          = permute(MVSscaling.AmplitudesPerSubjectMRI, [1 3 2]);
    case 'SingleSubjects'
        DataForPlotting.PlotType      = 'SingleSubjects';
        DataForPlotting.ScalingData   = zeros(size(MVSscaling.ScalingPerSubject,1),size(MVSscaling.ScalingPerSubject,2));
        DataForPlotting.AmplitudesData= zeros(size(MVSscaling.ScalingPerSubject,1),MVSscaling.Design.N_MRIs,size(MVSscaling.ScalingPerSubject,2));
        DataForPlotting.DataTitle     = cell(size(MVSscaling.ScalingPerSubject,2),1);
        
        %sort by group to make it easier to view
        for IndGroup = 1:MVSscaling.Design.NGroups
            DataForPlotting.NSubjPerGroup(IndGroup,1) = length(find(MVSscaling.Design.GroupNrPerSubj==IndGroup));
        end
        [B,IX] = sort(MVSscaling.Design.GroupNrPerSubj);
        
        for IndSubj = 1:length(IX)
            SubjInd = IX(IndSubj);
            DataForPlotting.DataTitle{IndSubj}          = ['S',num2str(SubjInd),'(Gr',num2str(MVSscaling.Design.GroupNrPerSubj(SubjInd)),')'];
            if(UseAbs)
                DataForPlotting.ScalingData(:,IndSubj)  = abs(MVSscaling.ScalingPerSubject(:,SubjInd));
            else
                DataForPlotting.ScalingData(:,IndSubj)  =     MVSscaling.ScalingPerSubject(:,SubjInd);
            end
            DataForPlotting.AmplitudesData(:,:,IndSubj) = squeeze(MVSscaling.AmplitudesPerSubjectMRI(:,SubjInd,:));
        end
end


end