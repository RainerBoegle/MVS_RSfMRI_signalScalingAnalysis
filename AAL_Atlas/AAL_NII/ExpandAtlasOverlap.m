function [OutputStruct]=ExpandAtlasOverlap(varargin)
% This function will take the overlap of a Mask and an Atlas (i.e. Mask.*AAL atlas or other) and
% check all voxels of the original Mask that have not been assigned for their neighbors' values
% and assign the most often found label in the neighborhood the label of that voxel and reiterating again.
%
% The extend of the neighbors searched can be changed but tests have shown
% that nearest neighbors give best results without changing assignments too
% much from expected results, especially if spacing is sufficient.
%
% The algorithm either stops if all voxels are assigned or if the maximum
% iterations, which is 10-times-the number of voxels (default) to be assign, 
% is reached or 20 consequtive iterations (default) didn't produce any new added voxels.
%
% In a second step, AAL Atlas areas are combined.
% The "Basic" combination is to search for areas with equal label while ignoring 
% "L"eft or "R"ight side, e.g. "Thalamus_L" & "Thalamus_R" are just labeled "Thalamus"
% and also numbered subparts are combined, e.g. Vermis1 ... Vermis4 is just considered as Vermis.
% NB: This leaves "Temporal_Pole_Mid" and "Temporal_Pole_Sup" and such separated!
%
% The "General" combination further combines the labels from "basic", that are at the top of the hierarchy, 
% i.e. "Temporal_Pole_Mid" and "Temporal_Pole_Sup" are just labeled as "Temporal"
%
%Usage:
%      [OutputStruct]=ExpandAtlasOverlap(ShowWarnings,MaxIterFactor,NCyclesNoChangeMax); Default: Don't show warnings and Maximum Iteration Factor == 10 and NCyclesNoChangeMax == 20
%      [OutputStruct]=ExpandAtlasOverlap(1);        %Show warnings & Iterate  10-times-NVoxels to assign if not 20 cycles are unchanged.  
%      [OutputStruct]=ExpandAtlasOverlap(1,100);    %Show warnings & Iterate 100-times-NVoxels to assign if not 20 cycles are unchanged. 
%      [OutputStruct]=ExpandAtlasOverlap(1,100,10); %Show warnings & Iterate 100-times-NVoxels to assign if not 10 cycles are unchanged. 
%
%
%NB: use [I1,I2,I3] = ind2sub(OutputStruct.Dim,OutputStruct.CurrentIndices); 
%    or  [I1,I2,I3] = ind2sub(OutputStruct.Dim,OutputStruct.IndsOverstepped);
%    to get 3D voxel-coordinates of the unassigned voxels.
%
%Date: 15.08.2014 (V1)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

if(nargin==1)
    ShowWarnings       = varargin{1}; %show warnings?
    FactorMaxIterations= 10.5;        %go over voxels for a maximum of 10.5 times
    NCyclesNoChangeMax = 20;          %How often may n voxels not assigned be unchanged before we give up.
else
    if(nargin==2)
        ShowWarnings       = varargin{1}; %show warnings?
        FactorMaxIterations= varargin{2}; %go over voxels for a maximum of "Input" times
        NCyclesNoChangeMax = 20;          %How often may n voxels not assigned be unchanged before we give up.
    else
        if(nargin==3)
            ShowWarnings       = varargin{1}; %show warnings?
            FactorMaxIterations= varargin{2}; %go over voxels for a maximum of "Input" times
            NCyclesNoChangeMax = varargin{3}; %How often may n voxels not assigned be unchanged before we give up.
        else
            disp('Using Defaults...');
            ShowWarnings       = 0;    %show warnings?
            FactorMaxIterations= 10.5; %go over voxels for a maximum of 10.5 times
            NCyclesNoChangeMax = 20;   %How often may n voxels not assigned be unchanged before we give up.
        end
    end
end
%if iteration is choosen, let's remove small areas after nearest neighbor analysis step, such that next orders can assign this area according to higher orders
NMin = 27; %Voxel plus nearest neighbors, if area at end of iteration with NearestNeighbors is smaller than that, then remove it, but keep searching for area nearby that might replace it.


%% get Overlap_MaskEffectsOfMRI_AAL NIFIT
Overlap_MaskEffectsOfMRI_AAL_Path = spm_select(1,'image','Select NIFTI Overlap of EffectOfMRI with AAL atlas...');
NII_Overlap = nifti(Overlap_MaskEffectsOfMRI_AAL_Path);
Overlap_dat = round(NII_Overlap.dat(:,:,:));
OverlapOrg_dat = Overlap_dat;

%% get Mask_EffectsOfMRI NIFTI
Mask_EffectsOfMRI_Path = spm_select(1,'image','Select NIFTI Mask of EffectOfMRI...');
NII_MaskEffectsMRI = nifti(Mask_EffectsOfMRI_Path);
MaskEffectsMRI_dat = NII_MaskEffectsMRI.dat(:,:,:);

%% determine indices where Mask_EffectsOfMRI==1 & Overlap_MaskEffectsOfMRI_AAL==0, i.e. the voxels that were missed
Indices_org = find(MaskEffectsMRI_dat~=0 & Overlap_dat==0);
Indices_Backup = Indices_org; %for check

%% get AAL labels for later checks
load('AAL_LabelsOnly.mat');
EffectOfMRI_AALsplit.LabelsAAL = Labels;
clear Labels

%% go over 3D voxel coordinates and search in the neighborhood if there is some structure nearby
%i.e. 
% i.) collect the values of the NEAREST neighborhood from the atlas (at first Overlap_MaskEffectsOfMRI_AAL) that are NONZERO and assign
%     the mode to the values if not empty.
% ii.)remove the assigned voxel from the list and update the atlas
% repeat from i.) until last voxel is assigned or nearest neighbors are all zero --> warning.

choice_search = questdlg({'How deep should the search for neighbors go?'; 'Nearest Neighbors only or manually define order?'; '[NB: increasing the order beyond 2 will take very long!]'},'How far should the search be extended?','NearestNeighbors','ManualNeighbors','IterateOrders','ManualNeighbors');
switch(choice_search)
    case 'NearestNeighbors'
        NOrd = 1;
    case {'ManualNeighbors','IterateOrders'}
        %ask user for order
        if(strcmp(choice_search,'ManualNeighbors'))
            prompt   ={['Order: ']};
        else
            prompt   ={['MaxOrder: ']};
        end
        dlg_title='NHood-Order';
        num_lines=1;
        defAns   ={'2'};
        
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        answer   = inputdlg(prompt,dlg_title,num_lines,defAns,options);
        
        if(strcmp(choice_search,'ManualNeighbors'))
            NOrd = eval(answer{1});
        else
            NOrd = 1:eval(answer{1});
            choice_remoutliers = questdlg('Remove outliers, i.e. unassigned voxels after NearestNeighbor iteration-step?','Remove Outliers?','Yes','No','Yes');
            switch(choice_remoutliers)
                case 'Yes'
                    RemOutliers = 1;
                otherwise
                        RemOutliers = 0;
            end
        end
end

disp(' ');
disp('Trying to expand map of assigned AAL areas...');
if(strcmp(choice_search,'IterateOrders'))
    NAllIterations = zeros(length(NOrd),1);
else
    NAllIterations = 0;
end
for CurrOrd = NOrd
    disp(['Currently searching in ',num2str(CurrOrd),'.Neighborhood...']);
    Dim                = size(Overlap_dat);
    CurrentIndices     = Indices_org;
    NIndsLeft          = length(Indices_org);
    NIndsLeftLast      = [];
    IndsOverstepped    = []; %collect voxel indices that could not be assigned, just in case it doesn't work out in the limit of maximum Iterations
    CurrentIteration   = 1; %just for keeping track
    StopNext           = 0; %stop on next iteration
    KeepGrowing        = 1;
    while(KeepGrowing)
        NIndsLeftLast  = [NIndsLeftLast; NIndsLeft]; %start with all
        if(~isempty(CurrentIndices))
            % get 3D voxel coordinates of the missing voxels
            [I1,I2,I3] = ind2sub(Dim,CurrentIndices(1));
            % collect nearest neighbors that are not zero
            CollectedVoxels = [];
            for Ind1 = max([I1-CurrOrd,1]):min([I1+CurrOrd,Dim(1)])
                for Ind2 = max([I2-CurrOrd,1]):min([I2+CurrOrd,Dim(2)])
                    for Ind3 = max([I3-CurrOrd,1]):min([I3+CurrOrd,Dim(3)])
                        if(Overlap_dat(Ind1,Ind2,Ind3)~=0)
                            CollectedVoxels = [CollectedVoxels; Overlap_dat(Ind1,Ind2,Ind3)];
                        end
                    end
                end
            end
            % decide based on the collected nearest neighbors
            if(~isempty(CollectedVoxels))
                Overlap_dat(I1,I2,I3) = mode(CollectedVoxels);
                if(any(CurrentIndices(1)==IndsOverstepped)) %this index was once overstepped, so remove it from that list
                    IndsOverstepped(find(CurrentIndices(1)==IndsOverstepped)) = [];
                end
                CurrentIndices(1) = []; %remove this voxel!
                NIndsLeft = NIndsLeft-1;
                
                if(length(unique(CollectedVoxels))>1) %several structures are in nearest neighborhood
                    NCollected = length(CollectedVoxels);
                    NMode      = length(find(mode(CollectedVoxels)==CollectedVoxels));
                    OtherPossibilities = zeros(length(unique(CollectedVoxels(CollectedVoxels~=mode(CollectedVoxels)))),1);
                    OtherPossibilities(:,1) = unique(CollectedVoxels(CollectedVoxels~=mode(CollectedVoxels)));
                    NOthers = length(CollectedVoxels(CollectedVoxels~=mode(CollectedVoxels)));
                    if(ShowWarnings)
                        disp(['WARNING! Assigning mode(CollectedVoxels)= ',num2str(mode(CollectedVoxels)),'(',num2str(NMode),'of',num2str(NCollected),'/27max) to (',num2str([I1,I2,I3]),'), but ',num2str(OtherPossibilities'),'(combined ',num2str(NOthers),'of',num2str(NCollected),'/27) is also possible.']);
                    end
                end
            else %no neighbors belonging to any structure in the atlas at current iteration
                if(ShowWarnings)
                    disp(['WARNING! Can not find any neighbors for (',num2str([I1,I2,I3]),')!!! Will move it to the end of the list and continue.']);
                end
                if(length(CurrentIndices)>=2)
                    tmpCurrInds = [CurrentIndices(2:end); CurrentIndices(1)]; %move this voxel to the end of list to be evaluated again
                    IndsOverstepped = [IndsOverstepped; CurrentIndices(1)]; %add this index to the list of indices that have been overstepped.
                else %special case of only one left
                    if(any(IndsOverstepped==CurrentIndices(1))) %only one left and has been overstepped before, check if StopNext has been set and stop if set
                        if(~StopNext)
                            StopNext = 0; %stop on next iteration
                        else
                            if(ShowWarnings)
                                disp('This has happened too often for my taste. Let''s end this now.');
                            end
                            KeepGrowing = 0;
                        end
                    end
                    tmpCurrInds     = CurrentIndices(1);
                    IndsOverstepped = CurrentIndices(1);
                end
                CurrentIndices = tmpCurrInds; %assign moved voxels
            end
            % keep going as long as there are voxels to assign and as long as
            % this hasn't gone for too long. i.e. FactorMaxIterations-times the number of voxels to assign
            if(CurrentIteration<=ceil(FactorMaxIterations.*length(Indices_org))) %don't keep on going if you have been around more than FactorMaxIterations-times
                CurrentIteration = CurrentIteration+1;
            else %max number of tries to assign all voxels reached, let's end this.
                disp('Maximum number of iterations reached. Stopping now.');
                KeepGrowing = 0;
            end
            
            [ModeMultiInds,NModeMultiInds] = mode(NIndsLeftLast);
            if(NModeMultiInds==NCyclesNoChangeMax) %for twenty times we did not change the number of assigned voxels! Let's stop.
                disp(['No change in number of voxels to be assigned for ',num2str(NCyclesNoChangeMax),' cycles. Stopping now.']);
                KeepGrowing = 0;
            end
        else %all assigned
            disp('Finished iterating and all missing Voxels have been assigned!');
            disp('Well done, let''s check the results.');
            if(strcmp(choice_search,'IterateOrders'))
                NAllIterations(CurrOrd) = CurrentIteration;
            else
                NAllIterations = CurrentIteration;
            end
            KeepGrowing = 0;
        end
    end    
    if(~isempty(CurrentIndices))
        if(strcmp(choice_search,'IterateOrders'))
            if(CurrOrd==1)
                if(RemOutliers)
                    for Ind = 1:length(CurrentIndices)
                        Indices_org(Indices_org==CurrentIndices(Ind)) = []; %remove those that are outliers in first round.
                    end
                end
                for Ind = 1:length(EffectOfMRI_AALsplit.LabelsAAL)
                    if(~isempty(find(Overlap_dat==Ind)))
                        if(length(find(Overlap_dat==Ind))<NMin) %if for Label number there are only less than minimum number, then remove those to give higher interations a chance to associate this with another area.
                            Overlap_dat(Overlap_dat==Ind) = 0;
                        end
                    else
                        continue;
                    end
                end
            end
            NAllIterations(CurrOrd) = CurrentIteration;
        else
            NAllIterations = CurrentIteration;
        end
        disp(['Stopped iterating after ',num2str(CurrentIteration),' Iterations.']);
        disp(['Still ',num2str(length(CurrentIndices')),' Voxels left unassigned!']);
        disp(['Indices of Voxels are (',num2str(unique(IndsOverstepped)'),')']);
        IndsOverstepped=unique(IndsOverstepped);
    end
    disp(' ');
end
disp(' ');
if(~isempty(CurrentIndices))
    disp(['Stopped iterating (NHood-Mode: "',choice_search,'" NOrdMax= ',num2str(max(NOrd)),') after ',num2str(sum(NAllIterations)),' Iterations (',num2str(NAllIterations'),').']);
    disp(['Still ',num2str(length(CurrentIndices')),' Voxels left unassigned!']);
    disp(['Indices of Voxels are (',num2str(unique(IndsOverstepped)'),')']);
    IndsOverstepped=unique(IndsOverstepped);
    disp(' ');
else
    disp('Finished iterating and all missing Voxels have been assigned!');
    disp(['NHood-Mode: "',choice_search,'" (NOrdMax= ',num2str(max(NOrd)),') Used ',num2str(sum(NAllIterations)),' Iterations in total (',num2str(NAllIterations'),').']);
    disp('Well done, let''s check the results.');
    disp(' ');
end


%% write out results fname = ['Expanded_',Overlap_MaskEffectsOfMRI_AAL];
[BasePath,OrgFName,Ext] = fileparts(Overlap_MaskEffectsOfMRI_AAL_Path);

Vo = spm_vol(Overlap_MaskEffectsOfMRI_AAL_Path);
switch(choice_search)
    case 'NearestNeighbors'
        Vo.fname = [BasePath,filesep,'Expanded_NN_',OrgFName,'.nii'];
        OutName  = ['Expanded_NN_',OrgFName,'.nii'];
    case {'ManualNeighbors','IterateOrders'}
        Vo.fname = [BasePath,filesep,'Expanded_',choice_search,'_',num2str(max(NOrd)),'N_',OrgFName,'.nii'];
        OutName  = ['Expanded_',choice_search,'_',num2str(max(NOrd)),'N_',OrgFName,'.nii'];
end

Vo = spm_write_vol(Vo, round(Overlap_dat));

disp(' ');
disp(['Results have been written out as "',OutName,'" to "',BasePath,'"']);
if(~isempty(CurrentIndices)) %some are missing so let's show them separately
    MissingVoxels = zeros(size(Overlap_dat));
    MissingVoxels(CurrentIndices) = 1;
    
    VoNotAssigned = Vo;
    switch(choice_search)
        case 'NearestNeighbors'
            VoNotAssigned.fname = [BasePath,filesep,'StillMissingVoxels_Expanded_NN_',OrgFName,'.nii'];
            OutName  = ['StillMissingVoxels_Expanded_NN_',OrgFName,'.nii'];
        case {'ManualNeighbors','IterateOrders'}
            VoNotAssigned.fname = [BasePath,filesep,'StillMissingVoxels_Expanded_',choice_search,'_',num2str(max(NOrd(:))),'N_',OrgFName,'.nii'];
            OutName  = ['StillMissingVoxels_Expanded_',choice_search,'_',num2str(max(NOrd(:))),'N_',OrgFName,'.nii'];
    end    
    
    VoNotAssigned = spm_write_vol(VoNotAssigned, MissingVoxels);
    disp(['The voxels that are still missing have been written out as "',OutName,'" to "',BasePath,'"']);
end


%% reduce number of AALlabels and display what is included
disp(' ');
% EffectOfMRI_AALsplit.LabelsRed = Labels;
LabelsCombinedBasic   = cell(size(EffectOfMRI_AALsplit.LabelsAAL));
LabelsCombinedGeneral = cell(size(EffectOfMRI_AALsplit.LabelsAAL));
for IndLabel = 1:length(EffectOfMRI_AALsplit.LabelsAAL)
    [start_idx, end_idx, extents, matches, tokens, names, splits] = regexp(EffectOfMRI_AALsplit.LabelsAAL{IndLabel},'_'); %split current label by removing underscore "_"
    IndsAcceptedSplits = []; %collect splits that are acceptable then put them back together
    for IndSplit = 1:length(splits)
        if(~(strcmp(splits{IndSplit},'L')||strcmp(splits{IndSplit},'R')||~isempty(regexp(splits{IndSplit},'\d','ONCE')))) %this is what we want to avoid
            IndsAcceptedSplits = [IndsAcceptedSplits;IndSplit];
        end
    end
    
    NewLabel = []; %put label back together without L/R and so on..., then reassign
    if(length(IndsAcceptedSplits)>1) %fill up including spaces
        for IndSplit = 1:length(IndsAcceptedSplits)
            if(IndSplit~=length(IndsAcceptedSplits))
                NewLabel = [NewLabel,splits{IndsAcceptedSplits(IndSplit)},' '];
            else
                NewLabel = [NewLabel,splits{IndsAcceptedSplits(IndSplit)}];
            end
        end
    else
        NewLabel = splits{IndsAcceptedSplits(1)};
    end
    
    %display and assign
%     disp(['AAL: "',EffectOfMRI_AALsplit.LabelsAAL{IndLabel},'" --> "',NewLabel,'"']);
    LabelsCombinedBasic{IndLabel}   = NewLabel;
    LabelsCombinedGeneral{IndLabel} = splits{1}; %Label highest in Hierarchy
end
% disp(' ');
% check Labels for the equal titles and collect indices such that we can reduce number of masks
MatchingLabelsBasic = cell(size(LabelsCombinedBasic));
for IndLabel = 1:length(LabelsCombinedBasic)
    MatchLabels = [];
    for IndLCheck = 1:length(LabelsCombinedBasic)
        if(strcmp(LabelsCombinedBasic{IndLabel},LabelsCombinedBasic{IndLCheck}))
            MatchLabels = [MatchLabels; IndLCheck];
        end
    end
    MatchingLabelsBasic{IndLabel} = sort(MatchLabels);
%     disp([num2str(IndLabel,'%03.0f'),'"',LabelsCombinedBasic{IndLabel},'" found in (',num2str(MatchingLabelsBasic{IndLabel}'),')']);
end
MatchingLabelsGeneral = cell(size(LabelsCombinedGeneral));
for IndLabel = 1:length(LabelsCombinedGeneral)
    MatchLabels = [];
    for IndLCheck = 1:length(LabelsCombinedGeneral)
        if(strcmp(LabelsCombinedGeneral{IndLabel},LabelsCombinedGeneral{IndLCheck}))
            MatchLabels = [MatchLabels; IndLCheck];
        end
    end
    MatchingLabelsGeneral{IndLabel} = sort(MatchLabels);
%     disp([num2str(IndLabel,'%03.0f'),'"',LabelsCombinedBasic{IndLabel},'" found in (',num2str(MatchingLabelsBasic{IndLabel}'),')']);
end

EffectOfMRI_AALsplit.LabelsCombinedBasic   = LabelsCombinedBasic;
EffectOfMRI_AALsplit.LabelsCombinedGeneral = LabelsCombinedGeneral;
EffectOfMRI_AALsplit.MatchingLabelsBasic   = MatchingLabelsBasic;
EffectOfMRI_AALsplit.MatchingLabelsGeneral = MatchingLabelsGeneral;

%% combine areas
% setup masks according to unique masks that are of the same reduced type,
% i.e. search for indices that belong to the same type of area and replace
% them in the mask. At the end we only have the minimal number of labels
% used.

%AAL split of EffectOfMRI
NII_tmp                 = nifti(Vo.fname);
tmp_AALsplit_Input      = round(NII_tmp.dat(:,:,:));
clear NII_tmp

% assign masks according to minimal set of areas
AALsplitLabelNums = unique(round(tmp_AALsplit_Input(:)));
if(any(AALsplitLabelNums==0))
    AALsplitLabelNums(AALsplitLabelNums==0) = [];
end

%% Basic
disp(['Before "Basic": ',num2str(unique(tmp_AALsplit_Input(tmp_AALsplit_Input~=0))')]);
AALsplit_CombineBasic = tmp_AALsplit_Input;
IndicesToCheck = AALsplitLabelNums;
while(~isempty(IndicesToCheck))
    %check & replace
    CurrentInds = MatchingLabelsBasic{IndicesToCheck(1)};
    for Ind = 1:length(CurrentInds)
        if(1) %(any(AALsplit_CombineBasic==CurrentInds(Ind)))
            AALsplit_CombineBasic(AALsplit_CombineBasic==CurrentInds(Ind)) = IndicesToCheck(1); %replace with first from list
        end
    end
    for Ind = 1:length(CurrentInds)
        %throw what has been matched already
        if(any(IndicesToCheck==CurrentInds(Ind)))
            IndicesToCheck(IndicesToCheck==CurrentInds(Ind)) = [];
        end
    end
end
disp(['After "Basic": ',num2str(unique(AALsplit_CombineBasic(AALsplit_CombineBasic~=0))')]);

%% General
disp(['Before "General": ',num2str(unique(tmp_AALsplit_Input(tmp_AALsplit_Input~=0))')]);
AALsplit_CombineGeneral = tmp_AALsplit_Input;
IndicesToCheck = AALsplitLabelNums;
while(~isempty(IndicesToCheck))
    %check & replace
    CurrentInds = MatchingLabelsGeneral{IndicesToCheck(1)};
    for Ind = 1:length(CurrentInds)
        if(1) %(any(AALsplit_CombineGeneral==CurrentInds(Ind)))
            %disp([num2str(CurrentInds(Ind)),' = ',num2str(IndicesToCheck(1))]);
            AALsplit_CombineGeneral(AALsplit_CombineGeneral==CurrentInds(Ind)) = IndicesToCheck(1); %replace with first from list
        end
    end
    %throw what has been matched already
    for Ind = 1:length(CurrentInds)
        if(any(IndicesToCheck==CurrentInds(Ind)))
            IndicesToCheck(IndicesToCheck==CurrentInds(Ind)) = [];
        end
    end
end
disp(['After "General": ',num2str(unique(AALsplit_CombineGeneral(AALsplit_CombineGeneral~=0))')]);

%% write out the data
[BasePath,OrgFName,Ext] = fileparts(Vo.fname);

VoBasic = spm_vol(Overlap_MaskEffectsOfMRI_AAL_Path);
VoBasic.fname = [BasePath,filesep,'CombineBasic_',OrgFName,'.nii'];
OutName  = ['CombineBasic_',OrgFName,'.nii'];
disp(' ');
disp(['Writing out "',OutName,'"...']);
VoBasic = spm_write_vol(VoBasic, round(AALsplit_CombineBasic));
disp('done.');

VoGeneral = spm_vol(Overlap_MaskEffectsOfMRI_AAL_Path);
VoGeneral.fname = [BasePath,filesep,'CombineGeneral_',OrgFName,'.nii'];
OutName  = ['CombineGeneral_',OrgFName,'.nii'];
disp(' ');
disp(['Writing out "',OutName,'"...']);
VoGeneral = spm_write_vol(VoGeneral, round(AALsplit_CombineGeneral));
disp('done.');

%% display the splits
%check for the unique ones and count what is in it
for Ind = 1:3
    switch(Ind)
        case 1
            Data  = Overlap_dat;
            Labels= EffectOfMRI_AALsplit.LabelsAAL;
            Title = 'Expanded Overlap "EffectOfMRI" with AAL Atlas';
        case 2
            Data  = AALsplit_CombineBasic;
            Labels= LabelsCombinedBasic;
            Title = 'Combine "BASIC" of expanded Overlap';
        case 3
            Data  = AALsplit_CombineGeneral;
            Labels= LabelsCombinedGeneral;
            Title = 'Combine "GENERAL" of expanded Overlap';
    end
    
    EffectOfMRI_AALsplit.Info(Ind).Title = Title;
    UniqueAALsplitLabelNums = unique(round(Data(:)));
    if(any(UniqueAALsplitLabelNums==0))
        UniqueAALsplitLabelNums(UniqueAALsplitLabelNums==0) = [];
    end
%     disp(['Ind ',num2str(Ind),': ',num2str(UniqueAALsplitLabelNums'),' "',Title,'".']);
    
    N_AALsplitMasks          = length(UniqueAALsplitLabelNums);
    EffectOfMRI_AALsplit.Info(Ind).AllMasksName   = cell(N_AALsplitMasks,1);
    EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks = zeros(N_AALsplitMasks,1);
    for IndSplitMasks = 1:N_AALsplitMasks
        EffectOfMRI_AALsplit.Info(Ind).AllMasksName{IndSplitMasks}   = Labels{UniqueAALsplitLabelNums(IndSplitMasks)};
        EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks(IndSplitMasks) = length(find(Data(:)==UniqueAALsplitLabelNums(IndSplitMasks)));
    end
    NAll = sum(EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks(:));
    
    %check mask sizes/voxel counts
    %NB some might have gone missing by expansion process before not being able to assign them due to nearest neighbors.
    InfoTxt{1} = ['Split of "EffectOfMRI": ',num2str(NAll),'Voxels'];
    InfoTxt{2} =  '------------------------------------------------------';
    InfoTxt{3} = ['All (',num2str(N_AALsplitMasks),') in order of Mask number'];
    for IndMask = 1:N_AALsplitMasks
        InfoTxt{3+IndMask} = [EffectOfMRI_AALsplit.Info(Ind).AllMasksName{IndMask},': ',num2str(EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks(IndMask)),'Voxels (',num2str((EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks(IndMask)/NAll)*100),'%)'];
    end
    InfoTxt{3+length(EffectOfMRI_AALsplit.Info(Ind).AllMasksName)+1} = '...............................................';
    InfoTxt{3+length(EffectOfMRI_AALsplit.Info(Ind).AllMasksName)+2} = ['All (',num2str(N_AALsplitMasks),') sorted by voxel number'];
    [Tmp_Sort,IndicesSort] = sort(EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks,'descend');
    for IndMask = 1:N_AALsplitMasks
        InfoTxt{3+length(EffectOfMRI_AALsplit.Info(Ind).AllMasksName)+2+IndMask} = [EffectOfMRI_AALsplit.Info(Ind).AllMasksName{IndicesSort(IndMask)},': ',num2str(EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks(IndicesSort(IndMask))),'Voxels (',num2str((EffectOfMRI_AALsplit.Info(Ind).NVoxelAllMasks(IndicesSort(IndMask))/NAll)*100),'%)'];
    end
    EffectOfMRI_AALsplit.Info(Ind).InfoTxt = InfoTxt;
    helpdlg(InfoTxt,Title);
    clear InfoTxt
end

%% create output struct
OutputStruct.choice_search   = choice_search;
OutputStruct.NOrd            = NOrd;
OutputStruct.Vo              = Vo;
OutputStruct.VoBasic         = VoBasic;
OutputStruct.VoGeneral       = VoGeneral;
OutputStruct.CurrentIndices  = CurrentIndices;
if(~isempty(CurrentIndices)) %some are missing so let's show them separately
   OutputStruct.VoNotAssigned= VoNotAssigned;
end
OutputStruct.NAllIterations  = NAllIterations;
OutputStruct.Dim             = Dim;
OutputStruct.IndsOverstepped = IndsOverstepped;
OutputStruct.Indices_org     = Indices_org;
if(length(Indices_org)~=length(Indices_Backup))
    OutputStruct.Indices_Backup = Indices_Backup;
end

OutputStruct.EffectOfMRI_AALsplit = EffectOfMRI_AALsplit;

%% save labels back
clear Labels
Labels(1).Labels = EffectOfMRI_AALsplit.LabelsAAL; 
Labels(2).Labels = LabelsCombinedBasic;
Labels(3).Labels = LabelsCombinedGeneral;

Labels(1).Name   = 'AAL-Labels';
Labels(2).Name   = 'Combine(Basic) AAL-Labels';
Labels(3).Name   = 'Combine(General) AAL-Labels';

save('AAL_Labels.mat','Labels');

end