function Suggestions = SuggestMergers(EffectOfMRI_AALsplit,IndSplit,SearchDist)
% This function looks through a given split "IndSplit" for a given distance
% in Voxels "SearchDist" and suggest mergers of areas that are in this
% search range distance.

%% get map of atlas
NII_tmp        = nifti(EffectOfMRI_AALsplit.MPaths{IndSplit});
AALsplit_Input = round(NII_tmp.dat(:,:,:));
IndsAreas      = unique(AALsplit_Input);
if(any(IndsAreas==0))
    IndsAreas(IndsAreas==0) = [];
end
Dim            = size(AALsplit_Input);
clear NII_tmp

Repeat = 1; %in case we are not successful
while(Repeat)
    %% setup Suggestions-struct
    Suggestions.IndSplit            = IndSplit; %keep track of which one we looked at.
    Suggestions.VoxelDistForSearch  = SearchDist; %how far to look (in voxels) for other areas than the current area
    %Old: (now will grow such that we always end with least number of mergers)
    % Suggestions.IndsMergers         = cell(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1);  %each one gets a suggestion, at least staying itself.
    % Suggestions.NVoxelMergers       = cell(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1);  %for each one suggestion, note the size.
    % Suggestions.NVoxelTotalMergers  = zeros(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1); %the number of voxels if merged, such that we can sort later
    % Suggestions.NMergers            = zeros(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1); %keep track of how many areas are merged.
    
    %% do search
    N_NeighborsToSearch = Suggestions.VoxelDistForSearch; %how far to look for other areas than the current area
    [SortedList,IndsAscending] = sort(EffectOfMRI_AALsplit.Split(IndSplit).NVoxelAllMasks); %'ascend'
    SortedAreaInds = IndsAreas(IndsAscending); %indices of search from smallest area to largest
    KeepGoing = 1; %init true for start
    CurrInd   = 1; %start at 1 of course
    while(KeepGoing)
        if(~isempty(SortedAreaInds)) %still something to do
            CurrentAreaInd = SortedAreaInds(1); %start here, if nothing else is found then at least this index is to be noted
            CollectedVoxels     = [CurrentAreaInd]; %init: collect nearest neighbors that are not zero & not the current index
            CollectedVoxelsLast = []; %init empty
            KeepSearching = 1; %init true
            while(KeepSearching) %keep growing not possible any more
                if(~isempty(CollectedVoxelsLast)) %search has been done let's compare and remove to search only new ones
                    CollectedVoxelsNew = CollectedVoxelsLast;
                    for Ind = 1:length(CollectedVoxels)
                        CollectedVoxelsNew(CollectedVoxelsNew==CollectedVoxels(Ind))=[]; %remove already collected --> if this doesn't leave any new voxels, then we stop.
                    end
                    CollectedVoxels = CollectedVoxelsLast; %update already collected.
                else %first time around --> init.
                    CollectedVoxelsLast = CollectedVoxels; %init
                    CollectedVoxelsNew  = CollectedVoxels; %init
                end
                
                if(~isempty(CollectedVoxelsNew)) %search if there are indices to look for
                    IndsCurrentSearch = []; %init
                    for Ind = 1:length(CollectedVoxelsNew)
                        IndsCurrentSearch = [IndsCurrentSearch; find(AALsplit_Input==CollectedVoxelsNew(Ind))]; %init indices for the current search
                    end
                    [I1,I2,I3] = ind2sub(Dim,IndsCurrentSearch); % get 3D voxel coordinates of the voxels belonging to current area
                    for IndSearch = 1:length(IndsCurrentSearch)
                        for Ind1 = max([I1(IndSearch)-N_NeighborsToSearch,1]):min([I1(IndSearch)+N_NeighborsToSearch,Dim(1)])
                            for Ind2 = max([I2(IndSearch)-N_NeighborsToSearch,1]):min([I2(IndSearch)+N_NeighborsToSearch,Dim(2)])
                                for Ind3 = max([I3(IndSearch)-N_NeighborsToSearch,1]):min([I3(IndSearch)+N_NeighborsToSearch,Dim(3)])
                                    if(((AALsplit_Input(Ind1,Ind2,Ind3)~=0) && (AALsplit_Input(Ind1,Ind2,Ind3)~=CurrentAreaInd)))
                                        CollectedVoxelsLast = [CollectedVoxelsLast; AALsplit_Input(Ind1,Ind2,Ind3)];
                                    end
                                end
                            end
                        end
                    end
                    CollectedVoxelsLast = unique(CollectedVoxelsLast); %who are the unique neighbors?
                else
                    KeepSearching = 0; %nothing new to search, so let's stop
                end
            end %can't grow any more, so the current search is over. Note results and remove indices that have been merged from list of possible mergers
            %Note results of this merge
            N_Members = zeros(length(CollectedVoxels),1); %number of members, i.e. voxels for the indices
            for IndAreas = 1:length(CollectedVoxels)
                N_Members(IndAreas) = length(find(AALsplit_Input==CollectedVoxels(IndAreas)));
            end
            [TmpSort,IndsAreasNDescending] = sort(N_Members,'descend');
            Suggestions.IndsMergers{CurrInd}        = CollectedVoxels(IndsAreasNDescending); %assign highes number as first, i.e. the label to be named frist and descending for merger
            Suggestions.NVoxelMergers{CurrInd}      = TmpSort;
            Suggestions.NVoxelTotalMergers(CurrInd) = sum(N_Members);
            Suggestions.NMergers(CurrInd)           = length(CollectedVoxels);
            
            for Ind = 1:length(CollectedVoxels)
                SortedAreaInds(SortedAreaInds==CollectedVoxels(Ind)) = []; %remove indices that have been merged and keep going
            end
            CurrInd = CurrInd+1;
        else
            KeepGoing = 0; %nothing left to start a search from, so let's stop completely, sort results and return them
        end
    end
    
    %% resort Suggestions
    [TmpSort,IndsSuggestionsDescending] = sort(Suggestions.NVoxelTotalMergers,'descend');
    SuggestionsTmp = Suggestions; %backup
    for Ind = 1:length(IndsSuggestionsDescending)
        Suggestions.IndsMergers{Ind}        = SuggestionsTmp.IndsMergers{IndsSuggestionsDescending(Ind)};   %each one gets a suggestion, at least staying itself.
        Suggestions.NVoxelMergers{Ind}      = SuggestionsTmp.NVoxelMergers{IndsSuggestionsDescending(Ind)}; %for each one suggestion, note the size.
        Suggestions.NVoxelTotalMergers(Ind) = SuggestionsTmp.NVoxelTotalMergers(IndsSuggestionsDescending(Ind)); %the number of voxels if merged, such that we can sort later
        Suggestions.NMergers(Ind)           = SuggestionsTmp.NMergers(IndsSuggestionsDescending(Ind));      %keep track of how many areas are merged.
    end
    
    %% done make a short report OR repeat
    NAreasBeforeMerge = length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName);
    NAreasNow         = length(Suggestions.IndsMergers);
    InfoTxt{1} = ['Merging of subparts of "EffectsOfMRI"(Split "',EffectOfMRI_AALsplit.SplitNames{Suggestions.IndSplit},'") using a search up to ',num2str(Suggestions.VoxelDistForSearch),'.NearestNeighbors'];
    if(NAreasBeforeMerge~=NAreasNow)
        InfoTxt{2} = ['Has reduced number of areas in the suggestion from ',num2str(NAreasBeforeMerge),' to ',num2str(NAreasNow),'.'];
        for IndMergers = 1:NAreasNow
            if(Suggestions.NMergers(IndMergers)~=1)
                InfoTxt{length(InfoTxt)+1} = ['Merger ',num2str(IndMergers),': Merging ',num2str(Suggestions.NMergers(IndMergers)),' subparts with ',num2str(Suggestions.NVoxelTotalMergers(IndMergers)),'Voxels in total. Merging indices (',num2str(Suggestions.IndsMergers{IndMergers}'),') with NVoxels (',num2str(Suggestions.NVoxelMergers{IndMergers}'),') each.'];
                AreasStr = [];
                for Ind = 1:length(Suggestions.IndsMergers{IndMergers})
                    if(Ind==1)
                        if(Ind==length(Suggestions.IndsMergers{IndMergers}))
                            AreasStr = [AreasStr,EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{Suggestions.IndsMergers{IndMergers}(Ind)}];
                        else
                            AreasStr = [AreasStr,EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{Suggestions.IndsMergers{IndMergers}(Ind)},'+'];
                        end
                    else
                        if(Ind~=length(Suggestions.IndsMergers{IndMergers}))
                            AreasStr = [AreasStr,EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{Suggestions.IndsMergers{IndMergers}(Ind)},'+'];
                        else
                            AreasStr = [AreasStr,EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{Suggestions.IndsMergers{IndMergers}(Ind)}];
                        end
                    end
                end
                InfoTxt{length(InfoTxt)+1} = ['Areas merged: "',AreasStr,'"'];
            else
                InfoTxt{length(InfoTxt)+1} = ['Merger ',num2str(IndMergers),': NOT A REAL MERGE! I.e. subpart "',EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{Suggestions.IndsMergers{IndMergers}},'" could not be connected to others and therefore stays like it was. (',num2str(Suggestions.NVoxelTotalMergers(IndMergers)),'Voxels in total.)'];
            end
            InfoTxt{length(InfoTxt)+1} = ' ';
            InfoTxt{length(InfoTxt)+1} = '       HINT:   KEEP THIS WINDOW OPEN   :HINT';
        end
        H = helpdlg(InfoTxt,'Merge suggestion report');
        uiwait(H,2); %wait 2 seconds such that user can take a look at it.
        Repeat = 0;
    else
        if(NAreasBeforeMerge<NAreasNow)
            InfoTxt{1} = ['ERROR ',num2str(NAreasBeforeMerge),'==NAreasBeforeMerge<NAreasNow==',num2str(NAreasNow),'!!! Check code, because this should be impossible.'];
            Repeat = 0;
        else
            InfoTxt{2} = '     !!!WAS NOT ABLE TO REDUCE NUMBER OF SUBPARTS!!!';
            InfoTxt{3} = ['DO YOU WANT TO REPEAT WITH A SEARCH EXTENDED TO ONE NEIGHBORHOOD FURTHER, I.E. N = ',num2str(Suggestions.VoxelDistForSearch+1),'?'];
            choice = questdlg(InfoTxt,'Repeat with LARGER search?','Yes','No','Yes');
            switch(choice)
                case 'Yes'
                    SearchDist = SearchDist+1;
                    Repeat = 1;
                    clear Suggestions
                case 'No'
                    Repeat = 0;
            end
        end
    end
end

%% DONE (finally) add infotxt
Suggestions.Info.NAreasBeforeMerge = NAreasBeforeMerge;
Suggestions.Info.NAreasAfterMerge  = NAreasNow;
Suggestions.Info.InfoTxt = InfoTxt;

end