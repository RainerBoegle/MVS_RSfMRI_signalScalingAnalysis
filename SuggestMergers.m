function Suggestions = SuggestMergers(EffectOfMRI_AALsplit,IndSplit,SearchDist,NAreasMax,NVoxelsMax,SuggestionsInit)
% This function looks through a given split "IndSplit" for a given distance
% in Voxels "SearchDist" and suggest mergers of areas that are in this
% search range distance.
%
%Date: 24.08.2014 (V3)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% get map of atlas
NII_tmp        = nifti(EffectOfMRI_AALsplit.MPaths{IndSplit});
AALsplit_Input = round(NII_tmp.dat(:,:,:));
Dim            = size(AALsplit_Input);
clear NII_tmp

%% work on Initial suggestions, i.e. remove them from AALsplit_Input and
%% save them first before trying to add further.
if(~isempty(SuggestionsInit))
    AreaIndsToCensor = [];
    for Ind = 1:length(SuggestionsInit.IndsMergers)
        Suggestions.IndsMergers{Ind}       = SuggestionsInit.IndsMergers{Ind};
        Suggestions.NVoxelMergers{Ind}     = SuggestionsInit.NVoxelMergers{Ind};
        Suggestions.NVoxelTotalMergers(Ind)= SuggestionsInit.NVoxelTotalMergers(Ind);
        Suggestions.NMergers(Ind)          = SuggestionsInit.NMergers(Ind);
        
        AreaIndsToCensor = [AreaIndsToCensor; SuggestionsInit.IndsMergers{Ind}];
    end
    %update
    for IndRemove = 1:length(AreaIndsToCensor)
        AALsplit_Input(AALsplit_Input==AreaIndsToCensor(IndRemove)) = 0; %remove those that are already in the suggestions(init).
    end
    IndsAreas      = unique(AALsplit_Input); %Indices of areas in the map
    if(any(IndsAreas==0))
        IndsAreas(IndsAreas==0) = [];
    end
    NVoxelAllMasks = zeros(length(IndsAreas),1); %number of voxels per areas/mask so we can sort later
    for Ind = 1:length(IndsAreas)
        NVoxelAllMasks(Ind) = length(find(AALsplit_Input==IndsAreas(Ind)));
    end
    StartInd = length(SuggestionsInit.IndsMergers) + 1; %start index for suggestions of mergers
else
    IndsAreas      = unique(AALsplit_Input); %Indices of areas in the map
    if(any(IndsAreas==0))
        IndsAreas(IndsAreas==0) = [];
    end
    NVoxelAllMasks = zeros(length(IndsAreas),1); %number of voxels per areas/mask so we can sort later
    for Ind = 1:length(IndsAreas)
        NVoxelAllMasks(Ind) = length(find(AALsplit_Input==IndsAreas(Ind)));
    end
    StartInd = 1; %start index for suggestions of mergers
end

%% start algorithm
NAreaTotal = length(IndsAreas);
NVoxelTotal= length(find(AALsplit_Input~=0));

Repeat = 1; %in case we are not successful
while(Repeat)
    %% setup Suggestions-struct
    Suggestions.IndSplit            = IndSplit; %keep track of which one we looked at.
    Suggestions.VoxelDistForSearch  = SearchDist; %how far to look (in voxels) for other areas than the current area
    Suggestions.NAreasMax           = NAreasMax;  %Maximum number of areas that are allowed in a merge
    Suggestions.NVoxelsMax          = NVoxelsMax; %Maximum number of voxels for a merge. NB: this is ignored if NMerge<=2 and (NVoxel<=NVoxelsMax*3 or NVoxel<=NVoxelTotal*0.11)
    %Old: (now will grow such that we always end with least number of mergers)
    % Suggestions.IndsMergers         = cell(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1);  %each one gets a suggestion, at least staying itself.
    % Suggestions.NVoxelMergers       = cell(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1);  %for each one suggestion, note the size.
    % Suggestions.NVoxelTotalMergers  = zeros(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1); %the number of voxels if merged, such that we can sort later
    % Suggestions.NMergers            = zeros(length(EffectOfMRI_AALsplit.Split(IndSplit).AllMasksName),1); %keep track of how many areas are merged.
    
    %% do search
    N_NeighborsToSearch = Suggestions.VoxelDistForSearch; %how far to look for other areas than the current area
    [SortedList,IndsAscending] = sort(NVoxelAllMasks); %'ascend'
    SortedAreaInds = IndsAreas(IndsAscending); %indices of search from smallest area to largest
    KeepGoing = 1; %init true for start
    CurrInd   = StartInd; %start at ind suggested by SuggestionsInit
    while(KeepGoing)
        if(~isempty(SortedAreaInds)) %still something to do
            CurrentAreaInd = SortedAreaInds(1); %start here, if nothing else is found then at least this index is to be noted
            CollectedVoxels     = [CurrentAreaInd]; %init: collect nearest neighbors that are not zero & not the current index
            CollectedVoxelsLast = []; %init empty
            KeepSearching = 1; %init true
            while(KeepSearching) %keep growing until not possible any more or break criteria reached.
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
                    if(length(CollectedVoxelsLast)<=2) %might be okay to keep merging more
                        NVoxel = zeros(length(CollectedVoxelsLast),1);
                        for IndArea = 1:length(CollectedVoxelsLast)
                            NVoxel(IndArea) = length(find(AALsplit_Input==CollectedVoxelsLast(IndArea)));
                        end
                        if((sum(NVoxel)>=2*NVoxelsMax) || (sum(NVoxel)>NVoxelTotal.*0.1)) %stop if this merge is already more than 10% or 3*NVoxelsMax
                            KeepSearching = 0; %let's stop because this will get too large!
                            if((sum(NVoxel)>=3*NVoxelsMax) || (sum(NVoxel)>NVoxelTotal.*0.2)) %don't even merge these.
                                disp(' ');
                                disp('Merger judged too big!');
                                disp(['Will keep "',EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{CollectedVoxels},'" separate.']);
                            else
                                CollectedVoxels = CollectedVoxelsLast; %allow the merging of these two.
                                disp(' ');
                                disp('Merging stopped! Only: ');
                                for Ind = 1:length(CollectedVoxels)
                                    disp(['"',EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{CollectedVoxels(Ind)},'"']);
                                end
                                disp('...merged');
                            end
                        end
                    else
                        if(length(CollectedVoxelsLast)<=NAreasMax) %might be okay to merge these
                            NVoxel = zeros(length(CollectedVoxelsLast),1);
                            for IndArea = 1:length(CollectedVoxelsLast)
                                NVoxel(IndArea) = length(find(AALsplit_Input==CollectedVoxelsLast(IndArea)));
                            end
                            if((sum(NVoxel)>NVoxelsMax) || (sum(NVoxel)>NVoxelTotal.*0.1)) %allow the merging of these but not more.
                                KeepSearching = 0; %let's stop because this will get too large!
                                CollectedVoxels = CollectedVoxelsLast; %allow the merging of these.
                                disp(' ');
                                disp('Merging stopped! Only: ');
                                for Ind = 1:length(CollectedVoxels)
                                    disp(['"',EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{CollectedVoxels(Ind)},'"']);
                                end
                                disp('...merged');
                            end
                        else
                            KeepSearching = 0; %let's stop because this will get too large!
                            NVoxel = zeros(length(CollectedVoxelsLast),1);
                            for IndArea = 1:length(CollectedVoxelsLast)
                                NVoxel(IndArea) = length(find(AALsplit_Input==CollectedVoxelsLast(IndArea)));
                            end
                            if((sum(NVoxel)<=NVoxelsMax) || (sum(NVoxel)<=NVoxelTotal.*0.1)) %allow the merging of these but not more.
                                CollectedVoxels = CollectedVoxelsLast; %allow the merging of these.
                            else
                                NVoxel2 = zeros(length(CollectedVoxels),1); %this is smaller than NVoxel because it is from last iteration
                                for IndArea = 1:length(CollectedVoxels)
                                    NVoxel2(IndArea) = length(find(AALsplit_Input==CollectedVoxels(IndArea)));
                                end
                                if((sum(NVoxel)/sum(NVoxel2))>1.15) %if this is only 15% larger since last time, then we let it pass
                                    disp('I let it slip this time... ;)');
                                    CollectedVoxels = CollectedVoxelsLast; %allow the merging of these.
                                else
                                    %leave the additions of last search off the list
                                    disp(' ');
                                    disp('The areas added in the last search did not only increase the number of mergers over maximum number of areas in merge, "NAreasMax",')
                                    disp('but also has now stepped over maximum number of voxels in the merge "NVoxelsMax".');
                                    disp('Therefore this is concluded as a significant increase and this last addition of areas will be rejected, keeping the interation before that.');
                                end
                            end
                            disp(' ');
                            disp('Merging stopped! Only: ');
                            for Ind = 1:length(CollectedVoxels)
                                disp(['"',EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{CollectedVoxels(Ind)},'"']);
                            end
                            disp('...merged.');
                        end
                    end 
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
                AALsplit_Input(AALsplit_Input==CollectedVoxels(Ind)) = 0;  %also remove indices from mask
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
        end
        InfoTxt{length(InfoTxt)+1} = ' ';
        InfoTxt{length(InfoTxt)+1} = '       HINT:   KEEP THIS WINDOW OPEN   :HINT';
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