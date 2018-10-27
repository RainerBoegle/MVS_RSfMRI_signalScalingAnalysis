function Suggestions = ManualMergers(EffectOfMRI_AALsplit,IndSplit)
% This function helps the user to make mergers of areas manually as a preparation
% for the automatic algorithm for suggesting mergers for the remaining areas.
% 
% Show possible areas and ask user to select multiple or a single area for
% a merger. If user chooses quit then stop with manual merging, return
% results and go over to setting up the automatic algorithm.
%
%Date: 24.08.2014 (V1)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% setup Suggestions-Struct
Suggestions.IndSplit = IndSplit;

%% get map of atlas
NII_tmp        = nifti(EffectOfMRI_AALsplit.MPaths{IndSplit});
AALsplit_Input = round(NII_tmp.dat(:,:,:));
IndsAreas      = unique(AALsplit_Input);
if(any(IndsAreas==0))
    IndsAreas(IndsAreas==0) = [];
end
Dim            = size(AALsplit_Input);
clear NII_tmp

%% get labels for all areas contained in the mask
AllAreasAvailable = cell(length(IndsAreas),1);
for Ind = 1:length(IndsAreas)
    AllAreasAvailable{Ind} = EffectOfMRI_AALsplit.Labels(Suggestions.IndSplit).Labels{IndsAreas(Ind)};
end

%% start showing the areas to user for selection of mergers
KeepMerging  = 1; %init true
IndMerge     = 1; %start with the first of course
AreaList     = AllAreasAvailable; %init with all areas
IndsAreasList= IndsAreas; %init with all indices corresponding to areas
while(KeepMerging)
    [Selection,ok] = listdlg('ListString',AreaList,'SelectionMode','multiple','Name','Select Mergers','OKString','Merge','CancelString','Quit');
    if(ok) %add a merge
        Suggestions.IndsMergers{IndMerge} = zeros(length(Selection),1);
        for Ind = 1:length(Selection)
            Suggestions.IndsMergers{IndMerge}(Ind) = IndsAreasList(Selection(Ind));
        end
        %update
        for Ind = 1:length(Suggestions.IndsMergers{IndMerge})
            IndsAreasList(IndsAreasList==Suggestions.IndsMergers{IndMerge}(Ind)) = []; %remove indices
        end
        if(~isempty(IndsAreasList))
            AreaList = cell(length(IndsAreasList),1);
            for Ind = 1:length(IndsAreasList)
                AreaList{Ind} = AllAreasAvailable{find(IndsAreas==IndsAreasList(Ind))}; %make list of remaining areas
            end
            IndMerge = IndMerge + 1;
        else
            KeepMerging = 0; %nothing to merge any more so let's stop.
        end
    else
        KeepMerging = 0; %stopped by user
    end
end %done with merging

%% write all information to suggestion-struct
for Ind = 1:length(Suggestions.IndsMergers)
    Indices = Suggestions.IndsMergers{Ind};
    Suggestions.NVoxelMergers{Ind} = zeros(length(Indices),1);
    for IndsIndices = 1:length(Indices)
        Suggestions.NVoxelMergers{Ind}(IndsIndices) = length(find(AALsplit_Input==Indices(IndsIndices)));
    end
    Suggestions.NVoxelTotalMergers(Ind,1) = sum(Suggestions.NVoxelMergers{Ind});
    Suggestions.NMergers(Ind,1)           = length(Indices);
end

end