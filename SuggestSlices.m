function [SuggestedSlices] = SuggestSlices(NIFTI_files,Mode,varargin)
% This function can suggest slices for display with slover.
% The suggested slices are combined from all input NIFTI_files.
% The mode of suggesting slices can be set in the following ways:
%Mode:
%    1. 'Simple'  --> Just pick all unique slices from lowest to highest slice
%                     in all nifit files that are not empty, i.e. contain values
%                     other than zero (Background). 
%                     If more then MaxNumSlices, then thin out by keeping
%                     only those that are more than minimal resolution apart.
%    2. 'Cluster' --> ASSUME that all input images are Cluster maps or masks.
%                     Pick slices that cover each cluster in each image as 
%                     best as possible, i.e. top & bottom and up to three slices
%                     with most voxels in them (sorted descending).
%                     NB: if only one cluster then we assume that is a [0,1] mask
%                     and take N slices, 
%                     top M and just above median + median + just below median + two below median if possible
%                               just above 3rdQrt + 3rdQrt + just below 3rdQrt if possible
%                               just above 1stQrt + 1stQrt + just below 1stQrt if possible 
%                     if possible.
%                     Then thin out till MaxNumSlices is reached. 
%    3. 'Stats'   --> ASSUME that all input images are Statistics maps.
%                     For each image find the maximum statistics value and remove
%                     the two neighboring slices and repeat with the remaining slices.
%                     Do the same starting with the minimum and ascending.
%                     Thin out each list for all images till MaxNumSlices is reached.
%
%V2.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V2.1: (13.02.2015): 'Cluster'-mode: if only one value per input, i.e. [0,1] mask --> pick more slices. V2.0(28.02.2015): initial implementation (changed version from original subfunction of various display tools functions)

%% Defaults
MaxNumSlices = 42;

%% extra inputs?
if(nargin>2)
    if(mod(length(varargin),2)~=0)
        error('Extra inputs have to be in pairs, i.e. in the form SuggestSlices(NIFTI_files,Mode,''Command-String'',Command-Value,...)');
    end
    CommStrings = varargin(1:2:end);
    CommVals    = varargin(2:2:end);
    for Ind = 1:length(CommStrings)
        switch(CommStrings{Ind})
            case 'Dir'
                if(~strcmpi(CommVals{Ind},'z'))
                    error('Sorry not supported yet!');
                end
            case 'MaxNumSlices'
                if(isnumeric(CommVals{Ind}))
                    MaxNumSlices = CommVals{Ind};
                else
                    error('MaxNumSlices Command-Value must be numeric!');
                end
            otherwise
                error(['Command-String "',CommStrings{Ind},'" unknown!']);
        end
    end
end

%% suggest slices dependent on mode
switch(Mode)
    case {'Simple','simple','Simply','simply'}
        try
            SuggestedSlices = getAllPossibleSlices(NIFTI_files,1); %get all & thin out neighbors by MinResZ
        catch CATCHsimple
            disp_catch(CATCHsimple,[mfilename,'>getAllPossibleSlices>Simple'],'CatchObjSimple');
            SuggestedSlices = [];
            return;
        end
        if(length(SuggestedSlices)>MaxNumSlices) %thin out because too many
            SuggestedSlices = ThinOutSuggestedSlices(SuggestedSlices,MaxNumSlices);
        end
    case {'Cluster','cluster','Clusters','clusters'}
        try
            [SuggestedSlices,Res_z] = SuggestSlicesForClusterImages(NIFTI_files);
            %% use min(Res_z) to sort out those slices that are too close to neighbors
            if(length(SuggestedSlices)>MaxNumSlices) %thin out because too many
                MinResZ  = min(Res_z);
                [RemoveIndices,SortedCoords_z] = RemoveTooCloseNeighbors(SuggestedSlices,MinResZ);
                if(~isempty(RemoveIndices))
                    SortedCoords_z(RemoveIndices) = []; %remove the ones marked for removal
                end
                SuggestedSlices = unique(SortedCoords_z);
            end
        catch CATCHcluster
            disp_catch(CATCHcluster,[mfilename,'>SuggestSlicesForClusterImages'],'CatchObjCluster');
            SuggestedSlices = [];
            return;
        end
        if(length(SuggestedSlices)>MaxNumSlices) %thin out because too many
            SuggestedSlices = ThinOutSuggestedSlices(SuggestedSlices,MaxNumSlices);
        end
    case {'Stats','stats','Statistics','statistics'}
        try
            [SuggestedSlices_max,Res_z] = SuggestSlicesForStatsImages(NIFTI_files,'max');
             SuggestedSlices_min        = SuggestSlicesForStatsImages(NIFTI_files,'min');
            SuggestedSlices = unique([SuggestedSlices_max; SuggestedSlices_min]);
            %% use min(Res_z) to sort out those slices that are too close to neighbors
            if(length(SuggestedSlices)>MaxNumSlices) %thin out because too many
                MinResZ  = min(Res_z);
                [RemoveIndices,SortedCoords_z] = RemoveTooCloseNeighbors(SuggestedSlices,MinResZ);
                if(~isempty(RemoveIndices))
                    SortedCoords_z(RemoveIndices) = []; %remove the ones marked for removal
                end
                SuggestedSlices = unique(SortedCoords_z);
            end
        catch CATCHstats
            disp_catch(CATCHstats,[mfilename,'>SuggestSlicesForStatsImages'],'CatchObjStats');
            SuggestedSlices = [];
            return;
        end
        if(length(SuggestedSlices)>MaxNumSlices) %thin out because too many
            SuggestedSlices = ThinOutSuggestedSlices(SuggestedSlices,MaxNumSlices);
        end
    otherwise
        error(['Mode-string "',Mode,'" unknown!']);
end


end


%% subfunctions
%% SuggestSlicesForStatsImages
function [SuggestedSlices,Res_z] = SuggestSlicesForStatsImages(NIFTI_files,StatsMode)
% ASSUME that all input images are Statistics maps.
% For each image find the maximum statistics value and remove
% the two neighboring slices and repeat with the remaining slices.
% Do the same starting with the minimum and ascending.
% Thin out each list for all images till MaxNumSlices is reached.

SuggestedSlices = [];
Res_z           = [];
for IndData = 1:length(NIFTI_files)
    %% get image data
    NII_tmp = nifti(NIFTI_files{IndData});
    Data    = NII_tmp.dat(:,:,:);
    v2m     = NII_tmp.mat;
    Res_z   = [Res_z; abs(v2m(3,3))];
    
    while(~isempty(find(Data(:)~=0)))
        %% get coordinates with non-zero values
        vox = getNonZeroVox(Data); %get voxel coordinates of non-zero voxels, i.e. where atlas or maps are.
        if(~isempty(vox))
            Coords_z = zeros(size(vox)); %init temporary coords
            for i=1:size(vox,1)
                Coords_z(i,1:3)=vox(i,:)*v2m(1:3,1:3) + v2m(1:3,4)'; %z-coordinates values
            end
            [Coords_z, UniqueInds] = unique(Coords_z(:,3));
            vox = vox(UniqueInds,:);
        else
            break; %for safety make sure that we break out if no voxels are found
        end
        
        %% get data check mode and find max/min stats values and add slice to list then remove surounding slices ie set to zero
        if(strcmp(StatsMode,'max'))
            %% get data per slice
            DataPerSlice = zeros(size(vox,1),1);
            for IndSlice = 1:size(vox,1)
                DataPerSlice(IndSlice) = max(max(Data(:,:,vox(IndSlice,3))));
            end
            %% find maximum stats value and add slice to list
            [MaxVal,MaxInds]= max(DataPerSlice);
            SuggestedSlices = [SuggestedSlices; Coords_z(MaxInds(1))];
            
            %% remove data around this max slice
            RemoveInds = [vox(MaxInds(1),3); min([size(Data,3),vox(MaxInds(1),3)+1]); max([1,vox(MaxInds(1),3)-1])];
            Data(:,:,RemoveInds) = 0;
        else
            %% get data per slice
            DataPerSlice = zeros(size(vox,1),1);
            for IndSlice = 1:size(vox,1)
                DataPerSlice(IndSlice) = min(min(Data(:,:,vox(IndSlice,3))));
            end
            %% find maximum stats value and add slice to list
            [MinVal,MinInds]= min(DataPerSlice);
            SuggestedSlices = [SuggestedSlices; Coords_z(MinInds(1))];
            
            %% remove data around this max slice
            RemoveInds = [vox(MinInds(1),3); min([size(Data,3),vox(MinInds(1),3)+1]); max([1,vox(MinInds(1),3)-1])];
            Data(:,:,RemoveInds) = 0;
        end
    end
end

end

%% SuggestSlicesForClusterImages
function [SuggestedSlices,Res_z] = SuggestSlicesForClusterImages(NIFTI_files)
% ASSUME that all input images are Cluster maps or masks.
% Pick cluster numbers and number of voxels per slice.
% Pick slices that cover each cluster in each image as 
% best as possible, i.e. top & bottom and up to three slices
% with most voxels in them (sorted descending).


SuggestedSlices = [];
Res_z           = [];
for IndData = 1:length(NIFTI_files)
    %% get image data
    NII_tmp = nifti(NIFTI_files{IndData});
    Data    = NII_tmp.dat(:,:,:);
    v2m     = NII_tmp.mat;
    Res_z   = [Res_z; abs(v2m(3,3))];
    
    %% get cluster numbers 
    ClusterInds= unique(Data(Data~=0)); %remove 0, i.e. background
    
    %% per z-slice determine number of slices per Cluster Ind
    if(~isempty(ClusterInds))
        for IndCl = 1:length(ClusterInds)
            [Coords_z,Res_z_tmp,vox] = getCoords(Data,v2m,ClusterInds(IndCl));
            clear Res_z_tmp
            NVoxSlices = zeros(size(vox,1),1);
            for IndSlice = 1:length(NVoxSlices)
                NVoxSlices(IndSlice) = length(find(Data(:,:,vox(IndSlice,3))==ClusterInds(IndCl)));
            end
            
            %% sort according to Nvoxels  
            [SortedNVoxSlices,SortInds] = sort(NVoxSlices,'descend');
            clear SortedNVoxSlices
            
            %% take top three per cluster ind & add top and bottom slice
            %% NB: if there is only one cluster in total then this is probably a [0,1] mask,
            %%     therefore take top four indices + just above median + median + just below median & two below median  & add top and bottom slice
            if(length(ClusterInds)~=1)
                if(length(SortInds)>=3)
                    TopN = SortInds(1:3);
                else
                    TopN = SortInds(:);
                end
            else
                %only one cluster in total --> probably [0,1] mask --> take N if possible
                %top M and just above median + median + just below median + two below median if possible
                %          just above 3rdQrt + 3rdQrt + just below 3rdQrt if possible
                %          just above 1stQrt + 1stQrt + just below 1stQrt if possible 
                %Note that redundances will be removed later, so this just maximized visibility.
                if(length(SortInds)>=10)
                    %       TopM eight   ;                      median+1        ;                median              ;                median-1              ;                median-2              ;                       3rdQrt+1         ;                       3rdQrt         ;                       3rdQrt-1         ;                       1stQrt+1         ;                       1stQrt         ;                       1stQrt-1          ;         
                    TopN = [SortInds(1:8); SortInds(round(length(SortInds)/2)+1); SortInds(round(length(SortInds)/2)); SortInds(round(length(SortInds)/2)-1); SortInds(round(length(SortInds)/2)-2); SortInds(round(length(SortInds)*3/4)+1); SortInds(round(length(SortInds)*3/4)); SortInds(round(length(SortInds)*3/4)-1); SortInds(round(length(SortInds)*1/4)+1); SortInds(round(length(SortInds)*1/4)); SortInds(round(length(SortInds)*1/4)-1)];
                elseif(length(SortInds)>=7)
                    %       TopM four    ;                      median+1        ;                median              ;                median-1              ;                median-2              ;                       3rdQrt+1         ;                       3rdQrt         ;                       3rdQrt-1         ;                       1stQrt+1         ;                       1stQrt         ;                       1stQrt-1          ;
                    TopN = [SortInds(1:4); SortInds(round(length(SortInds)/2)+1); SortInds(round(length(SortInds)/2)); SortInds(round(length(SortInds)/2)-1); SortInds(round(length(SortInds)/2)-2); SortInds(round(length(SortInds)*3/4)+1); SortInds(round(length(SortInds)*3/4)); SortInds(round(length(SortInds)*3/4)-1); SortInds(round(length(SortInds)*1/4)+1); SortInds(round(length(SortInds)*1/4)); SortInds(round(length(SortInds)*1/4)-1)];
                elseif(length(SortInds)>4)
                    %       TopM four    ;                      median+1        ;                median              ;                median-1              ;                       3rdQrt         ;                       3rdQrt-1         ;                       1stQrt+1         ;                       1stQrt         
                    TopN = [SortInds(1:4); SortInds(round(length(SortInds)/2)+1); SortInds(round(length(SortInds)/2)); SortInds(round(length(SortInds)/2)-1); SortInds(round(length(SortInds)*3/4)); SortInds(round(length(SortInds)*3/4)-1); SortInds(round(length(SortInds)*1/4)+1); SortInds(round(length(SortInds)*1/4))];
                else
                    TopN = SortInds(:);
                end
            end
            TopN = unique(TopN); %remove redundant
            if(length(TopN)>4)
                SuggestedSlices = [SuggestedSlices; Coords_z(TopN)]; 
            else
                SuggestedSlices = [SuggestedSlices; Coords_z(TopN); Coords_z(1); Coords_z(end)]; 
            end
        end
    else
        continue;
    end
end

%% get unique ones
SuggestedSlices = unique(SuggestedSlices);

end

%% getAllPossibleSlices
function [SuggestedSlices,MinResZ,Coords_z,Res_z] = getAllPossibleSlices(NIFTI_files,CheckThinOut)
% This function finds ALL possible Slices from the Lowest z-direction slice
% & highest z-direction slice and the slice step-width in all input images.
%
% SuggestedSlices = [1stImage  Lowest-z-direction-slice:Step-width-z-direction-slice:Highest-z-direction-slice,
%                    ...,
%                    LastImage Lowest-z-direction-slice:Step-width-z-direction-slice:Highest-z-direction-slice];
%

Coords_z = [];
Res_z    = [];
for Ind = 1:length(NIFTI_files)
    NII_tmp = nifti(NIFTI_files{Ind});
    Data    = NII_tmp.dat(:,:,:);
    v2m     = NII_tmp.mat;
    
    [Coords_z_tmp,Res_z_tmp] = getCoords(Data,v2m);
    Coords_z = [Coords_z; Coords_z_tmp];
    Res_z    = [Res_z; Res_z_tmp];
end
if(isempty(Coords_z))
    error('Could not extract slices containing data. INPUT VOLUMES SEEM EMPTY!');
end

%% use min(Res_z) to sort out those slices that are too close to neighbors
MinResZ  = min(Res_z);
if(CheckThinOut)
    [RemoveIndices,SortedCoords_z] = RemoveTooCloseNeighbors(Coords_z,MinResZ);
    if(~isempty(RemoveIndices))
        SortedCoords_z(RemoveIndices) = []; %remove the ones marked for removal
    end
    SuggestedSlices = unique(SortedCoords_z);
else
    SuggestedSlices = sort(Coords_z);
end

end

%% getCoords
function [Coords_z,Res_z,vox] = getCoords(Data,v2m,varargin)
% get nonZero coords
if(nargin==2)
    vox = getNonZeroVox(Data); %get voxel coordinates of non-zero voxels, i.e. where atlas or maps are.
else
    if(nargin==3)
        TestValue = varargin{1};
        vox = getNonZeroVox(Data,TestValue); %get voxel coordinates of non-zero voxels, i.e. where atlas or maps are.
    else
        error('getCoords>Inputs to getNonZeroVox have to be 2 at maximum!');
    end
end
if(~isempty(vox))
    Coords_z = zeros(size(vox)); %init temporary coords
    for i=1:size(vox,1)
        Coords_z(i,1:3)=vox(i,:)*v2m(1:3,1:3) + v2m(1:3,4)'; %z-coordinates values
    end
    Res_z = abs(v2m(3,3)); %avoid any odd errors by defining dimensions positive. This usually shouldn't happen because only SPM uses a negative dimension for X-direction.
    [Coords_z, UniqueInds] = unique(Coords_z(:,3));
    vox = vox(UniqueInds,:);
else
    Coords_z = [];
    Res_z    = [];
end

end

%% getNonZeroVox
function [vox] = getNonZeroVox(Data,varargin)
if(nargin==1)
    TestValue = [];
else
    if(nargin==2)
        TestValue = varargin{1};
    else
        error('Inputs to getNonZeroVox have to be 2 at maximum!');
    end
end
if(isempty(TestValue))
    [X,Y,Z] = ind2sub(size(Data),find(Data~=0));
else
    [X,Y,Z] = ind2sub(size(Data),find(Data==TestValue));
end
vox = zeros(length(Z),3);
vox(:,1) = X;
vox(:,2) = Y;
vox(:,3) = Z;

end

%% RemoveTooCloseNeighbors
function [RemoveIndices,SortedCoords_z] = RemoveTooCloseNeighbors(Coords_z,MinResZ)
SortedCoords_z = sort(Coords_z);
RemoveIndices = [];
for Ind = 2:length(SortedCoords_z)
    if(sqrt((SortedCoords_z(Ind-1)-SortedCoords_z(Ind)).^2)<MinResZ) %potentially a slice that should be marked for removal (otherwise continue with next one)
        if(~isempty(find(RemoveIndices==Ind-1))) %but if the one just one step before it is listed for removal, then we also need to check the distance to the one two steps before it.
            if((Ind-2)>0) %need to make sure the index exists at all
                if(sqrt((SortedCoords_z(Ind-2)-SortedCoords_z(Ind)).^2)<MinResZ) %remove it! (otherwise keep it, i.e. don't mark for removal)
                    RemoveIndices = [RemoveIndices, Ind];
                end
            end
        else
            %remove it because the one before has not been marked for
            %removal, i.e. this one should be good enough for display,
            %given that it is less than MinResZ away.
            RemoveIndices = [RemoveIndices, Ind];
        end
    end
end

end

%% ThinOutSuggestedSlices
function SuggestedSlices = ThinOutSuggestedSlices(SuggestedSlices,MaxNumSlices)
% Thin out slices by trying several distances for removing close neighbors
for TryOutDist = 2:5
    if(length(SuggestedSlices)>MaxNumSlices) %too many let's try differently
        RemoveIndices = RemoveTooCloseNeighbors(SuggestedSlices,TryOutDist); %try with TryOutDist-mm
        if(~isempty(RemoveIndices))
            SuggestedSlices(RemoveIndices) = []; %remove the ones marked for removal
        end
    else
        break; %stop, we are done.
    end
end
    
end