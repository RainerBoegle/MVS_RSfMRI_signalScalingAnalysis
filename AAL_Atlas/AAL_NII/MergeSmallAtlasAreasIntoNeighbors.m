%Use OutputStruct from ExpandAtlasOverlap
%Display what we have
%Set threshold & suggest 27Voxels
%Take Atlas numbers of areas
%Go over areas and search neighbors different than them and not zero.
%collect neighbor indices for all voxel of that area
%go over neighbor indices and find biggest neighbor area including how many have this one as a neighbor
%assign this neighbor indice and remove this indice from the list of areas
%to examine if it is there.

NB start with smallest, such that it might merge with next biggest of the others if this is the case.

try this for different cutoffs.


%% To Do check if small areas can reach neighbors and if yes check the neighbor number and whoever is larger decides.
%% This will hopefully work from a NearestNeighbar iteration that has preceeded this algorithm to remove areas that are unnecessary.
%% Need to make a list of areas in input and keep track which areas have been incorporated and display this.
%% Then we can let the minimum number & maximum extend grow and check the results.

%% get all voxels belonging to a number that is smaller than the cutoff
%% check which voxels lie at the "surface" by checking which have neighbors that are not of their own type
%% for these voxels check what is in their Nth (best would be nearest) neighborhood that is not of their own type.
%% keep track of what is found --> check sizes of adjacent areas that have been touched and decide according to max(NVoxelsTouchingNeighbor*NumberOfVoxelsNeighbor) for Neighbor Ind==X.
%% area gets same number as that neighbor.
%% next