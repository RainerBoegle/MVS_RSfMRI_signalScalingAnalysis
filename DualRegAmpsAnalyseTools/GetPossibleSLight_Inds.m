function [SLightIndsInMaskCell,SLightIndsInMaskCellSize] = GetPossibleSLight_Inds(NHood,StartIndsVol,AllowedIndsVol,VolDim,RestrictAlsoToStartInds)
% This function takes each (linear-)index from the total volume that is part of the mask 
% and generates the indices RELATIVE TO THE MASK(!) that are within the SearchLight and the mask.
% 
% The parameter NHood indicates how many neighborhoods around the center
% voxel are included in the searchlight (as long as this fits the mask).
%
% NB:
%    This function is called by "GenerateSLight.m".
%    If you want to "go backwards", i.e. look up a coordinate and show 
%    the corresponding searchlight in the brain try using "CheckOutLambdaAtSL.m"
%    to see how it is done.
%
%Usage:
%      [SLightIndsInMaskCell,SLightIndsInMaskCellSize] = GetPossibleSLight_Inds(NHood,StartIndsVol,AllowedIndsVol,VolDim,RestrictAlsoToStartInds);
%      [SLightIndsInMaskCell,SLightIndsInMaskCellSize] = GetPossibleSLight_Inds(NHood,StartIndsVol,AllowedIndsVol,VolDim,1); %crop searchlights to start indices & allowed indices, -could be whole-brain-mask&&mask-significant-regions-some-analysis
%      [SLightIndsInMaskCell,SLightIndsInMaskCellSize] = GetPossibleSLight_Inds(NHood,StartIndsVol,AllowedIndsVol,VolDim,0); %DEFAULT: crop searchlights to stay within AllowedIndsVol, -could be whole-brain mask
%      [SLightIndsInMaskCell,SLightIndsInMaskCellSize] = GetPossibleSLight_Inds(NHood,StartIndsVol,AllowedIndsVol,VolDim);   %DEFAULT: crop searchlights to stay within AllowedIndsVol, -could be whole-brain mask
%
%
%V1.1
%Date: V1.1(23.7.2015) (more comments) V1.0(06.02.2015) (initial implementation based on test script for analysis of scaling data.)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% init
SLightIndsInMaskCell     = cell( length(StartIndsVol),1);
SLightIndsInMaskCellSize = zeros(length(StartIndsVol),1);
NVoxNHood = (2*NHood+1)^3;
try
    if(RestrictAlsoToStartInds)
        disp('Will restrict searchlight using both masks.');
    else
        disp('Will restrict searchlight using the "AllowedIndsVol" mask only.');
    end
catch
    RestrictAlsoToStartInds = 0;
    disp('Will restrict searchlight using the "AllowedIndsVol" mask only.');
end

%% for each start Index generate subscript indices and then get all in the neightborhood and converte back, then take out those that are not on the StartInds list.
H_waitbar = waitbar(0,['Preparing SearchLights...']);
CheckCenterIndexIsFirst = 1; %init as correct and change otherwise the use this to inform user later.
for IndCenterVoxel = 1:length(StartIndsVol)
    CenterCheckIndex = []; %init empty --> may be obsolete later
    CenterCheckIndexRelative2WholeBrainMask = []; %init empty --> may be obsolete later
    [xCenter,yCenter,zCenter]=ind2sub(VolDim,StartIndsVol(IndCenterVoxel)); %find center voxel coord in 3D
    IndsSL = zeros(NVoxNHood,1); %init 
    Ind = 0; %init
    IndsX = max([1,(xCenter-NHood)]):min([VolDim(1),(xCenter+NHood)]); IndsX = [xCenter, IndsX(IndsX~=xCenter)]; %make sure we have the center of the searchlight as the first entry AND only entries within the volume to avoid errors. Later check if inside mask as well.
    IndsY = max([1,(yCenter-NHood)]):min([VolDim(2),(yCenter+NHood)]); IndsY = [yCenter, IndsY(IndsY~=yCenter)]; %make sure we have the center of the searchlight as the first entry AND only entries within the volume to avoid errors. Later check if inside mask as well.
    IndsZ = max([1,(zCenter-NHood)]):min([VolDim(3),(zCenter+NHood)]); IndsZ = [zCenter, IndsZ(IndsZ~=zCenter)]; %make sure we have the center of the searchlight as the first entry AND only entries within the volume to avoid errors. Later check if inside mask as well.
    for Indx = IndsX
        for Indy = IndsY
            for Indz = IndsZ
                Ind = Ind + 1;
                IndsSL(Ind) = sub2ind(VolDim,Indx,Indy,Indz); %find LINEAR-Index in NHood distance
                if(Ind==1)
                    if((Indx==xCenter)&&(Indx==xCenter)&&(Indx==xCenter))
                        CenterCheckIndex = IndsSL(Ind);
                        CenterCheckIndexRelative2WholeBrainMask = find(CenterCheckIndex==AllowedIndsVol);
                    else
                        disp('WARNING Index 1 is not the center of searchlight!');
                        CheckCenterIndexIsFirst = 0;
                    end
                end
                if(RestrictAlsoToStartInds)
                    if((isempty(find(IndsSL(Ind)==StartIndsVol)))||(isempty(find(IndsSL(Ind)==AllowedIndsVol)))) %NOT in the whole brain mask OR analysis mask; -because not on the list! --> mark for removal
                        IndsSL(Ind) = 0;
                    end
                else
                    if(isempty(find(IndsSL(Ind)==AllowedIndsVol))) %NOT in the whole brain mask list! --> mark for removal
                        IndsSL(Ind) = 0;
                    end
                end
            end
        end
    end
    IndsSL(IndsSL==0) = []; %remove those inds that are not on the list!
    IndsSLfinal = zeros(length(IndsSL),1);
    for IndsList = 1:length(IndsSL)
        IndsSLfinal(IndsList) = find(IndsSL(IndsList)==AllowedIndsVol); %assign such that they indicate the indices RELATIVE TO THE WHOLE BRAIN MASK!
    end
    if(~isempty(CenterCheckIndexRelative2WholeBrainMask))
        if(CenterCheckIndexRelative2WholeBrainMask~=IndsSLfinal(1))
            disp('WARNING center is not first entry in searchlight index list!');
            CheckCenterIndexIsFirst = 0;
            if(~any(IndsSLfinal==CenterCheckIndexRelative2WholeBrainMask))
                disp('WARNING!!!! Center of Searchlight has gone missing!!!!!!!!!');
            end
        end
    end
    SLightIndsInMaskCell{IndCenterVoxel}    =IndsSLfinal; %Done, next center.
    SLightIndsInMaskCellSize(IndCenterVoxel)=length(IndsSLfinal)/NVoxNHood; %save the size.
    H_waitbar = waitbar(IndCenterVoxel/length(StartIndsVol),H_waitbar,['Preparing SearchLights ',num2str(IndCenterVoxel*100/length(StartIndsVol)),'%done...']);
    disp(['SL',num2str(IndCenterVoxel,['%0',num2str(ceil(log10(length(StartIndsVol)))),'g']),'of',num2str(length(StartIndsVol)),' contains ',num2str(length(IndsSL),['%0',num2str(ceil(log10(NVoxNHood))),'g']),'of',num2str(NVoxNHood),'Voxels, i.e. ',num2str(length(IndsSL)*100/NVoxNHood,3),'% of a neighborhood order ',num2str(NHood),'.']);
end 
close(H_waitbar);

%% Done
disp(' ');
if(CheckCenterIndexIsFirst)
    disp('Done.');
else
    disp('Done (some strange behavior though; -check messages...).');
end

end


