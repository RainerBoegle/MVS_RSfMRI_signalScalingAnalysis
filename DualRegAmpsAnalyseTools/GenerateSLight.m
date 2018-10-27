function SLight = GenerateSLight(V_WholeBrainMask,V_SLmask,NHood,RestrictAlsoToStartInds,varargin)
% This function will generate a SearchLight collection of indices given a
% whole brain template mask and possibly a more specific analysis mask volume
% for defining the size and allowed regions.
%
% The whole brain mask can be taken from the MELODIC analysis, which will be
% very inclusive probably extending beyond the brain.
% The SLmask can be a mask that includes the brain more tightly (brain/gray matter/whatever).
% 
% The idea is that the center voxels for each searchlight come from the SLmask,
% but the searchlight can go a bit outside this mask, i.e. into the whole brain mask.
%
% NB: if RestrictAlsoToStartInds==1 (true) then the voxels of the searchlight are only from the SLmask,
%     i.e. in the case as before, they are more tightly around the brain/gray matter/whatever
%
%
%Usage:
%       SLight = GenerateSLight(V_WholeBrainMask,V_SLmask,NHood,RestrictAlsoToStartInds,varargin);
%       SLight = GenerateSLight(V_WholeBrainMask,V_SLmask,NHood,RestrictAlsoToStartInds,OutputDir); %automatically save to output directory.
%       SLight = GenerateSLight(V_WholeBrainMask,V_SLmask,NHood,RestrictAlsoToStartInds,OutputMatPath); %EXTENSION SENSITIVE --> has to be *.mat!!! %automatically save to output *.mat-file specified in OutputMatPath.
%
%Inputs:
%       V_WholeBrainMask        <--   spm-vol struct or path to nifti
%       V_SLmask                <--   spm-vol struct or path to nifti
%       NHood                   <--   Order of the neighborhood to be included for searchlight. I.e. NHood=1 --> 27Voxels (26 nearest neighbors around the center voxel)
%       RestrictAlsoToStartInds <--   if set to "1" (true) searchlights will be restricted to SLmask, i.e. can't expand into WholeBrainMask.
%       varargin                <--   either a directory to save searchlight definition in with generic name based on input mask and NHood, 
%                                     OR give output path with 'Name.mat' to save to specific location with your chosen name. NB: must be *.mat.
%
%V1.1
%Date: V1.1(23.07.2015): (more comments+change disp_catch) V1.0(08.02.2015): (initial implementation based on test script for analysis of scaling data.)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)


%% check inputs
try %V_WholeBrainMask
   if(~isstruct(V_WholeBrainMask))
       if(ischar(V_WholeBrainMask))%probably the path to the mask
           SLight.V_WholeBrainMask = spm_vol(V_WholeBrainMask);
           SLight.WholeBrainRaw    = SLight.V_WholeBrainMask.private.dat(:)~=0; %mask
       else
           if(isempty(V_WholeBrainMask)) %empty input --> user select
               SLight.V_WholeBrainMask = spm_vol(spm_select(1,'image','Select WHOLE-BRAIN MASK volume (nifti-file)...'));
               SLight.WholeBrainRaw    = SLight.V_WholeBrainMask.private.dat(:)~=0; %mask
           else
               error('V_WholeBrainMask input is not a volume-struct or path to a nifti-file.');
           end
       end
   else
       if(~isfield(V_WholeBrainMask,'dim')||~isfield(V_WholeBrainMask,'private'))
           if(~isfield(V_WholeBrainMask,'dim'))
               if(~isfield(V_WholeBrainMask,'private'))
                   error('V_WholeBrainMask is missing the important fields "dim" & "private"!!!');
               else
                   error('V_WholeBrainMask is missing the important field "dim"!');
               end
           else
               error('V_WholeBrainMask is missing the important field "private"!');
           end
       else
           SLight.V_WholeBrainMask = V_WholeBrainMask;
           SLight.WholeBrainRaw    = SLight.V_WholeBrainMask.private.dat(:)~=0; %mask
       end
   end
catch CATCH_inputs_V_WholeBrainMask
    disp_catch(CATCH_inputs_V_WholeBrainMask,[mfilename,'>>Inputs::V_WholeBrainMask'],'CATCH_inputs_V_WholeBrainMask');
    SLight.V_WholeBrainMask= spm_vol(spm_select(1,'image','Select WHOLE-BRAIN MASK volume (nifti-file)...'));
    SLight.WholeBrainRaw   = SLight.V_WholeBrainMask.private.dat(:)~=0; %mask
end
try %V_SLmask
   if(~isstruct(V_SLmask))
       if(ischar(V_SLmask))%probably the path to the mask
           SLight.V_SLmask = spm_vol(V_SLmask);
           SLight.SLmaskRaw= SLight.V_SLmask.private.dat(:)~=0; %mask
       else
           if(isempty(V_SLmask)) %empty input --> user select
               SLight.V_SLmask = spm_vol(spm_select(1,'image','Select mask volume (nifti-file)...'));
               SLight.SLmaskRaw= SLight.V_SLmask.private.dat(:)~=0; %mask
           else
               error('V_SLmask input is not a volume-struct or path to a nifti-file.');
           end
       end
   else
       if(~isfield(V_SLmask,'dim')||~isfield(V_SLmask,'private'))
           if(~isfield(V_SLmask,'dim'))
               if(~isfield(V_SLmask,'private'))
                   error('V_SLmask is missing the important fields "dim" & "private"!!!');
               else
                   error('V_SLmask is missing the important field "dim"!');
               end
           else
               error('V_SLmask is missing the important field "private"!');
           end
       else
           SLight.V_SLmask = V_SLmask;
           SLight.SLmaskRaw= V_SLmask.private.dat(:)~=0; %mask
       end
   end
catch CATCH_inputs_V_SLmask
    disp_catch(CATCH_inputs_V_SLmask,[mfilename,'>>Inputs::V_SLmask'],'CATCH_inputs_V_SLmask');
    SLight.V_SLmask = spm_vol(spm_select(1,'image','Select mask volume (nifti-file)...'));
    SLight.SLmaskRaw= SLight.V_SLmask.private.dat(:)~=0; %mask
end
%% check dims & orientation
if(sum((SLight.V_SLmask.dim(:)-SLight.V_WholeBrainMask.dim(:)).^2)~=0)
    error('Dimensions are not the same for whole brain mask and analysis mask!');
end
if(sum((SLight.V_WholeBrainMask.mat(:)-SLight.V_SLmask.mat(:)).^2)~=0)
    error('Voxel2World matrices are not the same for whole brain mask and analysis mask!');
end

%% NHood
try %NHood
    if(isempty(NHood)) %Empty input --> user select
        answer_SLight = inputdlg({'SearchLight extend NHood= '},'SLight def',1,{'2'});
        SLight.NHood  = eval(answer_SLight{1}); clear answer_SLight
    else
        if(NHood<0)
            error('NHood must be zero or positive integer!');
        else
            if(rem(NHood,1)~=0)
                error('NHood must be an integer!');
            else
                SLight.NHood = NHood;
            end
        end
    end
catch CATCH_inputs_NHood
    disp_catch(CATCH_inputs_NHood,[mfilename,'>>Inputs::NHood'],'CATCH_inputs_NHood');
    answer_SLight = inputdlg({'SearchLight Extend NHood= '},'SLight def',1,{'2'});
    SLight.NHood = eval(answer_SLight{1}); clear answer_SLight
end
SLight.NVoxNHood = (2*SLight.NHood+1)^3;

%% RestrictAlsoToStartInds
try
    if(~isempty(RestrictAlsoToStartInds))
        if(RestrictAlsoToStartInds)
            disp('Will restrict searchlight using both masks.');
        else
            disp('Will restrict searchlight using the "AllowedIndsVol" mask only.');
        end
    else
        if(strcmp('Yes',questdlg({'Restrict the searchlights also to the start indices as well as the whole brain mask?'; 'The default is only to restrict them to the whole brain but not additionally to the start indices/analysis mask.'},'How to crop searchlights?','Yes','No(Default)','No(Default)')))
            RestrictAlsoToStartInds = 1;
        else
            RestrictAlsoToStartInds = 0;
        end
    end
catch
    RestrictAlsoToStartInds = 0;
    disp('Will restrict searchlight using the "AllowedIndsVol" mask only.');
end
SLight.RestrictAlsoToStartInds = RestrictAlsoToStartInds;

%% OutputDir
if(nargin==5)
    OutputDir = varargin{1};
    if(isempty(OutputDir)) %empty --> user select
        OutputDir = spm_select(1,'dir','Select output directoy for saving...');
        OutputDir_IsDir = 1; %OutputDir is really a directory.
    else
        if(ischar(OutputDir))
            if(exist(OutputDir)~=7) %not a directory
                [basedir,fname,ext] = fileparts(OutputDir);
                if(isempty(ext)) %probably intended to be the OutputDir directory?
                    if(~exist(OutputDir))
                        mkdir(OutputDir); clear basedir fname ext
                        OutputDir_IsDir = 1; %OutputDir is really a directory.
                    end
                else %has an extension --> is a filepath
                    if(~strcmpi(ext,'.mat'))
                        error('If you want to save the searchlight definition into a specific file, then you need to use a *.mat-file!');
                    else
                        if(~exist(basedir))
                            mkdir(basedir); clear basedir fname ext
                        end
                        OutputDir_IsDir = 0; %OutputDir is NOT really a directory.
                    end
                end
            else
                OutputDir_IsDir = 1; %OutputDir is really a directory.
            end
        else
            error('OutputDir must be a string (path to save *.mat-file or directory for standard savefile).');
        end
    end
else
    if(nargin<=4)
        OutputDir = []; %don't save
    else
        error('wrong number of inputs!');
    end
end
    
   
%% get indices that are in mask & generate searchlights
SLight.StartIndsVol  = find(SLight.SLmaskRaw~=0);
SLight.AllowedIndsVol= find(SLight.WholeBrainRaw~=0);
if(SLight.NHood~=0)
    [SLight.SLightIndsInMaskCell,SLight.SLightIndsInMaskCellSize] = GetPossibleSLight_Inds(SLight.NHood,SLight.StartIndsVol,SLight.AllowedIndsVol,SLight.V_SLmask.dim,SLight.RestrictAlsoToStartInds);
else
    disp('NHood==0 --> NO SEARCHLIGHT, but single voxels!!!');
    SLight.SLightIndsInMaskCell = cell(length(SLight.StartIndsVol),1);
    for Ind = 1:length(SLight.SLightIndsInMaskCell)
        SLight.SLightIndsInMaskCell{Ind} = find(SLight.StartIndsVol(Ind)==SLight.AllowedIndsVol); %index relative to mask!
    end
    SLight.SLightIndsInMaskCellSize= ones(length(SLight.StartIndsVol),1);
end

%% save
if(~isempty(OutputDir))
    if(OutputDir_IsDir)
        [basedir,fnameMask]   = fileparts(SLight.V_SLmask.fname);
        [basedir,fnameWBMask] = fileparts(SLight.V_WholeBrainMask.fname);
        if(SLight.RestrictAlsoToStartInds)
            save([OutputDir,filesep,'SLdef_Crop_',regexprep(fnameMask,' ','_'),'_W_', regexprep(fnameWBMask,' ','_'),'_NHood',num2str(SLight.NHood),'.mat'],'SLight');
        else
            save([OutputDir,filesep,'SLdef_',     regexprep(fnameMask,' ','_'),'_In_',regexprep(fnameWBMask,' ','_'),'_NHood',num2str(SLight.NHood),'.mat'],'SLight');
        end
    else
        save(OutputDir,'SLight'); %OutputDir is a path for saving!
    end
end

end

%% NB: disp_catch is now a separate function
% %% disp_catch
% function [] = disp_catch(CATCHobj,varargin)
% if(nargin==2)
%     disp(['Error occurred in function "',mfilename,'>',varargin{1},'"...']);
% else
%     disp(['Error occurred in function "',mfilename,'"...']);
% end
% disp([CATCHobj.identifier,': ',CATCHobj.message]);
% 
% end