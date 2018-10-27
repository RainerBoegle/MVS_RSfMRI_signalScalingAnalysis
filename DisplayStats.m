function [varargout] = DisplayStats(varargin)
%DisplayStats.m
%
% This Script can be used to produce Slice Overlay plots from NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING STATISTICS IMAGES.
%
% The SPM toolbox "slover" originally created by Matthew Brett (http://imaging.mrc-cbu.cam.ac.uk/imaging/MatthewBrett)
% is used to produce the Overlays.
%
% The main purpose of this Script is to make the use of "slover" easier.
%
% It is expected that the user inputs one or more masks, i.e. volumes with
% integer valued only, that are then displayed in as many colors as there
% are inputs or unique values of the input, in case there is only one input.
%
% All parameters are fixed (more or less).
%
% USAGE:
%       H=DisplayStats(MaskFiles); %MaskFiles can be empty, not input then spm_select is used. 
%
%
%
%
%NB: A future extension will be made with the function AutoDisplayStats.m
%    That function will allow settings to be given via a struct, such that this function 
%    can be called in a loop for plotting a series of images automatically.
%
%V1.5
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.5: (16.02.2015): Additional outputs and inputs for more control. V1.0: (14.12.2014): initial implementation (changed version from DisplayMasksNIFTIsOnSlices.m)

%% my quick fix 28.07.2015
PlotSpecial = 0; %0==normal behavior
if(PlotSpecial)
    Threshold_pos =  2.75; %4.9;
    Threshold_neg = -1.4; %-4.9;
end

%% select image or images for display
if(nargin<1)
    NIFTI_files = spm_select([1 2],'image','Select NIFTI-file(s) for creation of Overlay ...');
    NIFTI_files = cellstr(NIFTI_files);
    WaitTime = 2; %s DEFAULT;
else
    if(nargin<=2)
        NIFTI_files = varargin{1};
        if(~iscell(NIFTI_files))
            NIFTI_files = cellstr(NIFTI_files);
        end
        if(nargin==2)
            WaitTime = varargin{2};
        elseif(nargin>2)
            error('Too many input arguments');
        elseif(nargin==1)
            WaitTime = 2; %s DEFAULT;
        end
    end
end

%% Slices selection & orientation & transparency
if(PlotSpecial)
    SliceIndices = [-40; -36; -32; -24; -20; -16;  -8;  -4;  -2;   0;   2;   4;  12;  16;  20;  24;  36;  44];
else
    SliceIndices = SuggestSlices(NIFTI_files,'Stats'); %'Simple'); %
end
SliceIndicesStr = '[';
for IndSlice = 1:length(SliceIndices)
    if(IndSlice==1)
        SliceIndicesStr = [SliceIndicesStr,num2str(SliceIndices(IndSlice))];
    else
        SliceIndicesStr = [SliceIndicesStr,',',num2str(SliceIndices(IndSlice))];
    end
end
SliceIndicesStr = [SliceIndicesStr,']'];
% check if too many
if(length(SliceIndices)<=42)
    def_slices   = []; %init
else
    def_slices   = SliceIndicesStr;
end


Orientation = 'axial';
if(~isempty(def_slices))
    %Dialogbox title and Text
    dlg_title = 'Input for Slice display';
    prompt = {'Enter Slice Numbers [Start_z:Step:Stop_z] in mm:'};
    num_lines = 1;
    
    %Defaults
    def = {def_slices}; %{'-32:2:72','axial'};
    
    %Options for Dialogbox (make it resizable)
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    
    Slice_check = 0; %check var for Slice parameters.
    while(~Slice_check)
        %open dialogbox:
        answer = inputdlg(prompt,dlg_title,num_lines,def,options);
        
        %write out the dialogbox answers to variables
        SliceIndicesStr     = answer{1};
        SliceIndices        = round(eval(SliceIndicesStr));
        
        %% check SliceIndices
        if(~isvector(SliceIndices))
            disp('Slice Indices not a vector! Enter a vector of integers please.');
            Slice_check = 0;
        else
            Slice_check = 1;
        end
        %% check Orientation
        switch lower(Orientation)
            case {'axial','coronal','sagittal'}
                Slice_check = 1;
                %disp(['Orientation: ',Orientation]);
            otherwise
                disp('Orientation is neither "axial", "coronal" or "sagittal"! Please enter a valid orientation.');
                Slice_check = 0;
        end
    end
    %clear variables for dialog boxes to be reused without pain. ;)
    clear Slice_check dlg_title prompt num_lines def options answer
end
disp(['Will produce overlay using ',Orientation,' slices ',SliceIndicesStr,'.']);

%% standard SPM Structural or select image?
structural_img = [fileparts(which('spm.m')),filesep,'canonical',filesep,'single_subj_T1.nii'];       

%% check if Structural image can be reached
if(~exist(structural_img))
    disp(['Structural image (',structural_img,') not found! Check paths.']);
    pause(2);
    if((evalin('base','exist(''prevsect'')')==1))
        structural_img = evalin('base','prevsect');
    else
        structural_img = spm_select(1,'image','Select structural image for background...');
    end
    if(strcmp(structural_img((length(structural_img)-1):end),',1'))
        structural_img = structural_img(1:(length(structural_img)-2));
    end
end
disp(['Using "',structural_img,'" as structural background img']);

%% check inputs and select colors
NInputs= length(NIFTI_files);
Colors = cell(length(NIFTI_files),1);
Ranges = cell(length(NIFTI_files),1);
for IndInput = 1:NInputs
    vol = spm_vol(NIFTI_files{IndInput});
    DataPos = vol.private.dat(find(vol.private.dat(:)>0))';
    DataNeg = vol.private.dat(find(vol.private.dat(:)<0))';
    
    if(~isempty(DataPos))
        RangePos = minmax(DataPos);
    else
        RangePos = [];
    end
    if(~isempty(DataNeg))
        RangeNeg = minmax(DataNeg);
    else
        RangeNeg = [];
    end
    Ranges{IndInput} = [RangePos, RangeNeg];
    
    if(length(Ranges{IndInput}(:))==4) %pos&neg
        PosColor = hot(128);
        PosColor = [repmat(PosColor(20,:),19,1);PosColor(20:end,:)];
        NegColor = winter(128);
    else
        if(length(Ranges{IndInput}(:))==2) %pos OR neg
            if(all(Ranges{IndInput}>0)) %pos
                PosColor = hot(128);
                PosColor = [repmat(PosColor(20,:),19,1);PosColor(20:end,:)];
                NegColor = [];
            else
                if(all(Ranges{IndInput}<0)) %neg
                    PosColor = [];
                    NegColor = winter(128);
                else
                    error(['Something is wrong with the limits! Ranges{IndInput}(:)==',num2str(Ranges{IndInput}(:)')]);
                end
            end
        else
            error(['Something is wrong with the limits! length(Ranges{IndInput}(:))==',num2str(length(Ranges{IndInput}(:)))]);
        end
    end
    Colors{IndInput} = {NegColor;PosColor};
end

%% Fill params for slover call: "ACTIVATION OVERLAY FROM NIFTI"
TransparencyOverlay = 1; %no transparency

params = cell(1,1);
params{1}.slices    = SliceIndices;   %Slice Indices of Slices to display
params{1}.transform = Orientation; %Slice Orientation

%structural image as background
params{1}.img(1).vol   = spm_vol(structural_img);  %get Structural Image
params{1}.img(1).cmap  = gray(256); %ColorMap of Structural Image
params{1}.img(1).range = minmax(params{1}.img(1).vol.private.dat(:)'); %Displayed Range of Structural Image i.e. all values
params{1}.img(1).prop  = 1; %Transparency setting for Structural Image
params{1}.img(1).type  = 'truecolour'; %Image which can be overlayed by other image.

%overlays
N = 1; %init
for IndInput = 1:NInputs
    Cols = Colors{IndInput};
    if(~isempty(Cols{1})&&~isempty(Cols{2})) %negative & positive
        params{1}.img(N+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+IndInput).cmap = Cols{1};  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+IndInput).range= [min(Ranges{IndInput}(:)),0]; %,max(Ranges{IndInput}(:))]; %range for colors
        if(PlotSpecial)
            params{1}.img(N+IndInput).func = ['i1(i1>=',num2str(Threshold_neg),')=NaN;']; %remove range above zeros
        else
            params{1}.img(N+IndInput).func = ['i1(i1>=0)=NaN;']; %remove range above zeros
        end
        params{1}.img(N+IndInput).hold = 0; %nearest neighbor interpolation
        params{1}.img(N+IndInput).prop = TransparencyOverlay; %Transparency setting for Overlay Image
        if(TransparencyOverlay<1)
            params{1}.img(N+IndInput).type = 'truecolour'; %ie change colors to show overlap
        else
            params{1}.img(N+IndInput).type = 'split';      %ie replace Structural below with its Value/Color
        end
        N = length(params{1}.img);
        params{1}.img(N+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+IndInput).cmap = Cols{2};  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+IndInput).range= [0,max(Ranges{IndInput}(:))*2/3]; %range for colors
        if(PlotSpecial)
            params{1}.img(N+IndInput).func = ['i1(i1<=',num2str(Threshold_pos),')=NaN;']; %remove range above zeros
        else
            params{1}.img(N+IndInput).func = ['i1(i1<=0)=NaN;']; %remove range above zeros
        end
        params{1}.img(N+IndInput).hold = 0; %nearest neighbor interpolation
        params{1}.img(N+IndInput).prop = TransparencyOverlay; %Transparency setting for Overlay Image
        if(TransparencyOverlay<1)
            params{1}.img(N+IndInput).type = 'truecolour'; %ie change colors to show overlap
        else
            params{1}.img(N+IndInput).type = 'split';      %ie replace Structural below with its Value/Color
        end
        N = length(params{1}.img);
    elseif(~isempty(Cols{1})) %negative
        params{1}.img(N+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+IndInput).cmap = Cols{1};  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+IndInput).range= [min(Ranges{IndInput}(:)),0]; %,max(Ranges{IndInput}(:))]; %range for colors
        if(PlotSpecial)
            params{1}.img(N+IndInput).func = ['i1(i1>=',num2str(Threshold_neg),')=NaN;']; %remove range above zeros
        else
            params{1}.img(N+IndInput).func = ['i1(i1>=0)=NaN;']; %remove range above zeros
        end
        params{1}.img(N+IndInput).hold = 0; %nearest neighbor interpolation
        params{1}.img(N+IndInput).prop = TransparencyOverlay; %Transparency setting for Overlay Image
        if(TransparencyOverlay<1)
            params{1}.img(N+IndInput).type = 'truecolour'; %ie change colors to show overlap
        else
            params{1}.img(N+IndInput).type = 'split';      %ie replace Structural below with its Value/Color
        end
        N = length(params{1}.img);
    elseif(~isempty(Cols{2})) %positive
        params{1}.img(N+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "stats<0"
        params{1}.img(N+IndInput).cmap = Cols{2};  %ColorMap of Overlay (here: "activation" or "deactivation" or both)
        params{1}.img(N+IndInput).range= [0,max(Ranges{IndInput}(:))]; %range for colors
        if(PlotSpecial)
            params{1}.img(N+IndInput).func = ['i1(i1<=',num2str(Threshold_pos),')=NaN;']; %remove range above zeros
        else
            params{1}.img(N+IndInput).func = ['i1(i1<=0)=NaN;']; %remove range above zeros
        end
        params{1}.img(N+IndInput).hold = 0; %nearest neighbor interpolation
        params{1}.img(N+IndInput).prop = TransparencyOverlay; %Transparency setting for Overlay Image
        if(TransparencyOverlay<1)
            params{1}.img(N+IndInput).type = 'truecolour'; %ie change colors to show overlap
        else
            params{1}.img(N+IndInput).type = 'split';      %ie replace Structural below with its Value/Color
        end
        N = length(params{1}.img);
    end
    
end
if(NInputs>1)
    params{1}.cbar = 1+[1:NInputs];    %Only display Colorbar for Overlay
else
    if(length(params{1}.img)>2)
        params{1}.cbar = 2:3;
    else
        params{1}.cbar = 2;
    end
end
if(PlotSpecial)
    params{1}.xslices = 6; %number of slices per row (max)
    params{1}.cbar = []; %no colorbar
end

%% make text
[BaseDir, fName, ext] = fileparts(NIFTI_files{1});
[BaseDir, TopDir1] = fileparts(BaseDir);
[BaseDir, TopDir2] = fileparts(BaseDir);

ResultsLast2Dirs = ['..',filesep,TopDir2,filesep,TopDir1,filesep];

params{1}.printfile    = strrep([fName,'Overlay'],' ','_');

text_annot   = cell(NInputs+1,1);
text_annot{1}= ['Structural-Image: ',structural_img];
for IndInput = 1:NInputs
    [BaseDir, fName, ext] = fileparts(NIFTI_files{IndInput});
    text_annot{1+IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',filesep,fName,ext];
end

clear BaseDir fName TopDir1 TopDir2 ResultsLast2Dirs % Step_cmap thresInd

%clear to avoid pain on rerun. ;)
clear check SliceIndices Orientation structural_img 

%% make overlay & ask for saveing or not
%% create slover object & "print" to graphics window
obj = cell(1,1);
if(isfield(params{1},'img'))
    %% add a new window and try to leave space for the text 
    params{1}.figure        = spm_figure(); %new figure
    params{1}.area.position = [0.005,0.005,0.99,0.95];
    params{1}.area.units    = 'normalized';
    params{1}.area.halign   = 'center';
    params{1}.area.valign   = 'middle';
    pause(1);
    
    %% Call slover to construct obj for paint
    obj{1} = slover(params{1});
    if(1)%for debugging
        obj{1}.printstr = [obj{1}.printstr,' -append'];
    end
    
    %% paint slices
%     spm_figure('GetWin','Graphics');
    paint(obj{1});
    drawnow;
    
    %% write annotations from SPM.mat & xSPM-Struct
    Position= [0.005,0.96,0.99,0.05]; %[0.005,0.95,0.99,0.05]
    % create annotations
    axes('Position',Position,'Visible','off');
    text(0,0.0,text_annot,'FontSize',10);
    
    clear Position
    pause(0.5);
    drawnow;
end

%% output figure handle and/or colors and/or background colors used?
if(nargout<=1)
    varargout = cell(1,1);
    varargout{1} = params{1}.figure;
else
    if(nargout>1)
        varargout = cell(nargout,1);
        varargout{1} = params{1}.figure; %the first one is still the handle the rest is empty, maybe future extension will change that.
        if(nargout>2)
            varargout{2} = colors;
            varargout{3} = bg;
        else
            if(nargout==2)
                varargout{2} = colors;
            end
        end
    end
end

%% done
disp(' ');
disp('Done.');
disp(' ');
end

%% subfunction
%% disp_catch
function [] = disp_catch(CATCHobj,varargin)
if(nargin==2)
    disp(['Error occurred in function "',mfilename,'>',varargin{1},'"...']);
else
    disp(['Error occurred in function "',mfilename,'"...']);
end
disp([CATCHobj.identifier,': ',CATCHobj.message]);

end

%% getPossibleSlices
function [SuggestedSlices] = getPossibleSlices(NIFTI_files)
% This function finds a possible Lowest z-direction slice & highest
% z-direction slice and the slice step-width.
%
% SuggestedSlices(1) = Lowest     z-direction slice
% SuggestedSlices(2) = Step-width z-direction slice
% SuggestedSlices(3) = Highest    z-direction slice
%

Data = cell(length(NIFTI_files),1);
v2m  = cell(length(NIFTI_files),1);
for Ind = 1:length(NIFTI_files)
    NII_tmp  = nifti(NIFTI_files{Ind});
    Data{Ind}= NII_tmp.dat(:,:,:);
    V_tmp    = spm_vol(NIFTI_files{Ind});
    v2m{Ind} = V_tmp.mat;
end
[Coords_z,Res_z] = getCoords(Data,v2m);

try
    SuggestedSlices(1) = floor(min(Coords_z(:)));
    SuggestedSlices(2) = floor(min(Res_z(:)));
    SuggestedSlices(3) = ceil(max(Coords_z(:)));
catch
    error('Could not extract slices containing data. INPUT VOLUME SEEMS EMPTY!');
end
end

%% getCoords
function [Coords_z,Res_z] = getCoords(Data,v2m)
Coords_z = [];
Res_z    = [];
for Ind = 1:length(Data)
    v2m_tmp    = v2m{Ind};  %transformation matrix from voxel coordinate to world(mm) coordinate
    vox        = getNonZeroVox(Data{Ind}); %get voxel coordinates of non-zero voxels, i.e. where atlas or maps are.
    Coords_tmp = zeros(size(vox)); %init temporary coords
    for i=1:size(vox,1)
        Coords_tmp(i,1:3)=vox(i,:)*v2m_tmp(1:3,1:3) + v2m_tmp(1:3,4)';
    end    
    Coords_z = [Coords_z; Coords_tmp(:,3)];%z-coordinates values
    Res_z    = [Res_z; abs(v2m_tmp(3,3))]; %avoid any odd errors by defining dimensions positive. This usually shouldn't happen because only SPM uses a negative dimension for X-direction.
end

end

%% getNonZeroVox
function [vox] = getNonZeroVox(Data)
[X,Y,Z] = ind2sub(size(Data),find(Data~=0));
vox = zeros(length(Z),3);
vox(:,1) = X;
vox(:,2) = Y;
vox(:,3) = Z;

end