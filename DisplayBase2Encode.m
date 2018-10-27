function [varargout] = DisplayBase2Encode(SliceSelType,varargin)
%DisplayClusters.m
%
% This Script can be used to produce Slice Overlay plots from NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING BASE2 ENCODED COMBINATIONS OF 3 MASKS.
% (Currently 3 is the limit, i.e. 1,2,4 --> overlaps are 3 5 6 & 7)
% THIS HAS TO BE ONE FILE!!!
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
%       H=DisplayBase2Encode(SliceSelType,Base2EncodeFile); %MaskFiles can be empty or not input then spm_select is used. 
%                                                  %SliceSelType can be empty or not input then the default 'Cluster' is used.
%
%
%SliceSelType can either be
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
%
%BUT IT ONLY MAKES SENSE TO USE 'Cluster' HERE.
%(See SuggestSlices.m for more help.)
%
%
%
%
%NB: A future extension will be made with the function AutoDisplayClusters.m
%    That function will allow settings to be given via a struct, such that this function 
%    can be called in a loop for plotting a series of images automatically.
%
%V1.6
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.6: 15.03.2015: Additional support of slice selection beyond 'Simple'. V1.5: (31.01.2015): Additional outputs and inputs for more control. V1.0: (14.12.2014): initial implementation (changed version from DisplayMasksNIFTIsOnSlices.m)

%% check basic inputs --> slice selection type
try
    if(ischar(SliceSelType))
        disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
    end
catch
    SliceSelType = 'Cluster'; %'Simple'
    disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
end

%% select image or images for display
if(isempty(varargin))
    NIFTI_files = spm_select(1,'image','Select NIFTI-file for creation of BASE2ENCODE Overlay ...');
    NIFTI_files = cellstr(NIFTI_files);
    WaitTime = 2; %s DEFAULT;
else
    if(length(varargin)<=2)
        NIFTI_files = varargin{1};
        if(~iscell(NIFTI_files))
            NIFTI_files = cellstr(NIFTI_files);
        end
        if(length(varargin)==2)
            WaitTime = varargin{2};
        elseif(length(varargin)>2)
            error('Too many input arguments');
        elseif(length(varargin)==1)
            WaitTime = 2; %s DEFAULT;
        end
    end
end

%% Slices selection & orientation & transparency
if(isnumeric(SliceSelType))
    SliceIndices = SliceSelType;
else
    SliceIndices = SuggestSlices(NIFTI_files,SliceSelType); %'Cluster'); %'Simple'); %
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
NInputs   = length(NIFTI_files);
DataTmp   = cell(NInputs,1); %the volumes to be analysed for unique values.
DataFName = cell(NInputs,1); %keep names for comparison
DataExtNum= ones(NInputs,1); %keep extension numbers for selection of volume in case of 4D file that matches with others
UniqueNonZeroVals = cell(NInputs,1); %the unique values in the volumes.
Is4D      = zeros(NInputs,1); %keep track of 4D files
disp('Creating overlay from: ');
for IndInput = 1:NInputs
    disp(NIFTI_files{IndInput});
    [BaseTmp,FNameTmp,ExtTmp] = fileparts(NIFTI_files{IndInput}); clear BaseTmp
    DataFName{IndInput} = FNameTmp;
    StartInd = regexp(ExtTmp,'\d');
    if(isempty(StartInd))
        DataExtNum(IndInput) = 1;
    else
        DataExtNum(IndInput) = str2num(ExtTmp(StartInd));
    end
    NIItmp = nifti(NIFTI_files{IndInput});
    Dim    = size(NIItmp.dat);
    NDim   = length(Dim);
    if(NDim>3)
        DataTmp{IndInput} = NIItmp.dat(:,:,:,:); 
        if(Dim(4)>1)
            Is4D(IndInput)    = 1;
        end
    else
        if(NDim==3)
            DataTmp{IndInput} = NIItmp.dat(:,:,:);
        else
            error(['size(dat)=',num2str(Dim),' NDim=',num2str(NDim)]);
        end
    end        
end
if(any(Is4D)) %there are 4D files let's just assign the correct volume according to spm-select ie DataExtNum
    Indices = find(Is4D);
    for IndInput = 1:length(Indices)
        Dat = DataTmp{Indices(IndInput)};
        DataTmp{Indices(IndInput)} = Dat(:,:,:,DataExtNum(Indices(IndInput)));
    end
end
AllUniqueNonZeroVals = [];
for IndInput = 1:NInputs
    uniqueTmp = unique(DataTmp{IndInput});
    if(any(uniqueTmp(:)==0)) %remove zero
        uniqueTmp(uniqueTmp==0) = [];
    end
    uniqueTmp(isnan(uniqueTmp)) = []; %remove NaNs
    if(~isempty(uniqueTmp))
        UniqueNonZeroVals{IndInput} = uniqueTmp(:);
        AllUniqueNonZeroVals = unique([AllUniqueNonZeroVals;uniqueTmp(:)]);
    else
        error(['Input ',num2str(IndInput),' "',DataFName{IndInput},'" has no NONZERO voxels! (',NIFTI_files{IndInput},')']);
    end
end
if(size(AllUniqueNonZeroVals,2)>1)
    AllUniqueNonZeroVals = AllUniqueNonZeroVals'; %transpose to get column vector
end

%decide if one color per input file or n-colors per n-unique values for a single input
if(NInputs>1)
    if(strcmp('NColors=NInputs',questdlg('As many colors as input files or as many as the sum of all unique values in the NInputs?','NColors=?','NColors=NInputs','NColors=sum(unique(Data(1:NInputs)))','NColors=NInputs')))
        NColors = NInputs;
    else
        NColors = 0;
        for IndInput = 1:NInputs
            NColors = NColors+length(UniqueNonZeroVals{IndInput}); %FUTURE EXTENSION WILL BE TO ALLOW MULTIPLE COLORS PER INPUT, BUT HAVE TO ASK USER THEN.
        end
    end
else
    if(NInputs==1)
        NColors = length(UniqueNonZeroVals{1});
    end
end


%% beautiful colors needed#111111!!!
%           1      2     3(1+2)    4     5(1+4)  6(2+4)  7(1+2+4)
colors = [0 1 0; 0 0 1;  0 1 1;  1 0 0;  1 1 0;  1 0 1;  1 1 1]; %future extension will be combining more cases, but that is currently not possible with just red blue green combination 


%% Fill params for slover call: "ACTIVATION OVERLAY FROM NIFTI"
if(NInputs>1)
    TransparencyOverlay = 1;%.9999; %a little bit transparent
else
    TransparencyOverlay = 1;  %no   transparency
end

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
for IndInput = 1:NInputs
    params{1}.img(1+IndInput).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "atlas"
    params{1}.img(1+IndInput).cmap = colors;  %ColorMap of Overlay (here: "activation")
    
    params{1}.img(1+IndInput).range= [.95 7.01]; %first to last value
    params{1}.img(1+IndInput).func = 'i1(i1<1)=NaN;'; %remove zeros
    
    params{1}.img(1+IndInput).hold = 0; %nearest neighbor interpolation
    params{1}.img(1+IndInput).prop = TransparencyOverlay; %Transparency setting for Overlay Image
    if(TransparencyOverlay<1)
        params{1}.img(1+IndInput).type = 'truecolour'; %ie change colors to show overlap
    else
        params{1}.img(1+IndInput).type = 'split';      %ie replace Structural below with its Value/Color
    end
end
params{1}.cbar = 2;
params{1}.xslices = 6; %number of slices per row (max)

%% make text
[BaseDir, fName, ext] = fileparts(NIFTI_files{1});
[BaseDir, TopDir1] = fileparts(BaseDir);
[BaseDir, TopDir2] = fileparts(BaseDir);

ResultsLast2Dirs = ['..',filesep,TopDir2,filesep,TopDir1,filesep];

params{1}.printfile    = strrep([fName,'Overlay'],' ','_');

text_annot   = cell(NInputs+2,1);
text_annot{1}= ['Structural-Image: ',structural_img];
for IndInput = 1:NInputs
    [BaseDir, fName, ext] = fileparts(NIFTI_files{IndInput});
    text_annot{1+IndInput} = ['Overlay-Image ',num2str(IndInput),':   ..',filesep,fName,ext];
end
text_annot{2+NInputs} = ['Areas are assigned by ',num2str(NColors),' distinguishable colors.'];

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
%     text(0,0.0,text_annot,'FontSize',10);
    
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

