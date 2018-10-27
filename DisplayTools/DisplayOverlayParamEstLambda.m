function varargout = DisplayOverlayParamEstLambda(SliceSelType,NIFTI_files,Colors,Limits)
% This function can be used to produce Slice Overlay plots from NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING PARAMETER-ESTIMATE ANALYSIS RESULTS
% OF MVSscaling ANALYSIS.
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
% SliceSelType can be empty or not input then the default is 'Cluster', see SuggestSlices.m for more help.
% NIFTI_files can be empty, not input then spm_select is used. 
% Colors determines the colors that are used for the overlay, ONE FOR EACH VOLUME IN NIIpaths. 
% Limits is just additional to colors for determining the correct numbers of inputs and limits.
% NB: multiple files will be plotted in separate windows 
%     However, all of them will have the same slices picked for them.
%
% Suggested Limits = [0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)];
% Suggested Colors:   D-blue;  blue;  cyan; D-green; green; D-yellow; yellow; orange;   red;      red+;       red++;        redLimit;  white;
%           Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1]; 
% Suggested SlicesForDisplay = [-42; -36; -26; -18; -12; -4; +2; +10; +16; +24; +30; +36; +42; +48];
%
%
% USAGE:
%                    varargout              = DisplayOverlayParamEstLambda(SliceSelType,NIIpaths,Colors,Limits); 
%       [slover_obj,H,Colors,Limits,Slices] = DisplayOverlayParamEstLambda(SliceSelType,NIIpaths,Colors,Limits); 
%
%
%V1.1
%Date: V1.1(18.02.2016) small edit to add current file name to figure. V1.0(27.07.2015) (based on version 1.5 of DisplayClusters.m & 1.0 DisplayCIniis.m.)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)


%% check basic inputs --> slice selection type
try
    if(ischar(SliceSelType))
        disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
    end
catch
    SliceSelType = 'Cluster'; %'Simple'
    disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
end

%% select images for display
try
    if(isempty(NIFTI_files))
        NIFTI_files = spm_select([1 Inf],'image','Select NIFTI-file(s) for creation of Overlay ...');
        NIFTI_files = cellstr(NIFTI_files);
    else
        if(ischar(NIFTI_files))
            NIFTI_files = cellstr(NIFTI_files);
        end
    end
catch
    NIFTI_files = spm_select([1 Inf],'image','Select NIFTI-file(s) for creation of Overlay ...');
    NIFTI_files = cellstr(NIFTI_files);
end
NInputs = length(NIFTI_files);

%% check colors together with Limits 
try
    if(isempty(Limits))
        [Limits,SortInds] = sort([0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)],'descend');
    end
catch
    [Limits,SortInds] = sort([0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)],'descend');
end

try
    if(isempty(Colors))
        %         D-blue;  blue;  cyan; D-green; green; D-yellow; yellow; orange;   red;      red+;       red++;        redLimit;  white;
        Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1]; %future extension in case we interpolate more.
        if(exist('SortInds','var'))
            Colors = Colors(SortInds,:);
        end
    end
catch
    %         D-blue;  blue;  cyan; D-green; green; D-yellow; yellow; orange;   red;      red+;       red++;        redLimit;  white;
    Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1]; %future extension in case we interpolate more.
    if(exist('SortInds','var'))
        Colors = Colors(SortInds,:);
    end
end
%don't show cbar
cbar = [];

%check
if(length(Limits)~=size(Colors,1))
    error('Number of Colors must be the same as Limits.');
end

WaitTime = 2; %s DEFAULT;

%% Slices selection & orientation & transparency
if(ischar(SliceSelType))
    SliceIndices = SuggestSlices(NIFTI_files,SliceSelType); %'Cluster'); %'Simple'); %
else
    if(isnumeric(SliceSelType))
        SliceIndices = SliceSelType;
    else
        error('SliceSelType is wrong format.');
    end
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

%% Fill params for slover call
TransparencyOverlay = 1;  %1==no   transparency

params = cell(NInputs,1); %init parameters for slover object
obj    = cell(NInputs,1); %init the slover object
disp('Displaying: ');
for IndInput = 1:NInputs
    params{IndInput}.slices    = SliceIndices;   %Slice Indices of Slices to display
    params{IndInput}.transform = Orientation; %Slice Orientation
    
    %structural image as background
    params{IndInput}.img(1).vol   = spm_vol(structural_img);  %get Structural Image
    params{IndInput}.img(1).cmap  = gray(256); %ColorMap of Structural Image
    params{IndInput}.img(1).range = minmax(params{IndInput}.img(1).vol.private.dat(:)'); %Displayed Range of Structural Image i.e. all values
    params{IndInput}.img(1).prop  = 1; %Transparency setting for Structural Image
    params{IndInput}.img(1).type  = 'truecolour'; %Image which can be overlayed by other image.
    
    
    %overlays
    params{IndInput}.img(2).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "atlas"
    params{IndInput}.img(2).cmap = Colors;  %ColorMap of Overlay (here: "activation")
    params{IndInput}.img(2).range= [1 length(Limits)]; %V1: [1 6]; %V2: [1 13]; %minmax(params{IndInput}.img(2).vol.private.dat(:)'); %Displayed Range i.e. all values %first to last value
    params{IndInput}.img(2).func = 'i1(i1<.0001)=NaN;'; %remove zeros
    
    %transparency & interpolation + colorbar
    params{IndInput}.img(2).hold = 0; %nearest neighbor interpolation
    params{IndInput}.img(2).prop = TransparencyOverlay; %Transparency setting for Overlay Image
    if(TransparencyOverlay<1)
        params{IndInput}.img(2).type = 'truecolour'; %ie change colors to show overlap
    else
        params{IndInput}.img(2).type = 'split';      %ie replace Structural below with its Value/Color
    end
    params{IndInput}.cbar = cbar; %(don't?) show colorbars
    params{IndInput}.xslices = 6; %number of slices per row (max)
    
    
    %% make text
    [BaseDir, fName, ext] = fileparts(NIFTI_files{IndInput});
    [BaseDir, TopDir1] = fileparts(BaseDir);
    [BaseDir, TopDir2] = fileparts(BaseDir);
    
    ResultsLast2Dirs = ['..',filesep,TopDir2,filesep,TopDir1,filesep];
    
    text_annot    = cell(3,1);
    text_annot{1} = ['Structural-Image: ',structural_img];
    text_annot{2} = ['Overlay-Image "',fName,ext,'".'];
    text_annot{3} = ['Overlay-Image from ',ResultsLast2Dirs];
    
    params{1}.printfile    = strrep([fName,'Overlay'],' ','_');
    
    clear BaseDir fName TopDir1 TopDir2 ResultsLast2Dirs % Step_cmap thresInd
    
    %% make overlay & ask for saveing or not
    %% create slover object & "print" to graphics window
    if(isfield(params{1},'img'))
        %% add a new window and try to leave space for the text
        params{IndInput}.figure        = spm_figure(); %new figure
        params{IndInput}.area.position = [0.005,0.005,0.99,0.95];
        params{IndInput}.area.units    = 'normalized';
        params{IndInput}.area.halign   = 'center';
        params{IndInput}.area.valign   = 'middle';
        pause(1);
        
        %% Call slover to construct obj for paint
        obj{IndInput} = slover(params{IndInput});
        if(1)%for debugging
            obj{IndInput}.printstr = [obj{IndInput}.printstr,' -append'];
        end
        
        %% paint slices
        %     spm_figure('GetWin','Graphics');
        paint(obj{IndInput});
        drawnow;
        
        %% write annotations from SPM.mat & xSPM-Struct
        Position= [0.005,0.96,0.99,0.05]; %[0.005,0.95,0.99,0.05]
        % create annotations
        axes('Position',Position,'Visible','off');
        text(0,0.0,text_annot,'FontSize',10);
        
        clear Position
        pause(0.5);
        drawnow;
    else
        disp(['For Input ',num2str(IndInput),' field "img" is missing.']);
    end
end
%% output figure handle and/or colors and/or background colors used?
H = cell(NInputs,1);
for IndInput = 1:NInputs
    H{IndInput} = params{IndInput}.figure;
end
varargout{1} = obj; %the slover object for later print & saving
varargout{2} = H; %the handle of the figure
varargout{3} = Colors; %colors of course
varargout{4} = Limits; %Limits of course
varargout{5} = SliceIndices; %Slices


%% done
disp(' ');
disp('Done.');
disp(' ');
end