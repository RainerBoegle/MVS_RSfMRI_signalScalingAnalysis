function [varargout] = DisplayLogDniis(SliceSelType,varargin)
%DisplayCIniis.m
%
% This Script can be used to produce Slice Overlay plots from NIFTI-files.
% THIS SCRIPT IS MAINLY USEFUL FOR PLOTTING Log2/10/whatever(d)-FILES INDICATING THE ACCEPTANCE OR
% REJECTION REGIONS FOR H0 (Lambda = 2sqrt(2)) OF MVSscaling ANALYSIS.
%
% The SPM toolbox "slover" originally created by Matthew Brett (http://imaging.mrc-cbu.cam.ac.uk/imaging/MatthewBrett)
% is used to produce the Overlays.
%
% The main purpose of this Script is to make the use of "slover" easier.
%
% All parameters are fixed (more or less).
%
% USAGE:
%       H=DisplayLogDniis(SliceSelType,LogDniifiles); %LogDniifiles can be empty, not input then spm_select is used. 
%                                              %NB: multiple files will be plotted in separate windows 
%                                                   However, all of them will have the same slices picked for them.
%                                              %SliceSelType can be empty or not input then the default is 'Cluster',
%                                               see SuggestSlices.m for more help.
%
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0: initial implementation (based on version 1.5 of DisplayClusters.m)

%% check basic inputs --> slice selection type
try
    if(ischar(SliceSelType))
        disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
    else
        if(isnumeric(SliceSelType))
            disp('Using preselected slices...');
        end
    end
catch
    SliceSelType = 'Cluster'; %'Simple'
    disp(['Will select Slices using Selection mode "',SliceSelType,'"']);
end

%% select image or images for display
if(isempty(varargin))
    NIFTI_files = spm_select([1 Inf],'image','Select NIFTI-file(s) for creation of Overlay ...');
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

%% check inputs and select colors
NInputs   = length(NIFTI_files);

%% my own color design --> different in the future maybe
% -Inf ... -Gpoint --> winter(64) ending in [0 1 0.5]
NegCols = winter(64);
NegCols = NegCols(end:-1:1,:);

% -Gpoint... +Gpoint --> bluegreenyellow(64) starting in [0 1 0.5] middle is [0 1 0] ending in [1 1 0]
Gpoint = 0.1; %0.25; %0.5; %bluegreenyellow range amplitude
bluegreenyellow         = repmat([0 1 0],64,1);
%bluegreenyellow(1:32,3) = linspace(0.5,0,32);
bluegreenyellow(32:-1:1,3)= (1-(linspace(1,0,32).^2))./2;
%bluegreenyellow(32:64,1)= linspace(0,1,33);
bluegreenyellow(32:64,1)= 1-(linspace(1,0,33).^2);

% +Gpoint... Inf   --> flipped autumn(64)
PosCols = autumn(64);
PosCols = PosCols(end:-1:1,:);

Colors = cell(3,1);
Colors{1} = NegCols; 
Colors{2} = bluegreenyellow; 
Colors{3} = PosCols;

%% determine ranges
for IndInput = 1:length(NIFTI_files)
    NII_tmp = nifti(NIFTI_files{IndInput});
    Data = NII_tmp.dat(:);
    Data = Data(~isnan(Data(:))); %remove NaNs
    Data = Data(~isinf(Data(:))); %remove Infinite parts
    NegData = Data(Data<0); %negative values
    PosData = Data(Data>0); %positive values
    
    %the winter part --> negative
    if(~isempty(NegData))
    if((min(NegData(:))<2*median(NegData(:)))||(min(NegData(:))<1.5*quantile(NegData(:),0.25))) %difference is too extreme
        disp(['Input ',num2str(IndInput),' negative data minimum (',num2str(min(NegData(:))),') is too far away from median (',num2str(median(NegData(:))),') or lower quartile (',num2str(quantile(NegData(:),0.25)),')! Defaulting to lower quartile.']);
        RangesOfInputs(:,1,IndInput) = [-Gpoint quantile(NegData(:),0.25)]; %limit by lower quartile
        if(-Gpoint<quantile(NegData(:),0.25)) %check
            disp(['Input ',num2str(IndInput),' negative data range is deficient! Defaulting to -2*Gpoint=',num2str(-2*Gpoint),'.']);
            RangesOfInputs(:,1,IndInput) = [-Gpoint -2*Gpoint];
        end
    else
        RangesOfInputs(:,1,IndInput) = [-Gpoint min(NegData(:))];
        if(-Gpoint<min(NegData(:))) %check
            disp(['Input ',num2str(IndInput),' negative data range is deficient! Defaulting to -2*Gpoint=',num2str(-2*Gpoint),'.']);
            RangesOfInputs(:,1,IndInput) = [-Gpoint -2*Gpoint];
        end
    end
    else
        disp('NO NEGATIVE DATA!');
        RangesOfInputs(:,1,IndInput) = [-Gpoint -2*Gpoint];
    end
    
    %the bluegreenyellow part --> region around zero 
    RangesOfInputs(:,2,IndInput) = [-Gpoint Gpoint]; %the bluegreenyellow part
    
    %the autumn part --> positive
    if((max(PosData(:))>2*median(PosData(:)))||(max(PosData(:))>1.5*quantile(PosData(:),0.75))) %difference is too extreme
        disp(['Input ',num2str(IndInput),' positive data maximum (',num2str(max(PosData(:))),') is too far away from median (',num2str(median(PosData(:))),') or upper quartile (',num2str(quantile(PosData(:),0.75)),')! Defaulting to upper quartile.']);
        RangesOfInputs(:,3,IndInput) = [Gpoint quantile(PosData(:),0.75)]; %limit by upper quartile
        if(Gpoint>quantile(PosData(:),0.75)) %check
            disp(['Input ',num2str(IndInput),' positive data range is deficient! Defaulting to +2*Gpoint=',num2str(+2*Gpoint),'.']);
            RangesOfInputs(:,3,IndInput) = [Gpoint 2*Gpoint];
        end
    else
        RangesOfInputs(:,3,IndInput) = [Gpoint max(PosData(:))];
        if(Gpoint>max(PosData(:))) %check
            disp(['Input ',num2str(IndInput),' positive data range is deficient! Defaulting to +2*Gpoint=',num2str(+2*Gpoint),'.']);
            RangesOfInputs(:,3,IndInput) = [Gpoint 2*Gpoint];
        end
    end
end
if(length(NIFTI_files)>1)
    RangesForOverlay = zeros(2,3); %init
    disp('Combining ranges of all inputs');
    disp('Negative: ');
    for IndInput = 1:length(NIFTI_files)
        disp(['Input ',num2str(IndInput),': [',num2str(RangesOfInputs(1,1,IndInput)),',',num2str(RangesOfInputs(2,1,IndInput)),']']);
    end
    RangesForOverlay(:,1) = min(squeeze(RangesOfInputs(:,1,:)),[],2); %negative range
    disp(['Combined: [',num2str(RangesForOverlay(1,1)),',',num2str(RangesForOverlay(2,1)),'].']);
    disp(' ');
    disp('NearZeroRange: ');
    for IndInput = 1:length(NIFTI_files)
        disp(['Input ',num2str(IndInput),': [',num2str(RangesOfInputs(1,2,IndInput)),',',num2str(RangesOfInputs(2,2,IndInput)),']']);
    end
    RangesForOverlay(:,2) = [-Gpoint Gpoint]; %the bluegreenyellow part
    disp(['Combined: [',num2str(RangesForOverlay(1,2)),',',num2str(RangesForOverlay(2,2)),'].']);
    disp(' ');
    disp('Positive: ');
    for IndInput = 1:length(NIFTI_files)
        disp(['Input ',num2str(IndInput),': [',num2str(RangesOfInputs(1,3,IndInput)),',',num2str(RangesOfInputs(2,3,IndInput)),']']);
    end
    RangesForOverlay(:,3) = max(squeeze(RangesOfInputs(:,3,:)),[],2); %negative range
    disp(['Combined: [',num2str(RangesForOverlay(1,3)),',',num2str(RangesForOverlay(2,3)),'].']);
    disp(' ');    
else
    RangesForOverlay = RangesOfInputs;
end
disp('Final ranges for overlays: ');
disp(['Negative:      [',num2str(RangesForOverlay(1,1)),',',num2str(RangesForOverlay(2,1)),'].']);
disp(['NearZeroRange: [',num2str(RangesForOverlay(1,2)),',',num2str(RangesForOverlay(2,2)),'].']);
disp(['Positive:      [',num2str(RangesForOverlay(1,3)),',',num2str(RangesForOverlay(2,3)),'].']);
disp(' ');

%don't show cbar
cbar = [2 3 4];


%% Show colorbar as a whole 
figure(81); imagesc(repmat((1:64)',1,10)); colormap([NegCols(end:-1:1,:);bluegreenyellow;PosCols]); title('Colors for plotting'); colorbar; axis('off');
H = helpdlg('These are the colors that have been generated.','Color-Generation Results');
if(isinf(WaitTime))
    uiwait(H);
else
    uiwait(H,WaitTime); 
end
try
    close(H);
end

%% Fill params for slover call: "ACTIVATION OVERLAY FROM NIFTI"
OpacityOverlay = 1;  %no   transparency

params = cell(NInputs,1); %init parameters for slover object
obj    = cell(NInputs,1); %init the slover object
disp('Displaying: ');
for IndInput = 1:NInputs
    disp(NIFTI_files{IndInput});
    params{IndInput}.slices    = SliceIndices;   %Slice Indices of Slices to display
    params{IndInput}.transform = Orientation; %Slice Orientation
    
    %structural image as background
    params{IndInput}.img(1).vol   = spm_vol(structural_img);  %get Structural Image
    params{IndInput}.img(1).cmap  = gray(256); %ColorMap of Structural Image
    params{IndInput}.img(1).range = minmax(params{IndInput}.img(1).vol.private.dat(:)'); %Displayed Range of Structural Image i.e. all values
    params{IndInput}.img(1).prop  = 1; %Transparency setting for Structural Image
    params{IndInput}.img(1).type  = 'truecolour'; %Image which can be overlayed by other image.
    
    %overlays
    for IndOverlay = 1:size(RangesForOverlay,2)
        params{IndInput}.img(1+IndOverlay).vol  = spm_vol(NIFTI_files{IndInput}); %get NIFTI assumed as "atlas"
        params{IndInput}.img(1+IndOverlay).cmap = Colors{IndOverlay};  %ColorMap of Overlay (here: "activation")
        
        params{IndInput}.img(1+IndOverlay).range= RangesForOverlay(:,IndOverlay); %Displayed Range i.e. all values %first to last value
        if(IndOverlay==1) %negative
            params{IndInput}.img(1+IndOverlay).func = ['i1(i1>=',num2str(-Gpoint),')=NaN;']; %remove upper & zero parts
        else
            if(IndOverlay==2) %zero part
                params{IndInput}.img(1+IndOverlay).func = ['i1(abs(i1)>',num2str(Gpoint),')=NaN;']; %remove upper & lower parts, but not zero parts
            else
                if(IndOverlay==3) %positive
                    params{IndInput}.img(1+IndOverlay).func = ['i1(i1<=',num2str(Gpoint),')=NaN;']; %remove lower & zero parts
                end
            end
        end
        
        %transparency & interpolation + colorbar
        params{IndInput}.img(1+IndOverlay).hold = 0; %nearest neighbor interpolation
        params{IndInput}.img(1+IndOverlay).prop = OpacityOverlay; %Transparency setting for Overlay Image
        if(OpacityOverlay<1)
            params{IndInput}.img(1+IndOverlay).type = 'truecolour'; %ie change colors to show overlap
        else
            params{IndInput}.img(1+IndOverlay).type = 'split';      %ie replace Structural below with its Value/Color
        end
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
    
    params{IndInput}.printfile    = strrep([fName,'Overlay'],' ','_');
    
    clear BaseDir fName TopDir1 TopDir2 ResultsLast2Dirs % Step_cmap thresInd
    
    %% make overlay & ask for saveing or not
    %% create slover object & "print" to graphics window
    if(isfield(params{IndInput},'img'))
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
%         Position= [0.005,0.96,0.99,0.05]; %[0.005,0.95,0.99,0.05]
%         % create annotations
%         axes('Position',Position,'Visible','off');
%         text(0,0.0,text_annot,'FontSize',10);
%         
%         clear Position
%         pause(0.5);
%         drawnow;
    else
        disp(['For Input ',num2str(IndInput),' field "img" is missing.']);        
    end
end

%% output figure handle and/or colors and/or background colors used?
Hfigs = cell(NInputs,1);
for IndInput = 1:NInputs
    Hfigs{IndInput} = params{IndInput}.figure;
end
varargout{1} = obj; %the slover object for later print & saving
varargout{2} = Hfigs; %the handle of the figure
varargout{3} = Colors; %colors of course

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

