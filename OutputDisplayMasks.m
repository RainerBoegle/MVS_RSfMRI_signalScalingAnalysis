%% Get MVSscaling struct
ScalingAnaMatPath = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
load(ScalingAnaMatPath);

%% Get DataForPlotting struct
DataForPlottingMatPath = spm_select(1,'mat','Select DataForPlotting.mat-file...');
load(DataForPlottingMatPath);

[BasePath, FName, ext] = fileparts(DataForPlottingMatPath);

%% select a mask that should be output
[MaskInds,OK] = listdlg('ListString',DataForPlotting.EffectOfMRI_split.AllMasksNames,'Name','Select mask','PromptString','Select masks that should be output: ','CancelString','Quit');
if(~OK)
    return;
end

%% make file name and get volume template struct for writing back out
OutDir = [BasePath,filesep,'DisplayMaskOutput'];
if(~exist(OutDir))
    mkdir(OutDir);
end
V_out = spm_vol(MVSscaling.Masks.MPaths{1});
if(V_out.dt(1)<16)
    V_out.dt(1) = 16; %not necessary but save
end

%% get masks and make into nifti
for IndMask = 1:length(MaskInds)
    MaskData = DataForPlotting.EffectOfMRI_split.AllMasks(:,MaskInds(IndMask));
    
    V_out.fname = [OutDir,filesep,regexprep(DataForPlotting.EffectOfMRI_split.AllMasksNames{MaskInds(IndMask)},' ','_'),'_',FName,'.nii'];
    
    %% write back out
    Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
    Y(MVSscaling.Masks.WholeBrainRaw~=0) = MaskData;
    Y = reshape(Y,V_out.dim);
    
    V_out = spm_write_vol(V_out,Y);
end

%% write out all together in one ROI map
V_outAll = V_out;
Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
for Ind= 1:length(DataForPlotting.EffectOfMRI_split.AllMasksNames)
    Y(MVSscaling.Masks.WholeBrainRaw~=0) = Y(MVSscaling.Masks.WholeBrainRaw~=0)+Ind.*DataForPlotting.EffectOfMRI_split.AllMasks(:,Ind);
end
Y = reshape(Y,V_out.dim);
V_outAll.fname = [OutDir,filesep,'AllMasks_',FName,'.nii'];

V_outAll = spm_write_vol(V_outAll,Y);

%% helpdlg for info
H = helpdlg(DataForPlotting.EffectOfMRI_split.AllMasksNames,'AllMasksNames');
H2= DisplayClusters(V_outAll.fname);

%% Done.
disp(' ');
disp('Done.');