function [ExpandAtlasOverlap_OutputStruct]=MakeAALOverlap(EffectsOfMRIMask_Path,OutputDir,TitleOfOverlap)
% This function creates the overlap with the AAL-Atlas given a mask 
% i.e. a (0/1)-image of the "EffectsOfMRI". 
%
% Furthermore it extends the overlap created if parts
% are missing and creates combined maps of areas for two reduced sets of Labels
% derived from the AAL-Atlas label-set.
%
% This includes a "basic combination":
%  - left & right parts of areas and multiple numbered labels are joined/combined,
%    i.e. Temporal_Mid_L & Temporal_Mid_R are joined into Temporal_Mid (bilateral)
%    and Vermis1 ... Vermis4 are joined into Vermis (total of all these parts.)
% and a "general combination":
%  - the joining as above and further joining by removing of secondary denoters of areas
%    e.g. like "Mid", "Inf", "Orb" and so on, i.e. Temporal_Mid & Temporal_Pole & Temporal_Inf
%    are joined into Temporal (total of all these parts).
%
%Usage:
%      [ExpandAtlasOverlap_OutputStruct]=MakeAALOverlap(EffectsOfMRIMask_Path,OutputDir);
%
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)

%% default path to AAL-Atlas
[FunctionDir] = fileparts(mfilename('fullpath'));
AAL_Atlas_Dir = [FunctionDir,filesep,'AAL_Atlas'];

%% check inputs, create matlabbatch and run batch
[matlabbatch,Output_fName] = FillBatch(EffectsOfMRIMask_Path,AAL_Atlas_Dir,OutputDir,TitleOfOverlap);
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch);

%% bugfix round results
NII_Overlap = nifti([OutputDir,filesep,Output_fName]); %get newly made overlap nii-file
Overlap_dat = round(NII_Overlap.dat(:,:,:)); %round to make it into index
Vo          = spm_vol([OutputDir,filesep,Output_fName]); %get volume handle for spm
if(Vo.dt(1)<16)
    Vo.dt(1) = 16; %not necessary but save
end
Vo          = spm_write_vol(Vo, round(Overlap_dat)); %write out rounded results to get indices written out.

%% Call expand overlap to create further masks
[ExpandAtlasOverlap_OutputStruct] = ExpandAtlasOverlap([OutputDir,filesep,Output_fName],EffectsOfMRIMask_Path,0,100,30); %no warning are displayed, 100-times the number of voxels is maximum iteration and if after 30 cycles nothing is changed we stop ahead of time.

end

%% subfunction fill matlabbatch
function [matlabbatch,Output_fName] = FillBatch(EffectsOfMRIMask_Path,AAL_Atlas_Dir,OutputDir,TitleOfOverlap)
% This function checks the inputs and fills the matlabbatch for spm processing.

%% check EffectsOfMRI file
if(~exist(EffectsOfMRIMask_Path))
    EffectsOfMRIMask_Path = spm_select(1,'image',['Select Mask-"',TitleOfOverlap,'" *.nii file...']);
end
[BaseDir,FName,Ext] = fileparts(EffectsOfMRIMask_Path);
Output_fName        = ['AAL_Overlap_',FName,Ext];

%% check AAL-Atlas exists
AAL_Atlas_Path      = [AAL_Atlas_Dir,filesep,'aal.nii'];
if(~exist(AAL_Atlas_Path))
    AAL_Atlas_Path = spm_select(1,'image','Select AAL-Atlas *.nii file...');
end

%% check if output directory exists otherwise make it
if(~exist(OutputDir,'dir'))
    mkdir(OutputDir);
end

%% fill matlabbatch
matlabbatch{1}.spm.util.imcalc.input          = {EffectsOfMRIMask_Path; AAL_Atlas_Path}; %Mask-"EffectsOfMRI" & AAL-Atlas
matlabbatch{1}.spm.util.imcalc.output         = Output_fName;
matlabbatch{1}.spm.util.imcalc.outdir         = {OutputDir};
matlabbatch{1}.spm.util.imcalc.expression     = 'i1.*i2'; %make sure orientation and resolution + size of the mask "EffectsOfMRI" are used.
matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;

end