%% warning
h = helpdlg('WARNING!!! The paths are hardcoded for a special case, i.e the DMN component. For others you need to change the paths accordingly.','WARNING!!!');
uiwait(h);

%% zData files
DataFolder = 'C:\Users\Reniar\Dropbox\CodesBackup\MVSfMRI_ScalingAnalysis\MVSscalingData_DMNmod\TestScalingDataSEARCHLIGHT_NHood1\SpecialCollectionOfResults\orgZVals_MedDifferent2sqrt2';
zAgg_path  = 'SL1_RANKSIGNzValsAggregate_MedDifferent2sqrt2.nii';
zMedS_path = 'SL1_RANKSIGNzValsMedS_MedDifferent2sqrt2.nii';
zMedSL_path= 'SL1_RANKSIGNzValsMedSL_MedDifferent2sqrt2.nii';

disp(['Loading zData from folder "',DataFolder,'"...']);
NII_zAgg  = nifti([DataFolder,filesep,zAgg_path]);
NII_zMedS = nifti([DataFolder,filesep,zMedS_path]);
NII_zMedSL= nifti([DataFolder,filesep,zMedSL_path]);

%% Brain Mask
BrainMask_path = ['C:\Users\Reniar\Dropbox\CodesBackup\MVSfMRI_ScalingAnalysis\MVSscalingData_DMNmod\TestScalingDataSEARCHLIGHT_NHood1\MedDifferent2sqrt2_V1',filesep,'SL1_UsableData_CropFINAL_NewWholeBrainMask_MELODIC_mask_MedDifferent2sqrt2.nii'];

NII_BrainMask = nifti(BrainMask_path);
BrainMask = NII_BrainMask.dat(:,:,:)>0;

%% Get zData of each that is in BrainMask and calculate d, log10d & log2d
CutOff = 1.96; %z-val cutoff

zData_Agg  = NII_zAgg.dat(:,:,:);
zData_MedS = NII_zMedS.dat(:,:,:);
zData_MedSL= NII_zMedSL.dat(:,:,:);

disp(['Treating zData "',zAgg_path,'"...']);
d_Agg = nan(size(zData_Agg)); dlog10_Agg = d_Agg; dlog2_Agg = d_Agg; %init to NaN
d_Agg(BrainMask>0)      = CutOff./abs(zData_Agg(BrainMask>0));
dlog10_Agg(BrainMask>0) = log10(d_Agg(BrainMask>0));
dlog2_Agg( BrainMask>0) = log2( d_Agg(BrainMask>0));

disp(['Treating zData "',zMedS_path,'"...']);
d_MedS = nan(size(zData_MedS)); dlog10_MedS = d_MedS; dlog2_MedS = d_MedS; %init to NaN
d_MedS(BrainMask>0)      = CutOff./abs(zData_MedS(BrainMask>0));
dlog10_MedS(BrainMask>0) = log10(d_MedS(BrainMask>0));
dlog2_MedS( BrainMask>0) = log2( d_MedS(BrainMask>0));

disp(['Treating zData "',zMedSL_path,'"...']);
d_MedSL = nan(size(zData_MedSL)); dlog10_MedSL = d_MedSL; dlog2_MedSL = d_MedSL; %init to NaN
d_MedSL(BrainMask>0)      = CutOff./abs(zData_MedSL(BrainMask>0));
dlog10_MedSL(BrainMask>0) = log10(d_MedSL(BrainMask>0));
dlog2_MedSL( BrainMask>0) = log2( d_MedSL(BrainMask>0));

%% write out to NIFTI
OutputFolder = 'C:\Users\Reniar\Dropbox\CodesBackup\MVSfMRI_ScalingAnalysis\MVSscalingData_DMNmod\TestScalingDataSEARCHLIGHT_NHood1\SpecialCollectionOfResults\2sqrt2_H0_AcceptanceRegions';

Vout = spm_vol([DataFolder,filesep,zAgg_path]);

disp(' ');
disp(['Writing dData to folder "',OutputFolder,'"...']);
disp( 'Writing Aggregate zData d, dlog10 and dlog2.');
%d_Agg
Vout.fname = [OutputFolder,filesep,'d_',zAgg_path];
spm_write_vol(Vout,d_Agg);
%dlog10_Agg
Vout.fname = [OutputFolder,filesep,'dlog10_',zAgg_path];
spm_write_vol(Vout,dlog10_Agg);
%dlog2_Agg
Vout.fname = [OutputFolder,filesep,'dlog2_',zAgg_path];
spm_write_vol(Vout,dlog2_Agg);

disp('Writing MedianSubject zData d, dlog10 and dlog2.');
%d_MedS
Vout.fname = [OutputFolder,filesep,'d_',zMedS_path];
spm_write_vol(Vout,d_MedS);
%dlog10_MedS
Vout.fname = [OutputFolder,filesep,'dlog10_',zMedS_path];
spm_write_vol(Vout,dlog10_MedS);
%dlog2_MedS
Vout.fname = [OutputFolder,filesep,'dlog2_',zMedS_path];
spm_write_vol(Vout,dlog2_MedS);

disp('Writing MedianSL zData d, dlog10 and dlog2.');
%d_MedSL
Vout.fname = [OutputFolder,filesep,'d_',zMedSL_path];
spm_write_vol(Vout,d_MedSL);
%dlog10_MedS
Vout.fname = [OutputFolder,filesep,'dlog10_',zMedSL_path];
spm_write_vol(Vout,dlog10_MedSL);
%dlog2_MedS
Vout.fname = [OutputFolder,filesep,'dlog2_',zMedSL_path];
spm_write_vol(Vout,dlog2_MedSL);

%% Done.
disp(' ');
disp('DONE.');