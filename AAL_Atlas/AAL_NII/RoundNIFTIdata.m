%% get NIFTI to be rounded
NIIorg_Path = spm_select(1,'image','Select NIFTI-file to be rounded...');

NIIorg = nifti(NIIorg_Path);
NIIorg_data = NIIorg.dat(:);
%% get Info from NIIorg
V_NIIorg = spm_vol(NIIorg_Path);

%% change file name for output
V_NIIrounded = V_NIIorg;
V_NIIdiff    = V_NIIorg;

[BasePath,FName,Ext] = fileparts(V_NIIorg.fname);
V_NIIrounded.fname = [BasePath,filesep,'Rounded_',FName,Ext];
V_NIIdiff.fname    = [BasePath,filesep,'Diff_',   FName,Ext];
%% go over entries and round them
Yrounded = zeros(V_NIIrounded.dim);
Yrounded = round(NIIorg_data);

Ydiff = zeros(V_NIIdiff.dim);
Ydiff = Yrounded(:) - NIIorg_data(:);

%% write out volumes
Vo_rounded = spm_write_vol(V_NIIrounded, reshape(Yrounded,V_NIIrounded.dim));

Vo_diff    = spm_write_vol(V_NIIdiff,    reshape(Ydiff   ,V_NIIdiff.dim));

%% done
disp('Done');