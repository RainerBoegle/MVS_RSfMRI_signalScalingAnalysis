%% get scaling analysis results
ScalingAnaMatPath = spm_select(1,'mat','Select MVS_ScalingPerSubject.mat-file...',[],pwd,'^MVS_ScalingPerSubject.mat',1);
load(ScalingAnaMatPath);

[BasePath, FName, ext] = fileparts(ScalingAnaMatPath);

%% get Info from Mask
V_mask = spm_vol(ScalingPerSubject.MaskPath);
NII_mask = nifti(ScalingPerSubject.MaskPath);
IndicesMask = find(NII_mask.dat(:)==1);

%% change file name for output
V_file = V_mask;
if(V_file.dt(1)<16)
    V_file.dt(1) = 16; %not necessary but save
end

%% write out each volume into tempdir
AnzahlStellen = ceil(log10(size(ScalingPerSubject.Data2DMasked,2)))+1;
FormatStr = ['%0',num2str(AnzahlStellen),'.0f'];

tempdir = [BasePath,filesep,'tempdir'];
if(~exist(tempdir,'dir'))
    mkdir(tempdir);
end

for IndFile = 1:size(ScalingPerSubject.Data2DMasked,2)
    Y = zeros(V_file.dim);
    for idx = 1:length(IndicesMask)
        Y(IndicesMask(idx)) = ScalingPerSubject.Data2DMasked(idx,IndFile);
    end
    
    V_file.fname = [tempdir,filesep,'Inputs_',FName,'_',num2str(IndFile,FormatStr),'.nii'];
    Vo(IndFile) = spm_write_vol(V_file, Y);
end

%% use fslmerge -t to merge into 4D vol
V4 = spm_file_merge(Vo,[BasePath,filesep,'Inputs4D_',FName,'.nii'],0);

%% delete tempdir & files
try
    if(exist(tempdir,'dir'))
        delete([tempdir,filesep,'*.*']);
        rmdir(tempdir);
    end
end

%% save Grouping in same folder as 4D NIFTI
Grouping = ScalingPerSubject.Design;
Grouping.Associated4DFile = [BasePath,filesep,'Inputs4D_',FName,'.nii'];
Grouping.InfoStr = {'Asians==1';'Caucasians==2'};
save([BasePath,filesep,'Grouping_Inputs4D_',FName,'.mat'],'Grouping');

%% done
disp('Done');