%% This is just a test script
%% Future version needs to improve on this
%% select filelist and find associated stdev data for full timeseries --> make a list
%% select Stdev of residuals together with dual-reg folder --> select in order and do processing

%% select filelist and then select Stdev of full timeseries
FileListPath = spm_select(1,'.filelist','Select filelist...');
FileList = importdata(FileListPath);

StdevFullTS_FileList = cellstr(spm_select(Inf,'image','Select Stdev_"FullTimeSeries" files...'));
rMSSDFullTS_FileList = cellstr(spm_select(Inf,'image','Select rMSSD_"FullTimeSeries" files...'));

%% make new filelist that brings this in order that .filelist has --> that is the order in the dr_stage2_ic????.nii files
NewFileListStdevFullTS  = cell(length(FileList),1);
NewFileList_rMSSDFullTS = cell(length(FileList),1);
for IndFile = 1:length(FileList)
    CurrStr = FileList{IndFile};
    CurrSubj= regexp(CurrStr,'Subj_p0\d\d','match');
    MRIStr = regexp(CurrStr,'/AERA/|/GE/','match');
    if(strcmp(MRIStr{1},'/AERA/'))
        CurrMRINum = 1;
    else
        CurrMRINum = 2;
    end
    CurrRun = regexp(CurrStr,'rsfMRI_\d','match');
    
    
    SearchStr = ['Stdev_',CurrSubj{1},'_MRI',num2str(CurrMRINum),'_',CurrRun{1},'_filtered_func_data.nii'];
    Idx = find(~cellfun(@isempty,strfind(StdevFullTS_FileList,SearchStr)));
    NewFileListStdevFullTS{IndFile} = StdevFullTS_FileList{Idx};
    NewFileList_rMSSDFullTS{IndFile}= rMSSDFullTS_FileList{Idx};
end

%% select dual reg dir & Stdev of residuals in order...
dual_reg_dir = spm_select(1,'dir','Select directory with DualRegAmps...');
MaskPath = spm_select(1,'image','Select mask...');

FileListStdevResiduals = cellstr(spm_select(Inf,'any','Select Stdev of Residuals timeseries IN ORDER MATCHING VOLUMES IN dr_stage2_ic????.nii FILES...'));
FileListrMSSDResiduals = cellstr(spm_select(Inf,'any','Select rMSSD of Residuals timeseries IN ORDER MATCHING VOLUMES IN dr_stage2_ic????.nii FILES...'));



%% do processing
DualRegStage2ICFiles     = cellstr(spm_select('List',dual_reg_dir,'^dr_stage2_ic\d\d\d\d.nii'));
DualRegStage2ICFilePaths = cell(length(DualRegStage2ICFiles),1);
for IndICFile = 1:length(DualRegStage2ICFiles)
    DualRegStage2ICFilePaths{IndICFile} = [dual_reg_dir,filesep,DualRegStage2ICFiles{IndICFile}];
    CalcAdjDualRegAmps([dual_reg_dir,filesep,DualRegStage2ICFiles{IndICFile}],NewFileListStdevFullTS, MaskPath,'AdjStdevInputs_');
    CalcAdjDualRegAmps([dual_reg_dir,filesep,DualRegStage2ICFiles{IndICFile}],NewFileList_rMSSDFullTS,MaskPath,'AdjrMSSDInputs_');
    CalcAdjDualRegAmps([dual_reg_dir,filesep,DualRegStage2ICFiles{IndICFile}],FileListStdevResiduals, MaskPath,'AdjStdevRes_');
    CalcAdjDualRegAmps([dual_reg_dir,filesep,DualRegStage2ICFiles{IndICFile}],FileListrMSSDResiduals, MaskPath,'AdjrMSSDRes_');
end

%% save settings
save([dual_reg_dir,filesep,'CollectionForAdjustmentDualRegAmps_',datestr(now,'ddmmmyyyy_HHMM'),'.mat'],'FileListPath','FileList','StdevFullTS_FileList','NewFileListStdevFullTS','NewFileList_rMSSDFullTS','dual_reg_dir','MaskPath','FileListStdevResiduals','DualRegStage2ICFilePaths');

disp('DONE.');