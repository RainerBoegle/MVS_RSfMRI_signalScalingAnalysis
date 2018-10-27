%% This is just a test script to check out some TSNR and rMSSD data
%% Future idea regarding MVS detection in fMRI data: Compare various summary stats like TSNR or Stdev over time or rMSSD for Gray Matter and White Matter, for Certain Areas that are MVS areas and for Areas that are not, maybe resting state networks as masks???
%% best would be just to take spheres at certain coordinates and sample data this way (ULTRA SIMPLE)
%% But the main thing probably has to be what is going on in certain MVS areas vs white matter compared to R-Networks/ROIs not MVS vs white matter...

DataDir = spm_select(1,'dir','Select directory with TSNR and rMSSD files...');

TSNR_files        = cellstr(spm_select('List',DataDir,'^TSNR_Subj_p.*.nii'));  %let's hope the order is always the same
Mean_files        = cellstr(spm_select('List',DataDir,'^Mean_Subj_p.*.nii'));  %let's hope the order is always the same
Stdev_files       = cellstr(spm_select('List',DataDir,'^Stdev_Subj_p.*.nii')); %let's hope the order is always the same
rMSSD_files       = cellstr(spm_select('List',DataDir,'^rMSSD_Subj_p.*.nii')); %let's hope the order is always the same
rMSSD2_files      = cellstr(spm_select('List',DataDir,'^rMSSD2nd_Subj_p.*.nii')); %let's hope the order is always the same
StdevMSSD_files   = cellstr(spm_select('List',DataDir,'^StdevMSSD_Subj_p.*.nii')); %let's hope the order is always the same
StdevMSSD2nd_files= cellstr(spm_select('List',DataDir,'^StdevMSSD2nd_Subj_p.*.nii')); %let's hope the order is always the same

MaskFilePath = spm_select('List',DataDir,'^mask.nii');

%% Get Subjects from data (will only work on current data)
SubjectsMRIsRunsCellStr = cellfun(@(x)x,regexp(TSNR_files,'Subj_p\d\d\d_MRI\d_rsfMRI_\d','match'));

[DataIdx,ok] = listdlg('ListString',SubjectsMRIsRunsCellStr,'SelectionMode','single','InitialValue',5,'Name','Select Data','PromptString','Pick one');
if(~ok)
    return;
end

%% load data
NII_mask         = nifti([DataDir,filesep,MaskFilePath]);
NII_TSNR         = nifti([DataDir,filesep,TSNR_files{        DataIdx}]);
NII_Mean         = nifti([DataDir,filesep,Mean_files{        DataIdx}]);
NII_Stdev        = nifti([DataDir,filesep,Stdev_files{       DataIdx}]);
NII_rMSSD        = nifti([DataDir,filesep,rMSSD_files{       DataIdx}]);
NII_rMSSD2       = nifti([DataDir,filesep,rMSSD2_files{      DataIdx}]);
NII_StdevMSSD    = nifti([DataDir,filesep,StdevMSSD_files{   DataIdx}]);
NII_StdevMSSD2nd = nifti([DataDir,filesep,StdevMSSD2nd_files{DataIdx}]);

Mask_data         = NII_mask.dat(:);
TSNR_data         = NII_TSNR.dat(:);         TSNR_data         = TSNR_data(        Mask_data~=0);
Mean_data         = NII_Mean.dat(:);         Mean_data         = Mean_data(        Mask_data~=0);
Stdev_data        = NII_Stdev.dat(:);        Stdev_data        = Stdev_data(       Mask_data~=0);
rMSSD_data        = NII_rMSSD.dat(:);        rMSSD_data        = rMSSD_data(       Mask_data~=0);
rMSSD2_data       = NII_rMSSD2.dat(:);       rMSSD2_data       = rMSSD2_data(      Mask_data~=0);
StdevMSSD_data    = NII_StdevMSSD.dat(:);    StdevMSSD_data    = StdevMSSD_data(   Mask_data~=0);
StdevMSSD2nd_data = NII_StdevMSSD2nd.dat(:); StdevMSSD2nd_data = StdevMSSD2nd_data(Mask_data~=0);

DataTypesCellStr  = {'TSNR_data'; 'Mean_data'; 'Stdev_data'; 'rMSSD_data'; 'rMSSD2_data'; 'StdevMSSD_data'; 'StdevMSSD2nd_data'};

%% ask user to select xAxis Data
[xAxisIdx,ok] = listdlg('ListString',DataTypesCellStr,'SelectionMode','single','PromptString','Pick a x-Axis');
if(~ok)
    return;
end

[yAxisIdx,ok] = listdlg('ListString',DataTypesCellStr,'SelectionMode','multiple','PromptString','Pick one or several y-Axis');
if(~ok)
    return;
end

%% pick data for plot
xAxisData = eval(DataTypesCellStr{xAxisIdx});

yAxisData = cell(length(yAxisIdx),1);
for IndData = 1:length(yAxisIdx)
    yAxisData{IndData} = eval(DataTypesCellStr{yAxisIdx(IndData)});
end

%% plot together
Cols = distinguishable_colors(length(yAxisIdx),{'w','k','y','c'});
MarkerStr = {'+';'o';'*';'.';'x';'s';'d';'^';'v';'>';'<';'p';'h'};

FigNumBase = 420;
TitleStr = cell(2,1);
figure(FigNumBase); clf; hold on
for IndData = 1:length(yAxisIdx)
    plot(xAxisData,yAxisData{IndData},MarkerStr{IndData},'Color',Cols(IndData,:)); hold on
    if(IndData==1)
        TitleStr{1} = DataTypesCellStr{yAxisIdx(IndData)};
    else
        TitleStr{1} = [TitleStr{1},' & ',DataTypesCellStr{yAxisIdx(IndData)}];
    end
end
TitleStr{2} = ['over ',DataTypesCellStr{xAxisIdx}];
title(TitleStr,'Interpreter','none');
xlabel(DataTypesCellStr{xAxisIdx},'Interpreter','none');
legend(DataTypesCellStr(yAxisIdx),'Interpreter','none');

%% plot each one alone.
for IndData = 1:length(yAxisIdx)
    figure(FigNumBase+IndData); clf;  
    plot(xAxisData,yAxisData{IndData},MarkerStr{IndData},'Color',Cols(IndData,:)); hold on
    title([DataTypesCellStr{yAxisIdx(IndData)},' over ',DataTypesCellStr{xAxisIdx}],'Interpreter','none');
    xlabel(DataTypesCellStr{xAxisIdx},'Interpreter','none');
    ylabel(DataTypesCellStr{yAxisIdx(IndData)},'Interpreter','none');
end

%% plot special interest: rMSSD/sqrt(2)*Stdev
figure(811); clf;
plot(rMSSD_data,sqrt(2).*Stdev_data,'kx'); title('sqrt(2).*Stdev_data over rMSSD_data','Interpreter','none');

figure(812); clf;
plot(Mean_data,rMSSD_data./Stdev_data,'kx'); title('rMSSD_data./Stdev_data over Mean_data','Interpreter','none');

%% plot special interest: Mean/rMSSD over TSNR (is different estimate of variation over time making a difference in TSNR estimate from rMSSD?)
figure(81); clf;
plot(TSNR_data,Mean_data./rMSSD_data,'kx'); title('Mean_data./rMSSD_data over TSNR_data(=Mean_data./Stdev_data)','Interpreter','none');


%% boxplot of TSNR per for selected data
[BoxPlotIdx,ok] = listdlg('ListString',SubjectsMRIsRunsCellStr,'SelectionMode','multiple','InitialValue',5,'Name','Select Data','PromptString','Pick some data for boxplot of TSNR');
if(~ok)
    return;
end

BPdata = cell(length(BoxPlotIdx),1);
Labels = cell(length(BoxPlotIdx),1);
for IndData = 1:length(BoxPlotIdx)
    NII_TSNR  = nifti([DataDir,filesep,TSNR_files{BoxPlotIdx(IndData)}]);
    TSNR_data = NII_TSNR.dat(:);
    TSNR_data = TSNR_data(Mask_data~=0);
    
    BPdata{IndData} = TSNR_data;
    Labels{IndData} = SubjectsMRIsRunsCellStr{BoxPlotIdx(IndData)};
end

BPdataVec = zeros(length(BoxPlotIdx)*length(TSNR_data),1);
Grouping  = zeros(length(BoxPlotIdx)*length(TSNR_data),1);
for IndData = 1:length(BoxPlotIdx)
    BPdataVec((IndData-1)*length(TSNR_data)+(1:length(TSNR_data))) = BPdata{IndData};
    Grouping( (IndData-1)*length(TSNR_data)+(1:length(TSNR_data))) = IndData.*ones(size(BPdata{IndData}));
end

figure(42); clf;
boxplot(BPdataVec,Grouping,'notch','on','labels',Labels); title('TSNR')