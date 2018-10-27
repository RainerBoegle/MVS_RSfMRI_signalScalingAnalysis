%% get MVSscaling data
MVSscalingDataPath = spm_select(1,'MVSscaling_','Select MVSscaling (fMRI-)data...');
load(MVSscalingDataPath);

%% get Searchlight definition
SLdefFilePath = spm_select(1,'SLdef_','Select Searchlight definition *.mat file...');
load(SLdefFilePath);

%% select output base-directory
OutputBaseDir = spm_select(1,'dir','Select Base-Output Directory...');

%% Select coordinates of searchlight centers to be extracted
Def_Org     = {'[2,38,24;2,42,24;-2,42,24;2,38,20;-2,38,24;2,42,20;-2,42,20;-2,38,20;-2,22,24;-2,24,20;-2,26,20;-2,28,20;+2,22,24;+2,24,20;+2,26,20;+2,28,20;2,30,-8;2,32,-8;42,-6,-8;-42,-6,-8;-40,-18,0;+40,-18,0;-2,-54,-20;+2,-54,-20;-2,-48,-20;+2,-48,-20;-2,-48,-16;+2,-48,-16;-2,-62,-36;+2,-62,-36;-16,-34,-36;+16,-34,-36;-10,-74,16;-10,-74,18;-10,-74,20;-12,-70,12;46,-30,18;-20,-62,-36;2,-54,-18;-2,-54,-18;2,-54,-14;-2,-54,-14;2,-54,-12;-2,-54,-12;0,-22,-45;0,-62,-42;-2,-62,-42;2,-62,-42;-12,-72,10;-12,-72,12;-12,-72,14;-12,-72,16;-2,27,16;-2,22,16;-2,22,18;-2,22,20;-2,22,22;-2,22,24]'};
Def_CingAnt1 = {'[-2,42,24;2,42,24;2,38,24;0,40,24;-2,42,20;2,42,20;2,38,20;0,40,20]'};
Def_CingAnt1b= {'[2,38,20;2,38,24;2,42,24;-2,42,24]'};
Def_CingAnt1c= {'[2,38,20;2,38,24;-2,42,24]'};
Def_CingAnt1d= {'[2,38,20;-2,42,24]'};
Def_CingAnt2 = {'[-2,24,20;-2,26,20;-2,22,20;-2,22,18;-2,24,24;-2,26,24;-2,22,24]'};
Def_CingAnt2b= {'[-2,24,20;-2,26,20;-2,22,20;-2,22,18]'};
% Def_CingAntSelect= {'[-2,24,20;-2,26,20;-2,22,20;-2,22,18;2,38,20]'};
% Def_CingAntSelect2= {'[-2,26,20;-2,22,20;2,34,20]'};
Def_CingAntSelect2b= {'[-2,26,20;2,34,20]'};
 
Def_Vermis1Full= {'[2,-54,-12;-2,-54,-12;2,-54,-14;-2,-54,-14;-2,-54,-18;-2,-54,-20;2,-50,-12;-2,-50,-12;2,-50,-14;-2,-50,-14;-2,-50,-18;-2,-50,-20;2,-58,-12;-2,-58,-12;2,-58,-14;-2,-58,-14;-2,-58,-18;-2,-58,-20]'};
% Def_Vermis1Select={'[-2,-50,-14;-2,-50,-12;-2,-54,-18]'};
Def_Vermis1Select2={'[-2,-50,-12;-2,-54,-16]'};

Def_Vermis2Full= {'[-2,-62,-42;0,-62,-42;2,-62,-42;-2,-62,-36;0,-62,-36;2,-62,-36;-2,-62,-46;0,-62,-46;2,-62,-46;-2,-62,-38;0,-62,-38;2,-62,-38]'};
Def_Vermis2Select= {'[-2,-62,-42;-2,-62,-36]'};

Def_CalcarineAll = {'[-10,-74,12;-14,-74,12; -10,-74,16;-14,-74,16;  -10,-74,20]'};
Def_CalcarineSelect = {'[-10,-74,16;-14,-74,16]'};

%                   CingAntSelect2b      Vermis1Select2         Vermis2Select          CalcarineSelect
Def_AllSelect = {'[-2,26,20;2,34,20; -2,-50,-12;-2,-54,-16; -2,-62,-42;-2,-62,-36; -10,-74,16;-14,-74,16]'};
answer_SLcentersMM = inputdlg({'Enter mm-coordinates of searchlight(s) to extract: '},'Which searchlights?',1,Def_AllSelect); %combine index?

CoordsToSearch = eval(answer_SLcentersMM{1});
invMat = inv(SLight.V_SLmask.mat); %World2Voxel-Matrix

VoxelsToSearch = zeros(size(CoordsToSearch));
for IndVox = 1:size(CoordsToSearch,1)
    VoxelsToSearch(IndVox,:) = round(CoordsToSearch(IndVox,:)*invMat(1:3,1:3)+invMat(1:3,4)');
end

%% mark selected voxels in a mask
SelVoxMask = zeros(SLight.V_SLmask.dim);
for IndVox = 1:size(VoxelsToSearch,1)
    SelVoxMask(VoxelsToSearch(IndVox,1),VoxelsToSearch(IndVox,2),VoxelsToSearch(IndVox,3)) = 1;
end

Vout = SLight.V_SLmask;
if(Vout.dt(1)<16)
    Vout.dt(1) = 16; %not necessary but save
end
StartDateTime = datestr(now,'yymmmdd_HHMM');
Vout.fname    = [OutputBaseDir,filesep,'SelectedVoxelsMask_',StartDateTime,'.nii'];

Vout = spm_write_vol(Vout,SelVoxMask);

%% mark each individual and all combined searchlights
%NB: if this mask creation works then we can be sure to extract scaling data using SL_LinIndsVoxelsToSearch in the same way...
LinIndVoxelsToSearch = zeros(size(VoxelsToSearch,1),1);
IndicesRelativeToWholeBrainMask = cell(size(VoxelsToSearch,1),1);
SL_LinIndsVoxelsToSearch        = cell(size(VoxelsToSearch,1),1);
AggregateScalingdataPerSLnSubj  = cell(size(VoxelsToSearch,1),1);
AllAggregateScalingDataPerSubj  = []; %to be filled...
CollectedIndicesRelativeToWholeBrainMask = []; %to be filled...
AllData = zeros(Vout.dim); %prep empty set in size of mask
for IndVox = 1:size(VoxelsToSearch,1)
    LinIndVoxelsToSearch(IndVox)     = sub2ind(SLight.V_SLmask.dim,VoxelsToSearch(IndVox,1),VoxelsToSearch(IndVox,2),VoxelsToSearch(IndVox,3)); %find center voxel coord in linear index
    IndicesRelativeToWholeBrainMask{IndVox}  = SLight.SLightIndsInMaskCell{find(SLight.StartIndsVol==LinIndVoxelsToSearch(IndVox))}; %these are the ones that we need to get the scaling data our of MVSscaling.ScalingPerSubject
    CollectedIndicesRelativeToWholeBrainMask = [CollectedIndicesRelativeToWholeBrainMask; IndicesRelativeToWholeBrainMask{IndVox}(:)];
    SL_LinIndsVoxelsToSearch{IndVox} = SLight.AllowedIndsVol(IndicesRelativeToWholeBrainMask{IndVox}); %these are the whole volume indices to create a mask.
    
    AggregateScalingdataPerSLnSubj{IndVox} = abs(MVSscaling.ScalingPerSubject(IndicesRelativeToWholeBrainMask{IndVox},:));
    
    Data = zeros(Vout.dim); %create empty data in size of the whole volume
    Data(SL_LinIndsVoxelsToSearch{IndVox}) = 1; %assign data/mask at the corresponding place in the mask
    Vout.fname = [OutputBaseDir,filesep,'SelectedVoxelsMask_',StartDateTime,'_SLofVox',num2str(IndVox),'.nii'];
    
    Vout = spm_write_vol(Vout,Data); %write out.
    AllData = AllData+Data; %add to all data
end
Vout.fname = [OutputBaseDir,filesep,'SelectedVoxelsMask_',StartDateTime,'_CombinedSLofAllVoxels.nii'];
Vout = spm_write_vol(Vout,AllData); %write out.

CollectedIndicesRelativeToWholeBrainMask = unique(CollectedIndicesRelativeToWholeBrainMask); %only keep each one without double occurrence
AllAggregateScalingDataPerSubj = abs(MVSscaling.ScalingPerSubject(CollectedIndicesRelativeToWholeBrainMask,:));

%% display the data per voxel/SL and then everything combined
UseYLim = 1;
UseSpecialIndicators = 1; %use 80% & 95% lines to show where the data lies
if(UseSpecialIndicators)
    YLimits = [0 ceil(quantile(AllAggregateScalingDataPerSubj(:),0.95))+1]; %[0 12];
else
    YLimits = [0 12];
end

figure(10); clf;
boxplot([nanmedian(AllAggregateScalingDataPerSubj,2);AllAggregateScalingDataPerSubj(:);nanmedian(AllAggregateScalingDataPerSubj,1)'],[ones(length(nanmedian(AllAggregateScalingDataPerSubj,2)),1);2.*ones(length(AllAggregateScalingDataPerSubj(:)),1);3.*ones(length(nanmedian(AllAggregateScalingDataPerSubj,1)'),1)],'labels',{'MedS';'Aggregate';'MedROI'},'notch','on'); title('All ScalingData of all SLs combined.'); hold on;
plot(1:3,(2*sqrt(2)).*ones(3,1),'r-'); hold on;
plot(1:3,sqrt(2).*ones(3,1),'r--'); hold on;
if(UseSpecialIndicators)
    plot(0.5:1.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,2),0.8).*ones(2,1),'g--'); hold on; %80% of the MedianSubject (VoxelDistribution) data lie below this line
    plot(0.5:1.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,2),0.95).*ones(2,1),'g-'); hold on; %95% of the MedianSubject (VoxelDistribution) data lie below this line
    plot(0.5:1.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,2),0.2).*ones(2,1),'g--'); hold on; %80% of the MedianSubject (VoxelDistribution) data lie below this line
    plot(0.5:1.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,2),0.05).*ones(2,1),'g-'); hold on; %95% of the MedianSubject (VoxelDistribution) data lie below this line

    plot(1.5:2.5,quantile(AllAggregateScalingDataPerSubj(:),0.8).*ones(2,1),'b--'); hold on; %80% of the aggregate data lie below this line
    plot(1.5:2.5,quantile(AllAggregateScalingDataPerSubj(:),0.95).*ones(2,1),'b-'); hold on; %95% of the aggregate data lie below this line
    plot(1.5:2.5,quantile(AllAggregateScalingDataPerSubj(:),0.2).*ones(2,1),'b--'); hold on; %80% of the aggregate data lie below this line
    plot(1.5:2.5,quantile(AllAggregateScalingDataPerSubj(:),0.05).*ones(2,1),'b-'); hold on; %95% of the aggregate data lie below this line

    plot(2.5:3.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,1)',0.8).*ones(2,1),'m--'); hold on; %80% of the MedianROI (subject distribution) data lie below this line
    plot(2.5:3.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,1)',0.95).*ones(2,1),'m-'); hold on; %95% of the MedianROI (subject distribution) data lie below this line
    plot(2.5:3.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,1)',0.2).*ones(2,1),'m--'); hold on; %80% of the MedianROI (subject distribution) data lie below this line
    plot(2.5:3.5,quantile(nanmedian(AllAggregateScalingDataPerSubj,1)',0.05).*ones(2,1),'m-'); hold on; %95% of the MedianROI (subject distribution) data lie below this line
end

if(UseYLim)
    ylim(YLimits);
end

for IndSL = 1:length(AggregateScalingdataPerSLnSubj)
    CurrData = AggregateScalingdataPerSLnSubj{IndSL};
    figure(10+IndSL); clf;
    boxplot([nanmedian(CurrData,2);CurrData(:);nanmedian(CurrData,1)'],[ones(length(nanmedian(CurrData,2)),1);2.*ones(length(CurrData(:)),1);3.*ones(length(nanmedian(CurrData,1)'),1)],'labels',{'MedS';'Aggregate';'MedSL'},'notch','on'); title(['All ScalingData of SL ',num2str(IndSL),' Center==(',num2str(CoordsToSearch(IndSL,1)),',',num2str(CoordsToSearch(IndSL,2)),',',num2str(CoordsToSearch(IndSL,3)),').']); hold on;
    plot(1:3,(2*sqrt(2)).*ones(3,1),'r-'); hold on;
    plot(1:3,sqrt(2).*ones(3,1),'r--'); hold on;
    if(UseSpecialIndicators)
        plot(0.5:1.5,quantile(nanmedian(CurrData,2),0.8).*ones(2,1),'g--'); hold on; %80% of the MedianSubject (VoxelDistribution) data lie below this line
        plot(0.5:1.5,quantile(nanmedian(CurrData,2),0.95).*ones(2,1),'g-'); hold on; %95% of the MedianSubject (VoxelDistribution) data lie below this line
        plot(0.5:1.5,quantile(nanmedian(CurrData,2),0.2).*ones(2,1),'g--'); hold on; %80% of the MedianSubject (VoxelDistribution) data lie below this line
        plot(0.5:1.5,quantile(nanmedian(CurrData,2),0.05).*ones(2,1),'g-'); hold on; %95% of the MedianSubject (VoxelDistribution) data lie below this line
        
        plot(1.5:2.5,quantile(CurrData(:),0.8).*ones(2,1),'b--'); hold on; %80% of the aggregate data lie below this line
        plot(1.5:2.5,quantile(CurrData(:),0.95).*ones(2,1),'b-'); hold on; %95% of the aggregate data lie below this line
        plot(1.5:2.5,quantile(CurrData(:),0.2).*ones(2,1),'b--'); hold on; %80% of the aggregate data lie below this line
        plot(1.5:2.5,quantile(CurrData(:),0.05).*ones(2,1),'b-'); hold on; %95% of the aggregate data lie below this line
    
        plot(2.5:3.5,quantile(nanmedian(CurrData,1)',0.8).*ones(2,1),'m--'); hold on; %80% of the MedianROI data lie below this line
        plot(2.5:3.5,quantile(nanmedian(CurrData,1)',0.95).*ones(2,1),'m-'); hold on; %95% of the MedianROI data lie below this line
        plot(2.5:3.5,quantile(nanmedian(CurrData,1)',0.2).*ones(2,1),'m--'); hold on; %80% of the MedianROI data lie below this line
        plot(2.5:3.5,quantile(nanmedian(CurrData,1)',0.05).*ones(2,1),'m-'); hold on; %95% of the MedianROI data lie below this line
    end
    
    if(UseYLim)
        ylim(YLimits);
    end
end

figure(20); clf;
subplot(1,3,1); hist(nanmedian(AllAggregateScalingDataPerSubj,2),[0:0.13:max(nanmedian(AllAggregateScalingDataPerSubj,2))]); title('histogram of MedianSubject (VoxelDistribution) data of complete ROI');
subplot(1,3,2); hist(AllAggregateScalingDataPerSubj(:),[0:0.01:ceil(quantile(AllAggregateScalingDataPerSubj(:),0.8)),(ceil(quantile(AllAggregateScalingDataPerSubj(:),0.8))+0.1):0.1:(ceil(quantile(AllAggregateScalingDataPerSubj(:),0.95))-0.1),ceil(quantile(AllAggregateScalingDataPerSubj(:),0.95)):(max(AllAggregateScalingDataPerSubj(:))-ceil(quantile(AllAggregateScalingDataPerSubj(:),0.95)))/20:max(AllAggregateScalingDataPerSubj(:))]); title('my special histogram of aggregate data of complete ROI');
subplot(1,3,3); hist(nanmedian(AllAggregateScalingDataPerSubj,1)',[0:0.2:max(nanmedian(AllAggregateScalingDataPerSubj,1)')]); title('histogram of MedianROI (subject distribution) data of complete ROI');