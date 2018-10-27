%% Get MVSscaling struct
ScalingAnaMatPath = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
load(ScalingAnaMatPath);

[BasePath, FName, ext] = fileparts(ScalingAnaMatPath);
MVSscalingTest.InputPath = ScalingAnaMatPath;

%% make magnitude scaling data
if(0)
    ScalingPerSubject = abs(MVSscaling.ScalingPerSubject);
else
    ScalingPerSubject =     MVSscaling.ScalingPerSubject;
end

%% SLight smooth?
%to do

%% get overlap mask
OverlapMask_Path = spm_select(1,'image','Select DMN and EffectsOfMRI overlap mask...');
V_OverlapMask = spm_vol(OverlapMask_Path);
OverlapMaskData = V_OverlapMask.private.dat(:);

if(length(OverlapMaskData)~=length(MVSscaling.Masks.WholeBrainRaw))
    error('Masks do not have the same dimension!');
end

%% select some of the data
Data{1} = ScalingPerSubject(OverlapMaskData(MVSscaling.Masks.WholeBrainRaw~=0)~=0,:); %Total "interesting" voxels
Data{2} = ScalingPerSubject((OverlapMaskData(MVSscaling.Masks.WholeBrainRaw~=0)>0)&(OverlapMaskData(MVSscaling.Masks.WholeBrainRaw~=0)<2),:); %Total "EffectsOfMRI" voxels
Data{3} = ScalingPerSubject((OverlapMaskData(MVSscaling.Masks.WholeBrainRaw~=0)>1)&(OverlapMaskData(MVSscaling.Masks.WholeBrainRaw~=0)<3),:); %Total "DMN-WO-EffectsOfMRI" voxels

%% check it out
figure(42); clf;
for IndSubj = 1:size(ScalingPerSubject,2)
    ax(1)=subplot(1,3,1); hist(Data{1}(abs(Data{1}(:,IndSubj))<6,IndSubj),[-6*sqrt(2):0.01:6*sqrt(2)]); title(['Distribution "interesting"         voxels Subj ',num2str(IndSubj)]); hold on;
    subplot(1,3,1); hist(sqrt(2).*ones(10,1),[-6*sqrt(2):0.01:6*sqrt(2)]); hold on
    subplot(1,3,1); hist(2*sqrt(2).*ones(15,1),[-6*sqrt(2):0.01:6*sqrt(2)]); hold on
    ax(2)=subplot(1,3,2); hist(Data{2}(abs(Data{2}(:,IndSubj))<6,IndSubj),[-6*sqrt(2):0.01:6*sqrt(2)]); title(['Distribution "EffectsOfMRI"        voxels Subj ',num2str(IndSubj)]); hold on
    subplot(1,3,2); hist(sqrt(2).*ones(10,1),[-6*sqrt(2):0.01:6*sqrt(2)]); hold on
    subplot(1,3,2); hist(2*sqrt(2).*ones(15,1),[-6*sqrt(2):0.01:6*sqrt(2)]); hold on
    ax(3)=subplot(1,3,3); hist(Data{3}(abs(Data{3}(:,IndSubj))<6,IndSubj),[-6*sqrt(2):0.01:6*sqrt(2)]); title(['Distribution "DMN-WO-EffectsOfMRI" voxels Subj ',num2str(IndSubj)]); hold on
    subplot(1,3,3); hist(sqrt(2).*ones(10,1),[-6*sqrt(2):0.01:6*sqrt(2)]); hold on
    subplot(1,3,3); hist(2*sqrt(2).*ones(15,1),[-6*sqrt(2):0.01:6*sqrt(2)]); hold on
    linkaxes([ax(1),ax(2),ax(3)],'x');
    h = helpdlg(['Subj ',num2str(IndSubj)],'wait?');
    uiwait(h);
end