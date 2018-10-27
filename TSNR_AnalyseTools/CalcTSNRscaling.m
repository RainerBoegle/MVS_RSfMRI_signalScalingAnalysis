function TSNRscaling = CalcTSNRscaling(TSNR_Design,MaskPath,SLight)
% This function calculates the scaling of TSNR values as the fraction of the TSNR values for the 2ndMRI to the 1stMRI,
% just like MVSscaling analysis of DualReg Amplitudes.
% in the first case only the raw data is used, -after averaging over the Runs per MRI,
% whereas a searchlight is used in the second case to average across space and Runs.
%
%Usage:
%       TSNRscaling = CalcTSNRscaling(TSNR_Design,MaskPath,SLight);
%
%
%V1.0
%Date: V1.0(18.09.2015): initial implementation based on test script.
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%% remember
% TSNRscaling.Design  = TSNR_Design;
% TSNRscaling.Design.filelist
% TSNRscaling.Design.SubjNrMRInrRunNr
% TSNRscaling.Design.UniqueSubjNrs
% TSNRscaling.Design.UniqueMRInrs
% TSNRscaling.Design.NSubjs
% TSNRscaling.Design.N_MRIs
%
% TSNRscaling.Mask.MaskNII_FilePath= MaskPath;
% ...
% TSNRscaling.SLight  = SLight;
% TSNRscaling.AmplitudesPerSubjectMRI = cell(2,1); %AmplitudesPerSubjectMRI{1or2} will be (NVoxelWholeBrain,NSubj,NMRI); -i.e. average over Runs OR in case of TSNR{2} nanmedian over space (searchlight) and Runs
% TSNRscaling.ScalingPerSubject       = cell(2,1); %ScalingPerSubject{1or2}       will be (NVoxelWholeBrain,NSubj);      -i.e. fraction of MRI==2 divided by MRI==1. 

%% create "design" from FileList
TSNRscaling.Design = TSNR_Design;

%% get mask
TSNRscaling.Mask.MaskNII_FilePath  = MaskPath;
TSNRscaling.Mask.NII               = nifti(TSNRscaling.Mask.MaskNII_FilePath);
TSNRscaling.Mask.Mask3D            = TSNRscaling.Mask.NII.dat(:,:,:)>0; %make sure it is 0,1 only
TSNRscaling.Mask.MaskLinearIndices = find(TSNRscaling.Mask.Mask3D(:)==1);
TSNRscaling.Mask.NVoxelsInMask     = length(TSNRscaling.Mask.MaskLinearIndices);
TSNRscaling.Mask.mat               = TSNRscaling.Mask.NII.mat;
TSNRscaling.Mask.dim               = TSNRscaling.Mask.NII.dat.dim;

%% get TSNR data
disp('Loading all TSNR files...');
TSNRscaling.InputData.NII      = cell(length(TSNRscaling.Design.filelist),1); %init
TSNRscaling.InputData.NII{1}   = nifti(TSNRscaling.Design.filelist{1});
TSNRscaling.InputData.mat      = TSNRscaling.InputData.NII{1}.mat;
TSNRscaling.InputData.dim      = TSNRscaling.InputData.NII{1}.dat.dim;
TSNRscaling.InputData.AllData4D= zeros([TSNRscaling.InputData.dim,length(TSNRscaling.Design.filelist)]);
for IndFile = 1:length(TSNRscaling.Design.filelist)
    TSNRscaling.InputData.NII{IndFile}             = nifti(TSNRscaling.Design.filelist{IndFile});
    TSNRscaling.InputData.AllData4D(:,:,:,IndFile) = TSNRscaling.InputData.NII{IndFile}.dat(:,:,:);
end
TSNRscaling.InputData.AllData2D      = reshape(TSNRscaling.InputData.AllData4D,[],size(TSNRscaling.InputData.AllData4D,4));

%% apply mask
TSNRscaling.InputData.DataInMask = TSNRscaling.InputData.AllData2D(TSNRscaling.Mask.MaskLinearIndices,:);

%% add searchlight definition
TSNRscaling.SLight = SLight;

%% init fields for amplitudes and scaling data
TSNRscaling.ApproachesInfo = {'VoxWise','i.e. raw TSNR values per Voxel are used for scaling value calculation.'; ...
                              'SLavg'  ,'i.e. searchlight average of TSNR values (spatial average using nanmedian) per subject and calculation of scaling values from this average per subject.'; ...
                             };
TSNRscaling.AmplitudesPerSubjectMRI = cell(size(TSNRscaling.ApproachesInfo,1),1);
TSNRscaling.ScalingPerSubject       = cell(size(TSNRscaling.ApproachesInfo,1),1);

%% do averaging of amplitudes over runs per subject and MRI, and then create Scaling data.
%% do this for all voxels and per searchlight in two runs. (this means the first is fast the second takes longer...)
for IndApproach = 1:size(TSNRscaling.ApproachesInfo,1)
    switch(TSNRscaling.ApproachesInfo{IndApproach,1})
        case {'VoxWise','SLavg'}
            TSNRscaling.AmplitudesPerSubjectMRI{IndApproach} = zeros(TSNRscaling.Mask.NVoxelsInMask,TSNRscaling.Design.NSubjs,TSNRscaling.Design.N_MRIs);
            TSNRscaling.ScalingPerSubject{      IndApproach} = NaN(  TSNRscaling.Mask.NVoxelsInMask,TSNRscaling.Design.NSubjs); 
            %NB: Keep the voxels relative to the total mask such that later processing does not have to be changed alot for update of method to V1.1.
            %    For "SLavg" use searchlight to do averaging and project results back to the center index of the searchlight (as should be) instead of using searchlight index for this.
            %    This might seem strange but will make further development and adjustment of functions easier while keeping the backward compatibility.
        otherwise
            error(['Approach "',TSNRscaling.ApproachesInfo{IndApproach,1},'" unknown!']);
    end
    
    disp( '------------------------------------------');
    disp( 'Creating TSNR-fMRI scaling data per subject ...');
    disp(['Approach: ',TSNRscaling.ApproachesInfo{IndApproach,1},': ',TSNRscaling.ApproachesInfo{IndApproach,2}]);
    if(strcmp(TSNRscaling.ApproachesInfo{IndApproach,1},'VoxWise'))
        fprintf('Subj ');
    end
    for IndSubj = 1:TSNRscaling.Design.NSubjs        
        switch(TSNRscaling.ApproachesInfo{IndApproach,1}) %this might make it slower than it could be, but probably doesn't matter much as processing of each subject is the bottle neck...
            case 'VoxWise'
                fprintf('%02d...',IndSubj);
                %average over runs per MRI
                for IndMRI = 1:TSNRscaling.Design.N_MRIs
                    MRIindex = find(TSNRscaling.Design.SubjNrMRInrRunNr(:,1)==TSNRscaling.Design.UniqueSubjNrs(IndSubj) & TSNRscaling.Design.SubjNrMRInrRunNr(:,2)==TSNRscaling.Design.UniqueMRInrs(IndMRI));
                    TSNRscaling.AmplitudesPerSubjectMRI{IndApproach}(:,IndSubj,IndMRI) = mean(TSNRscaling.InputData.DataInMask(:,MRIindex),2);
                    clear MRIindex
                end
            case 'SLavg'
                fprintf('Subj %02d:',IndSubj);
                %average over runs per MRI AND searchlight, i.e. remove mean over RUNs per MRI just select the data and then use this in the nanmedian for averaging
                for IndMRI = 1:TSNRscaling.Design.N_MRIs
                    fprintf('MRI %02d:',IndMRI);
                    MRIindex = find(TSNRscaling.Design.SubjNrMRInrRunNr(:,1)==TSNRscaling.Design.UniqueSubjNrs(IndSubj) & TSNRscaling.Design.SubjNrMRInrRunNr(:,2)==TSNRscaling.Design.UniqueMRInrs(IndMRI)); %this is actually the RUNs that are extracted here per subject and MRI/Scanner
                    %old: AllAmplitudes = mean(MVSscaling.InputData.DataInMask(:,MRIindex),2);
                    AllAmplitudes = TSNRscaling.InputData.DataInMask(:,MRIindex); %new give the whole data of both runs and later use SLaverage (nanmedian) to get the average over RUNs AND space (searchlight)
                    
                    NextPercentDone = 0; %init %for indicating the progress 1%-stepwise
                    reverseStr = '';
                    for IndSLCenterVox = 1:length(TSNRscaling.SLight.SLightIndsInMaskCell)
                        CurrInds = TSNRscaling.SLight.SLightIndsInMaskCell{IndSLCenterVox};
                        CenterInd= TSNRscaling.SLight.SLightIndsInMaskCell{IndSLCenterVox}(1);
                        CurrData = AllAmplitudes(CurrInds,:); %data in searchlight for all RUNs of this MRI/Scanner
                        TSNRscaling.AmplitudesPerSubjectMRI{IndApproach}(CenterInd,IndSubj,IndMRI) = nanmedian(CurrData(:)); %average over space and 
                        
                        %progress report
                        if(NextPercentDone==0)%start
                            %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
                            msg = sprintf('SLavg(Amps): %03.0f percent done.', floor(IndSLCenterVox*100/length(TSNRscaling.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
                            fprintf([reverseStr, msg]);
                            reverseStr = repmat(sprintf('\b'), 1, length(msg));
                            NextPercentDone = 1;
                        else
                            if(floor(IndSLCenterVox*100/length(TSNRscaling.SLight.SLightIndsInMaskCell))==NextPercentDone)
                                %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
                                msg = sprintf('SLavg(Amps): %03.0f percent done.', floor(IndSLCenterVox*100/length(TSNRscaling.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
                                fprintf([reverseStr, msg]);
                                reverseStr = repmat(sprintf('\b'), 1, length(msg));
                                NextPercentDone = NextPercentDone+1;
                            end
                        end
                        clear CurrInds CenterInd
                    end
                    clear MRIindex AllAmplitudes 
                end
        end
        %create scaling
        TSNRscaling.ScalingPerSubject{IndApproach}(:,IndSubj) = abs(TSNRscaling.AmplitudesPerSubjectMRI{IndApproach}(:,IndSubj,2)./TSNRscaling.AmplitudesPerSubjectMRI{IndApproach}(:,IndSubj,1));
        if(strcmp(TSNRscaling.ApproachesInfo{IndApproach,1},'SLavg'))
            fprintf('\n');
        end
    end
    fprintf('DONE.\n');
end
fprintf('ALL DONE.\n');

end