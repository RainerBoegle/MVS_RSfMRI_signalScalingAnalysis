function MVSscaling = CreateAmplitudeScalingData(DualRegStage2_FilePath,Mask_FilePath,FileListInfo,SLight)
% This function extracts the dual regression stage 2 amplitude/weights from a 4D file using a mask
% and then averages the amplitudes over runs per subject and MRI using the information extracted from the FileList.
% These average amplitudes are then used to calculate the scaling values Lambda.
% (This is as described in Boegle et al., 2015.)
%NB:
% Since V1.1 the above method for calculating the amplitudes and scaling values is the first entry
% in the fields .AmplitudesPerSubjectMRI & .ScalingPerSubject, (e.g. .ScalingPerSubject{1}(IndVoxel,IndSubj))
% that are now cells of size=2-x-1.
% The second entry (e.g. .ScalingPerSubject{2}(IndCenterVox(IndSL),IndSubj)) contains the amplitude values per searchlight that were created from averaging
% (using nanmedian) the amplitudes values in the respective searchlight per subject to get the amplitude values for calculating the scaling values.
% However, the scaling value calculation is still the same, i.e. the fraction of the (spatially searchlight averaged) amplitude values.
% Otherwise if in both calculations searchlight are used this would lead to a double smoothing that is problematic in regard to localization.
% This is an attempt to reduce noise in the amplitudes and thus scaling values via spatial average (nanmedian).
%
%
%
%USAGE:
%       MVSscaling = CreateAmplitudeScalingData(DualRegStage2_FilePath,Mask_FilePath,FileListInfo,SLight);
%
%
%V1.1
%Date: V1.1(03&04.09.2015): adjustments for allowing searchlight averaging of amplitudes before scaling analysis in addition to original method. + extra careful check if searchlight handling gives the right indices, -so far it seems right and this is the 4th check that I did, -so I should start trusting it. V1.0(22.July.2015): based on earlier script, here we skip the ROI analysis part, because it gets replaced by searchlights later.
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)

%% get IC number from DualRegStage2_FilePath if the form of the file is as expected, give a warning otherwise.
[BaseDir,FName,ext] = fileparts(DualRegStage2_FilePath);
AsExpected  = 0; %init expected form as NOT AS EXPECTED
DispWarning = 0; %init warning as not necessary
if(~strcmp(ext,'.nii'))
    DispWarning = 1;
else
    AsExpected = AsExpected+1; %extension is correct
end
if(~strcmp(cellfun(@(x)x,regexp(regexp(FName,'dr_stage2_ic\d\d\d\d','match'),'dr_stage2_ic','match')),'dr_stage2_ic'))
    DispWarning = 1;
else
    AsExpected = AsExpected+1; %beginning of filename is correct
end
if(isempty(cellfun(@(x)x,regexp(regexp(FName,'dr_stage2_ic\d\d\d\d','match'),'\d\d\d\d','match'))))
    DispWarning = 1;
else
    ICnumExtract = cellfun(@(x)x,regexp(regexp(FName,'dr_stage2_ic\d\d\d\d','match'),'\d\d\d\d','match'));
    if(~ischar(ICnumExtract))
        if(iscell(ICnumExtract))
            ICnumExtract = ICnumExtract{1};
            ICnumOrg     = str2num(ICnumExtract);
            AsExpected   = AsExpected+1; %numbering of ICs is (probably) as expected
        end
    else
        ICnumOrg   = str2num(ICnumExtract);
        AsExpected = AsExpected+1; %numbering of ICs is (probably) as expected
    end
end
    
if(DispWarning)
    warning('Dual regression Stage2 file is expected to be a NII-file with filename "dr_stage2_ic####.nii"')
end

%% assign ICnum
if(AsExpected==3)
    MVSscaling.ICnum = ICnumOrg+1;
else
    %This will not produce errors in computation but outputs to NIFTI will be stopped and user has
    %to do this by hand, such that correct outputs are (hopefully) produced by the user in that case.
    MVSscaling.ICnum = 'IC number could not be determined from input DRstage2-file, check "MVSscaling.InputData.InputsNII_Path" to see what was input.';
    warning(MVSscaling.ICnum);
end

%% create "design" from FileList
MVSscaling.Design = FileListInfo;

%% get mask
MVSscaling.Mask.MaskNII_FilePath  = Mask_FilePath;
MVSscaling.Mask.NII               = nifti(MVSscaling.Mask.MaskNII_FilePath);
MVSscaling.Mask.Mask3D            = MVSscaling.Mask.NII.dat(:,:,:)>0; %make sure it is 0,1 only
MVSscaling.Mask.MaskLinearIndices = find(MVSscaling.Mask.Mask3D(:)==1);
MVSscaling.Mask.NVoxelsInMask     = length(MVSscaling.Mask.MaskLinearIndices);
MVSscaling.Mask.mat               = MVSscaling.Mask.NII.mat;
MVSscaling.Mask.dim               = MVSscaling.Mask.NII.dat.dim;

%% get dualreg weights
MVSscaling.InputData.InputsNII_Path = DualRegStage2_FilePath;
MVSscaling.InputData.NII            = nifti(MVSscaling.InputData.InputsNII_Path);
MVSscaling.InputData.AllData4D      = MVSscaling.InputData.NII.dat(:,:,:,:);
MVSscaling.InputData.AllData2D      = reshape(MVSscaling.InputData.AllData4D,[],size(MVSscaling.InputData.AllData4D,4));
MVSscaling.InputData.mat            = MVSscaling.InputData.NII.mat;
MVSscaling.InputData.dim            = MVSscaling.InputData.NII.dat.dim;

%% apply mask
MVSscaling.InputData.DataInMask = MVSscaling.InputData.AllData2D(MVSscaling.Mask.MaskLinearIndices,:);

%% add searchlight definition
MVSscaling.SLight = SLight;

%% init fields for amplitudes and scaling data
MVSscaling.ApproachesInfo = {'VoxWise','i.e. raw amplitudes per Voxel are used for scaling value calculation.'; ...
                             'SLavg'  ,'i.e. searchlight average of amplitudes (spatial average using nanmedian) per subject and calculation of scaling values from this average per subject.'; ...
                            };
MVSscaling.AmplitudesPerSubjectMRI = cell(size(MVSscaling.ApproachesInfo,1),1);
MVSscaling.ScalingPerSubject       = cell(size(MVSscaling.ApproachesInfo,1),1);

%% do averaging of amplitudes over runs per subject and MRI, and then create Scaling data.
%% do this for all voxels and per searchlight in two runs. (this means the first is fast the second takes longer...)
for IndApproach = 1:size(MVSscaling.ApproachesInfo,1)
    switch(MVSscaling.ApproachesInfo{IndApproach,1})
        case {'VoxWise','SLavg'}
            MVSscaling.AmplitudesPerSubjectMRI{IndApproach} = zeros(MVSscaling.Mask.NVoxelsInMask,MVSscaling.Design.NSubjs,MVSscaling.Design.N_MRIs);
            MVSscaling.ScalingPerSubject{      IndApproach} = NaN(  MVSscaling.Mask.NVoxelsInMask,MVSscaling.Design.NSubjs); 
            %NB: Keep the voxels relative to the total mask such that later processing does not have to be changed alot for update of method to V1.1.
            %    For "SLavg" use searchlight to do averaging and project results back to the center index of the searchlight (as should be) instead of using searchlight index for this.
            %    This might seem strange but will make further development and adjustment of functions easier while keeping the backward compatibility.
        otherwise
            error(['Approach "',MVSscaling.ApproachesInfo{IndApproach,1},'" unknown!']);
    end
    
    disp( '------------------------------------------');
    disp( 'Creating fMRI scaling data per subject ...');
    disp(['Approach: ',MVSscaling.ApproachesInfo{IndApproach,1},': ',MVSscaling.ApproachesInfo{IndApproach,2}]);
    if(strcmp(MVSscaling.ApproachesInfo{IndApproach,1},'VoxWise'))
        fprintf('Subj ');
    end
    for IndSubj = 1:MVSscaling.Design.NSubjs        
        switch(MVSscaling.ApproachesInfo{IndApproach,1}) %this might make it slower than it could be, but probably doesn't matter much as processing of each subject is the bottle neck...
            case 'VoxWise'
                fprintf('%02d...',IndSubj);
                %average over runs per MRI
                for IndMRI = 1:MVSscaling.Design.N_MRIs
                    MRIindex = find(MVSscaling.Design.SubjNrScannerNrRunNr(:,1)==MVSscaling.Design.UniqueSubjNrs(IndSubj) & MVSscaling.Design.SubjNrScannerNrRunNr(:,2)==MVSscaling.Design.UniqueScannerNrs(IndMRI));
                    MVSscaling.AmplitudesPerSubjectMRI{IndApproach}(:,IndSubj,IndMRI) = mean(MVSscaling.InputData.DataInMask(:,MRIindex),2);
                    clear MRIindex
                end
            case 'SLavg'
                fprintf('Subj %02d:',IndSubj);
                %average over runs per MRI AND searchlight, i.e. remove mean over RUNs per MRI just select the data and then use this in the nanmedian for averaging
                for IndMRI = 1:MVSscaling.Design.N_MRIs
                    fprintf('MRI %02d:',IndMRI);
                    MRIindex = find(MVSscaling.Design.SubjNrScannerNrRunNr(:,1)==MVSscaling.Design.UniqueSubjNrs(IndSubj) & MVSscaling.Design.SubjNrScannerNrRunNr(:,2)==MVSscaling.Design.UniqueScannerNrs(IndMRI)); %this is actually the RUNs that are extracted here per subject and MRI/Scanner
                    %old: AllAmplitudes = mean(MVSscaling.InputData.DataInMask(:,MRIindex),2);
                    AllAmplitudes = MVSscaling.InputData.DataInMask(:,MRIindex); %new give the whole data of both runs and later use SLaverage (nanmedian) to get the average over RUNs AND space (searchlight)
                    
                    NextPercentDone = 0; %init %for indicating the progress 1%-stepwise
                    reverseStr = '';
                    for IndSLCenterVox = 1:length(MVSscaling.SLight.SLightIndsInMaskCell)
                        CurrInds = MVSscaling.SLight.SLightIndsInMaskCell{IndSLCenterVox};
                        CenterInd= MVSscaling.SLight.SLightIndsInMaskCell{IndSLCenterVox}(1);
                        CurrData = AllAmplitudes(CurrInds,:); %data in searchlight for all RUNs of this MRI/Scanner
                        MVSscaling.AmplitudesPerSubjectMRI{IndApproach}(CenterInd,IndSubj,IndMRI) = nanmedian(CurrData(:)); %average over space and 
                        
                        %progress report
                        if(NextPercentDone==0)%start
                            %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
                            msg = sprintf('SLavg(Amps): %03.0f percent done.', floor(IndSLCenterVox*100/length(MVSscaling.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
                            fprintf([reverseStr, msg]);
                            reverseStr = repmat(sprintf('\b'), 1, length(msg));
                            NextPercentDone = 1;
                        else
                            if(floor(IndSLCenterVox*100/length(MVSscaling.SLight.SLightIndsInMaskCell))==NextPercentDone)
                                %disp([ICinfoStr,'Parameter estimate & testing of Lambda distribution ',num2str(floor(IndSLCenterVox*100/length(MVSscalingTest.SLight.SLightIndsInMaskCell))),'% done.']);
                                msg = sprintf('SLavg(Amps): %03.0f percent done.', floor(IndSLCenterVox*100/length(MVSscaling.SLight.SLightIndsInMaskCell))); %Don't forget this semicolon
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
        MVSscaling.ScalingPerSubject{IndApproach}(:,IndSubj) = abs(MVSscaling.AmplitudesPerSubjectMRI{IndApproach}(:,IndSubj,2)./MVSscaling.AmplitudesPerSubjectMRI{IndApproach}(:,IndSubj,1));
        if(strcmp(MVSscaling.ApproachesInfo{IndApproach,1},'SLavg'))
            fprintf('\n');
        end
    end
    fprintf('DONE.\n');
end
fprintf('ALL DONE.\n');

end