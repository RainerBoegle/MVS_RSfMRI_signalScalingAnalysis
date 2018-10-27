function [ScalingPerSubject,SLightSmooth] = SLsmooth(ScalingPerSubject,SLightSmooth,Mode)
% This function can smooth the scaling data given a searchlight definition
% per subject!
%
%
%V1.1
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment(19.February.2015): Outlier detection for the center of the searchlight and replace or remove possible now.

%% check inputs
%% SLightSmooth
if(isempty(SLightSmooth))
    SLightSmoothDefPath = spm_select(1,'mat','Select SearchLight definition *.mat file...');
    if(isempty(SLightSmoothDefPath))
        answer_SLightSmooth = inputdlg({'SearchLight extend NHood= '},'SLight def',1,{'2'});
        SLightSmooth = GenerateSLight([],[],eval(answer_SLightSmooth{1}),[],spm_select(1,'dir','Select directory to save searchlight definition *.mat'));
    else
        load(SLightSmoothDefPath);
        if(exist('SLight','var'))
            SLightSmooth = SLight;
        else
            error(['Could not retrieve SLightSmooth from "',SLightSmoothDefPath,'"']);
        end
    end
else
    if(ischar(SLightSmooth)) %assume it is a path
        SLightSmoothDefPath = SLightSmooth; clear SLightSmooth
        load(SLightSmoothDefPath);
        if(exist('SLight','var'))
            SLightSmooth = SLight;
        else
            error(['Could not retrieve SLightSmooth from "',SLightSmoothDefPath,'"']);
        end
    else
        if(~isstruct(SLightSmooth))
            error('SLightSmooth must be a struct or string, i.e. path to the *.mat-file containing the structure. (Use empty input to prompt load. Quit load to promt generation of SLightSmooth struct.)');
        end
    end
end
%% ScalingPerSubject
if(length(SLightSmooth.AllowedIndsVol)>size(ScalingPerSubject,1))
    error('"SLightSmooth" might be incompatible with input Scaling Data "ScalingPerSubject". "length(SLightSmooth.AllowedIndsVol)>size(ScalingPerSubject,1)"!!!');
end
%% Mode
try
    if(~isempty(Mode))
        switch(Mode)
            case 'IfOutlierAssignMedianSL'
                AssignIfOutlier = 1;
                %% outlier w definition
                w = 1.5;
                SLightSmooth.Outliers.w = w;
                SLightSmooth.Outliers.Assign = 'median';
            case 'AssignMedianSL'
                AssignIfOutlier = 0;
            case 'IfOutlierAssignNaN'
                AssignIfOutlier = 2;
                %% outlier w definition
                w = 1.5;
                SLightSmooth.Outliers.w = w;
                SLightSmooth.Outliers.Assign = 'NaN';
            otherwise
                Mode = 'AssignMedianSL';
                AssignIfOutlier = 0;
        end
    else
        Mode = 'AssignMedianSL';
        AssignIfOutlier = 0;
    end
catch
    Mode = 'AssignMedianSL';
    AssignIfOutlier = 0;
end
SLightSmooth.Mode = Mode;


%% do median smoothing
H_waitbar = waitbar(0,['SearchLight-SMOOTHING(',SLightSmooth.Mode,') centered at every voxel starting...']);
for IndSLCenterVox = 1:length(SLightSmooth.SLightIndsInMaskCell)
    CurrInds = SLightSmooth.SLightIndsInMaskCell{IndSLCenterVox};
    
    if(AssignIfOutlier~=0)
        Data = ScalingPerSubject(CurrInds,:);
        Qrt  = quantile(Data,[.25 .75],1);
        CheckValHigh = Qrt(2,:)+w*(Qrt(2,:)-Qrt(1,:));
        CheckValLow  = Qrt(1,:)-w*(Qrt(2,:)-Qrt(1,:));
        for IndSubj = 1:size(ScalingPerSubject,2)
            if(ScalingPerSubject(IndSLCenterVox,IndSubj)>=CheckValHigh(IndSubj)) %those with center voxel higher than outlier limit
                if(AssignIfOutlier==1)
                    ScalingPerSubject(IndSLCenterVox,IndSubj) = nanmedian(ScalingPerSubject(CurrInds,IndSubj),1);
                else
                    ScalingPerSubject(IndSLCenterVox,IndSubj) = NaN(class(ScalingPerSubject));
                end
            else
                if(ScalingPerSubject(IndSLCenterVox,IndSubj)<=CheckValLow(IndSubj)) %those with center voxel lower  than outlier limit
                    if(AssignIfOutlier==1)
                        ScalingPerSubject(IndSLCenterVox,IndSubj) = nanmedian(ScalingPerSubject(CurrInds,IndSubj),1);
                    else
                        ScalingPerSubject(IndSLCenterVox,IndSubj) = NaN(class(ScalingPerSubject));
                    end
                end
            end
        end
    else
        ScalingPerSubject(IndSLCenterVox,:) = nanmedian(ScalingPerSubject(CurrInds,:),1);
    end
    try
        H_waitbar = waitbar(IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell),H_waitbar,['SearchLight-SMOOTHING(',SLightSmooth.Mode,') centered at every voxel ',num2str(100*IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell)),'%done...']);
    catch
        H_waitbar = waitbar(IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell),['SearchLight-SMOOTHING(',SLightSmooth.Mode,') centered at every voxel ',num2str(100*IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell)),'%done...']);
    end
end
try
    H_waitbar = waitbar(IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell),H_waitbar,['SearchLight-SMOOTHING(',SLightSmooth.Mode,') centered at every voxel ',num2str(100*IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell)),'%done...']);
catch
    H_waitbar = waitbar(IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell),['SearchLight-SMOOTHING(',SLightSmooth.Mode,') centered at every voxel ',num2str(100*IndSLCenterVox/length(SLightSmooth.SLightIndsInMaskCell)),'%done...']);
end
uiwait(H_waitbar,1);
close(H_waitbar);


end