function CovarsStruct = GetCovarsViaSubjNrs(CovarsMatFilePath,FileListInfo)
% This function determines the covars from the covars *.mat file, matching subject numbers 
% and writes them into the structure "CovarsStruct".
%
%Usage:
%       CovarsStruct = GetCovarsViaSubjNrs(CovarsMatFilePath,FileListInfo);
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment(22.July.2015): initial implementation based on test script.


%% get covars *.mat file
EyeMovementDataRaw = load(CovarsMatFilePath); %gets EyeMovementDataRaw.ScalingFactorsLambda; EyeMovementDataRaw.SPVdata; EyeMovementDataRaw.InfoStruct;

%% get Lambda, SPV in head coil and angles per subject from EyeMovementDataRaw Info struct into RawCovars
RawCovars.IDsubjEyeMovementData = zeros(size(EyeMovementDataRaw.InfoStruct.ID_Subj)); %number assign to subjects from eye movement data --> assign from ID_Subj
for IndSubj = 1:length(EyeMovementDataRaw.InfoStruct.ID_Subj)
    RawCovars.IDsubjEyeMovementData(IndSubj) = str2num(EyeMovementDataRaw.InfoStruct.ID_Subj{IndSubj}(2:end));
end

RawCovars.ScalingFactorsLambda = EyeMovementDataRaw.ScalingFactorsLambda;
RawCovars.SPVcorrData = EyeMovementDataRaw.SPVdata(:,[2,6]);
RawCovars.HeadAngleData = zeros(length(EyeMovementDataRaw.InfoStruct.ID_Subj),2);
for IndSubj = 1:length(EyeMovementDataRaw.InfoStruct.ID_Subj)
    RawCovars.HeadAngleData(IndSubj,1) = EyeMovementDataRaw.InfoStruct.Angles{IndSubj,1}(2); %MRI 1
    RawCovars.HeadAngleData(IndSubj,2) = EyeMovementDataRaw.InfoStruct.Angles{IndSubj,2}(2); %MRI 2
end

%% match subject numbers using FileListInfo.UniqueSubjNrs
CovarsStruct.CovarName{1} = 'lambda';
CovarsStruct.CovarName{2} = 'SPV-MRI1';
CovarsStruct.CovarName{3} = 'SPV-MRI2';
CovarsStruct.CovarName{4} = 'HeadAngle-MRI1';
CovarsStruct.CovarName{5} = 'HeadAngle-MRI2';
CovarsStruct.Regressors   = zeros(length(FileListInfo.UniqueSubjNrs),length(CovarsStruct.CovarName));
% assign subjects' data according to coding of IDs
for Ind = 1:length(FileListInfo.UniqueSubjNrs)
    CorrespondingSubjInd = find(RawCovars.IDsubjEyeMovementData==FileListInfo.UniqueSubjNrs(Ind)); %the index that fits the subject number relative to eye movement data for the fMRI data that we use here.
    
    CovarsStruct.Regressors(Ind,1) = RawCovars.ScalingFactorsLambda(CorrespondingSubjInd); %this should be the correct scaling data for the correlation analysis
    CovarsStruct.Regressors(Ind,2) = RawCovars.SPVcorrData(CorrespondingSubjInd,1); %this should be the correct SPV data (MRI 1) for the correlation analysis
    CovarsStruct.Regressors(Ind,3) = RawCovars.SPVcorrData(CorrespondingSubjInd,2); %this should be the correct SPV data (MRI 2) for the correlation analysis
    CovarsStruct.Regressors(Ind,4) = RawCovars.HeadAngleData(CorrespondingSubjInd,1); %this should be the correct HeadAngle data (MRI 1) for the correlation analysis
    CovarsStruct.Regressors(Ind,5) = RawCovars.HeadAngleData(CorrespondingSubjInd,2); %this should be the correct HeadAngle data (MRI 2) for the correlation analysis
end


%% get also EyeMovementDataRaw & FileListInfo data and assign output
CovarsStruct.FileListInfo       = FileListInfo;
CovarsStruct.EyeMovementDataRaw = EyeMovementDataRaw;

end