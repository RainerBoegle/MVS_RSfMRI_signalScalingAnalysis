function [VolsCI,VolsPVals,VolsPValsINVERSE,VolsZVals,VolsPVals_sign,VolsPVals_signINVERSE,VolsZVals_sign,VolsDataQuality,VolsSkewKurt,VolsCorr,VolsAmpCorr]=WriteNII_SLResults(MVSscalingTest,OutputDir)
% This function allows to write out the results of a searchlight analysis
% of scaling data.
%
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment(09.February.2015): initial implementation based on test script.

%% check inputs
try
    if(ischar(MVSscalingTest))
        load(MVSscalingTest); %assume it is a path
    else
        if(~isstruct(MVSscalingTest))
            error('MVSscalingTest is not a struct!');
        end
    end
catch
    load(spm_select(1,'mat','Select MVSscalingTest *.mat-file...',[],pwd,'^ResultsSL_',1));
end
try
    if(~isempty(OutputDir))
        if(~exist(OutputDir))
            mkdir(OutputDir);
        end
    else
        OutputDir = spm_select(1,'dir','Select output directory for nifti...');
    end
catch
    OutputDir = spm_select(1,'dir','Select output directory for nifti...');
end
    
%% smooth used?
try
    if(MVSscalingTest.UseSLsmooth)
        SmoothStr = 's';
    else
        SmoothStr = '';
    end
catch
    SmoothStr = '';
end

%% signrank stats is a bit special therefore ask the user based on the statistics type, what should be done.
%% p<alpha or (1-p)<alpha? --> inform about statistics type!
disp(' ');
UsePVals = 1; %init
PStr = ''; %init empty
if(strcmp('1-p',questdlg('Use p-values or 1-p values?','p or 1-p?','1-p','p','1-p')))
    UsePVals = 0;
    PStr = 'Not'; %1-p
    disp('-------');
    disp('(p=1-p)');
end

%% FWE correction? (on number of input SL or whole brain voxels?)
MultCorrStr = ''; %init empty
if(strcmp('Yes',questdlg('Apply multiple testing correction?','MultTestCorr?','Yes','No','No')))
    ans_alpha = inputdlg({'alpha= '},'alpha?',1,{'0.05'});
    try
        alpha = eval(ans_alpha{1});
    catch
        alpha = 0.05;
    end
    UseThres = 1;
    UseFWE = 1; %init
    UseNVoxWholeBrain = 0; %init
    
    if(strcmp('none',questdlg('Use "FWE" correction or "none"?','FWE or none?','FWE','none','none'))) %maybe FDR as well?
        UseFWE = 0;
        UseNVoxWholeBrain = 0;
        MultCorrStr = 'NONE'; 
    else
        MultCorrStr = 'FWE'; 
        if(strcmp('N-WholeBrain',questdlg('Use all voxels of whole brain for FWE correction or the number of all searchlights?','N?','N-WholeBrain','N-SLights','N-SLights')))
            UseNVoxWholeBrain = 1;
        end
    end
else
    alpha = [];
    UseThres = 0;
    UseFWE = 0;
    UseNVoxWholeBrain = 0;
end

%% prep filenames
[basedir,fnameMask]   = fileparts(MVSscalingTest.SLight.V_SLmask.fname);
[basedir,fnameWBMask] = fileparts(MVSscalingTest.SLight.V_WholeBrainMask.fname);

%% prep results fields
FieldsToOutput = fieldnames(MVSscalingTest.Results);

%% do CIoverlap 
VolsCI = cell(length(FieldsToOutput),1); %init
for IndFields = 1:length(FieldsToOutput)
    if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
        Data = MVSscalingTest.Results.(FieldsToOutput{IndFields}).CIoverlap;
        if(all(Data(:)==0))
            disp([FieldsToOutput{IndFields},'-CIoverlap contains only zeros! Skipping...']);
            VolsCI{IndFields} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CI',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CI',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsCI{IndFields} = V_out;
        end
    end
end

%% in case of higher or lower check with median?
if(strcmp(MVSscalingTest.StatisticsSettings.ChoiceTest,'Different'))
    CheckMedian = 0; %don't check with expected median
else
    ChoiceCheckMedian = questdlg('Check expected median with sample median OR with (appropriate) Confidence Interval?','Check with what?','Median','CI','Do not check!','Do not check!');
    if(strcmp('Median',ChoiceCheckMedian))
        CheckMedian = 1; %check median with expected median
    else
        if(strcmp('CI',ChoiceCheckMedian))
            switch(MVSscalingTest.StatisticsSettings.ChoiceTest)
                case 'Higher'
                    CheckMedian = 4; %check low notch with expected median
                otherwise
                    CheckMedian = 5; %check high notch with expected median
            end
        else
            CheckMedian = 0; %don't check with expected median
        end
    end
end

%% signrank test
disp(' ');
disp('signrank test');
%% output pvals & INVERSE
VolsPVals       = cell(length(FieldsToOutput),1); %init
VolsPValsINVERSE= cell(length(FieldsToOutput),1); %init
for IndFields = 1:length(FieldsToOutput)
    if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
        if(UsePVals)
            Data =   MVSscalingTest.Results.(FieldsToOutput{IndFields}).p;
        else
            Data = 1-MVSscalingTest.Results.(FieldsToOutput{IndFields}).p;
        end
        if(CheckMedian~=0)
            MedCheckData = MVSscalingTest.Results.(FieldsToOutput{IndFields}).MedianQrtCI(:,CheckMedian);
            switch(MVSscalingTest.StatisticsSettings.ChoiceTest)
                case 'Higher'
                    MedCheckData_tmp = MedCheckData;
                    MedCheckData = zeros(size(MedCheckData_tmp)); %reinit
                    MedCheckData(MedCheckData_tmp>MVSscalingTest.StatisticsSettings.ExpectedMedian) = 1;
                otherwise %Lower
                    MedCheckData_tmp = MedCheckData;
                    MedCheckData = zeros(size(MedCheckData_tmp)); %reinit
                    MedCheckData(MedCheckData_tmp<MVSscalingTest.StatisticsSettings.ExpectedMedian) = 1;
            end
        else
            MedCheckData = ones(size(MVSscalingTest.Results.(FieldsToOutput{IndFields}).p));
        end
        if(UseThres)
            if(CheckMedian~=0)
                switch(MVSscalingTest.StatisticsSettings.ChoiceTest)
                    case 'Higher'
                        if(CheckMedian==1)
                            ShortTitle = ['PT',MultCorrStr,PStr,'hiMed'];
                        else
                            ShortTitle = ['PT',MultCorrStr,PStr,'hiCI'];
                        end
                    otherwise
                        if(CheckMedian==1)
                            ShortTitle = ['PT',MultCorrStr,PStr,'loMed'];
                        else
                            ShortTitle = ['PT',MultCorrStr,PStr,'loCI'];
                        end
                end
            else
                ShortTitle = ['PThres',MultCorrStr,PStr];
            end
            if(~UsePVals)
                ShortTitle = ['1',regexprep(ShortTitle,'PThres','PTr')];
            end
            if(UseNVoxWholeBrain)
                ShortTitle = [ShortTitle,'WB'];
            end
            if(UseFWE)
                if(UseNVoxWholeBrain)
                    alpha_test = alpha/length(find(MVSscalingTest.SLight.WholeBrainRaw~=0));
                else
                    alpha_test = alpha/length(find(MVSscalingTest.SLight.SLmaskRaw~=0));
                end
            else
                alpha_test = alpha;
            end
            Data_tmp = Data;
            Data = zeros(size(Data_tmp)); %reinit
            Data((Data_tmp<alpha_test)&MedCheckData) = 1;
            if(any((Data_tmp<alpha_test)&MedCheckData))
                disp(['For "',FieldsToOutput{IndFields},'": ',num2str(length(find((Data_tmp<alpha_test)&MedCheckData))),'-Voxels are significant at alpha= ',num2str(alpha_test)]);
            end
        else
            ShortTitle = ['PVals',MultCorrStr,PStr];
        end
        if(UseThres&&all(Data(:)==0))
            if(exist('Data_tmp','var'))
                disp([FieldsToOutput{IndFields},'-Pvals not significant! (min(p)=',num2str(min(Data_tmp(:))),') Skipping...']);
            else
                if(min(Data(:))~=0)
                    disp([FieldsToOutput{IndFields},'-Pvals not significant! (min(p)=',num2str(min(Data(:))),') Skipping...']);
                else
                    disp([FieldsToOutput{IndFields},'-Pvals not significant! Skipping...']);
                end
            end
            VolsPVals{IndFields} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsPVals{IndFields} = V_out;
            
            %% INVERSE
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_INVERSE_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_INVERSE_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = ~Data; %INVERSE Mask, i.e. NOT(Data)
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsPValsINVERSE{IndFields} = V_out;
            
        end
    end
end

%% output zvals
VolsZVals = cell(length(FieldsToOutput),1); %init
for IndFields = 1:length(FieldsToOutput)
    if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
        Data = MVSscalingTest.Results.(FieldsToOutput{IndFields}).zval;
        if(all(Data(:)==0))
            disp([FieldsToOutput{IndFields},'-zvals contains only zeros! Skipping...']);
            VolsZVals{IndFields} = [];
        else
            if(any(isinf(Data(:))))
                Data(isinf(Data(:))) = max(Data(~isinf(Data(:)))).*2; %maybe nonsense but worth a try.
                disp([FieldsToOutput{IndFields},'-zvals contain Infinite-Values! Adjusting...']);
            end
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_zVals',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_zVals',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsZVals{IndFields} = V_out;
        end
    end
end

%% signtest
disp(' ');
disp('signtest');
%% signtest - output pvals & INVERSE
VolsPVals_sign       = cell(length(FieldsToOutput),1); %init
VolsPVals_signINVERSE= cell(length(FieldsToOutput),1); %init
for IndFields = 1:length(FieldsToOutput)
    if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
        if(UsePVals)
            Data =   MVSscalingTest.Results.(FieldsToOutput{IndFields}).p_sign;
        else
            Data = 1-MVSscalingTest.Results.(FieldsToOutput{IndFields}).p_sign;
        end
        if(CheckMedian~=0)
            MedCheckData = MVSscalingTest.Results.(FieldsToOutput{IndFields}).MedianQrtCI(:,CheckMedian);
            switch(MVSscalingTest.StatisticsSettings.ChoiceTest)
                case 'Higher'
                    MedCheckData_tmp = MedCheckData;
                    MedCheckData = zeros(size(MedCheckData_tmp)); %reinit
                    MedCheckData(MedCheckData_tmp>MVSscalingTest.StatisticsSettings.ExpectedMedian) = 1;
                otherwise %Lower
                    MedCheckData_tmp = MedCheckData;
                    MedCheckData = zeros(size(MedCheckData_tmp)); %reinit
                    MedCheckData(MedCheckData_tmp<MVSscalingTest.StatisticsSettings.ExpectedMedian) = 1;
            end
        else
            MedCheckData = ones(size(MVSscalingTest.Results.(FieldsToOutput{IndFields}).p_sign));
        end
        if(UseThres)
            if(CheckMedian~=0)
                switch(MVSscalingTest.StatisticsSettings.ChoiceTest)
                    case 'Higher'
                        if(CheckMedian==1)
                            ShortTitle = ['PT',MultCorrStr,PStr,'hiMed'];
                        else
                            ShortTitle = ['PT',MultCorrStr,PStr,'hiCI'];
                        end
                    otherwise
                        if(CheckMedian==1)
                            ShortTitle = ['PT',MultCorrStr,PStr,'loMed'];
                        else
                            ShortTitle = ['PT',MultCorrStr,PStr,'loCI'];
                        end
                end
            else
                ShortTitle = ['PThres',MultCorrStr,PStr];
            end
            if(~UsePVals)
                ShortTitle = ['1',regexprep(ShortTitle,'PThres','PTr')];
            end
            if(UseNVoxWholeBrain)
                ShortTitle = [ShortTitle,'WB'];
            end
            if(UseFWE)
                if(UseNVoxWholeBrain)
                    alpha_test = alpha/length(find(MVSscalingTest.SLight.WholeBrainRaw~=0));
                else
                    alpha_test = alpha/length(find(MVSscalingTest.SLight.SLmaskRaw~=0));
                end
            else
                alpha_test = alpha;
            end
            Data_tmp = Data;
            Data = zeros(size(Data_tmp)); %reinit
            Data((Data_tmp<alpha_test)&MedCheckData) = 1;
            if(any((Data_tmp<alpha_test)&MedCheckData))
                disp(['For "',FieldsToOutput{IndFields},'": ',num2str(length(find((Data_tmp<alpha_test)&MedCheckData))),'-Voxels are significant at alpha= ',num2str(alpha_test)]);
            end
        else
            ShortTitle = ['PVals',MultCorrStr,PStr];
        end
        if(UseThres&&all(Data(:)==0))
            disp([FieldsToOutput{IndFields},'-Pvals not significant! Skipping...']);
            VolsPVals_sign{IndFields} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SIGNTEST_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SIGNTEST_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsPVals_sign{IndFields} = V_out;
            
            %% INVERSE
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SIGNTESTINVERSE_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SIGNTESTINVERSE_',ShortTitle,'_',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = ~Data; %INVERSE Mask, i.e. NOT(Data)
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsPVals_signINVERSE{IndFields} = V_out;
            
        end
    end
end

%% signtest - output zvals
VolsZVals_sign = cell(length(FieldsToOutput),1); %init
for IndFields = 1:length(FieldsToOutput)
    if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
        Data = MVSscalingTest.Results.(FieldsToOutput{IndFields}).zval_sign;
        if(all(Data(:)==0))
            disp([FieldsToOutput{IndFields},'-zvals_sign contains only zeros! Skipping...']);
            VolsZVals_sign{IndFields} = [];
        else
            if(any(isinf(Data(:))))
                if(all(isinf(Data(:))))
                    disp([FieldsToOutput{IndFields},'-zvals_sign contains ONLY Infinite-Values! Skipping...']);
                    VolsZVals_sign{IndFields} = Inf;
                    continue;
                else
                    disp([FieldsToOutput{IndFields},'-zvals_sign contain Infinite-Values! Adjusting...']);
                    Data(isinf(Data(:))) = max(Data(~isinf(Data(:)))).*2; %maybe nonsense but worth a try.
                end
            end
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SIGNTESTzVals',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SIGNTESTzVals',regexprep(regexprep(FieldsToOutput{IndFields},'MedianSLight','MedSL'),'MedianSubj','MedS'),'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',PStr,MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsZVals_sign{IndFields} = V_out;
        end
    end
end


%% DATA QUALITY
%% DataQuality
% MVSscalingTest.Results.DataQuality.UsableData     = ones( length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as all data usable
% MVSscalingTest.Results.DataQuality.SomeMissingData= zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as no missing data
% MVSscalingTest.Results.DataQuality.OutliersGeneral= zeros(length(MVSscalingTest.SLight.SLightIndsInMaskCell),1); %init as no outliers generally
% MVSscalingTest.Results.DataQuality.OutliersCenter 
VolsDataQuality = cell(4,1);
if(isfield(MVSscalingTest.Results,'DataQuality'))
    %% UsableData
    if(isfield(MVSscalingTest.Results.DataQuality,'UsableData'))
        Data = MVSscalingTest.Results.DataQuality.UsableData;
        if(all(Data(:)==0))
            disp(['DataQuality.UsableData contains only zeros! Skipping...']);
            VolsDataQuality{1} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_UsableData_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_UsableData_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsDataQuality{1} = V_out;
        end
    end
    
    %% SomeMissingData
    if(isfield(MVSscalingTest.Results.DataQuality,'SomeMissingData'))
        Data = MVSscalingTest.Results.DataQuality.SomeMissingData;
        if(all(Data(:)==0))
            disp(['DataQuality.SomeMissingData contains only zeros! Skipping...']);
            VolsDataQuality{2} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SomeMissingData_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SomeMissingData_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsDataQuality{2} = V_out;
        end
    end
    
    %% OutliersGeneral
    if(isfield(MVSscalingTest.Results.DataQuality,'OutliersGeneral'))
        Data = MVSscalingTest.Results.DataQuality.OutliersGeneral;
        if(all(Data(:)==0))
            disp(['DataQuality.OutliersGeneral contains only zeros! Skipping...']);
            VolsDataQuality{3} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_OutliersGeneral_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_OutliersGeneral_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsDataQuality{3} = V_out;
        end
    end
    
    %% OutliersGeneral
    if(isfield(MVSscalingTest.Results.DataQuality,'OutliersCenter'))
        Data = MVSscalingTest.Results.DataQuality.OutliersCenter;
        if(all(Data(:)==0))
            disp(['DataQuality.OutliersCenter contains only zeros! Skipping...']);
            VolsDataQuality{4} = [];
        else
            V_out= MVSscalingTest.SLight.V_SLmask;
            if(V_out.dt(1)<16)
                V_out.dt(1) = 16; %not necessary but save
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_OutliersCenter_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_OutliersCenter_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
            VolsDataQuality{4} = V_out;
        end
    end
else
    disp('DataQuality was not checked in this analysis.');
end

%% Skewness and Kurtosis
VolsSkewKurt = cell(length(FieldsToOutput),2);
if(isfield(MVSscalingTest.Results.Aggregate,'SkewKurt'))
    for IndFields = 1:length(FieldsToOutput)
        if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
            for IndSkewKurt = 1:2
                if(IndSkewKurt==1)
                    TitleSkewKurt = 'SKEWNESS';
                else
                    TitleSkewKurt = 'KURTOSIS';
                end
                Data = MVSscalingTest.Results.(FieldsToOutput{IndFields}).SkewKurt(:,IndSkewKurt);
                
                if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                    V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_',TitleSkewKurt,'_',FieldsToOutput{IndFields},'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
                else
                    V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_',TitleSkewKurt,'_',FieldsToOutput{IndFields},'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
                end
                Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
                Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
                Y = reshape(Y,V_out.dim);
                
                try
                    V_out = spm_write_vol(V_out,Y);
                catch
                    V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                    V_out = spm_write_vol(V_out,Y);
                end
                VolsSkewKurt{IndFields,IndSkewKurt} = V_out;
            end
        end
    end
else
    VolsSkewKurt = [];
    disp('Skewness and Kurtosis was not checked in this analysis.');
end

%% sign test
if(strcmp(MVSscalingTest.StatisticsSettings.ChoiceTest,'Different'))
    for IndFields = 1:length(FieldsToOutput)
        if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
            Data = ~MVSscalingTest.Results.(FieldsToOutput{IndFields}).h_sign;
            if(all(Data==0))
                disp(['No SignH0 simple voxels found for ',FieldsToOutput{IndFields},'.']);
                continue;
            else
                disp(['Writing ',num2str(length(find(Data))),' simple significant (Sign==H0) voxels found for ',FieldsToOutput{IndFields},'.']);
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SignH0_simple_',FieldsToOutput{IndFields},'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SignH0_simple_',FieldsToOutput{IndFields},'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
%             VolsDataQuality{4} = V_out;
        end
    end
    %% minimal correction
    for IndFields = 1:length(FieldsToOutput)
        if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
            Data = ~MVSscalingTest.Results.(FieldsToOutput{IndFields}).h_sign;
            Data((1-MVSscalingTest.Results.(FieldsToOutput{IndFields}).p_sign(Data==1))>(alpha/length(find(Data==1)))) = 0; %not significant using multiple comparison correction
            if(all(Data==0))
                disp(['No SignH0 simple voxels found for ',FieldsToOutput{IndFields},'.']);
                continue;
            else
                disp(['Writing ',num2str(length(find(Data))),' simple significant (Sign==H0) voxels found for ',FieldsToOutput{IndFields},'.']);
            end
            if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SignH0_minmultcorr_',FieldsToOutput{IndFields},'_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            else
                V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_SignH0_minmultcorr_',FieldsToOutput{IndFields},'_All', regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'_Med',MVSscalingTest.StatisticsSettings.ChoiceTest,regexprep(MVSscalingTest.StatisticsSettings.answer_SetupStats{1},'*|\(|\)',''),'.nii'];
            end
            Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
            Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
            Y = reshape(Y,V_out.dim);
            try
                V_out = spm_write_vol(V_out,Y);
            catch
                V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
                V_out = spm_write_vol(V_out,Y);
            end
%             VolsDataQuality{4} = V_out;
        end
    end
end

%% correlation results
VolsCorr = cell(8,1);
if(isfield(MVSscalingTest.Results.MedianSLight,'Corr'))
    V_out= MVSscalingTest.SLight.V_SLmask;
    
    %pearson corr==1
    P_Data = MVSscalingTest.Results.MedianSLight.Corr.P(:,1);
    R_Data = MVSscalingTest.Results.MedianSLight.Corr.Rho(:,1);
    Pperm  = MVSscalingTest.Results.MedianSLight.Corr.Pperm(:,1);
    Pperm2 = MVSscalingTest.Results.MedianSLight.Corr.Pperm2(:,1);
    
    %output p as -log10
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_log10P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_log10P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = -log10(P_Data);
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{1} = V_out;
    %output 1-p 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_1P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_1P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-P_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{2} = V_out;
    
    %output Rho 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_R_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_R_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = R_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{3} = V_out;
    %output Rho as z-vals using Fisher z-transformation
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_ZtrafoR_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_ZtrafoR_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 0.5.*log((1+R_Data)./(1-R_Data));
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{4} = V_out;
    
    %output 1-p PERM-Test
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_1Pperm_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_1Pperm_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm;  %same dir as original Rho only considered
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{5} = V_out;
    
    %output 1-p PERM-Test ABS
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_1PpermABS_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_CorrLambda_1PpermABS_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm2; %both ways
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{6} = V_out;
    
    %spearman rank corr==2
    P_Data = MVSscalingTest.Results.MedianSLight.Corr.P(:,2);
    R_Data = MVSscalingTest.Results.MedianSLight.Corr.Rho(:,2);
    Pperm  = MVSscalingTest.Results.MedianSLight.Corr.Pperm(:,2);
    Pperm2 = MVSscalingTest.Results.MedianSLight.Corr.Pperm2(:,2);
    
    %output p as -log10
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_log10P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_log10P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = -log10(P_Data);
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{7} = V_out;
    
    %output 1-p 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_1P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_1P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-P_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{8} = V_out;
    
    %output Rho 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_R_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_R_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = R_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{9} = V_out;
    
    %output Rho as z-vals using Fisher z-transformation
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_ZtrafoR_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_ZtrafoR_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 0.5.*log((1+R_Data)./(1-R_Data));
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{10} = V_out;
    
    %output 1-p PERM-Test
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_1Pperm_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_1Pperm_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm;  %same dir as original Rho only considered
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{11} = V_out;
    
    %output 1-p PERM-Test ABS
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_1PpermABS_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_RankCorrLambda_1PpermABS_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm2; %both ways
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{12} = V_out;
end

%% AmpCorr & AbsAmpCorr
VolsAmpCorr = {};
if(isfield(MVSscalingTest.Results.MedianSLight,'AmpCorr')||isfield(MVSscalingTest.Results.MedianSLight,'AbsAmpCorr'))
    V_out= MVSscalingTest.SLight.V_SLmask;
    
    %% AmpCorr
    %pearson corr==1
    P_Data = MVSscalingTest.Results.MedianSLight.AmpCorr.P(:,1);
    R_Data = MVSscalingTest.Results.MedianSLight.AmpCorr.Rho(:,1);
    Pperm  = MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm(:,1);
    Pperm2 = MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm2(:,1);
    
    %output 1-p 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_1P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_1P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-P_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output Rho 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_R_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_R_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = R_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{end+1} = V_out;
    
    %output 1-p PERM-test 1 (incl. direction)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_1Pperm_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_1Pperm_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output 1-p PERM-test 2 (ABS rho)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_1PpermABS_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpCorrSPV_1PpermABS_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm2;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %% spearman corr==2
    P_Data = MVSscalingTest.Results.MedianSLight.AmpCorr.P(:,2);
    R_Data = MVSscalingTest.Results.MedianSLight.AmpCorr.Rho(:,2);
    Pperm  = MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm(:,2);
    Pperm2 = MVSscalingTest.Results.MedianSLight.AmpCorr.Pperm2(:,2);
    
    %output 1-p 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_1P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_1P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-P_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output Rho 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_R_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_R_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = R_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{end+1} = V_out;
    
    %output 1-p PERM-test 1 (incl. direction)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_1Pperm_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_1Pperm_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output 1-p PERM-test 2 (ABS rho)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_1PpermABS_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AmpRankCorrSPV_1PpermABS_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm2;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    
    
    %% ABS AmpCorr
    %% ABS AmpCorr pearson corr==1
    P_Data = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.P(:,1);
    R_Data = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Rho(:,1);
    Pperm  = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm(:,1);
    Pperm2 = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm2(:,1);
    
    %output 1-p 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_1P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_1P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-P_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output Rho 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_R_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_R_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = R_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{end+1} = V_out;
    
    %output 1-p PERM-test 1 (incl. direction)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_1Pperm_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_1Pperm_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output 1-p PERM-test 2 (ABS rho)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_1PpermABS_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpCorrSPV_1PpermABS_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm2;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %% ABS AmpCorr spearman corr==2
    P_Data = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.P(:,2);
    R_Data = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Rho(:,2);
    Pperm  = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm(:,2);
    Pperm2 = MVSscalingTest.Results.MedianSLight.AbsAmpCorr.Pperm2(:,2);
    
    %output 1-p 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_1P_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_1P_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-P_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output Rho 
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_R_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_R_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = R_Data;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsCorr{end+1} = V_out;
    
    %output 1-p PERM-test 1 (incl. direction)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_1Pperm_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_1Pperm_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
    
    %output 1-p PERM-test 2 (ABS rho)
    if(MVSscalingTest.SLight.RestrictAlsoToStartInds)
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_1PpermABS_Crop',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    else
        V_out.fname = [OutputDir,filesep,SmoothStr,'SL',num2str(MVSscalingTest.SLight.NHood),'_AbsAmpRankCorrSPV_1PpermABS_All',regexprep(fnameMask,' ','_'),'_',regexprep(fnameWBMask,' ','_'),'.nii'];
    end
    Y = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = 1-Pperm2;
    Y = reshape(Y,V_out.dim);
    try
        V_out = spm_write_vol(V_out,Y);
    catch
        V_out.fname = regexprep(regexprep(regexprep(V_out.fname,'MedianSLight','MedSL'),'MedianSubj','MedS'),'Aggregate','Agg');
        V_out = spm_write_vol(V_out,Y);
    end
    VolsAmpCorr{end+1} = V_out;
end

%% Done.
disp('Done with NIFTI output.');
disp(' ');
end