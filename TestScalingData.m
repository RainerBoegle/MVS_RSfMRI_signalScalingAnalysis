%% Get MVSscaling struct
ScalingAnaMatPath = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
load(ScalingAnaMatPath);

[BasePath, FName, ext] = fileparts(ScalingAnaMatPath);
MVSscalingTest.InputPath = ScalingAnaMatPath;


%% create directory for output
OutDir = [BasePath,filesep,'TestScalingData'];
if(~exist(OutDir))
    mkdir(OutDir);
end

%% define Searchlight (NHood and such)
%todo

%% smooth data by applying median of searchlight on each subject and write median to center voxel?
%to do

%% define test statistics (need to be robust so I will test the median somehow???)
answer_expectedMedian = inputdlg({'ExpectedMedian= '},'Median=?',1,{'2*sqrt(2)'});
ExpectedMedian = eval(answer_expectedMedian{1});
ChoiceTest = questdlg({'Do you want to test if the data''s median is:'; ' '; ['A. "different" from expected median(',answer_expectedMedian{1},')']; ['B. "higher" than the expected median(',answer_expectedMedian{1},')']; 'OR ';  ['C. "lower" than the expected median(',answer_expectedMedian{1},')']},'Which test?','Different','Higher','Lower','Higher');
switch(ChoiceTest)
    case 'Different'
        TailType = 'both';
    case 'Higher'
        TailType = 'right';
    case 'Lower'
        TailType = 'left';
end
MVSscalingTest.StatisticsSettings.answer_expectedMedian = answer_expectedMedian;
MVSscalingTest.StatisticsSettings.ExpectedMedian        = ExpectedMedian;
MVSscalingTest.StatisticsSettings.ChoiceTest            = ChoiceTest;
MVSscalingTest.StatisticsSettings.TailType              = TailType;

%% test each voxel for effect
MVSscalingTest.p         = ones( size(MVSscaling.ScalingPerSubject,1),1); %init not significant
MVSscalingTest.h         = zeros(size(MVSscaling.ScalingPerSubject,1),1); %init not significant
MVSscalingTest.signedrank= zeros(size(MVSscaling.ScalingPerSubject,1),1); %init not significant
MVSscalingTest.zval      = zeros(size(MVSscaling.ScalingPerSubject,1),1); %init not significant

H_waitbar = waitbar(0,['Voxel-wise test starting...']);
for IndVox = 1:size(MVSscaling.ScalingPerSubject,1)
    Data = MVSscaling.ScalingPerSubject(IndVox,:);
    try
        [MVSscalingTest.p(IndVox),MVSscalingTest.h(IndVox),stats] = signrank(Data(:),ExpectedMedian,'tail',TailType);
        MVSscalingTest.signedrank(IndVox) = stats.signedrank;
        if(isfield(stats,'zval'))
            MVSscalingTest.zval(IndVox) = stats.zval;
        end
        
    catch
        MedianData = median(Data(:));
        [MVSscalingTest.p(IndVox),MVSscalingTest.h(IndVox),stats] = signrank(Data(:),ExpectedMedian);
        switch(ChoiceTest)
            case 'Higher'
                if(MedianData<ExpectedMedian);
                   MVSscalingTest.p(IndVox)         = 1;
                   MVSscalingTest.h(IndVox)         = 0;
                   MVSscalingTest.signedrank(IndVox)= Inf;
                   MVSscalingTest.zval(IndVox)= 0;
                else
                    if(isfield(stats,'zval'))
                        MVSscalingTest.zval(IndVox) = MVSscalingTest.h(IndVox).*abs(stats.zval);
                    end
                end
            case 'Lower'
                if(MedianData>ExpectedMedian);
                   MVSscalingTest.p(IndVox)         = 1;
                   MVSscalingTest.h(IndVox)         = 0;
                   MVSscalingTest.signedrank(IndVox)= Inf;
                   MVSscalingTest.zval(IndVox)      = 0;
                else
                    if(isfield(stats,'zval'))
                        MVSscalingTest.zval(IndVox) = MVSscalingTest.h(IndVox).*abs(stats.zval);
                    end
                end
            otherwise
                MVSscalingTest.signedrank(IndVox)= stats.signedrank;
                if(isfield(stats,'zval'))
                    MVSscalingTest.zval(IndVox)  = MVSscalingTest.h(IndVox).*abs(stats.zval);
                end
        end
    end
    try
        H_waitbar = waitbar(IndVox/size(MVSscaling.ScalingPerSubject,1),H_waitbar,['Voxel-wise test ',num2str(IndVox*100/size(MVSscaling.ScalingPerSubject,1)),'% done...']);
    catch
        H_waitbar = waitbar(IndVox/size(MVSscaling.ScalingPerSubject,1),['Voxel-wise test ',num2str(IndVox*100/size(MVSscaling.ScalingPerSubject,1)),'% done...']);
    end
end
try
    H_waitbar = waitbar(IndVox/size(MVSscaling.ScalingPerSubject,1),H_waitbar,'Voxel-wise test is finished!');
catch
    H_waitbar = waitbar(IndVox/size(MVSscaling.ScalingPerSubject,1),'Voxel-wise test is finished!');
end
uiwait(H_waitbar,2);
close(H_waitbar);

%% also try testing the aggregated data in each searchlight? I.e. totally different approach where all data from all voxels in searchlight and from all subjects is aggregated and then tested for effect.


%% save the results
save([OutDir,filesep,'ResultsScalingTest_Median',ChoiceTest,regexprep(answer_expectedMedian{1},'*',''),'.mat'],'MVSscalingTest');

%% write out a mask
V_out = spm_vol(MVSscaling.Masks.MPaths{1});
if(V_out.dt(1)<16)
    V_out.dt(1) = 16; %not necessary but save
end
V_out.fname = [OutDir,filesep,'ResultsScalingTest_MaskSignif_Median',ChoiceTest,regexprep(answer_expectedMedian{1},'*',''),'.nii'];
Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
Y(MVSscaling.Masks.WholeBrainRaw~=0) = MVSscalingTest.h;
Y = reshape(Y,V_out.dim);

V_out = spm_write_vol(V_out,Y);

%% save inverse of mask?
if(strcmp('Yes',questdlg('Save INVERSE of mask?','Save InverseMask?','Yes','No','No')))
    V_out.fname = [OutDir,filesep,'ResultsScalingTest_INVERSEMaskSignif_Median',ChoiceTest,regexprep(answer_expectedMedian{1},'*',''),'.nii'];
    Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
    Y(MVSscaling.Masks.WholeBrainRaw~=0) = ~MVSscalingTest.h;
    Y = reshape(Y,V_out.dim);

    V_out = spm_write_vol(V_out,Y);
end

%% write out zvals
V_out.fname = [OutDir,filesep,'ResultsScalingTest_zVals_Median',ChoiceTest,regexprep(answer_expectedMedian{1},'*',''),'.nii'];
Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
Y(MVSscaling.Masks.WholeBrainRaw~=0) = MVSscalingTest.zval;
Y = reshape(Y,V_out.dim);

V_out = spm_write_vol(V_out,Y);


%% Done.
disp(' ');
disp('Done.');