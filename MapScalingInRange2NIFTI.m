%% get scaling analysis results
ScalingAnaMatPath = spm_select(1,'mat','Select MVSfMRIscaling-Results.mat-file...',[],pwd,'^MVSscaling_',1);
load(ScalingAnaMatPath);

[BasePath, FName, ext] = fileparts(ScalingAnaMatPath);

%% define the median value and the range that is accepted.
answer_median_n_range = inputdlg({'Median Value: ','Range in % around Median: '},'median & %-range?',1,{'2*sqrt(2)','10'});
MedianValue = eval(answer_median_n_range{1});
Percentage  = eval(answer_median_n_range{2});
Range = zeros(2,1);
Range(1) = MedianValue-(MedianValue*Percentage/100);
Range(2) = MedianValue+(MedianValue*Percentage/100);

%% get median subject scaling values
MedSubjScaling = median(abs(MVSscaling.ScalingPerSubject),2);

MedSubjScaling(~(MedSubjScaling>=Range(1)&MedSubjScaling<=Range(2)))= 0; %set all else to zero.

%% make file name and get volume template struct for writing back out
OutDir = [BasePath,filesep,'Median_pm_Percentage_NIIs'];
if(~exist(OutDir))
    mkdir(OutDir);
end
fname_out = [FName,'Median',num2str(MedianValue),'pm',num2str(Percentage),'Percent.nii'];
V_out = spm_vol(MVSscaling.Masks.MPaths{1});

%% change file name for output
V_out.fname = [OutDir,filesep,fname_out];
if(V_out.dt(1)<16)
    V_out.dt(1) = 16; %not necessary but save
end

%% write back out
Y = zeros(size(MVSscaling.Masks.WholeBrainRaw));
Y(MVSscaling.Masks.WholeBrainRaw~=0) = MedSubjScaling;
Y = reshape(Y,V_out.dim);

V_out = spm_write_vol(V_out,Y);