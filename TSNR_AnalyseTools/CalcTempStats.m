function [TSNR,Mu,Stdev,Mask2D,V_out,OutputNIIpath,rMSSD,StdevMSSD,SucDiff,rMSSD2,StdevMSSD2,SucDiff2] = CalcTempStats(InputFuncNIIpath,MaskPath,OutputDir,Options)
% This function calculates the temporal SNR, i.e. mean(timecourse)/std(timecourse) of a functional MRI timeseries.
% NB: TSNR is mostly useful for resting-state fMRI, for other designs it can be misleading if not enough dynamics are contained in the timeseries. 
%     (See SNR & CNR definition paper by Welvaert PLOSone 2013 and also Murphy NI 2007 and most importantly Triantafyllou NI 2005!)
% The calculated mean and stdeviation will also be output.
% AND also 1./Stdev as "InvStdev_..."
%
%Usage:
%       [TSNR,Mu,Stdev] = CalcTempStats(InputFuncNIIpath,MaskPath,OutputDir,Options); 
%                           %InputFuncNIIpath: is the location of the 4D-NIFTI file of the fMRI timeseries. (NOT SPM-STYLE!!! I.e. not a list of files but one 4D NIFTI-file, given as a single string location!)
%                           %MaskPath:         is the location of the 3D-NIFTI file defining a [0,1]-mask. This can be empty or omitted, then all voxels are used in the calculation, i.e. no brain vs background distinction.
%                           %OutputDir:        is the directory for the 3D-NIFTI files that contain the temporal SNR at the voxels indicated in mask, other voxels will be NaN (Not a Number).
%                           %                  as well as the mean and the standard deviation images as 3D-NIFTI files.
%                           %                  The filename of the input (InputFuncNIIpath) will be used and the prefix "TSNR_" or "Mean_" or "Stdev_" will be added.
%                           %                  If it is empty or missing then no output NIFTI-file will be generated, -and only the return values of this function are returned. 
%                           %Options:          is a cell containing a option-string,value pair, right now only e.g. {'detrend',[100 3]}.
%                                              This means that detrending is done before the Stdev is calculated and the highpass filter has a cutoff of 100seconds and assuming a TR of 3seconds in this case.
%                                              [NEEDS SPM!!!]
%
%       [TSNR,Mu,Stdev] = CalcTempStats(InputFuncNIIpath,[],pwd); %use all voxels (no mask) and output NIFTI-file containing TSNR, mean & stdev in the current directory & NO detrending.
%       [TSNR,Mu,Stdev] = CalcTempStats(InputFuncNIIpath,[],pwd,{'detrend',[100 3]}); %as above & detrending with cutoff of 100seconds for highpass filter and assuming TR of 3seconds [NEEDS SPM!!!].
%
%
%V1.3
%Date: V1.3(04.11.2015): also allow option for detrending. V1.2(19.09.2015): also output 1./Stdev as "InvStdev_...". V1.1(18.09.2015): add mean image and stdev image output. V1.0(08.09.2015): initial implementation based on test script.
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)

%Idea for future version: instead of mask path or empty allow also single numeric value that indicates the threshold above which voxels are considered to be in the mask. Probably best to make this some kind of "percentage of the total mean of the mean timeseries"???

%% Check inputs
%InputFuncNIIpath
V_InputFuncData = spm_vol(InputFuncNIIpath); %mainly need the private part of the struct that contains the nifti object.
if(isstruct(V_InputFuncData))
    if(length(V_InputFuncData)>1)
        if(V_InputFuncData(1).private.dat.dim(4)<=1)
            error(['input functional data is not a 4D-NIFTI file. (dim(4)=',num2str(V_InputFuncData(1).private.dat.dim(4)),')']);
        end
    else
        error('input functional data is not a 4D-NIFTI file.');
    end
else
    error('input functional data is not a NIFTI file.');
end

%MaskPath
try
    if(isempty(MaskPath))
        disp('No Mask: using all voxels...');
        Mask2D = ones(V_InputFuncData(1).private.dat.dim(1)*V_InputFuncData(1).private.dat.dim(2)*V_InputFuncData(1).private.dat.dim(3),1);
    elseif(iscellstr(MaskPath))
        MaskPath_tmp = MaskPath{1};  clear MaskPath
        MaskPath     = MaskPath_tmp; clear MaskPath_tmp
        V_Mask = spm_vol(MaskPath);
        Mask2D = V_Mask.private.dat(:);
        Mask2D = Mask2D~=0; %make into [0,1]-mask if not is already.
    elseif(ischar(MaskPath))
        V_Mask = spm_vol(MaskPath);
        Mask2D = V_Mask.private.dat(:);
        Mask2D = Mask2D~=0; %make into [0,1]-mask if not is already.
    else
        error('MaskPath must be a string/char pointing to the location of the mask that should be used. (or empty)');
    end
catch CATCH_MaskPath
    disp_catch(CATCH_MaskPath,[mfilename,'>CalcTSNR>MaskPath'],'CATCH_MaskPath');
    disp('No Mask: using all voxels...');
    Mask2D = ones(V_InputFuncData(1).private.dat.dim(1)*V_InputFuncData(1).private.dat.dim(2)*V_InputFuncData(1).private.dat.dim(3),1);
end

%OutputNIIpath from OutputDir
try
    if(isempty(OutputDir))
        OutputNIIpath = [];
        disp('No output directory: No NIFTI-file will be created.');
    else
        if(iscell(OutputDir))
            OutputDir_tmp = OutputDir{1};  clear OutputDir
            OutputDir     = OutputDir_tmp; clear OutputDir_tmp
        else
            if(~ischar(OutputDir))
                error('OutputDir must be a string/char pointing to the location where the output should be written to. (or empty)');
            end
        end
        disp(['Will try to output files to directory "',OutputDir,'".']);
        if(exist(OutputDir,'dir')) %is it a directory? --> append filename of input with TSNR_
            [tmp,fName] = fileparts(V_InputFuncData(1).fname); clear tmp
            disp(['Will use the filename of the input 4D-NIFTI with the prefix "TSNR_", "Mean_" and "Stdev_" (& "InvStdev_") to create output file "TSNR_',fName,'.nii" etc...']);
            OutputNIIpath    = cell(8,1);
            OutputNIIpath{1} = [OutputDir,filesep,'TSNR_',        fName,'.nii'];
            OutputNIIpath{2} = [OutputDir,filesep,'Mean_',        fName,'.nii'];
            OutputNIIpath{3} = [OutputDir,filesep,'Stdev_',       fName,'.nii'];
            OutputNIIpath{4} = [OutputDir,filesep,'InvStdev_',    fName,'.nii'];
            OutputNIIpath{5} = [OutputDir,filesep,'rMSSD_',       fName,'.nii'];
            OutputNIIpath{6} = [OutputDir,filesep,'StdevMSSD_',   fName,'.nii'];
            OutputNIIpath{7} = [OutputDir,filesep,'rMSSD2nd_',    fName,'.nii'];
            OutputNIIpath{8} = [OutputDir,filesep,'StdevMSSD2nd_',fName,'.nii'];
        else
            error('OutputDir must be a string/char pointing to the location where the output should be written to. (or empty)');
        end
    end
catch CATCH_OutputNIIpath
    disp_catch(CATCH_OutputNIIpath,[mfilename,'>CalcTSNR>OutputNIIpath'],'CATCH_OutputNIIpath');
    disp('No output path or directory entered: No NIFTI-file will be created.');
    OutputNIIpath = [];
end

%Options: NB: will need to be extended in the future
try
    if(isempty(Options))
        UseDetrend = 0;
        disp('No detrending.');
    elseif(iscell(Options))
        if(strcmp(Options{1},'detrend')||strcmp(Options{1},'detrend2'))
            CutOff = Options{2}(1);
            TR     = Options{2}(2);
            if(strcmp(Options{1},'detrend'))
                UseDetrend = 1;
                disp(['Will detrend with CutOff= ',num2str(CutOff),'seconds and assumed TR= ',num2str(TR),'seconds.']);
            else
                UseDetrend = 2;
                disp(['Will detrend METHOD-2 with CutOff= ',num2str(CutOff),'seconds and assumed TR= ',num2str(TR),'seconds.']);
            end
        else
            UseDetrend = 0;
            disp(['UNKNOWN OPTION "',Options{1},'" --> Default to no detrending.']);
        end
    else
        UseDetrend = 0;
        disp('UNKNOWN OPTION --> Default to no detrending.');
    end
catch CATCH_Options
    disp_catch(CATCH_Options,[mfilename,'>CalcTSNR>Options'],'CATCH_Options');
    disp('No detrending.');
    UseDetrend = 0;
end

%% do processing
Data4D = V_InputFuncData(1).private.dat(:,:,:,:);
Data2D = reshape(Data4D,[],size(Data4D,4));
Data2DinMask = Data2D(Mask2D~=0,:);


%% detrend?
if(UseDetrend==0)
    Mu = mean(Data2DinMask,2);
else
    K.RT     = TR;
    K.row    = 1:size(Data2DinMask,2);
    K.HParam = CutOff;
    
    if(UseDetrend==1)
        disp('Detrending...');
        Mu           = mean(Data2DinMask,2);
        Data2DinMask = spm_filter(K,Data2DinMask')';
        Data2DinMask = Data2DinMask - repmat(mean(Data2DinMask,2),1,size(Data2DinMask,2)); %remove residual mean
        Data2DinMask = Data2DinMask + repmat(Mu,1,size(Data2DinMask,2)); %now it should be like detrend with Mean added in again, as FSL usually does it and we should get similar behavior for FSL & SPM preprocessed data.
    else
        disp('Detrend2ing...');
        Data2DinMask = spm_filter(K,Data2DinMask')';
        Mu           = mean(Data2DinMask,2); %use residual mean after detrend as mean of images
    end
end
    
%% continue processing TSNR and so on
Stdev  =  std(Data2DinMask,0,2);
TSNR   = Mu./Stdev;
%MSSD
SucDiff   = Data2DinMask(:,2:end) - Data2DinMask(:,1:end-1); %successive difference of time series
MSSD      = mean(SucDiff.^2,2); %mean squared successive difference of time series
StdevMSSD = std(SucDiff.^2,0,2);%std  squared successive difference of time series
rMSSD     = sqrt(MSSD); %sqrt of MSSD to have it similar to Stdev!
%2ndIteration for MSSD ie differences of differences
SucDiff2  = SucDiff(:,2:end) - SucDiff(:,1:end-1); %2ndIteration: successive difference of successive difference!
MSSD2     = mean(SucDiff2.^2,2); %2ndIteration: mean squared successive difference of time series
StdevMSSD2= std(SucDiff2.^2,0,2);%2ndIteration: std  squared successive difference of time series
rMSSD2    = sqrt(MSSD2); %2ndIteration: sqrt of MSSD to have it similar to Stdev!

%% output?
if(~isempty(OutputNIIpath))
    V_out = cell(8,1);
    for IndFile = 1:length(OutputNIIpath)
        [OutputNIIdir,OutputNIIfilename,ext] = fileparts(OutputNIIpath{IndFile});
        disp(['Writing NIFTI-file "',OutputNIIfilename,ext,'" to directory "',OutputNIIdir,'".']);
        V_out{IndFile} = V_InputFuncData(1);
        V_out{IndFile}.fname = OutputNIIpath{IndFile}; %assign filename and directory
        if(V_out{IndFile}.dt(1)<16)
            V_out{IndFile}.dt(1) = 16; %just to be sure the encoding is good enough.
        end
        if(V_out{IndFile}.n(1)~=1 || V_out{IndFile}.n(2)~=1) %do not encode it as a 4D-NIFTI!!!
            V_out{IndFile}.n(1) = 1;
            V_out{IndFile}.n(2) = 1;
        end
        V_out{IndFile} = rmfield(V_out{IndFile},'private'); %prevent encoding as a 4D-NIFTI!!!
        
        Y = zeros(V_out{IndFile}.dim);
        switch(IndFile)
            case 1 %TSNR
                Y(Mask2D~=0) = TSNR;
            case 2 %Mean
                Y(Mask2D~=0) = Mu;
            case 3 %Stdev
                Y(Mask2D~=0) = Stdev;
            case 4 %1./Stdev
                Y(Mask2D~=0) = 1./Stdev;
            case 5 %rMSSD
                Y(Mask2D~=0) = rMSSD;
            case 6 %StdevMSSD
                Y(Mask2D~=0) = StdevMSSD;
            case 7 %rMSSD2
                Y(Mask2D~=0) = rMSSD2;
            case 8 %StdevMSSD2
                Y(Mask2D~=0) = StdevMSSD2;
            otherwise
                error('unknown file type for output! Must only be TSNR, Mean or Stdev.');
        end
        V_out{IndFile} = spm_write_vol(V_out{IndFile},Y); %write out.
    end
else
    V_out = [];
end

end