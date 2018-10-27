function [OutputNIIpaths,Limits,Colors,OutputNII_Vols,OutputDir] = MakeOverlayParamEstLambda(InputStruct,Limits,Colors,OutputDir,SlicesForDisplay)
% This function can be used to produce the overlay of the parameter estimates of median Lambda values
% and corresponding CIwidth, given limits and colors. 
% The data in InputStruct (either MVSscalingTest-struct or TSNRscalingTest-struct) is "bined" using the limits,
% e.g. from the highest limit to the lowest and all data gets marked in the overlay with the colors according to the "bin".
%
% NB: since V1.1 this will create two kinds of outputs for the two kinds of data creation, see comment for V1.1 below.
%
% Suggested Limits = [0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)];
% Suggested Colors:   D-blue;  blue;  cyan; D-green; green; D-yellow; yellow; orange;   red;      red+;       red++;        redLimit;  white;
%           Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1]; 
% Suggested SlicesForDisplay = [-42; -36; -26; -18; -12; -4; +2; +10; +16; +24; +30; +36; +42; +48];
%   
%Usage:
%      [OutputNIIpaths,Limits,Colors,OutputNII_Vols,OutputDir] = MakeOverlayParamEstLambda(MVSscalingTest,Limits,Colors,OutputDir,SlicesForDisplay);
%      [OutputNIIpaths,Limits,Colors,OutputNII_Vols,OutputDir] = MakeOverlayParamEstLambda(MVSscalingTest,[0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)],[0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1],[],[-42; -36; -26; -18; -12; -4; +2; +10; +16; +24; +30; +36; +42; +48]);
%      [OutputNIIpaths,Limits,Colors,OutputNII_Vols,OutputDir] = MakeOverlayParamEstLambda(MVSscalingTest);  %default all other inputs and ask user to select if necessary eg. output directory
%      [OutputNIIpaths,Limits,Colors,OutputNII_Vols,OutputDir] = MakeOverlayParamEstLambda(TSNRscalingTest); %same as above but with TSNR struct. %default all other inputs and ask user to select if necessary eg. output directory
%
%V2.0
%Date: V2.0(18.09.2015): The same as before in terms of functionality, but now it can check if "MVSscalingTest"-struct or "TSNRscalingTest"-struct was input and process data correspondingly. V1.1(04.09.2015): now includes changes necessary to output the two approaches to data preparation (voxel-wise or searchlight-average of amplitudes before scaling values calculation per subject). V1.0(27.07.2015) (initial implementation based on test script for analysis of scaling data.)
%Author: Rainer.Boegle (Rainer.Boegle@googlemail.com)


%% check Colors & Limits & SlicesForDisplay
try
    if(isempty(Limits))
        Limits = [0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)];
    end
catch
    Limits = [0; 1/8*sqrt(2); 1/4*sqrt(2); 1/2*sqrt(2); 3/4*sqrt(2); sqrt(2); 5/4*sqrt(2); 3/2*sqrt(2); 7/4*sqrt(2); 2*sqrt(2); 9/4*sqrt(2); 10/4*sqrt(2); 3*sqrt(2)];
end
try
    if(isempty(Colors))
        Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1]; 
    end
catch
    Colors = [0 0 .5; 0 0 1; 0 1 1;  0 .5 0; 0 1 0;   .5 1 0;  1 1 0; 1 .5 0; 1 0 0; 1 1/4 1/4; 1 7/10 7/10; 1 7.5/10 7.5/10;  1 1 1]; 
end
try
    if(isempty(SlicesForDisplay))
        SlicesForDisplay = [-42; -36; -26; -18; -12; -4; +2; +10; +16; +24; +30; +36; +42; +48];
    else
        if(any(isnan(SlicesForDisplay)))
            DisplayResults = 0; %Don't display now.
        else
            DisplayResults = 1; %Display when done.
        end
    end
catch
    SlicesForDisplay = [-42; -36; -26; -18; -12; -4; +2; +10; +16; +24; +30; +36; +42; +48];
    DisplayResults   = 1; %Display when done.
end

if(size(Colors,2)==3)
    if(size(Colors,1)==length(Limits))
        [LimTest,SortInds] = sort(Limits,'descend');
        if(~(all(LimTest(:)==Limits)&&isequal(SortInds,1:length(Limits))))
            Limits = LimTest;
            Colors = Colors(SortInds,:);
        end
    else
        error('Need the same number of colors as limits!');
    end
else
    error('Colors must be a Nx3 [r,g,b] color matrix, with N==number of limits/bins.');
end

%% check output directory
try
    if(isempty(OutputDir))
        OutputDir = spm_select(1,'dir','Select output directory...');
    else
        if(iscell(OutputDir))
                OutputDir = OutputDir{1}; 
        end
        if(ischar(OutputDir))
            if(~exist(OutputDir,'dir'))
                mkdir(OutputDir);
            end
        else
            OutputDir = spm_select(1,'dir','Select output directory (input was not usable)...');
        end
    end
catch
    OutputDir = spm_select(1,'dir','Select output directory...');
end

%% idInputStruct
%to do: make a subfunction idInputStruct(InputStruct) that determines which input was given, MVSscalingTest-struct or TSNRscalingTest-struct
%       need to set NHood, ApproachesInfo, Medians & CIwidths in a new way AND new variable DataTypeStr that indicates if it is MVSscaling from Amplitudes or TSNRscaling
DataTypeStr     = idInputStruct(InputStruct);
SLight          = InputStruct.SLight;
ApproachesInfo  = InputStruct.ApproachesInfo;
MedianQrtCIwidth= InputStruct.ParamEstLambdaTest.MedianQrtCIwidth;


%% set indices of data kind
if(SLight.NHood>0)
    IndsData = 1:3;
else
    IndsData = 2;
end

%% init outputs
OutputNIIpaths = cell(size(ApproachesInfo,1),3,2,2); %init
OutputNII_Vols = cell(size(ApproachesInfo,1),3,2,2); %init

%% loop over approaches
for IndApproach = 1:size(ApproachesInfo,1)
    disp(['Data is from approach: ',ApproachesInfo{IndApproach,1},': ',ApproachesInfo{IndApproach,2}]);
    %% get data
    Medians = squeeze(MedianQrtCIwidth{IndApproach}(:,1,:)); %2ndDim: MedS, Agg & MedSL
    CIwidths= squeeze(MedianQrtCIwidth{IndApproach}(:,4,:)); %2ndDim: MedS, Agg & MedSL
    
    %% Do assignment
    MedianBins = NaN(size(Medians)); %init as NaN such that overlay will be transparent
    CIwidthBins= NaN(size(Medians)); %init as NaN such that overlay will be transparent
    for IndSLcenter = 1:length(SLight.SLightIndsInMaskCell)
        for IndLim = 1:length(Limits)
            for IndData = IndsData
                %Medians
                if(IndLim==1)
                    if(Medians(IndSLcenter,IndData)>=Limits(IndLim))
                        MedianBins(IndSLcenter,IndData)   = IndLim;
                    end
                else
                    if((Medians(IndSLcenter,IndData)>=Limits(IndLim))&&(Medians(IndSLcenter,IndData)<Limits(IndLim-1)))
                        MedianBins(IndSLcenter,IndData)   = IndLim;
                    end
                end
                
                
                %CIwidths
                if(IndLim==1)
                    if(CIwidths(IndSLcenter,IndData)>=Limits(IndLim))
                        CIwidthBins(IndSLcenter,IndData)   = IndLim;
                    end
                else
                    if((CIwidths(IndSLcenter,IndData)>=Limits(IndLim))&&(CIwidths(IndSLcenter,IndData)<Limits(IndLim-1)))
                        CIwidthBins(IndSLcenter,IndData)   = IndLim;
                    end
                end
            end
        end
    end
    
    %% write to nifti
    DataType = {'Median';'CIwidth'};
    DataKind = {'MedS';'Agg';'MedSL'};
    for IndBinRaw = 1:2
        for IndDataType = 1:length(DataType) %Median & CIwidth
            for IndDataKind = 1:length(DataKind) %MedS, Agg & MedSL
                if(SLight.NHood==0)
                    if(IndDataKind~=2)
                        continue;
                    end
                end
                if(IndBinRaw==1)
                    FName = [ApproachesInfo{IndApproach,1},'_',DataType{IndDataType},'Lambda_',DataKind{IndDataKind},'_',num2str(length(Limits)),'Bins','_',DataTypeStr];
                else
                    FName = [ApproachesInfo{IndApproach,1},'_',DataType{IndDataType},'LambdaRAW_',DataKind{IndDataKind},'_',DataTypeStr];
                end
                
                %V = struct([]); %init
                
                V = SLight.V_SLmask;
                if(V.dt(1)<16)
                    V.dt(1) = 16; %not necessary but save
                end
                V.fname=[OutputDir,filesep,FName,'.nii'];
                OutputNIIpaths{IndApproach,IndDataKind,IndDataType,IndBinRaw}=V.fname;
                
                if(IndDataType==1)
                    if(IndBinRaw==1)
                        CurrData = MedianBins( :,IndDataKind); %Binned
                    else
                        CurrData = Medians(    :,IndDataKind); %RAW
                    end
                else
                    if(IndBinRaw==1)
                        CurrData = CIwidthBins(:,IndDataKind); %Binned
                    else
                        CurrData = CIwidths(   :,IndDataKind); %RAW
                    end
                end
                Y = zeros(size(SLight.SLmaskRaw));
                Y(SLight.SLmaskRaw~=0) = CurrData;
                Y = reshape(Y,V.dim);
                spm_write_vol(V,Y);
                
                OutputNII_Vols{IndApproach,IndDataKind,IndDataType,IndBinRaw} = V;
            end
        end
    end
    
    %% Display Overlay
    %% check SlicesForDisplay
    try
        if(~isempty(SlicesForDisplay)&&~any(isnan(SlicesForDisplay)))
            disp(['Using SlicesForDisplay: ',num2str(SlicesForDisplay(:)')]);
        end
    catch
        SlicesForDisplay = []; %use automatic slice selection later
    end
    
    %% Display overlays
    if(DisplayResults)
        for IndDataType = 1:length(DataType) %DataType = {'Median';'CIwidth'};
            %     for IndDataKind = 1:length(DataKind) %DataKind = {'MedS';'Agg';'MedSL'};
            H = helpdlg(['Displaying "',ApproachesInfo{IndApproach,1},'_',DataType{IndDataType},'Lambda".'],'Current Display');
            
            if(SLight.NHood>0)
                if(~isempty(SlicesForDisplay))
                    DisplayOverlayParamEstLambda(SlicesForDisplay,squeeze(OutputNIIpaths(IndApproach,:,IndDataType,1)),Colors,Limits);
                else
                    DisplayOverlayParamEstLambda(       'Cluster',squeeze(OutputNIIpaths(IndApproach,:,IndDataType,1)),Colors,Limits);
                end
            else
                if(~isempty(SlicesForDisplay))
                    DisplayOverlayParamEstLambda(SlicesForDisplay,squeeze(OutputNIIpaths(IndApproach,2,IndDataType,1)),Colors,Limits);
                else
                    DisplayOverlayParamEstLambda(       'Cluster',squeeze(OutputNIIpaths(IndApproach,2,IndDataType,1)),Colors,Limits);
                end
            end
            uiwait(H);
            %     end
        end
    end
    
end

end


%% subfunction idInputStruct
function DataTypeStr = idInputStruct(InputStruct)
% This function identifies which kind of struct is input, either MVSscalingTest or TSNRscalingTest.

if(isfield(InputStruct,'Design')) %either MVSscalingTest or TSNRscalingTest
    if(isfield(InputStruct.Design,'SubjNrScannerNrRunNr')) %Must be MVSscalingTest
        DataTypeStr = ['IC',num2str(InputStruct.ICnum,'%02g')];
    elseif(isfield(InputStruct.Design,'SubjNrMRInrRunNr')) %Must be TSNRscalingTest
        DataTypeStr = 'TSNR';
    else
        error('InputStruct unknown!');
    end
else
    error('InputStruct does not have a design field! Can not be identified for output of data.');
end

end
    
