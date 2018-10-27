%% Select mask
MaskPath = spm_select(1,'image','Select mask image...');
Vol_Mask = spm_vol(MaskPath);
MaskData = Vol_Mask.private.dat(:,:,:);
UniqueVals = unique(MaskData(:));
UV2 = UniqueVals(:); UV2(UV2==1) = []; UV2(UV2==0) = [];
if(length(UniqueVals)~=2||~isempty(UV2))
    disp('making mask a 0/1-mask via ">0".');
    MaskData = MaskData>0;
end

%% Select image to be masked with NaNs
InputImagePath = spm_select(Inf,'image','Select Image(s) to be masked...');
InputImagePath = cellstr(InputImagePath);

%% Select Output directory
OutDir = spm_select(1,'dir','Select output directory...');

%% do the masking
for IndInput = 1:length(InputImagePath)
    Vol = spm_vol(InputImagePath{IndInput});
    
    Data = Vol.private.dat(:,:,:);
    Data(~MaskData) = nan; %set everything not in mask to NaN.
    
    Vout = Vol;
    [tmp,FName,ext] = fileparts(Vol.fname);
    Vout.fname = [OutDir,filesep,FName,ext];
    
    if(Vout.dt(1)<16)
        Vout.dt(1)=16;
    end
    
    spm_write_vol(Vout,Data);
    clear Data
end