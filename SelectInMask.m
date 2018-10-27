%% Select mask file
MaskPath = spm_select(1,'image','Select Mask *.nii file...');

V_Mask = spm_vol(MaskPath);

%% enter coordinate regions for finding values
answer_regselect = inputdlg({'x[mm]= [start end]';'y[mm]= [start end]';'z[mm]= [start end]'},'region select',1,{':',':',':'});

%% determine voxel coordinates

x = []; %init
y = []; %init
z = []; %init
invMat = inv(V_Mask.mat);
for IndAns = 1:3
    if(strcmp(answer_regselect{IndAns},':'))
        switch(IndAns)
            case 1
                x = 1:V_Mask.dim(IndAns);
            case 2
                y = 1:V_Mask.dim(IndAns);
            case 3
                z = 1:V_Mask.dim(IndAns);
        end
    else
        MMcoords = eval(answer_regselect{IndAns});
        switch(IndAns)
            case 1
                VoxStartStop = sort([floor((MMcoords(1)*invMat(1,1))+invMat(1,4)); ceil((MMcoords(2)*invMat(1,1))+invMat(1,4))]);
                Inds = VoxStartStop(1):VoxStartStop(2);
                x = Inds;
            case 2
                VoxStartStop = sort([floor((MMcoords(1)*invMat(2,2))+invMat(2,4)); ceil((MMcoords(2)*invMat(2,2))+invMat(2,4))]);
                Inds = VoxStartStop(1):VoxStartStop(2);
                y = Inds;
            case 3
                VoxStartStop = sort([floor((MMcoords(1)*invMat(3,3))+invMat(3,4)); ceil((MMcoords(2)*invMat(3,3))+invMat(3,4))]);
                Inds = VoxStartStop(1):VoxStartStop(2);
                z = Inds;
        end
    end
end

%% make new mask
NewMask = zeros(V_Mask.dim);

OrgMask = V_Mask.private.dat(:,:,:);
NewMask(x,y,z) = OrgMask(x,y,z);
% NewMask(:,46:51,:) = round(OrgMask(:,46:51,:));

%% write out
[OutDir,FName,ext] = fileparts(MaskPath);
if(length(ext)>length('.nii')||length(ext)>length('.hdr'))
    ext = ext(1:4);
end

Vout = V_Mask;
Vout = rmfield(Vout,'private');
if(Vout.dt(1)<16)
    Vout.dt(1) = 16;
end
Vout.fname = [OutDir,filesep,'RegSelect_',datestr(now,'yymmmdd_HHMM'),'_',FName,ext];
spm_write_vol(Vout,NewMask);

%% done.
disp('DONE.');
disp(' ');

% %% Select Data files 
% DataFilesPath = spm_select(Inf,'image','Select data files...');
% DataFilesPath = cellstr(DataFilesPath);
% 
% for IndDataFile = 1:length(DataFilesPath)
%     V_Data = spm_vol(DataFilesPath{IndDataFile});
%     
%     OrgData = V_Data.private.dat(:,:,:);
%     NewData = OrgData;
%     NewData(NewMask>0&NewData>-.5) = abs(NewData(NewMask>0&NewData>-.5));
%     
%     VOut = V_Data;
%     [OutDir,FName,ext] = fileparts(VOut.fname);
%     if(length(ext)>length('.nii')||length(ext)>length('.hdr'))
%         ext = ext(1:4);
%     end
%     VOut.fname = [OutDir,filesep,'Edit_',FName,ext];
%     spm_write_vol(VOut,NewData);
%     disp(['Writing Data-File ',num2str(IndDataFile)]);
%     clear VOut
% end
% 
% %% done.
% disp('DONE.');
% disp(' ');
