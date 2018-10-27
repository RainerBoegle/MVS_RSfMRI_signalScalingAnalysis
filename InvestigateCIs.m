function [] = InvestigateCIs(MeanAssignDir,varargin)
% Assign values for confidence intervals to map out possible areas of
% certainty.
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment V1.0(15.03.2015): initial implementation 

try
    MeanAssignDir
catch
    MeanAssignDir = 1; %upward
end

%% get data & prep output
DataPath  = spm_select(1,'mat','Select data...');
OutputDir = spm_select(1,'dir','Select Output directory...');
load(DataPath); % load('ResultsSL1_CropFINAL_NewWholeBrainMask_W_MELODIC_mask.mat');

FieldsToOutput = fieldnames(MVSscalingTest.Results);

%% limits and data assignments
LimitsDataAssignment.DataPath  = DataPath;
LimitsDataAssignment.OutputDir = OutputDir;

if(0) %1st-Version
    if(nargin>1)
        CyanLimit = varargin{1};
        RedLimit  = varargin{2};
    else
        CyanLimit = 2; %close range around   sqrt(2)
        RedLimit  = 4; %close range around 2*sqrt(2)
    end
    try
        disp(['Using CyanLimit= ',num2str(CyanLimit)]);
        disp(['Using RedLimit = ',num2str(RedLimit)]);
    catch
        CyanLimit = 2; %close range around   sqrt(2)
        RedLimit  = 4; %close range around 2*sqrt(2)
        disp(['Using CyanLimit= ',num2str(CyanLimit)]);
        disp(['Using RedLimit = ',num2str(RedLimit)]);
    end
    InfoStr = ['CyanLimit_',num2str(CyanLimit),'_RedLimit_',num2str(RedLimit)];
    
    LimitsDataAssignment.Limits = cell(6,3,4); %ColorName & CI limits(CI1&CI2) & limits
    LimitsDataAssignment.Assign = 1:6;
    
    %NB: "["==">=";   "]"=="<=";       "(" ==">";   ")" == "<"
    LimitsDataAssignment.Limits{1,1,1} = 'Blue';
    LimitsBlue1                        = {'[',0,sqrt(2),')'};
    LimitsBlue2                        = {'[',0,sqrt(2),')'};
    LimitsDataAssignment.Limits{1,2,1} = LimitsBlue1{1};
    LimitsDataAssignment.Limits{1,2,2} = LimitsBlue1{2};
    LimitsDataAssignment.Limits{1,2,3} = LimitsBlue1{3};
    LimitsDataAssignment.Limits{1,2,4} = LimitsBlue1{4};
    LimitsDataAssignment.Limits{1,3,1} = LimitsBlue2{1};
    LimitsDataAssignment.Limits{1,3,2} = LimitsBlue2{2};
    LimitsDataAssignment.Limits{1,3,3} = LimitsBlue2{3};
    LimitsDataAssignment.Limits{1,3,4} = LimitsBlue2{4};
    
    LimitsDataAssignment.Limits{2,1,1} = 'Cyan';
    LimitsCyan1                        = {'(',     0 ,  sqrt(2),']'};
    LimitsCyan2                        = {'[',sqrt(2),CyanLimit,')'};
    LimitsDataAssignment.Limits{2,2,1} = LimitsCyan1{1};
    LimitsDataAssignment.Limits{2,2,2} = LimitsCyan1{2};
    LimitsDataAssignment.Limits{2,2,3} = LimitsCyan1{3};
    LimitsDataAssignment.Limits{2,2,4} = LimitsCyan1{4};
    LimitsDataAssignment.Limits{2,3,1} = LimitsCyan2{1};
    LimitsDataAssignment.Limits{2,3,2} = LimitsCyan2{2};
    LimitsDataAssignment.Limits{2,3,3} = LimitsCyan2{3};
    LimitsDataAssignment.Limits{2,3,4} = LimitsCyan2{4};
    
    LimitsDataAssignment.Limits{3,1,1} = 'Green';
    LimitsGreen1                       = {'(',0,  sqrt(2),']'};
    LimitsGreen2                       = {'[',2,2*sqrt(2),')'};
    LimitsDataAssignment.Limits{3,2,1} = LimitsGreen1{1};
    LimitsDataAssignment.Limits{3,2,2} = LimitsGreen1{2};
    LimitsDataAssignment.Limits{3,2,3} = LimitsGreen1{3};
    LimitsDataAssignment.Limits{3,2,4} = LimitsGreen1{4};
    LimitsDataAssignment.Limits{3,3,1} = LimitsGreen2{1};
    LimitsDataAssignment.Limits{3,3,2} = LimitsGreen2{2};
    LimitsDataAssignment.Limits{3,3,3} = LimitsGreen2{3};
    LimitsDataAssignment.Limits{3,3,4} = LimitsGreen2{4};
    
    LimitsDataAssignment.Limits{4,1,1} = 'Yellow';
    LimitsYellow1                      = {'(',sqrt(2),2*sqrt(2),')'};
    LimitsYellow2                      = {'(',sqrt(2),2*sqrt(2),']'};
    LimitsDataAssignment.Limits{4,2,1} = LimitsYellow1{1};
    LimitsDataAssignment.Limits{4,2,2} = LimitsYellow1{2};
    LimitsDataAssignment.Limits{4,2,3} = LimitsYellow1{3};
    LimitsDataAssignment.Limits{4,2,4} = LimitsYellow1{4};
    LimitsDataAssignment.Limits{4,3,1} = LimitsYellow2{1};
    LimitsDataAssignment.Limits{4,3,2} = LimitsYellow2{2};
    LimitsDataAssignment.Limits{4,3,3} = LimitsYellow2{3};
    LimitsDataAssignment.Limits{4,3,4} = LimitsYellow2{4};
    
    LimitsDataAssignment.Limits{5,1,1} = 'Red';
    LimitsRed1                         = {'(',  sqrt(2),2*sqrt(2),')'};
    LimitsRed2                         = {'(',2*sqrt(2),RedLimit ,']'};
    LimitsDataAssignment.Limits{5,2,1} = LimitsRed1{1};
    LimitsDataAssignment.Limits{5,2,2} = LimitsRed1{2};
    LimitsDataAssignment.Limits{5,2,3} = LimitsRed1{3};
    LimitsDataAssignment.Limits{5,2,4} = LimitsRed1{4};
    LimitsDataAssignment.Limits{5,3,1} = LimitsRed2{1};
    LimitsDataAssignment.Limits{5,3,2} = LimitsRed2{2};
    LimitsDataAssignment.Limits{5,3,3} = LimitsRed2{3};
    LimitsDataAssignment.Limits{5,3,4} = LimitsRed2{4};
    
    LimitsDataAssignment.Limits{6,1,1} = 'White';
    LimitsWhite1                       = {'(',2*sqrt(2),Inf,')'};
    LimitsWhite2                       = {'(',2*sqrt(2),Inf,')'};
    LimitsDataAssignment.Limits{6,2,1} = LimitsWhite1{1};
    LimitsDataAssignment.Limits{6,2,2} = LimitsWhite1{2};
    LimitsDataAssignment.Limits{6,2,3} = LimitsWhite1{3};
    LimitsDataAssignment.Limits{6,2,4} = LimitsWhite1{4};
    LimitsDataAssignment.Limits{6,3,1} = LimitsWhite2{1};
    LimitsDataAssignment.Limits{6,3,2} = LimitsWhite2{2};
    LimitsDataAssignment.Limits{6,3,3} = LimitsWhite2{3};
    LimitsDataAssignment.Limits{6,3,4} = LimitsWhite2{4};
    
    %% median assignment direction
    if(MeanAssignDir==1)
        MedAssignIndices = 1:length(LimitsDataAssignment.Assign); %top to bottom
    else
        MedAssignIndices = length(LimitsDataAssignment.Assign):-1:1; %top to bottom
    end
else
    LimitsDataAssignment.Limits = cell(13,3,4); %ColorName & CI limits(CI1&CI2) & limits
    LimitsDataAssignment.Assign = 1:13;
    
    if(MeanAssignDir==1)
        InfoStr = ['InterleavedQuartersSqrt2_',num2str(size(LimitsDataAssignment.Limits,1)),'_MedAssignUp'];
    else
        InfoStr = ['InterleavedQuartersSqrt2_',num2str(size(LimitsDataAssignment.Limits,1)),'_MedAssignDown'];
    end
    
    %NB: "["==">=";   "]"=="<=";       "(" ==">";   ")" == "<"
    LimitsDataAssignment.Limits{1,1,1} = 'Dark-Blue';
    LimitsDBlue1                       = {'[',0,1/4*sqrt(2),')'};
    LimitsDBlue2                       = {'[',0,1/4*sqrt(2),')'};
    LimitsDataAssignment.Limits{1,2,1} = LimitsDBlue1{1};
    LimitsDataAssignment.Limits{1,2,2} = LimitsDBlue1{2};
    LimitsDataAssignment.Limits{1,2,3} = LimitsDBlue1{3};
    LimitsDataAssignment.Limits{1,2,4} = LimitsDBlue1{4};
    LimitsDataAssignment.Limits{1,3,1} = LimitsDBlue2{1};
    LimitsDataAssignment.Limits{1,3,2} = LimitsDBlue2{2};
    LimitsDataAssignment.Limits{1,3,3} = LimitsDBlue2{3};
    LimitsDataAssignment.Limits{1,3,4} = LimitsDBlue2{4};
    
    LimitsDataAssignment.Limits{2,1,1} = 'Blue';
    LimitsBlue1                        = {'[',0,1/2*sqrt(2),')'};
    LimitsBlue2                        = {'[',0,1/2*sqrt(2),')'};
    LimitsDataAssignment.Limits{2,2,1} = LimitsBlue1{1};
    LimitsDataAssignment.Limits{2,2,2} = LimitsBlue1{2};
    LimitsDataAssignment.Limits{2,2,3} = LimitsBlue1{3};
    LimitsDataAssignment.Limits{2,2,4} = LimitsBlue1{4};
    LimitsDataAssignment.Limits{2,3,1} = LimitsBlue2{1};
    LimitsDataAssignment.Limits{2,3,2} = LimitsBlue2{2};
    LimitsDataAssignment.Limits{2,3,3} = LimitsBlue2{3};
    LimitsDataAssignment.Limits{2,3,4} = LimitsBlue2{4};
    
    LimitsDataAssignment.Limits{3,1,1} = 'Cyan';
    LimitsCyan1                        = {'[',1/4*sqrt(2),3/4*sqrt(2),')'};
    LimitsCyan2                        = {'(',1/4*sqrt(2),3/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{3,2,1} = LimitsCyan1{1};
    LimitsDataAssignment.Limits{3,2,2} = LimitsCyan1{2};
    LimitsDataAssignment.Limits{3,2,3} = LimitsCyan1{3};
    LimitsDataAssignment.Limits{3,2,4} = LimitsCyan1{4};
    LimitsDataAssignment.Limits{3,3,1} = LimitsCyan2{1};
    LimitsDataAssignment.Limits{3,3,2} = LimitsCyan2{2};
    LimitsDataAssignment.Limits{3,3,3} = LimitsCyan2{3};
    LimitsDataAssignment.Limits{3,3,4} = LimitsCyan2{4};
    
    LimitsDataAssignment.Limits{4,1,1} = 'Dark-Green';
    LimitsDGreen1                      = {'[',1/2*sqrt(2),sqrt(2),')'};
    LimitsDGreen2                      = {'(',1/2*sqrt(2),sqrt(2),']'};
    LimitsDataAssignment.Limits{4,2,1} = LimitsDGreen1{1};
    LimitsDataAssignment.Limits{4,2,2} = LimitsDGreen1{2};
    LimitsDataAssignment.Limits{4,2,3} = LimitsDGreen1{3};
    LimitsDataAssignment.Limits{4,2,4} = LimitsDGreen1{4};
    LimitsDataAssignment.Limits{4,3,1} = LimitsDGreen2{1};
    LimitsDataAssignment.Limits{4,3,2} = LimitsDGreen2{2};
    LimitsDataAssignment.Limits{4,3,3} = LimitsDGreen2{3};
    LimitsDataAssignment.Limits{4,3,4} = LimitsDGreen2{4};
    
    LimitsDataAssignment.Limits{5,1,1} = 'Green';
    LimitsGreen1                       = {'[',3/4*sqrt(2),5/4*sqrt(2),')'};
    LimitsGreen2                       = {'(',3/4*sqrt(2),5/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{5,2,1} = LimitsGreen1{1};
    LimitsDataAssignment.Limits{5,2,2} = LimitsGreen1{2};
    LimitsDataAssignment.Limits{5,2,3} = LimitsGreen1{3};
    LimitsDataAssignment.Limits{5,2,4} = LimitsGreen1{4};
    LimitsDataAssignment.Limits{5,3,1} = LimitsGreen2{1};
    LimitsDataAssignment.Limits{5,3,2} = LimitsGreen2{2};
    LimitsDataAssignment.Limits{5,3,3} = LimitsGreen2{3};
    LimitsDataAssignment.Limits{5,3,4} = LimitsGreen2{4};
    
    LimitsDataAssignment.Limits{6,1,1} = 'Dark-Yellow';
    LimitsDYellow1                     = {'[',sqrt(2),3/2*sqrt(2),')'};
    LimitsDYellow2                     = {'(',sqrt(2),3/2*sqrt(2),']'};
    LimitsDataAssignment.Limits{6,2,1} = LimitsDYellow1{1};
    LimitsDataAssignment.Limits{6,2,2} = LimitsDYellow1{2};
    LimitsDataAssignment.Limits{6,2,3} = LimitsDYellow1{3};
    LimitsDataAssignment.Limits{6,2,4} = LimitsDYellow1{4};
    LimitsDataAssignment.Limits{6,3,1} = LimitsDYellow2{1};
    LimitsDataAssignment.Limits{6,3,2} = LimitsDYellow2{2};
    LimitsDataAssignment.Limits{6,3,3} = LimitsDYellow2{3};
    LimitsDataAssignment.Limits{6,3,4} = LimitsDYellow2{4};
    
    LimitsDataAssignment.Limits{7,1,1} = 'Yellow';
    LimitsYellow1                      = {'[',5/4*sqrt(2),7/4*sqrt(2),')'};
    LimitsYellow2                      = {'(',5/4*sqrt(2),7/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{7,2,1} = LimitsYellow1{1};
    LimitsDataAssignment.Limits{7,2,2} = LimitsYellow1{2};
    LimitsDataAssignment.Limits{7,2,3} = LimitsYellow1{3};
    LimitsDataAssignment.Limits{7,2,4} = LimitsYellow1{4};
    LimitsDataAssignment.Limits{7,3,1} = LimitsYellow2{1};
    LimitsDataAssignment.Limits{7,3,2} = LimitsYellow2{2};
    LimitsDataAssignment.Limits{7,3,3} = LimitsYellow2{3};
    LimitsDataAssignment.Limits{7,3,4} = LimitsYellow2{4};
    
    LimitsDataAssignment.Limits{8,1,1} = 'Dark-Red'; %orange? [1 .5 0]
    LimitsDRed1                        = {'[',3/2*sqrt(2),2*sqrt(2),')'};
    LimitsDRed2                        = {'(',3/2*sqrt(2),2*sqrt(2),']'};
    LimitsDataAssignment.Limits{8,2,1} = LimitsDRed1{1};
    LimitsDataAssignment.Limits{8,2,2} = LimitsDRed1{2};
    LimitsDataAssignment.Limits{8,2,3} = LimitsDRed1{3};
    LimitsDataAssignment.Limits{8,2,4} = LimitsDRed1{4};
    LimitsDataAssignment.Limits{8,3,1} = LimitsDRed2{1};
    LimitsDataAssignment.Limits{8,3,2} = LimitsDRed2{2};
    LimitsDataAssignment.Limits{8,3,3} = LimitsDRed2{3};
    LimitsDataAssignment.Limits{8,3,4} = LimitsDRed2{4};
    
    LimitsDataAssignment.Limits{9,1,1} = 'Red';
    LimitsRed1                         = {'[',7/4*sqrt(2),9/4*sqrt(2),')'};
    LimitsRed2                         = {'(',7/4*sqrt(2),9/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{9,2,1} = LimitsRed1{1};
    LimitsDataAssignment.Limits{9,2,2} = LimitsRed1{2};
    LimitsDataAssignment.Limits{9,2,3} = LimitsRed1{3};
    LimitsDataAssignment.Limits{9,2,4} = LimitsRed1{4};
    LimitsDataAssignment.Limits{9,3,1} = LimitsRed2{1};
    LimitsDataAssignment.Limits{9,3,2} = LimitsRed2{2};
    LimitsDataAssignment.Limits{9,3,3} = LimitsRed2{3};
    LimitsDataAssignment.Limits{9,3,4} = LimitsRed2{4};
    
    LimitsDataAssignment.Limits{10,1,1} = 'Red+'; %rgb=[1 1/4 1/4]
    LimitsRedP1                         = {'[',2*sqrt(2),10/4*sqrt(2),')'};
    LimitsRedP2                         = {'(',2*sqrt(2),10/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{10,2,1} = LimitsRedP1{1};
    LimitsDataAssignment.Limits{10,2,2} = LimitsRedP1{2};
    LimitsDataAssignment.Limits{10,2,3} = LimitsRedP1{3};
    LimitsDataAssignment.Limits{10,2,4} = LimitsRedP1{4};
    LimitsDataAssignment.Limits{10,3,1} = LimitsRedP2{1};
    LimitsDataAssignment.Limits{10,3,2} = LimitsRedP2{2};
    LimitsDataAssignment.Limits{10,3,3} = LimitsRedP2{3};
    LimitsDataAssignment.Limits{10,3,4} = LimitsRedP2{4};
    
    LimitsDataAssignment.Limits{11,1,1} = 'Red++'; %rgb=[1 3/5 3/5]
    LimitsRedPP1                        = {'[',9/4*sqrt(2),11/4*sqrt(2),')'};
    LimitsRedPP2                        = {'(',9/4*sqrt(2),11/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{11,2,1} = LimitsRedPP1{1};
    LimitsDataAssignment.Limits{11,2,2} = LimitsRedPP1{2};
    LimitsDataAssignment.Limits{11,2,3} = LimitsRedPP1{3};
    LimitsDataAssignment.Limits{11,2,4} = LimitsRedPP1{4};
    LimitsDataAssignment.Limits{11,3,1} = LimitsRedPP2{1};
    LimitsDataAssignment.Limits{11,3,2} = LimitsRedPP2{2};
    LimitsDataAssignment.Limits{11,3,3} = LimitsRedPP2{3};
    LimitsDataAssignment.Limits{11,3,4} = LimitsRedPP2{4};
    
    LimitsDataAssignment.Limits{12,1,1} = 'RedLimit'; %rgb=[1 4/5 4/5]
    LimitsRedLimit1                     = {'[',10/4*sqrt(2),12/4*sqrt(2),')'};
    LimitsRedLimit2                     = {'(',10/4*sqrt(2),12/4*sqrt(2),']'};
    LimitsDataAssignment.Limits{12,2,1} = LimitsRedLimit1{1};
    LimitsDataAssignment.Limits{12,2,2} = LimitsRedLimit1{2};
    LimitsDataAssignment.Limits{12,2,3} = LimitsRedLimit1{3};
    LimitsDataAssignment.Limits{12,2,4} = LimitsRedLimit1{4};
    LimitsDataAssignment.Limits{12,3,1} = LimitsRedLimit2{1};
    LimitsDataAssignment.Limits{12,3,2} = LimitsRedLimit2{2};
    LimitsDataAssignment.Limits{12,3,3} = LimitsRedLimit2{3};
    LimitsDataAssignment.Limits{12,3,4} = LimitsRedLimit2{4};
    
    LimitsDataAssignment.Limits{13,1,1} = 'White';
    LimitsWhite1                        = {'[',12/4*sqrt(2),Inf,')'};
    LimitsWhite2                        = {'(',12/4*sqrt(2),Inf,')'};
    LimitsDataAssignment.Limits{13,2,1} = LimitsWhite1{1};
    LimitsDataAssignment.Limits{13,2,2} = LimitsWhite1{2};
    LimitsDataAssignment.Limits{13,2,3} = LimitsWhite1{3};
    LimitsDataAssignment.Limits{13,2,4} = LimitsWhite1{4};
    LimitsDataAssignment.Limits{13,3,1} = LimitsWhite2{1};
    LimitsDataAssignment.Limits{13,3,2} = LimitsWhite2{2};
    LimitsDataAssignment.Limits{13,3,3} = LimitsWhite2{3};
    LimitsDataAssignment.Limits{13,3,4} = LimitsWhite2{4};
    
    %% median assignment direction
    if(MeanAssignDir==1)
        MedAssignIndices = 1:length(LimitsDataAssignment.Assign); %top to bottom
    else
        MedAssignIndices = length(LimitsDataAssignment.Assign):-1:1; %top to bottom
    end
end

LimitsDataAssignment.CImapName     = cell(3,1);
LimitsDataAssignment.CImap         = cell(3,1);
LimitsDataAssignment.CIWidthBinsMap= cell(3,1);
LimitsDataAssignment.CIWidthMap    = cell(3,1);
LimitsDataAssignment.MedianRaw     = cell(3,1);

%% CI(1) = MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,1)-1.57*(MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,3)-MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,2))/sqrt(length(squeeze(median(Data,2))));
%% CI(2) = MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,1)+1.57*(MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,3)-MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,2))/sqrt(length(squeeze(median(Data,2))));
% MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,4)  
% MVSscalingTest.Results.Aggregate.MedianQrtCI(   IndSLCenterVox,5)  
% MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,4)   
% MVSscalingTest.Results.MedianSubj.MedianQrtCI(  IndSLCenterVox,5)  
% MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,4)  
% MVSscalingTest.Results.MedianSLight.MedianQrtCI(IndSLCenterVox,5)

%% CIoverlap 
Nmap = 0;
for IndFields = 1:length(FieldsToOutput)
    if(~strcmp(FieldsToOutput{IndFields},'Type')&&~strcmp(FieldsToOutput{IndFields},'DataQuality'))
        Nmap = Nmap+1; %prepare for next
        LimitsDataAssignment.CImapName{Nmap} = FieldsToOutput{IndFields};
        
        %init data
        Median    = squeeze(MVSscalingTest.Results.(FieldsToOutput{IndFields}).MedianQrtCI(:,1));
        CI(:,1:2) = squeeze(MVSscalingTest.Results.(FieldsToOutput{IndFields}).MedianQrtCI(:,4:5));
        CImap         = zeros(length(Median),1);
        CIWidthBinsMap= zeros(length(Median),1);
        
        % get raw data
        LimitsDataAssignment.CIWidthMap{Nmap} = CI(:,2)-CI(:,1);
        LimitsDataAssignment.MedianRaw{Nmap}  = Median;
        
        % go over searchlights to to assignments
        for IndSL = 1:length(Median)
            %do checks CImap
            CImapAssignments = []; %init
            for IndCheck = 1:length(LimitsDataAssignment.Assign)
                if(CheckLimit(CI(IndSL,1),squeeze(LimitsDataAssignment.Limits(IndCheck,2,:)))&&CheckLimit(CI(IndSL,2),squeeze(LimitsDataAssignment.Limits(IndCheck,3,:))))
                    CImapAssignments = [CImapAssignments; LimitsDataAssignment.Assign(IndCheck)];
                end
            end
            if(~isempty(CImapAssignments))
                if(length(CImapAssignments)>1)
                    disp(['Treating "',FieldsToOutput{IndFields},'" IndSL= ',num2str(IndSL,['%0',num2str(ceil(log10(length(Median)))),'g']),'of',num2str(length(Median)),': Possible assignments= [',num2str(CImapAssignments'),']']);
                    for IndCheck = MedAssignIndices
                        if(CheckLimit(Median(IndSL),squeeze(LimitsDataAssignment.Limits(IndCheck,2,:)))||CheckLimit(Median(IndSL),squeeze(LimitsDataAssignment.Limits(IndCheck,3,:))))
                            CImap(IndSL) = LimitsDataAssignment.Assign(IndCheck);
                            break;
                        end
                    end
                    disp(['Treating "',FieldsToOutput{IndFields},'" IndSL= ',num2str(IndSL,['%0',num2str(ceil(log10(length(Median)))),'g']),'of',num2str(length(Median)),': Decide by MEDIAN... ',num2str(CImap(IndSL))]);
                else
                    disp(['Treating "',FieldsToOutput{IndFields},'" IndSL= ',num2str(IndSL,['%0',num2str(ceil(log10(length(Median)))),'g']),'of',num2str(length(Median)),': DONE.']);
                    CImap(IndSL) = CImapAssignments;
                end
            else
                disp(['Treating "',FieldsToOutput{IndFields},'" IndSL= ',num2str(IndSL,['%0',num2str(ceil(log10(length(Median)))),'g']),'of',num2str(length(Median)),' is EMPTY!']);
                for IndCheck = MedAssignIndices
                    if(CheckLimit(Median(IndSL),squeeze(LimitsDataAssignment.Limits(IndCheck,2,:)))||CheckLimit(Median(IndSL),squeeze(LimitsDataAssignment.Limits(IndCheck,3,:))))
                        CImap(IndSL) = LimitsDataAssignment.Assign(IndCheck);
                        break;
                    end
                end
                disp(['Treating "',FieldsToOutput{IndFields},'" IndSL= ',num2str(IndSL,['%0',num2str(ceil(log10(length(Median)))),'g']),'of',num2str(length(Median)),': Decide by MEDIAN... ',num2str(CImap(IndSL))]);
            end
            
            %do checks CIWidthBinsMap
            disp(['Treating "',FieldsToOutput{IndFields},'" IndSL= ',num2str(IndSL,['%0',num2str(ceil(log10(length(Median)))),'g']),'of',num2str(length(Median)),'. [CIWidthBins]']);
            for IndCheck = MedAssignIndices
                if(CheckLimit(LimitsDataAssignment.CIWidthMap{Nmap}(IndSL),squeeze(LimitsDataAssignment.Limits(IndCheck,2,:)))||CheckLimit(LimitsDataAssignment.CIWidthMap{Nmap}(IndSL),squeeze(LimitsDataAssignment.Limits(IndCheck,3,:))))
                    CIWidthBinsMap(IndSL) = LimitsDataAssignment.Assign(IndCheck);
                    break;
                end
            end       
        end
        LimitsDataAssignment.CImap{Nmap}          = CImap;
        LimitsDataAssignment.CIWidthBinsMap{Nmap} = CIWidthBinsMap;
    end
end
disp(' ');
disp('done.');

%% write nifti
%% assignments CImap
NIFTIfiles = cell(length(LimitsDataAssignment.CImapName),1);
for IndMap = 1:length(LimitsDataAssignment.CImapName)
    V_out= MVSscalingTest.SLight.V_SLmask;
    if(V_out.dt(1)<16)
        V_out.dt(1) = 16; %not necessary but save
    end
    V_out.fname = [OutputDir,filesep,'CIassign_',InfoStr,'_',LimitsDataAssignment.CImapName{IndMap},'_',datestr(now,'yyyymmmdd_HHMM'),'.nii'];
    [BaseDir,fname,ext] = fileparts(V_out.fname);
    Y     = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Data  = LimitsDataAssignment.CImap{IndMap};
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
    Y     = reshape(Y,V_out.dim);
    V_out = spm_write_vol(V_out,Y);
    disp(['Writing out "',fname,ext,'" to "',BaseDir,'".']);
    
    LimitsDataAssignment.Vout{IndMap} = V_out;
    NIFTIfiles{IndMap} = V_out.fname; %filepaths for later plot
end

%% assignments CIWidthBinsMap
NIFTIfiles_CIWidthBinsMap = cell(length(LimitsDataAssignment.CImapName),1);
for IndMap = 1:length(LimitsDataAssignment.CImapName)
    V_out= MVSscalingTest.SLight.V_SLmask;
    if(V_out.dt(1)<16)
        V_out.dt(1) = 16; %not necessary but save
    end
    V_out.fname = [OutputDir,filesep,'CIWidthBinsAssign_',InfoStr,'_',LimitsDataAssignment.CImapName{IndMap},'_',datestr(now,'yyyymmmdd_HHMM'),'.nii'];
    [BaseDir,fname,ext] = fileparts(V_out.fname);
    Y     = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    Data  = LimitsDataAssignment.CIWidthBinsMap{IndMap};
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = Data;
    Y     = reshape(Y,V_out.dim);
    V_out = spm_write_vol(V_out,Y);
    disp(['Writing out "',fname,ext,'" to "',BaseDir,'".']);
    
    LimitsDataAssignment.Vout_CIWidthBinsMap{IndMap} = V_out;
    NIFTIfiles_CIWidthBinsMap{IndMap} = V_out.fname; %filepaths for later plot
end

%% raw stuff
NIFTIfilesCIwidth   = cell(length(LimitsDataAssignment.CImapName),1);
NIFTIfilesMedianRaw = cell(length(LimitsDataAssignment.CImapName),1);
for IndMap = 1:length(LimitsDataAssignment.CImapName)
    %CIwidths
    V_CIwidth= MVSscalingTest.SLight.V_SLmask;
    if(V_CIwidth.dt(1)<16)
        V_CIwidth.dt(1) = 16; %not necessary but save
    end
    V_CIwidth.fname = [OutputDir,filesep,'CIwidth_',LimitsDataAssignment.CImapName{IndMap},'_',datestr(now,'yyyymmmdd_HHMM'),'.nii'];
    [BaseDir,fname,ext] = fileparts(V_CIwidth.fname);
    Y     = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    DataCIWidth  = LimitsDataAssignment.CIWidthMap{IndMap};
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = DataCIWidth;
    Y     = reshape(Y,V_CIwidth.dim);
    V_CIwidth = spm_write_vol(V_CIwidth,Y);
    disp(['Writing out "',fname,ext,'" to "',BaseDir,'".']);
    
    LimitsDataAssignment.Vout_CIwidth{IndMap} = V_CIwidth;
    NIFTIfilesCIwidth{IndMap}                 = V_CIwidth.fname; %filepaths for later plot
    
    %MedianRaw
    V_MedianRaw= MVSscalingTest.SLight.V_SLmask;
    if(V_MedianRaw.dt(1)<16)
        V_MedianRaw.dt(1) = 16; %not necessary but save
    end
    V_MedianRaw.fname = [OutputDir,filesep,'MedianRaw_',LimitsDataAssignment.CImapName{IndMap},'_',datestr(now,'yyyymmmdd_HHMM'),'.nii'];
    [BaseDir,fname,ext] = fileparts(V_MedianRaw.fname);
    Y     = zeros(size(MVSscalingTest.SLight.SLmaskRaw));
    DataMedianRaw= LimitsDataAssignment.MedianRaw{IndMap};
    Y(MVSscalingTest.SLight.SLmaskRaw~=0) = DataMedianRaw;
    Y     = reshape(Y,V_MedianRaw.dim);
    V_MedianRaw = spm_write_vol(V_MedianRaw,Y);
    disp(['Writing out "',fname,ext,'" to "',BaseDir,'".']);
    
    LimitsDataAssignment.Vout_MedianRaw{IndMap} = V_MedianRaw;
    NIFTIfilesMedianRaw{IndMap}                 = V_MedianRaw.fname; %filepaths for later plot
end

%% save
disp('saving data to mat-file...');
save([OutputDir,filesep,'CIassign_',InfoStr,'_',datestr(now,'yyyymmmdd_HHMM'),'.mat'],'LimitsDataAssignment');

%% display distribution
figure(42); clf; boxplot(cell2mat(LimitsDataAssignment.CImap(:)'),'labels',LimitsDataAssignment.CImapName); title('Number of voxels & distribution')
DisplayCIniis('Cluster',NIFTIfiles); h=helpdlg('CI assignments','CI assignments plot'); uiwait(h);
DisplayCIniis('Cluster',NIFTIfiles_CIWidthBinsMap); h=helpdlg('CI width-bins assignments','CI width assignments plot'); uiwait(h);
DisplayCIniis('Cluster',NIFTIfilesCIwidth); h=helpdlg('CI width','CI width plot'); uiwait(h);
DisplayCIniis('Cluster',NIFTIfilesMedianRaw); h=helpdlg('Median raw','Median raw plot'); uiwait(h);

%% Done
disp(' ');
disp('DONE.');
disp(' ');

end

%% subfunction
function TFval = CheckLimit(Value,LimitsInfo)
% check if Value lies with in Limits.
%NB: "["==">=";   "]"=="<=";       "(" ==">";   ")" == "<"
switch(LimitsInfo{1})
    case '['
        CompareStr1 = '>=';
    case ']'
        CompareStr1 = '<=';
    case '('
        CompareStr1 = '>';
    case ')'
        CompareStr1 = '<';
    otherwise
        error(['unknown compare string "',LimitsInfo{1},'".']);
end
switch(LimitsInfo{4})
    case '['
        CompareStr2 = '>=';
    case ']'
        CompareStr2 = '<=';
    case '('
        CompareStr2 = '>';
    case ')'
        CompareStr2 = '<';
    otherwise
        error(['unknown compare string "',LimitsInfo{4},'".']);
end

TFval = eval([num2str(Value),CompareStr1,num2str(LimitsInfo{2}),'&&',num2str(Value),CompareStr2,num2str(LimitsInfo{3})]);

end
