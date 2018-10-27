function FileListInfo = Determine_SubjNrScannerRun_from_FileList(FileListPath)
% This function analyses the file-list ".filelist" such that subject nr, scanner nr and run nr are
% determined.
%
% Currently the search strings are hard coded.
%
%Usage:
%       FileListInfo = Determine_SubjNrScannerRun_from_FileList(); %all inputs are asked for and search strings are hard coded. (Will change in future version)
%       FileListInfo = Determine_SubjNrScannerRun_from_FileList(FileListPath); %FileListPath points to ".filelist" location and search strings are hard coded. (Will change in future version)
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment(21.July.2015): initial implementation based on test script.

%% get .filelist
try
    if(~exist('FileListPath','var'))
        FileListPath = spm_select(1,'.filelist','Select .filelist to determine SubjNr, Scanner & Run...');
    end
catch CATCH_FileListPath
    disp_catch(CATCH_FileListPath,[mfilename,'>>Inputs::FileListPath'],'CATCH_FileListPath');
    FileListPath = spm_select(1,'.filelist','Select .filelist to determine SubjNr, Scanner & Run...');
end

%% load .filelist
FileList = importdata(FileListPath);

%% key strings to look out for
SubjStr    = 'Subj_p\d\d\d'; %only one possibility for subject
ScannerStr = {'AERA','GE'}; %two possibility for scanner
RunStr     = {'rsfMRI_1','rsfMRI_2'}; %possibility for runs

%% apply regexp (NB: for loops are to unfold cell of cells into cellstr of length(FileList))
SubjNrs_Strings  = regexp(FileList, SubjStr, 'match');
for Ind = 1:length(SubjNrs_Strings)
    SubjNum = SubjNrs_Strings{Ind};
    count = 0; Err = 0; %try maximally 10 time then set error to Err=1 for allerting user
    while(iscell(SubjNum))
        SubjNrs_Strings{Ind} = SubjNum{1};
        SubjNum = SubjNrs_Strings{Ind};
        count = count+1;
        if(count>10)
            Err = 1;
            break;
        end
    end
    if(Err)
        error(['problem with SubjNum assignment at Ind ',num2str(Ind)]);
    end
end
        
Scanners_Strings = regexp(FileList, [ScannerStr{1},'|',ScannerStr{2}], 'match');
for Ind = 1:length(Scanners_Strings)
    ScannerNum = Scanners_Strings{Ind};
    count = 0; Err = 0; %try maximally 10 time then set error to Err=1 for allerting user
    while(iscell(ScannerNum))
        Scanners_Strings{Ind} = ScannerNum{1};
        ScannerNum = Scanners_Strings{Ind};
        count = count+1;
        if(count>10)
            Err = 1;
            break;
        end
    end
    if(Err)
        error(['problem with ScannerNum assignment at Ind ',num2str(Ind)]);
    end
end

Runs_Strings     = regexp(FileList, [RunStr{1},'|',RunStr{2}], 'match');
for Ind = 1:length(Runs_Strings)
    RunNum = Runs_Strings{Ind};
    count = 0; Err = 0; %try maximally 10 time then set error to Err=1 for allerting user
    while(iscell(RunNum))
        Runs_Strings{Ind} = RunNum{1};
        RunNum = Runs_Strings{Ind};
        count = count+1;
        if(count>10)
            Err = 1;
            break;
        end
    end
    if(Err)
        error(['problem with RunNum assignment at Ind ',num2str(Ind)]);
    end
end

%% make assigment matrix

SubjNrScannerNrRunNr = zeros(length(FileList),3);
for IndFile = 1:length(FileList)
    SubjNum = regexp(SubjNrs_Strings{IndFile}, '\d*', 'match');
    if(iscell(SubjNum))
        SubjNrScannerNrRunNr(IndFile,1) = eval(SubjNum{1});
    else
        if(ischar(SubjNum))
            SubjNrScannerNrRunNr(IndFile,1) = eval(SubjNum);
        end
    end
    if(strcmp(Scanners_Strings{IndFile},ScannerStr{1}))
        SubjNrScannerNrRunNr(IndFile,2) = 1;
    else
        SubjNrScannerNrRunNr(IndFile,2) = 2;
    end
    RunNum = regexp(Runs_Strings{IndFile}, '\d*', 'match');
    if(iscell(RunNum))
        SubjNrScannerNrRunNr(IndFile,3) = eval(RunNum{1});
    else
        if(ischar(RunNum))
            SubjNrScannerNrRunNr(IndFile,3) = eval(RunNum{1});
        end
    end
end

UniqueSubjNrs    = unique(SubjNrScannerNrRunNr(:,1));
UniqueScannerNrs = unique(SubjNrScannerNrRunNr(:,2));
UniqueRunNrs     = unique(SubjNrScannerNrRunNr(:,2));

%% output
FileListInfo.FileListPath         = FileListPath;
FileListInfo.FileList             = FileList;
FileListInfo.SubjNrScannerNrRunNr = SubjNrScannerNrRunNr;
FileListInfo.UniqueSubjNrs        = UniqueSubjNrs;
FileListInfo.UniqueScannerNrs     = UniqueScannerNrs;
FileListInfo.UniqueRunNrs         = UniqueRunNrs;
FileListInfo.NSubjs               = length(UniqueSubjNrs);
FileListInfo.N_MRIs               = length(UniqueScannerNrs);


end