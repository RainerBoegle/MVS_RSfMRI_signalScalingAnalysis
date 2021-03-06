%% get data
load('AAL_Labels.mat');

%% pick one 
FirstChoice = cell(length(Labels),1);
for Ind = 1:length(Labels)
    FirstChoice{Ind} = Labels(Ind).Name;
end
[SelectionInd_Split,ok] = listdlg('ListString',FirstChoice,'SelectionMode','single');
if(~ok)
    return;
else
    [SelectionInd_Label,ok] = listdlg('ListString',Labels(SelectionInd_Split).Labels,'SelectionMode','single');
    if(~ok)
        return;
    else
        [SelectionInd_Split2,ok] = listdlg('ListString',FirstChoice,'SelectionMode','single'); %the one that should be matched.
        if(~ok)
            return;
        end
    end
end

% MatchStr = regexptranslate('wildcard',regexprep(Labels(SelectionInd_Split).Labels{SelectionInd_Label},'_','|'));
% MatchStrOrg = Labels(SelectionInd_Split).Labels{SelectionInd_Label};
% Indices     = strfind(regexprep(MatchStrOrg,'_',' '),' ');
% MatchStrTmp    = cell(length(Indices),1);
% StartInd    = 1;
% MatchStr    = [];
% for Ind = 1:length(Indices)
%     MatchStrTmp{Ind} = MatchStrOrg(StartInd:(Indices(Ind)-1));
%     if(Ind<length(Indices))
%         StartInd = Indices(Ind)+1;
%     end
%     if(Ind>1)
%         MatchStr = [MatchStr,'|',MatchStrTmp{Ind}];
%     else
%         MatchStr = [MatchStr,MatchStrTmp{Ind}];
%     end
% end
% % MatchStr = ['[',MatchStr,']'];
MatchStrOrg = Labels(SelectionInd_Split).Labels{SelectionInd_Label};
MatchStr{1} = regexprep(MatchStrOrg,'_',' '); %first, look for FULL STRING WITHOUT UNDERSCORE BUT WHITESPACE
Indices     = strfind(MatchStr{1},' ');
if(~isempty(Indices))
    MatchStr{2} = MatchStrOrg(1:(Indices(1)-1)); %second, look for FIRST PART OF STRING
end

FoundOneAtLeast = 0;
ReportCellStr   = {'Report: '};
for Ind = 1:length(Labels(SelectionInd_Split2).Labels)
    MatchInd = 1;
    [startIndex] = regexpi(Labels(SelectionInd_Split2).Labels{Ind}, MatchStr{1});
    if(length(MatchStr)>1)
        if(isempty(startIndex)) %full string not found so let's try the shorter version
            MatchInd = 2;
            [startIndex] = regexpi(Labels(SelectionInd_Split2).Labels{Ind}, MatchStr{2});
        end
    end
%     startIndex = strfind(Labels(SelectionInd_Split2).Labels{Ind}, MatchStr);
    if(~isempty(startIndex))
        FoundOneAtLeast = FoundOneAtLeast+1;
        ReportCellStr{1+FoundOneAtLeast} = [MatchStr{1},'(',num2str(MatchInd),') =?= ',Labels(SelectionInd_Split2).Labels{Ind}];
    end
end
if(~FoundOneAtLeast)
    ReportCellStr{2} = ['No match found for ',MatchStrOrg];
end

H = helpdlg(ReportCellStr,'Report');
uiwait(H);
clear MatchStr

disp('Done');
disp(' ');