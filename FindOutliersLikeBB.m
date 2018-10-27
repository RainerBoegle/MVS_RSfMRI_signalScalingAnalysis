function [DataNew,IsOutlier,TestThisFunction] = FindOutliersLikeBB(Data,ReplaceValue,w,verbose)
% This function finds the outliers in the data vector "Data" according to
% the method introduced by John Tukey for creating box&whiskers-plots.
% (if matrix or array then FindOutliersLikeBB operates on data extracted from rows - BE CAREFUL!)
%
% I.e. all values below lLimit = Qrt1-w*IQR and above uLimit = Qrt3+w*IQR
% will be marked as outliers and replaced according to the command string
% or value (if~ischar(ReplaceValue)) in "ReplaceValue".
% NB: if "ReplaceValue" is empty or not given, then the default is to assign NaNs to the outliers
%
% Qrt1 is the frist quartile, i.e. the point below which 25% and therefore above which 75% of the data lie.
% Qrt3 is the third quartile, i.e. the point below which 75% and therefore above which 25% of the data lie.
% IQR  is the interquartile range, i.e. the length of the distance between the 3rd & 1st Quartile, IQR = Qrt3-Qrt1.
% w is a constant which has the default value of "1.5" just as John Tukey defined it at first.
% NB: If the data is normally distributed (mean==median==mode & symmetric [&kurtosis==3/excess==0])
%     then the range lLimit = Qrt1-w*IQR to uLimit = Qrt3+w*IQR should cover ~99% of the data,
%     leaving the remaining 1%(if present) as outliers.
%
%     This has worked well in many applications so far, therefore we set w == 1.5 per default.
%
% NB2: You can change w by inputing it as the third argument into "FindOutliersLikeBB". 
%      If you want you can even define as a two-element vector, if you wish, to get 
%                    lLimit = Qrt1-w(1)*IQR to uLimit = Qrt3+w(2)*IQR
%
%
% The outputs are
%       DataNew   --> same size as Data but with outliers replaced by ReplaceValue
%       IsOutlier --> NOT a logical vector (but kinda like it) of size Data
%                       with +1 in the places where the high outliers are, i.e those above of the upper limit 
%                     and
%                       with -1 (minus ones) in the places where the low outliers are, i.e those below the lower limit
%                     otherwise it is zero.
%
%
%Usage:
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data,ReplaceValue,w);
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data); %replace with NaN and use w = 1.5.
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data,'NaN',1.1); %replace with NaN and use w = 1.1.
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data,'median',1.1); %replace with nanmedian of data and use w = 1.1.
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data,nanmedian(Data,1),1.1); %equal to the version before: replace with nanmedian of data and use w = 1.1.
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data,'Inf',2); %replace with Inf(infinite) and use w = 2.
%      [DataNew,IsOutlier] = FindOutliersLikeBB(Data,81,2); %replace with 81 (love this number) and use w = 2.
%                            FindOutliersLikeBB(); %Display example including help.
%
%
%V1.0
%Author: Rainer Boegle (Rainer.Boegle@googlemail.com)
%Comment(19.February.2015): initial implementation.

%% Test function output
TestThisFunction = @TestFindOutliersLikeBB;
if(nargin==0) %just output the test function
    DataNew  = [];
    IsOutlier= [];
    if(nargout==0) %Also no outputs? --> run the test
        TestFindOutliersLikeBB(1); %run test and show helpdlg
    end
    return;
end

%% check inputs
%% verbose?
try
    if(verbose==1)
        disp('Will output messages...');
    end
catch
    verbose = 0;
end
%% data
if(~isnumeric(Data))
    error('Data has to be numeric.');
end
%% ReplaceValue
try
    if(ischar(ReplaceValue))
        switch(lower(ReplaceValue))
            case 'nan'
                ReplaceValue = NaN;
            case 'inf'
                ReplaceValue = Inf;
            case {'median','nanmedian'}
                ReplaceValue = nanmedian(Data,1);
            otherwise
                ReplaceValue = NaN;
                disp(['Unknown ReplaceValue Control-string "',ReplaceValue,'".']);
                disp( 'Will default to "NaN".');
        end
    else
        if(isnumeric(ReplaceValue))
            if(length(ReplaceValue)~=1) %not 1? --> error if not size 0f data when using repmat
                if(~isequal(size(Data),size(ReplaceValue))) %not equal size?
                    if(~isequal(size(Data),size(repmat(ReplaceValue,[size(Data,1), size(squeeze(ReplaceValue))])))) %not equal even when using repmat to increase?
                        error('ReplaceValue either has to be length 1 or size according to Data or size of reduced Data, i.e. can be same size as nanmedian(Data,1) would be.');
                    else
                        %size it up then using repmat
                        ReplaceValue = repmat(ReplaceValue,size(Data,1),1);
                    end
                end    
            end
        else
            error('ReplaceValue has to be a command string or numeric value.');
        end
    end         
catch
    ReplaceValue = NaN;
    if(verbose)
        disp('ReplaceValue will default to "NaN".');
    end
end
%% w
try
    if(isnumeric(w))
        if(length(w)==1)
            w = ones(2,1).*w;
        else
            if(length(w)~=2)
                error('"w" has to be a numeric value (or vector of length 2).');
            end
        end
    else
        error('"w" has to be a numeric value.');
    end
catch
    w = 1.5.*ones(2,1);
    if(verbose)
        disp('"w" will default to "1.5" (Tukey''s value leading to an expected coverage of 99% the of data, if data is normally distributed).');
    end
end

%% init outputs
DataNew  = Data;  %init all original data
IsOutlier= zeros(size(Data));%init none are outliers

%% apply
Qrt1     = repmat(quantile(Data,.25,1),size(Data,1),1);
Qrt3     = repmat(quantile(Data,.75,1),size(Data,1),1);
IQR      = Qrt3-Qrt1;
uLimit   = Qrt3+w(2).*IQR;
lLimit   = Qrt1-w(1).*IQR;
IsOutlier(Data>=uLimit) = +1;
IsOutlier(Data<=lLimit) = -1;
if(length(ReplaceValue)==1) %one for all
    DataNew(IsOutlier~=0) = ReplaceValue;
else
    %must be one value per reduced data matrix
    DataNew(IsOutlier~=0) = ReplaceValue(IsOutlier~=0);
end

end


%% TestThisFunction
function [Data,IsOutlier,IsOutlierALL] = TestFindOutliersLikeBB(varargin)
%Usage(TestFindOutliersLikeBB): 
%     [Data,IsOutlier,IsOutlierALL] = TestFindOutliersLikeBB();
%NB:
%      Data = randn(100,4);
%      [~,IsOutlier]    = FindOutliersLikeBB(Data);
%      [~,IsOutlierALL] = FindOutliersLikeBB(Data(:));
if(nargin==1)
    ShowHelp = varargin{1};
else
    ShowHelp = 0; %Don't show.
end
if(ShowHelp==1)
    h=helpdlg({'No Inputs and No Outputs therefore displaying Test-Function.'; ' '; 'If you want the pointer/handle of the test function,'; 'then call the function like this:'; ' '; '"[~,~,TestFindOutliersLikeBB]=FindOutliersLikeBB();"'; 'You can then run the test (which will follow now) via "[Data,IsOutlier,IsOutlierALL] = TestFindOutliersLikeBB();"'; 'NB:'; 'Data = randn(100,4);'; '[~,IsOutlier]    = FindOutliersLikeBB(Data);'; '[~,IsOutlierALL] = FindOutliersLikeBB(Data(:));'; ' '; 'And you can get the help for this function via:'; '"help FindOutliersLikeBB>TestFindOutliersLikeBB"'; ' '; ' '; 'The test you will see is based on data generated using the command randn(100,4), i.e. 100 samples in each of four columns drawn from a gaussian distribution.'; 'Therefore with default w=1.5 there should be only very few outliers detected, -on average about 1-4 for the total data and maybe 2-5 per column.'},'TestFunction');
    uiwait(h);
end
Data = randn(100,4);
[DataNew, IsOutlier]  = FindOutliersLikeBB(Data);
[DataNew2,IsOutlierALL] = FindOutliersLikeBB(Data(:));
IsOutlierALL = reshape(IsOutlierALL,size(Data));
clear DataNew DataNew2

if(ShowHelp~=-1) %make it possible that only the data is generated without display.
    figure(81); clf;
    subplot(2,2,1); imagesc(Data.*abs(IsOutlier),   [-(max(abs(Data(:)))+(max(abs(Data(:)))/100)),max(abs(Data(:)))+(max(abs(Data(:)))/100)]); title('Considering the columns individually'); colorbar;
    subplot(2,2,3); imagesc(Data.*abs(IsOutlierALL),[-(max(abs(Data(:)))+(max(abs(Data(:)))/100)),max(abs(Data(:)))+(max(abs(Data(:)))/100)]); title('Considering ALL the data together & reshape back for comparisons.'); colorbar;
    subplot(2,2,2); boxplot(Data); title('Considering the columns individually');
    subplot(2,2,4); boxplot(Data(:)); title('Considering ALL the data together.');
end

end