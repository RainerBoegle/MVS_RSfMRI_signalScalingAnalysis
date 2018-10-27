% Computes the appropriate threshold to use over the pvalues
% resulting from multiple tests of a null hypothesis,
% in order to have a particular False Discovery Rate (FDR).
%
% This is the proportion of false positives (incorrect rejections
% of the null hypothesis) among those tests for which the null
% hypothesis is rejected. In voxel terms, this is the expected proportion
% of voxels declared active (using the threshold) that is not
% active. 
%
% Note that "expected" here means that, were one to replicate the
% experiment many times, the average FDR over those replications
% would be no bigger than the value for which the threshold was
% determined. For ANY PARTICULAR data analysis, the actual FDR
% might be LARGER than that value.
%
% For more details, see "Controlling FDR in Neuroimaging", by
% Genovese, Lazar and Nichols (in NeuroImage).
%
% Nichols' FDR web page says, apropos of this:
%
% Note on the preprint: The preprint reflected the current results at the time of submission. 
% Since then, the independence FDR result of 1995 (BH 1995; in the above paper, "c(V)=1" on pg 9) 
% has been shown to be more general (BY 2001). 
% In particular, a technical condition known as positive regression dependency is sufficient to use the 1995 result.
% Examples of this condition include multivariate normal data where the covariance matrix has all positive elements;
% this seems a reasonable assumption for imaging, and hence the code above uses the c(V)=1 result. 
%
% so one should use the thresholdID output (which uses c(v)=1)
%
% In:
% - a vector of pvalues or a matrix whose columns are vectors of pvalues
% - a false discovery rate (e.g. 0.01, 0.05)
%
% Out:
% - the two thresholds described in the paper
%   - thresholdID - threshold based on independence or positive dependence
%   - thresholdN  - nonparametric pvalue threshold
%   - FDRstruct   - a structure/object that contains the important parameters of the computation and pointers to useful functions,
%   - TestFun     - pointer to a useful test function
%
% Dependencies:
%
% History
% - 2009 November 04 - created                         - fpereira@princeton.edu
% - 2015 February 17 - commented&adjusted              - Rainer.Boegle@googlemail.com
% - 2015 March    12 - 1. not significant return NaN   - Rainer.Boegle@googlemail.com
%                      2. if some p==0 check & adjust thresholds such that threshold is nonzero.
%
%NB(shortened and changed from Wikipedia on FDR):
%    The false discovery rate concept is less conservative, but arguably more appropriate for 
%    IDENTIFYING THE IMPORTANT FEW FROM THE TRIVIAL MANY EFFECTS that are tested.
%
%    I.e. when testing at every voxel, we do not really expect that all voxels should be
%    treated equally likely and that we should minimize the probability of getting EVEN ONE
%    falsely discovered effect/voxel due to multiple testing as done when using FWE-correction.
%
%	 What we want instead is to say that we have maximally Q%(FDR=0.01 --> 1%) falsely discovered
%    i.e. declared "active" voxels under all "ACTIVE" voxels.
%    ON AVERAGE, IF WE REPEATED EVERYTHING IN THE SAME WAY, HOWEVER IN SOME
%    INDIVIDUAL CASE THAT YOU MIGHT BE LOOKING AT IT MIGHT BE MORE.
%
%NB(one more thing... -final): 
%           A significant voxel at p-value (e.g. 0.05) should be considered as "one that is worthy of a second look".
%           "Let's have a closer look at these by testing something else at these locations or see what other literature
%           says about these or publish it for others to replicate this."
%
%           Generally for hypothesis testing at every voxel (or searchlight...) we have to consider that we might find
%           a large enough effect just by chance.
%           
%           What we want is not to restrict the probability that there are no such cases where discovery is just by chance
%           (false positive) for any one of the tests --> that would be FWER,
%           but we want to restrict the probability/the rate/the fraction of false discoveries in the number of discoveries.
%           (The number of discoveries here is itself a small fraction of the tests,
%            i.e. false positives are still very few of the total number of tests.)
%
%           We first rank all voxels according to their p-values from 1 to N.
%           We start out with a p-threshold at Bonferroni-level (p=q/N-tests) and increase the p-threshold
%           with a slope of q/C(V) from the most significant (ranked as number 1) to the least significant test (ranked as number N)
%           and see if there are some below this line. 
%           We then pick the largest and reject the null hypothesis for all that are more significant than this p-value.
%
%           Even if the most significant voxel does not get significant at Bonferroni-level, some "later voxel", 
%           i.e. a less significant one, might get significant, because of the "increased weighting". 
%           The "increased weighting" signifies that "some voxels have been tested already"
%           and therefore, there are less tests (less than N) still to perform.
%
%           Note that even though this does not seem to make sense, when one is used to view results from the
%           perspective of FWE-type testing, it makes sense overall when the perspective is moved to the proportion
%           of false discoveries.
%
%Usage:   
%      [thresholdID,thresholdN,FDRstruct,TestFun] = computePthresForFDR(pVals,qFDR);
%                                                   computePthresForFDR(); %for a demo
%
%Try the plotting function in the FDRstruct to see how your p-values are behaved.
%      FDRstruct.plot();
%       

function [thresholdID,thresholdN,FDRstruct,TestFun] = computePthresForFDR( varargin )
  %% no input no output --> run rest
  if(nargin==0&&nargout==0)
      testThis();
      thresholdID = [];
      thresholdN  = [];
      FDRstruct   = [];
      return;
  end
  
  %% assign test function
  TestFun = @testThis;
  
  %% check inputs and do calculations if inputs are given otherwise just output the test function and empty outputs
  l = length(varargin);
  if l < 1
    fprintf('syntax: computeFDR(<pvalue vector or matrix with pvalue vector columns>,<fdr>\n');
    fprintf('please see code comment header for more info\n');
    thresholdID = [];
    thresholdN  = [];
    FDRstruct   = [];
    return;
  end

  pValues = varargin{1}; [n,ncols] = size(pValues);
  fdr     = varargin{2};

  if n == 1; pValues = pValues'; n = ncols; ncols = 1; end %flip inputs if only one row --> dangerous, because this could lead to an error later! But that would be because the user made a mistake with the input being wrong shape and because it is the user's fault, they don't care. ;)
  
  
  [pValuesSorted,sortedIndices] = sort(pValues,1); %from smallest(most significant) to largest(least significant) p-value

  constantID = 1; %this means c(V)=1 --> see assumptions above concerning NORMALITY!
  constantN  = sum(1./(1:n)); %NON-PARAMETRIC assumption, i.e. weigh them with less and less likelyhood from most significant(smallest p-value) to least significant(largest p-value). --> this includes therefore the overall number of voxels "n", but is less stringent than Bonferroni (FWE) correction, which is just weighting all with 1/n
  
  ratioID = fdr / (n*constantID); %This would be a Bonferroni correction level, but note that later "pos"-vector will correct for the number of tests already applied, so it will be less stringent than Bonferroni correction.
  ratioN  = fdr / (n*constantN);  %This would be MORE THAN a Bonferroni correction level, but note that later "pos"-vector will correct for the number of tests already applied, so it will be less stringent than Bonferroni correction.

  pos = repmat((1:n)',1,ncols); %weighting for most to least significant --> this will be used to "increase" the significance of the voxels that are tested as we go along the listing of most significant to least significant.
  maskID = (pValuesSorted ./ pos) <= ratioID; %those that are below the ratio (assumed independence or positive dependence) are significant. NB: ratioID is at the Bonferroni-level, BUT division by "pos" includes a weighting that is not present in Bonferroni-correction (OR can be considered to be "1" for each of them). I.e. the most significant get's compared to a level like for FWE, because we used a factor of pos(1)=1 (pos(1)=="most significant"), but the least significant (and analog inbetween) is weighted with a factor of "n"(or 1/n to be precise), THIS MEANS that the (k>1)'th-significant voxel is thresholded with a corrected threshold that includes the assumption that the ones before were significant (because sorted!) and therefore this voxel does NOT need to be corrected for the total number of tests but the number of tests that are left to do. 
  [maxrowID] = max(maskID .* pos,[],1); %gonna be voxel "j" somewhere between 1:n (if any are significant at all).
  maskN = (pValuesSorted ./ pos) <= ratioN; %analog to above, BUT ratio is adjusted such that threshold is more stringent, i.e. dropping the independence or positive dependence (for multivariate-normal data) assumption.
  [maxrowN] = max(maskN .* pos,[],1); %gonna be voxel "h" somewhere between 1:n (if any are significant at all).

  thresholdID = NaN(1,ncols); %init as not significant in the sense that none of the data p-values if sufficiently low enough.
  thresholdN  = NaN(1,ncols); %init as not significant in the sense that none of the data p-values if sufficiently low enough.

  for iv = 1:ncols
    % if none survive the threshold is NaN
    if maxrowID(iv); thresholdID(iv) = pValuesSorted(maxrowID(iv),iv); end %assign significance value for ID case. (independent or positive dependent)
    if maxrowN(iv);  thresholdN(iv)  = pValuesSorted(maxrowN(iv),iv);  end %assign significance value for NON-Parametric case. (independent or positive dependent)
    if(thresholdID(iv)==0) %special case --> pthres should be FWE
        thresholdID(iv) = fdr/n;
    end
    if(thresholdN(iv)==0) %special case --> pthres should be FWE
        thresholdN(iv) = fdr/n;
    end
  end

%% collect some of the data in here for plotting if wished for
FDRstruct.plot          = @plot_computePthresForFDR; %function handle/pointer for plotting
FDRstruct.help          = @help_computePthresForFDR; %function handle/pointer for help
FDRstruct.fdr           = fdr;
FDRstruct.pFWEexpected  = fdr/n;
FDRstruct.pValuesSorted = pValuesSorted;
FDRstruct.constantN     = constantN;
FDRstruct.ratioID       = ratioID;  
FDRstruct.ratioN        = ratioN;
FDRstruct.Rank          = pos;
FDRstruct.maxrowID      = maxrowID;
FDRstruct.maxrowN       = maxrowN;
FDRstruct.thresholdID   = thresholdID;
FDRstruct.thresholdN    = thresholdN;
FDRstruct.Info.InfoTxt  = {'This struct contains all the important data generated when using "computePthresForFDR".'; 'If you want to plot the distribution of SORTED p-values and relationships used for generating the thresholds at FDR then use "FDRstruct.plot();".'};
FDRstruct.Info.InfoTitle= 'Info for FDRstruct';


%% nested function "plot_computePthresForFDR" for plotting the data   
function [h_fig] = plot_computePthresForFDR()
    h_fig = cell(size(FDRstruct.pValuesSorted,2),1);
    for IndData = 1:size(FDRstruct.pValuesSorted,2)
        %get data
        Ntotal  = size(FDRstruct.pValuesSorted,1);
        pVals   = FDRstruct.pValuesSorted(:,IndData);
        Rank    = FDRstruct.Rank(:,IndData);
        RankByN = FDRstruct.Rank(:,IndData)./size(FDRstruct.Rank,1);
        SlopeID = FDRstruct.ratioID*Ntotal;
        SlopeN  = FDRstruct.ratioN*Ntotal;
        P_ID    = FDRstruct.thresholdID(IndData);
        P_NP    = FDRstruct.thresholdN(IndData);
        %% do an extra check in case p-values contain zeros
        %% inform user
        if(any(pVals==0)||any(pVals<=FDRstruct.pFWEexpected))
            if(FDRstruct.pFWEexpected==P_ID||FDRstruct.pFWEexpected==P_NP)
                h=helpdlg('Note that, although some tests are significant, that the structure of the data does not allow a threshold that is better than FWE.','FWE still the best possible.');
                uiwait(h);
            end
        end
        %% adjust plotting
        if(any(pVals==0)) %this will produce errors in the plot --> make them nonzero
            if(any(pVals(pVals~=0)<eps('double'))) %make sure we do not get any funny behavior
                Odiff = abs(log10(min(pVals(pVals~=0))))-abs(log10(eps('double'))); %order of magnitude difference
                if(Odiff<=1)
                    pVals(pVals==0) = eps('double')/100; %enough separation to be suggestive enough
                    DataInfo = '[0=eps/100]';
                else
                    Odiff = ceil(Odiff)+1;
                    pVals(pVals==0) = eps('double')/10.^Odiff; %enough separation to be suggestive enough
                    DataInfo = ['[0=eps/',num2str(10.^Odiff),']'];
                end
            else
                pVals(pVals==0) = eps('double')/10; %enough separation to be suggestive enough
                DataInfo = '[0=eps/10]';
            end
        else
            DataInfo = [];
        end
        %% inform user
        if(~isempty(DataInfo))
            h=helpdlg('Some of the p-values are zero, i.e. extremely significant (or below precision of the used test). In this case the plotting will be adjusted by replacing all p-values which are zero with the "eps" function in relation to the next non-zero p-values, by considering the orders of magnitude separation. This is purely for display convenience and is not applied to the data.','Plotting adjusted.');
            uiwait(h);
        end
        
        %% plot & loglog-plot
        h_fig{IndData} = figure(); clf;
        if(size(FDRstruct.pValuesSorted,2)==1)
            TitleStr = 'sorted P-Values and threshold determination curves';
        else
            TitleStr = {'sorted P-Values and threshold determination curves'; ['Column ',num2str(IndData)]};
        end
        subplot(1,2,1);
        plot(RankByN,pVals,'kx','LineWidth',2); title(TitleStr); hold on
        plot(RankByN,RankByN.*SlopeID,'r-');
        plot(RankByN,RankByN.*SlopeN,'r--');
        plot(RankByN,ones(size(RankByN)).*FDRstruct.pFWEexpected,'m-');
        xlabel({'Rank/N'; ['(N=',num2str(Ntotal),')']});
        ylabel('P-Value');
        legend({['sorted P-Values(Data',DataInfo,')'],'ID: (Rank/N)*Q','NP: (Rank/N)*Q/C(N)','Expected pFWE(Q)'},'Location','NorthWest');
        plot(RankByN,pVals,'-','Color',[.5 .5 .5],'LineWidth',2); %linear trendline
        if(P_ID==0&&P_NP==0&&all(FDRstruct.pValuesSorted(:,IndData)~=0))
            text(.5,max(pVals(:))*.5,{upper(['not significant! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[1 0 1],'EdgeColor','black');
        else
            if(P_ID==0&&all(FDRstruct.pValuesSorted(:,IndData)~=0)) %should be not possible
                plot(RankByN(1:FDRstruct.maxrowN(IndData)),P_NP.*ones(length(1:FDRstruct.maxrowN(IndData)),1),'--','Color',gray(1));
                plot(RankByN(FDRstruct.maxrowN(IndData)).*ones(10,1),P_NP.*linspace(1,0,10)','--','Color',gray(1));
                text(RankByN(FDRstruct.maxrowN(IndData)),P_NP,'\leftarrow pNP');
                text(.5,max(pVals(:))*.5,{upper(['only significant for the NP-case! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[.75 1 0],'EdgeColor','black');
            else
                if(P_NP==0&&all(FDRstruct.pValuesSorted(:,IndData)~=0))
                    plot(RankByN(1:FDRstruct.maxrowID(IndData)),P_ID.*ones(length(1:FDRstruct.maxrowID(IndData)),1),'-','Color',gray(1));
                    plot(RankByN(FDRstruct.maxrowID(IndData)).*ones(10,1),P_ID.*linspace(1,0,10)','-','Color',gray(1));
                    text(RankByN(FDRstruct.maxrowID(IndData)),P_ID,'\leftarrow pID');
                    text(.5,max(pVals(:))*.5,{upper(['only significant for the ID-case! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[.25 1 0],'EdgeColor','black');
                else
                    plot(RankByN(1:FDRstruct.maxrowID(IndData)),P_ID.*ones(length(1:FDRstruct.maxrowID(IndData)),1),'-','Color',gray(1));
                    plot(RankByN(FDRstruct.maxrowID(IndData)).*ones(10,1),P_ID.*linspace(1,0,10)','-','Color',gray(1));
                    text(RankByN(FDRstruct.maxrowID(IndData)),P_ID,'\leftarrow pID');
                    
                    plot(RankByN(1:FDRstruct.maxrowN(IndData)),P_NP.*ones(length(1:FDRstruct.maxrowN(IndData)),1),'--','Color',gray(1));
                    plot(RankByN(FDRstruct.maxrowN(IndData)).*ones(10,1),P_NP.*linspace(1,0,10)','--','Color',gray(1));
                    text(RankByN(FDRstruct.maxrowN(IndData)),P_NP,'\leftarrow pNP');
                    text(.5,max(pVals(:))*.5,{upper(['both ID-&NP-case are significant! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[0 1 0],'EdgeColor','black');
                end
            end
        end
        hold off
        
        subplot(1,2,2);
        loglog(RankByN,pVals,'kx','LineWidth',2); title(TitleStr); hold on
        loglog(RankByN,RankByN.*SlopeID,'r-');
        loglog(RankByN,RankByN.*SlopeN,'r--');
        loglog(RankByN,ones(size(RankByN)).*FDRstruct.pFWEexpected,'m-');
        xlabel({'log(Rank/N)'; ['(N=',num2str(Ntotal),')']});
        ylabel('log(P-Value)');
        legend({['sorted P-Values(Data',DataInfo,')'],'ID: (Rank/N)*Q','NP: (Rank/N)*Q/C(N)','Expected pFWE(Q)'},'Location','NorthWest');
        loglog(RankByN,pVals,'-','Color',[.5 .5 .5],'LineWidth',2); %gray linear trend line
        if(P_ID==0&&P_NP==0)
            text(log10(.5),log10(max(pVals(:))*.5),{upper(['not significant! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[1 0 1],'EdgeColor','black');
        else
            if(P_ID==0) %should be not possible
                loglog(RankByN(1:FDRstruct.maxrowN(IndData)),P_NP.*ones(length(1:FDRstruct.maxrowN(IndData)),1),'--','Color',gray(1));
                loglog(RankByN(FDRstruct.maxrowN(IndData)).*ones(10,1),linspace(P_NP,pVals(1)*0.90,10)','--','Color',gray(1));
                text(log10(RankByN(FDRstruct.maxrowN(IndData))),log10(P_NP),'\leftarrow pNP');
                text(log10(.5),log10(max(pVals(:))*.5),{upper(['only significant for the NP-case! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[.75 1 0],'EdgeColor','black');
            else
                if(P_NP==0)
                    loglog(RankByN(1:FDRstruct.maxrowID(IndData)),P_ID.*ones(length(1:FDRstruct.maxrowID(IndData)),1),'-','Color',gray(1));
                    loglog(RankByN(FDRstruct.maxrowID(IndData)).*ones(10,1),linspace(P_ID,pVals(1)*0.90,10)','-','Color',gray(1));
                    text(log10(RankByN(FDRstruct.maxrowID(IndData))),log10(P_ID),'\leftarrow pID');
                    text(log10(.5),log10(max(pVals(:))*.5),{upper(['only significant for the ID-case! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[.25 1 0],'EdgeColor','black');
                else
                    loglog(RankByN(1:FDRstruct.maxrowID(IndData)),P_ID.*ones(length(1:FDRstruct.maxrowID(IndData)),1),'-','Color',gray(1));
                    loglog(RankByN(FDRstruct.maxrowID(IndData)).*ones(10,1),linspace(P_ID,pVals(1)*0.90,10)','-','Color',gray(1));
                    text(log10(RankByN(FDRstruct.maxrowID(IndData))),log10(P_ID),'\leftarrow pID');
                    
                    loglog(RankByN(1:FDRstruct.maxrowN(IndData)),P_NP.*ones(length(1:FDRstruct.maxrowN(IndData)),1),'--','Color',gray(1));
                    loglog(RankByN(FDRstruct.maxrowN(IndData)).*ones(10,1),linspace(P_NP,pVals(1)*0.90,10)','--','Color',gray(1));
                    text(log10(RankByN(FDRstruct.maxrowN(IndData))),log10(P_NP),'\leftarrow pNP');
                    text(log10(.5),log10(max(pVals(:))*.5),{upper(['both ID-&NP-case are significant! (q=',num2str(FDRstruct.fdr),')']); ['(pID=',num2str(P_ID),')']; ['(pNP=',num2str(P_NP),')']; ['(pFWE=',num2str(FDRstruct.pFWEexpected),')']},'FontUnits','normalized','FontWeight','bold','FontSize',.05,'HorizontalAlignment','center','BackgroundColor',[0 1 0],'EdgeColor','black');
                end
            end
        end
        hold off
        %% plotyy for stats determination -todo
%         [ax,h1,h2] = plotyy(x1,y1,x2,y2)
%         hold(ax(1))
%         plot(ax(1),x,y)
%         hold(ax(2))
%         plot(ax(2),x,y)        
    end
% %also
% [ax h1 h2]=plotyy([0,1],[0,1],[0,1],[0,1]); axes(ax(1)); cla; axes(ax(2)); cla;
% 
% axes(ax(1)); %make axis 1 the current axis hold(ax(1); now populate the first axis
% 
% axes(ax(2)); %make axis 2 the current axis hold(ax(2)); and finally populate the second axis.

end

%% nested function "help_computePthresForFDR" for output of help
function [] = help_computePthresForFDR()
    h_help = helpdlg(FDRstruct.Info.InfoTxt,FDRstruct.Info.InfoTitle);
    uiwait(h_help);
end

end

  
function [data,pvalues,thresholdID,thresholdNP,FDRstruct] = testThis(varargin)
% [data,pvalues,thresholdID,thresholdNP,FDRstruct] = testThis(varargin);
% [data,pvalues,thresholdID,thresholdNP,FDRstruct] = testThis(FinalMean,Sigma);
% [data,pvalues,thresholdID,thresholdNP,FDRstruct] = testThis(      0.9, 1.29); %two distributions for 100 voxels one with mean 0 & sigma 1 and one with increasing mean per voxel from 0 to 0.9(FinalMean) with sigma 1.29(Sigma)
% [data,pvalues,thresholdID,thresholdNP,FDRstruct] = testThis(); %two distributions for 100 voxels one with mean 0 & sigma 1 and one with increasing mean per voxel from 0 to 1(DEFAULT for FinalMean) with sigma 1.5(DEFAULT for Sigma)

if(nargin==1)
    FinalMean = varargin{1};
    Sigma     = 1;
else
    if(nargin==2)
        FinalMean = varargin{1};
        Sigma     = varargin{2};
    else
        FinalMean = 1;
        Sigma     = 1.75;
    end
end

  % 100 voxels, progressively more active
  means{1} = zeros(1,100);
  means{2} = linspace(0,FinalMean,100);
  sigma{1} = eye(length(means{1}));
  sigma{2} = Sigma.*sigma{1};
  
  for c=1:2
    data{c} = mvnrnd(means{c},sigma{c},32);
  end

  pvalues = zeros(1,100);
  for v = 1:100
    [h,p] = ttest2(data{1}(:,v),data{2}(:,v),0.05,-1);
    pvalues(v) = p;
  end

  [thresholdID,thresholdNP,FDRstruct] = computePthresForFDR( pvalues, 0.05 );
  FDRstruct.plot();
end

