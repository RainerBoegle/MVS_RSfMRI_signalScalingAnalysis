%% produce three distributions 
%% one skewed distribution i.e. violates signrank assumption of symmetry
%% and one symmetric i.e. good for signrank --> expect higher power
%% and a normal one with sigma of 1 around mean of 3
if(strcmp('Right',questdlg('For Data1 use skewness to the right or left for first simulated data?','Data1 skewness?','Right','Left','Right')))
    Data1 = [(-2:0.2:2)'; (-1.9:0.2:2)'; repmat(3,40,1);  10;  13]; %not symmetric & outliers
else
    Data1 = [ (4:0.2:8)';  (4.1:0.2:8)'; repmat(3,40,1); -10; -13]; %not symmetric & outliers
end

Data2 = 3+randi([-10,10],length(Data1),1);
Data3 = 3+randn(length(Data1),1);


%% display
figure(1); clf;
subplot(3,2,1); hist(Data1,length(Data1)); title(['Data1 NOT symmetric (median=',num2str(median(Data1),3),' & skewness=',num2str(skewness(Data1),3),' & kurtosis=',num2str(kurtosis(Data1),3),')']); 
subplot(3,2,2); boxplot(Data1);            title(['Data1 NOT symmetric (median=',num2str(median(Data1),3),' & skewness=',num2str(skewness(Data1),3),' & kurtosis=',num2str(kurtosis(Data1),3),')']); 
subplot(3,2,3); hist(Data2,length(Data2)); title(['Data2 symmetric (median=',num2str(median(Data2),3),' & skewness=',num2str(skewness(Data2),3),' & kurtosis=',num2str(kurtosis(Data2),3),')']); 
subplot(3,2,4); boxplot(Data2);            title(['Data2 symmetric (median=',num2str(median(Data2),3),' & skewness=',num2str(skewness(Data2),3),' & kurtosis=',num2str(kurtosis(Data2),3),')']); 
subplot(3,2,5); hist(Data3,length(Data3)); title(['Data3 normal (median=',num2str(median(Data3),3),' & skewness=',num2str(skewness(Data3),3),' & kurtosis=',num2str(kurtosis(Data3),3),')']); 
subplot(3,2,6); boxplot(Data3);            title(['Data3 normal (median=',num2str(median(Data3),3),' & skewness=',num2str(skewness(Data3),3),' & kurtosis=',num2str(kurtosis(Data3),3),')']); 


%% test for median values between data extremes and plot p-values
%% do each 10-times
MedianTestValues = min([Data1(:); Data2(:)]):0.2:max([Data1(:); Data2(:)]);
Nrepeat = 1;

AllTestDir = {'both';'left';'right'};
[SelInd,OK]= listdlg('ListString',AllTestDir,'SelectionMode','single','Name','Test-Tail direction','PromptString','Select test-tail direction');
if(OK)
    TestDir = AllTestDir{SelInd};
else
    TestDir = 'right'; %'both';
end

p_st     = ones( length(MedianTestValues),Nrepeat,3); %init not significant
h_st     = zeros(length(MedianTestValues),Nrepeat,3); %init not significant
zvals_st = zeros(length(MedianTestValues),Nrepeat,3); %init not significant

p_srt     = ones( length(MedianTestValues),Nrepeat,3); %init not significant
h_srt     = zeros(length(MedianTestValues),Nrepeat,3); %init not significant
zvals_srt = zeros(length(MedianTestValues),Nrepeat,3); %init not significant

for IndRepeat = 1:Nrepeat
    for IndTest = 1:length(MedianTestValues)
        [p1_st, h1_st, stats1_st] = signtest(Data1,MedianTestValues(IndTest),'tail',TestDir);
        [p2_st, h2_st, stats2_st] = signtest(Data2,MedianTestValues(IndTest),'tail',TestDir);
        [p3_st, h3_st, stats3_st] = signtest(Data3,MedianTestValues(IndTest),'tail',TestDir);
        
        [p1_srt, h1_srt, stats1_srt] = signrank(Data1,MedianTestValues(IndTest),'tail',TestDir);
        [p2_srt, h2_srt, stats2_srt] = signrank(Data2,MedianTestValues(IndTest),'tail',TestDir);
        [p3_srt, h3_srt, stats3_srt] = signrank(Data3,MedianTestValues(IndTest),'tail',TestDir);
        
        %% Data1
        p_st(IndTest,IndRepeat,1)    = p1_st;
        h_st(IndTest,IndRepeat,1)    = h1_st;
        zvals_st(IndTest,IndRepeat,1)= stats1_st.zval;
        
        p_srt(IndTest,IndRepeat,1)    = p1_srt;
        h_srt(IndTest,IndRepeat,1)    = h1_srt;
        zvals_srt(IndTest,IndRepeat,1)= stats1_srt.zval;
        
        %% Data2
        p_st(IndTest,IndRepeat,2)    = p2_st;
        h_st(IndTest,IndRepeat,2)    = h2_st;
        zvals_st(IndTest,IndRepeat,2)= stats2_st.zval;
        
        p_srt(IndTest,IndRepeat,2)    = p2_srt;
        h_srt(IndTest,IndRepeat,2)    = h2_srt;
        zvals_srt(IndTest,IndRepeat,2)= stats2_srt.zval;
        
        %% Data3
        p_st(IndTest,IndRepeat,3)    = p3_st;
        h_st(IndTest,IndRepeat,3)    = h3_st;
        zvals_st(IndTest,IndRepeat,3)= stats3_st.zval;
        
        p_srt(IndTest,IndRepeat,3)    = p3_srt;
        h_srt(IndTest,IndRepeat,3)    = h3_srt;
        zvals_srt(IndTest,IndRepeat,3)= stats3_srt.zval;
    end
end

%% plot p-values
figure(2); clf;
for IndRepeat = 1:Nrepeat
    ax1  = subplot(2,3,1); boxplot(Data1,'orientation','horizontal'); title({['Data1 NOT symmetric Test "',TestDir,'"']; ['(median=',num2str(median(Data1),3),' & skewness=',num2str(skewness(Data1),3),' & kurtosis=',num2str(kurtosis(Data1),3),')']}); xlim([-15 15]); %xlim([min(MedianTestValues(:)) max(MedianTestValues(:))])
    ax11 = subplot(2,3,4); plot(MedianTestValues,p_st( :,IndRepeat,1),'-');  title({['Data1 NOT symmetric Test "',TestDir,'"']; ['(median=',num2str(median(Data1),3),' & skewness=',num2str(skewness(Data1),3),' & kurtosis=',num2str(kurtosis(Data1),3),')']});  hold on  
           subplot(2,3,4); plot(MedianTestValues,p_srt(:,IndRepeat,1),'--'); hold on
end 
plot(repmat(median(Data1),21,1),(10.^(-20:0))','k-','LineWidth',2); hold on
plot(MedianTestValues,0.05.*ones(size(MedianTestValues))','r-','LineWidth',1.5);   xlim([-15 15]); hold on
legend(ax11,{'signtest','signrank-test','median','p=0.05'});
for IndRepeat = 1:Nrepeat
    ax2  = subplot(2,3,2); boxplot(Data2,'orientation','horizontal'); title({['Data2 symmetric Test "',TestDir,'"']; ['(median=',num2str(median(Data2),3),' & skewness=',num2str(skewness(Data2),3),' & kurtosis=',num2str(kurtosis(Data2),3),')']}); xlim([-15 15]); %xlim([min(MedianTestValues(:)) max(MedianTestValues(:))])
    ax21 = subplot(2,3,5); plot(MedianTestValues,p_st( :,IndRepeat,2),'-');  title({['Data1 NOT symmetric Test "',TestDir,'"']; ['(median=',num2str(median(Data2),3),' & skewness=',num2str(skewness(Data2),3),' & kurtosis=',num2str(kurtosis(Data2),3),')']}); hold on
           subplot(2,3,5); plot(MedianTestValues,p_srt(:,IndRepeat,2),'--'); hold on
end
plot(repmat(median(Data2),21,1),(10.^(-20:0))','k-','LineWidth',2);  hold on
plot(MedianTestValues,0.05.*ones(size(MedianTestValues))','r-','LineWidth',1.5);   xlim([-15 15]); hold on
legend(ax21,{'signtest','signrank-test','median','p=0.05'});
for IndRepeat = 1:Nrepeat
    ax3  = subplot(2,3,3); boxplot(Data3,'orientation','horizontal'); title({['Data3 normal Test "',TestDir,'"']; ['(median=',num2str(median(Data3),3),' & skewness=',num2str(skewness(Data3),3),' & kurtosis=',num2str(kurtosis(Data3),3),')']}); xlim([-15 15]); %xlim([min(MedianTestValues(:)) max(MedianTestValues(:))])
    ax31 = subplot(2,3,6); plot(MedianTestValues,p_st( :,IndRepeat,3),'-');  title({['Data1 NOT symmetric Test "',TestDir,'"']; ['(median=',num2str(median(Data3),3),' & skewness=',num2str(skewness(Data3),3),' & kurtosis=',num2str(kurtosis(Data3),3),')']}); hold on  
           subplot(2,3,6); plot(MedianTestValues,p_srt(:,IndRepeat,3),'--'); hold on
end
plot(repmat(median(Data3),21,1),(10.^(-20:0))','k-','LineWidth',2);  hold on
plot(MedianTestValues,0.05.*ones(size(MedianTestValues))','r-','LineWidth',1.5);   xlim([-15 15]); hold on
legend(ax31,{'signtest','signrank-test','median','p=0.05'});
linkaxes([ax11 ax21 ax31],'y'); set([ax11 ax21 ax31],'YScale','log');
linkaxes([ax2  ax1  ax3 ax11 ax21 ax31],'x');
