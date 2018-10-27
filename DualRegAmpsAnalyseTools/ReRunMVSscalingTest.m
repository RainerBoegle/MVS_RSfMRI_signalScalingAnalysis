function MVSscalingTest = ReRunMVSscalingTest(DataDirectory)

%load the data and rerun the tests

%% which ICs are prepared?
ICSelection = the ics that are detected


%% check if permtest is set up AND ask if it should be run
DoPermTest = ?
ContinueWithPermTest = ?


%% median test (signed rank test) and parameter estimates
disp(' ');
for IndIC = 1:length(ICSelection)
    %% assign
    disp(['Loading "MVSscalingTest"-struct for IC ',num2str(ICSelection(IndIC),'%02g'),' and testing Lambda values.']);
    load([DataDirectory,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat']);
    
    %% do the basic testing, i.e. PARAMETER ESTIMATION & SIGNED-RANK TEST OF LAMBDA, WITHOUT PERMUTATION TEST and write to MVSscalingTest
    MVSscalingTest = RunMVSscalingTest(MVSscalingTest);
    
    %% save MVSscalingTest to MVSscalingTest.mat
    disp('saving...');
    save([DataDirectory,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat'],'MVSscalingTest');
    clear MVSscaling MVSscalingTest
end
disp('Done with testing median(Lambda) and creating all "MVSscalingTest_IC#.mat" files.');

%% do permtest if wanted...
if(DoPermTest)
    if(~ContinueWithPermTest)
        h=helpdlg({'Permutation test will not be done now, but it is always possible to do it later based on saved setup and data.'; ' '; '(Never put off till tomorrow what you can do on the day after tomorrow just as well. -Mark Twain)'},'Do PermTest later...');
        uiwait(h);
    else
        %do permtest
        try
            clear MVSscalingTest
        end
        disp(' ');
        disp('Performing permutation test for correlation of covariates with Lambda values.');
        for IndIC = 1:length(ICSelection)
            %% load
            disp(['Loading "MVSscalingTest"-struct for IC ',num2str(ICSelection(IndIC),'%02g'),' for PermTest.']);
            load([DataDirectory,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat']);
            
            %% do permtest
            MVSscalingTest = RunMVScorrPermTest(MVSscalingTest);
            
            %% save MVSscalingTest to MVSscalingTest.mat
            disp('saving...');
            save([DataDirectory,'MVSscalingTest_IC',num2str(ICSelection(IndIC),'%02g'),'.mat'],'MVSscalingTest');
            clear MVSscaling MVSscalingTest
        end
    end
end

end