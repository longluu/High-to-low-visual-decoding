%%%%%%%%%%% Compare difference of repulsion between same and different previous trial %%%%%%%%%%%%%%%
%%%% Note that here we compute the mean repulsion for each subject and then
%%%% take the median and bootstrapped error bar across subjects. See
%%%% mainExp_CtrAnalysis_trial.m for a different approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
session = 1;
experimentNumber = 1;
markerSize = 8;
lineWidth = 0.5;
nBootstrap = 10000;
alpha = 0.05;
excludeIncorrectTrial = 2; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
                           
meanDiffEst_Cond1_same = NaN(1, length(subjectAll));
meanDiffEst_Cond1_diff = NaN(1, length(subjectAll));
meanDiffEst_Cond1 = NaN(1, length(subjectAll));

meanDiffEst_Cond2_same_1back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_1back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same_2back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_2back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same_3back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_3back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same_4back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_4back = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same_2back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_2back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same_3back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_3back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same_4back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff_4back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond2 = NaN(1, length(subjectAll));

meanDiffEst_Cond3_same_1back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_1back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same_2back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_2back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same_3back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_3back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same_4back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_4back = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same_2back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_2back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same_3back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_3back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same_4back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff_4back_cum = NaN(1, length(subjectAll));
meanDiffEst_Cond3 = NaN(1, length(subjectAll));

%% Compute the mean angle difference for same and difference groups
for ss = 1 : length(subjectAll)
    subject = subjectAll{ss};
    
    %% Experiment 1
    % First line
    experimentName = 'HighToLow_1lineShow1_49';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    % Load data
    load(dataFile)

    % Data array 
    %  Column 1: orientation line 1
    %  Column 2: orientation line 2
    %  Column 3: subject's estimate line 1
    %  Column 4: subject's estimate line 2
    %  Column 5: subject's reaction time line 1
    %  Column 6: subject's reaction time line 2
    estimateStim1_collapse = dataResponse(:, 3);

    % Second line
    experimentName = 'HighToLow_1lineShow1_54';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    % Load data
    load(dataFile)

    % Data array 
    %  Column 1: orientation line 1
    %  Column 2: orientation line 2
    %  Column 3: subject's estimate line 1
    %  Column 4: subject's estimate line 2
    %  Column 5: subject's reaction time line 1
    %  Column 6: subject's reaction time line 2
    estimateStim2_collapse = dataResponse(:, 4); 

    % Bootstrap the mean difference, correlation and percent correct
    experiment_condition = 1;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse);

    if excludeIncorrectTrial == 2            
        meanDiffEst_Cond1(1, ss) = sampleStat_all(1);
    else
        meanDiffEst_Cond1(1, ss) = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
    end
    
    %% Experiment 2
    experimentName = 'HighToLow';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    % Load data
    load(dataFile)

    % Data array 
    %  Column 1: orientation line 1
    %  Column 2: orientation line 2
    %  Column 3: subject's estimate line 1
    %  Column 4: subject's estimate line 2
    %  Column 5: subject's reaction time line 1
    %  Column 6: subject's reaction time line 2
    orientationLine1 = dataResponse(:, 1);
    orientationLine2 = dataResponse(:, 2);
    estimateLine1 = dataResponse(:, 3);
    estimateLine2 = dataResponse(:, 4);

    % Extract the stimulus orientation
    stimOrientation = unique(orientationLine1);

    % Collapse subjects' estimates across presentation order
    indStim1Est1 = orientationLine1 == stimOrientation(1);
    indStim1Est2 = orientationLine2 == stimOrientation(1);
    indStim2Est1 = orientationLine1 == stimOrientation(2);
    indStim2Est2 = orientationLine2 == stimOrientation(2);

    % Split subjects' estimates based on previous 1-back trial (same or different) 
    n_back = 1;
    trialOrder = params.trialOrder;
    [meanDiffEst_Cond2(ss), meanDiffEst_Cond2_same_1back(ss), meanDiffEst_Cond2_diff_1back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 2-back trial (same or different) 
    n_back = 2;
    [~, meanDiffEst_Cond2_same_2back(ss), meanDiffEst_Cond2_diff_2back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2); 
    
    % Split subjects' estimates based on previous 3-back trial (same or different) 
    n_back = 3;
    [~, meanDiffEst_Cond2_same_3back(ss), meanDiffEst_Cond2_diff_3back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   

    % Split subjects' estimates based on previous 4-back trial (same or different) 
    n_back = 4;
    [~, meanDiffEst_Cond2_same_4back(ss), meanDiffEst_Cond2_diff_4back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   

    % Split subjects' estimates based on previous 2-back trial cumulative (same or different) 
    n_back = 12;
    [~, meanDiffEst_Cond2_same_2back_cum(ss), meanDiffEst_Cond2_diff_2back_cum(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 123;
    [~, meanDiffEst_Cond2_same_3back_cum(ss), meanDiffEst_Cond2_diff_3back_cum(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 4-back trial cumulative (same or different) 
    n_back = 1234;
    [~, meanDiffEst_Cond2_same_4back_cum(ss), meanDiffEst_Cond2_diff_4back_cum(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    %% Experiment 3
    experimentName = 'HighToLow_separate';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    % Load data
    load(dataFile)

    % Data array 
    %  Column 1: orientation line 1
    %  Column 2: orientation line 2
    %  Column 3: subject's estimate line 1
    %  Column 4: subject's estimate line 2
    %  Column 5: subject's reaction time line 1
    %  Column 6: subject's reaction time line 2
    orientationLine1 = dataResponse(:, 1);
    orientationLine2 = dataResponse(:, 2);
    estimateLine1 = dataResponse(:, 3);
    estimateLine2 = dataResponse(:, 4);

    % Extract the stimulus orientation
    stimOrientation = unique(orientationLine1);

    % Collapse subjects' estimates across presentation order
    indStim1Est1 = orientationLine1 == stimOrientation(1);
    indStim1Est2 = orientationLine2 == stimOrientation(1);
    indStim2Est1 = orientationLine1 == stimOrientation(2);
    indStim2Est2 = orientationLine2 == stimOrientation(2);
    
    % Split subjects' estimates based on previous 1-back trial (same or different) 
    n_back = 1;
    trialOrder = params.trialOrder;
    [meanDiffEst_Cond3(ss), meanDiffEst_Cond3_same_1back(ss), meanDiffEst_Cond3_diff_1back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   
    
    % Split subjects' estimates based on previous 2-back trial (same or different) 
    n_back = 2;
    [~, meanDiffEst_Cond3_same_2back(ss), meanDiffEst_Cond3_diff_2back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2); 
    
    % Split subjects' estimates based on previous 3-back trial (same or different) 
    n_back = 3;
    [~, meanDiffEst_Cond3_same_3back(ss), meanDiffEst_Cond3_diff_3back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
    % Split subjects' estimates based on previous 4-back trial (same or different) 
    n_back = 4;
    [~, meanDiffEst_Cond3_same_4back(ss), meanDiffEst_Cond3_diff_4back(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
    % Split subjects' estimates based on previous 2-back trial cumulative (same or different) 
    n_back = 12;
    [~, meanDiffEst_Cond3_same_2back_cum(ss), meanDiffEst_Cond3_diff_2back_cum(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 123;
    [~, meanDiffEst_Cond3_same_3back_cum(ss), meanDiffEst_Cond3_diff_3back_cum(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 4-back trial cumulative (same or different) 
    n_back = 1234;
    [~, meanDiffEst_Cond3_same_4back_cum(ss), meanDiffEst_Cond3_diff_4back_cum(ss)] = split_diffEst_mean(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
end

%% Plot the results for 1-back
fig_title = "1-back";
[p_crossTrial_Cond2_1back, p_crossTrial_Cond3_1back] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_1back, meanDiffEst_Cond2_diff_1back, meanDiffEst_Cond3_same_1back, meanDiffEst_Cond3_diff_1back, fig_title, markerSize, nBootstrap);

%% Plot the results for 2-back
fig_title = "2-back";
[p_crossTrial_Cond2_2back, p_crossTrial_Cond3_2back] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_2back, meanDiffEst_Cond2_diff_2back, meanDiffEst_Cond3_same_2back, meanDiffEst_Cond3_diff_2back, fig_title, markerSize, nBootstrap);

%% Plot the results for 3-back
fig_title = "3-back";
[p_crossTrial_Cond2_3back, p_crossTrial_Cond3_3back] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_3back, meanDiffEst_Cond2_diff_3back, meanDiffEst_Cond3_same_3back, meanDiffEst_Cond3_diff_3back, fig_title, markerSize, nBootstrap);

%% Plot the results for 4-back
fig_title = "4-back";
[p_crossTrial_Cond2_4back, p_crossTrial_Cond3_4back] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_4back, meanDiffEst_Cond2_diff_4back, meanDiffEst_Cond3_same_4back, meanDiffEst_Cond3_diff_4back, fig_title, markerSize, nBootstrap);

%% Plot the results for all back
p_crossTrial_1back = [p_crossTrial_Cond2_1back p_crossTrial_Cond3_1back];
p_crossTrial_2back = [p_crossTrial_Cond2_2back p_crossTrial_Cond3_2back];
p_crossTrial_3back = [p_crossTrial_Cond2_3back p_crossTrial_Cond3_3back];
p_crossTrial_4back = [p_crossTrial_Cond2_4back p_crossTrial_Cond3_4back];

mean_p_1back = nanmedian(p_crossTrial_1back);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_1back);
sem_p_1back = nanstd(bootstat);
mean_p_2back = nanmedian(p_crossTrial_2back);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_2back);
sem_p_2back = nanstd(bootstat);
mean_p_3back = nanmedian(p_crossTrial_3back);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_3back);
sem_p_3back = nanstd(bootstat);
mean_p_4back = nanmedian(p_crossTrial_4back);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_4back);
sem_p_4back = nanstd(bootstat);

figure
hold on
plot([0 5], [1 1], '--k')
errorBarGraph([mean_p_1back; mean_p_2back; mean_p_3back; mean_p_4back], ...
    [mean_p_1back-sem_p_1back; mean_p_2back-sem_p_3back; mean_p_3back-sem_p_2back; mean_p_4back-sem_p_4back], ...
    [mean_p_1back+sem_p_1back; mean_p_2back+sem_p_3back; mean_p_3back+sem_p_2back; mean_p_4back+sem_p_4back], colorIndex)
box off
title('Collapse condition 2 and 3')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
xlabel('n-back')

%% Plot the results for 2-back cum
fig_title = "2-back-cumulative";
[p_crossTrial_Cond2_2back_cum, p_crossTrial_Cond3_2back_cum] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_2back_cum, meanDiffEst_Cond2_diff_2back_cum, meanDiffEst_Cond3_same_2back_cum, meanDiffEst_Cond3_diff_2back_cum, fig_title, markerSize, nBootstrap);

%% Plot the results for 3-back cum
fig_title = "3-back-cumulative";
[p_crossTrial_Cond2_3back_cum, p_crossTrial_Cond3_3back_cum] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_3back_cum, meanDiffEst_Cond2_diff_3back_cum, meanDiffEst_Cond3_same_3back_cum, meanDiffEst_Cond3_diff_3back_cum, fig_title, markerSize, nBootstrap);

%% Plot the results for 4-back cum
fig_title = "4-back-cumulative";
[p_crossTrial_Cond2_4back_cum, p_crossTrial_Cond3_4back_cum] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, ...
    meanDiffEst_Cond2_same_4back_cum, meanDiffEst_Cond2_diff_4back_cum, meanDiffEst_Cond3_same_4back_cum, meanDiffEst_Cond3_diff_4back_cum, fig_title, markerSize, nBootstrap);

%% Plot the results for all back cumulative trials
p_crossTrial_2back_cum = [p_crossTrial_Cond2_2back_cum p_crossTrial_Cond3_2back_cum];
p_crossTrial_3back_cum = [p_crossTrial_Cond2_3back_cum p_crossTrial_Cond3_3back_cum];
p_crossTrial_4back_cum = [p_crossTrial_Cond2_4back_cum p_crossTrial_Cond3_4back_cum];

mean_p_2back_cum = nanmedian(p_crossTrial_2back_cum);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_2back_cum);
sem_p_2back_cum = nanstd(bootstat);
mean_p_3back_cum = nanmedian(p_crossTrial_3back_cum);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_3back_cum);
sem_p_3back_cum = nanstd(bootstat);
mean_p_4back_cum = nanmedian(p_crossTrial_4back_cum);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_4back_cum);
sem_p_4back_cum = nanstd(bootstat);

figure
hold on
plot([0 5], [1 1], '--k')
errorBarGraph([mean_p_1back; mean_p_2back_cum; mean_p_3back_cum; mean_p_4back_cum], ...
    [mean_p_1back-sem_p_1back; mean_p_2back_cum-sem_p_2back_cum; mean_p_3back_cum-sem_p_3back_cum; mean_p_4back_cum-sem_p_4back_cum], ...
    [mean_p_1back+sem_p_1back; mean_p_2back_cum+sem_p_2back_cum; mean_p_3back_cum+sem_p_3back_cum; mean_p_4back_cum+sem_p_4back_cum], colorIndex)
box off
title('Collapse condition 2 and 3')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
xlabel('Cumulative n-back')
