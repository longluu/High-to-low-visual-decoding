%%%%%%%%%%% Compare difference of repulsion between same and different previous trial %%%%%%%%%%%%%%%
%%%% Average across subjects
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
session = 1;
experimentNumber = 1;
markerSize = 8;
lineWidth = 0.5;
nBootstrap = 10000;
alpha = 0.05;
excludeIncorrectTrial = 0; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
                                                      
diffEst_Cond2_same_1back = cell(1, length(subjectAll));
diffEst_Cond2_diff_1back = cell(1, length(subjectAll));
diffEst_Cond2_same_2back = cell(1, length(subjectAll));
diffEst_Cond2_diff_2back = cell(1, length(subjectAll));
diffEst_Cond2_same_3back = cell(1, length(subjectAll));
diffEst_Cond2_diff_3back = cell(1, length(subjectAll));
diffEst_Cond2_same_4back = cell(1, length(subjectAll));
diffEst_Cond2_diff_4back = cell(1, length(subjectAll));
diffEst_Cond2_same_2back_cum = cell(1, length(subjectAll));
diffEst_Cond2_diff_2back_cum = cell(1, length(subjectAll));
diffEst_Cond2_same_3back_cum = cell(1, length(subjectAll));
diffEst_Cond2_diff_3back_cum = cell(1, length(subjectAll));
diffEst_Cond2_same_4back_cum = cell(1, length(subjectAll));
diffEst_Cond2_diff_4back_cum = cell(1, length(subjectAll));
diffEst_Cond2_same_5back_cum = cell(1, length(subjectAll));
diffEst_Cond2_diff_5back_cum = cell(1, length(subjectAll));
diffEst_Cond2 = cell(1, length(subjectAll));

diffEst_Cond3_same_1back = cell(1, length(subjectAll));
diffEst_Cond3_diff_1back = cell(1, length(subjectAll));
diffEst_Cond3_same_2back = cell(1, length(subjectAll));
diffEst_Cond3_diff_2back = cell(1, length(subjectAll));
diffEst_Cond3_same_3back = cell(1, length(subjectAll));
diffEst_Cond3_diff_3back = cell(1, length(subjectAll));
diffEst_Cond3_same_4back = cell(1, length(subjectAll));
diffEst_Cond3_diff_4back = cell(1, length(subjectAll));
diffEst_Cond3_same_2back_cum = cell(1, length(subjectAll));
diffEst_Cond3_diff_2back_cum = cell(1, length(subjectAll));
diffEst_Cond3_same_3back_cum = cell(1, length(subjectAll));
diffEst_Cond3_diff_3back_cum = cell(1, length(subjectAll));
diffEst_Cond3_same_4back_cum = cell(1, length(subjectAll));
diffEst_Cond3_diff_4back_cum = cell(1, length(subjectAll));
diffEst_Cond3_same_5back_cum = cell(1, length(subjectAll));
diffEst_Cond3_diff_5back_cum = cell(1, length(subjectAll));

diffEst_Cond3 = cell(1, length(subjectAll));

%% Compute the mean difference for cond 1 across all subjects
estStim1 = NaN(length(subjectAll), 50);
estStim2 = NaN(length(subjectAll), 50);

% Extract data from all subjects
for ss = 1 : length(subjectAll)
    subject = subjectAll{ss};
    
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
    estStim1(ss, :) = dataResponse(:, 3);

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
    estStim2(ss, :) = dataResponse(:, 4); 
end

% Bootstrap the mean difference
experiment_condition = 1;
[confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                experiment_condition, alpha, estStim1(:), estStim2(:));

if excludeIncorrectTrial == 2            
    mean_diffEst_Cond1 = sampleStat_all(1);
else
    mean_diffEst_Cond1 = mean(estStim2(:)) - mean(estStim1(:));
end

%% Compute the mean angle difference for same and different groups for cond 2 and 3
for ss = 1 : length(subjectAll)
    subject = subjectAll{ss};    
    
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
    [diffEst_Cond2{ss}, diffEst_Cond2_same_1back{ss}, diffEst_Cond2_diff_1back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 2-back trial (same or different) 
    n_back = 2;
    [~, diffEst_Cond2_same_2back{ss}, diffEst_Cond2_diff_2back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2); 
    
    % Split subjects' estimates based on previous 3-back trial (same or different) 
    n_back = 3;
    [~, diffEst_Cond2_same_3back{ss}, diffEst_Cond2_diff_3back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   

    % Split subjects' estimates based on previous 4-back trial (same or different) 
    n_back = 4;
    [~, diffEst_Cond2_same_4back{ss}, diffEst_Cond2_diff_4back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   

    % Split subjects' estimates based on previous 2-back trial cumulative (same or different) 
    n_back = 12;
    [~, diffEst_Cond2_same_2back_cum{ss}, diffEst_Cond2_diff_2back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 123;
    [~, diffEst_Cond2_same_3back_cum{ss}, diffEst_Cond2_diff_3back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 1234;
    [~, diffEst_Cond2_same_4back_cum{ss}, diffEst_Cond2_diff_4back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 12345;
    [~, diffEst_Cond2_same_5back_cum{ss}, diffEst_Cond2_diff_5back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
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
    [diffEst_Cond3{ss}, diffEst_Cond3_same_1back{ss}, diffEst_Cond3_diff_1back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   
    
    % Split subjects' estimates based on previous 2-back trial (same or different) 
    n_back = 2;
    [~, diffEst_Cond3_same_2back{ss}, diffEst_Cond3_diff_2back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2); 
    
    % Split subjects' estimates based on previous 3-back trial (same or different) 
    n_back = 3;
    [~, diffEst_Cond3_same_3back{ss}, diffEst_Cond3_diff_3back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
    % Split subjects' estimates based on previous 4-back trial (same or different) 
    n_back = 4;
    [~, diffEst_Cond3_same_4back{ss}, diffEst_Cond3_diff_4back{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
    % Split subjects' estimates based on previous 2-back trial cumulative (same or different) 
    n_back = 12;
    [~, diffEst_Cond3_same_2back_cum{ss}, diffEst_Cond3_diff_2back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 123;
    [~, diffEst_Cond3_same_3back_cum{ss}, diffEst_Cond3_diff_3back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 1234;
    [~, diffEst_Cond3_same_4back_cum{ss}, diffEst_Cond3_diff_4back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 3-back trial cumulative (same or different) 
    n_back = 12345;
    [~, diffEst_Cond3_same_5back_cum{ss}, diffEst_Cond3_diff_5back_cum{ss}] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
end

%% Plot the results for 1-back
trial_type = '1-back';
[fr_mean_1back, fr_se_1back] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_1back, diffEst_Cond2_diff_1back, diffEst_Cond3_same_1back, diffEst_Cond3_diff_1back, trial_type);

%% Plot the results for 2-back
trial_type = '2-back';
[fr_mean_2back, fr_se_2back] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_2back, diffEst_Cond2_diff_2back, diffEst_Cond3_same_2back, diffEst_Cond3_diff_2back, trial_type);

%% Plot the results for 3-back
trial_type = '3-back';
[fr_mean_3back, fr_se_3back] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_3back, diffEst_Cond2_diff_3back, diffEst_Cond3_same_3back, diffEst_Cond3_diff_3back, trial_type);

%% Plot the results for 4-back
trial_type = '4-back';
[fr_mean_4back, fr_se_4back] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_4back, diffEst_Cond2_diff_4back, diffEst_Cond3_same_4back, diffEst_Cond3_diff_4back, trial_type);

%% Plot the results for all back
figure
hold on
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
plot([0 5], [1 1], '--k')
errorBarGraph([fr_mean_1back; fr_mean_2back; fr_mean_3back; fr_mean_4back], ...
    [fr_mean_1back-fr_se_1back; fr_mean_2back-fr_se_2back; fr_mean_3back-fr_se_3back; fr_mean_4back-fr_se_4back], ...
    [fr_mean_1back+fr_se_1back; fr_mean_2back+fr_se_2back; fr_mean_3back+fr_se_3back; fr_mean_4back+fr_se_4back], colorIndex)
box off
title('Collapse condition 2 and 3')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
xlabel('n-back')

%% Plot the results for 2-back cum
trial_type = '2-back-cum';
[fr_mean_2back_cum, fr_se_2back_cum] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_2back_cum, diffEst_Cond2_diff_2back_cum, diffEst_Cond3_same_2back_cum, diffEst_Cond3_diff_2back_cum, trial_type);

%% Plot the results for 3-back cum
trial_type = '3-back-cum';
[fr_mean_3back_cum, fr_se_3back_cum] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_3back_cum, diffEst_Cond2_diff_3back_cum, diffEst_Cond3_same_3back_cum, diffEst_Cond3_diff_3back_cum, trial_type);

%% Plot the results for 4-back cum
trial_type = '4-back-cum';
[fr_mean_4back_cum, fr_se_4back_cum] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_4back_cum, diffEst_Cond2_diff_4back_cum, diffEst_Cond3_same_4back_cum, diffEst_Cond3_diff_4back_cum, trial_type);

%% Plot the results for 5-back cum
trial_type = '5-back-cum';
[fr_mean_5back_cum, fr_se_5back_cum] = ...
    plot_ctrAnalysis_trial(mean_diffEst_Cond1, diffEst_Cond2_same_5back_cum, diffEst_Cond2_diff_5back_cum, diffEst_Cond3_same_5back_cum, diffEst_Cond3_diff_5back_cum, trial_type);

%% Plot the results for all back cumulative trials
figure
hold on
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
plot([0 5], [1 1], '--k')
errorBarGraph([fr_mean_1back; fr_mean_2back_cum; fr_mean_3back_cum; fr_mean_4back_cum; fr_mean_5back_cum], ...
    [fr_mean_1back-fr_se_1back; fr_mean_2back_cum-fr_se_2back_cum; fr_mean_3back_cum-fr_se_3back_cum; fr_mean_4back_cum-fr_se_4back_cum; fr_mean_5back_cum-fr_se_5back_cum], ...
    [fr_mean_1back+fr_se_1back; fr_mean_2back_cum+fr_se_2back_cum; fr_mean_3back_cum+fr_se_3back_cum; fr_mean_4back_cum+fr_se_4back_cum; fr_mean_5back_cum+fr_se_5back_cum], colorIndex)
box off
title('Collapse condition 2 and 3')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
xlabel('n-back-cum')
