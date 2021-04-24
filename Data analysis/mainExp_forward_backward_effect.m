%%%%%%%%%%% Compare forward vs. backward effect (i.e. whether the estimates of the same line is different when it's reported first vs. second) %%%%%%%%%%%%%%%
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
session = 1;
experimentNumber = 1;
nBootstrap = 10000;
markerSize = 8;
lineWidth = 0.5;
alpha = 0.05;
excludeIncorrectTrial = 0; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
bias_Cond2_stim1 = NaN(2, length(subjectAll));
bias_Cond2_stim2 = NaN(2, length(subjectAll));
bias_Cond3_stim1 = NaN(2, length(subjectAll));
bias_Cond3_stim2 = NaN(2, length(subjectAll));


%% Summary statistics
for ss = 1 : length(subjectAll)
    subject = subjectAll{ss};
    
    %% Experiment 1
    % First line
    experimentName = 'HighToLow_1lineShow1_49';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];
    load(dataFile)

    estimateStim1_collapse = dataResponse(:, 3);

    % Second line
    experimentName = 'HighToLow_1lineShow1_54';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];
    load(dataFile)

    estimateStim2_collapse = dataResponse(:, 4);
             
    est_stim1_cond1 = nanmean(estimateStim1_collapse);
    est_stim2_cond1 = nanmean(estimateStim2_collapse);
    
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
    estimateLine1(estimateLine1<0) = estimateLine1(estimateLine1<0)+180;
    estimateLine2(estimateLine2<0) = estimateLine2(estimateLine2<0)+180;
    
    % Extract the stimulus orientation
    stimOrientation = unique(orientationLine1);

    % Get the index for 49 deg stim (1st and 2nd report) and 54 deg (1st and 2nd report)
    % indStim1Est1 (49 deg, 1st), indStim1Est2 (49 deg, 2nd), indStim2Est1 (54 deg, 1st), indStim2Est2 (54 deg, 2nd)
    indStim1Est1 = orientationLine1 == stimOrientation(1);
    indStim1Est2 = orientationLine2 == stimOrientation(1);
    indStim2Est1 = orientationLine1 == stimOrientation(2);
    indStim2Est2 = orientationLine2 == stimOrientation(2);
    
    % Bias split by stimulus orientation and report order
    bias_Cond2_stim1(1, ss) = nanmean(estimateLine1(indStim1Est1)) - est_stim1_cond1;
    bias_Cond2_stim1(2, ss) = nanmean(estimateLine2(indStim1Est2)) - est_stim1_cond1;
    bias_Cond2_stim2(1, ss) = nanmean(estimateLine1(indStim2Est1)) - est_stim2_cond1;
    bias_Cond2_stim2(2, ss) = nanmean(estimateLine2(indStim2Est2)) - est_stim2_cond1;
    
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

    % Bias split by stimulus orientation and report order
    bias_Cond3_stim1(1, ss) = nanmean(estimateLine1(indStim1Est1)) - est_stim1_cond1;
    bias_Cond3_stim1(2, ss) = nanmean(estimateLine2(indStim1Est2)) - est_stim1_cond1;
    bias_Cond3_stim2(1, ss) = nanmean(estimateLine1(indStim2Est1)) - est_stim2_cond1;
    bias_Cond3_stim2(2, ss) = nanmean(estimateLine2(indStim2Est2)) - est_stim2_cond1;
    
end

%% Plot the results
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
hold on
subplot(1, 2, 1)
hold on
minPlot = min([bias_Cond2_stim1(:); bias_Cond2_stim2(:); bias_Cond3_stim1(:); bias_Cond3_stim2(:)]) - 2;
maxPlot = max([bias_Cond2_stim1(:); bias_Cond2_stim2(:); bias_Cond3_stim1(:); bias_Cond3_stim2(:)]) + 2;
legend_handle = NaN(1, length(subjectAll));
for ii = 1 : length(subjectAll)
    legend_handle(ii) = plot(bias_Cond2_stim1(1, ii), bias_Cond2_stim1(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    plot(bias_Cond2_stim2(1, ii), bias_Cond2_stim2(2, ii), 'd', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
title('Condition 2')
axis([minPlot, maxPlot, minPlot, maxPlot])
axis square
xlabel('Backward aftereffect (deg)')
ylabel('Forward aftereffect (deg)')

subplot(1, 2, 2)
hold on
legend_handle = NaN(1, length(subjectAll));
for ii = 1 : length(subjectAll)
    legend_handle(ii) = plot(bias_Cond3_stim1(1, ii), bias_Cond3_stim1(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    plot(bias_Cond3_stim2(1, ii), bias_Cond3_stim2(2, ii), 'd', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
title('Condition 3')
axis([minPlot, maxPlot, minPlot, maxPlot])
axis square
xlabel('Backward aftereffect (deg)')
ylabel('Forward aftereffect (deg)')

%% Perform statistical tests
p_effect_cond2 = signrank([bias_Cond2_stim1(1, :) bias_Cond2_stim2(1, :)], [bias_Cond2_stim1(2, :) bias_Cond2_stim2(2, :)],...
                                    'method', 'exact');
fprintf('p-value backward vs. forward effect cond 2: %8.6f \n', p_effect_cond2)                     
p_effect_cond3 = signrank([bias_Cond3_stim1(1, :) bias_Cond3_stim2(1, :)], [bias_Cond3_stim1(2, :) bias_Cond3_stim2(2, :)],...
                                    'method', 'exact');
fprintf('p-value backward vs. forward effect cond 3: %8.6f \n', p_effect_cond3)    
