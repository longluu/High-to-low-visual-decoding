%%%%%%%%%%% Analyze data from the experiment with custom bootstrap function (separate) %%%%%%%%%%%%%%%
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
nBootstrap = 10000;
markerSize = 8;
lineWidth = 0.5;
alpha = 0.05;
plotIndividualSubject = 0;
excludeIncorrectTrial = 2; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
meanDiffEst_Cond1 = NaN(4, length(subjectAll));
meanDiffEst_Cond2 = NaN(4, length(subjectAll));
meanDiffEst_Cond3 = NaN(4, length(subjectAll));
corr_Cond1 = NaN(4, length(subjectAll));
corr_Cond2 = NaN(4, length(subjectAll));
corr_Cond3 = NaN(4, length(subjectAll));
percentCorrect_Cond1 = NaN(4, length(subjectAll));
percentCorrect_Cond2 = NaN(4, length(subjectAll));
percentCorrect_Cond3 = NaN(4, length(subjectAll));
std_Cond2_split = NaN(2, length(subjectAll));
std_Cond3_split = NaN(2, length(subjectAll));
std_Cond1 = NaN(1, length(subjectAll));
std_Cond2 = NaN(1, length(subjectAll));
std_Cond3 = NaN(1, length(subjectAll));

%% Summary statistics
for ss = 1 : length(subjectAll)
    subject = subjectAll{ss};
    
    %% Experiment 1
    if useSplit1line
        %% First line
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

        %% Second line
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
         
    else
        experimentName = 'HighToLow_1lineShow1';
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
        estimateStim1 = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
        estimateStim2 = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
        estimateStim1_collapse = estimateStim1(~isnan(estimateStim1));
        estimateStim2_collapse = estimateStim2(~isnan(estimateStim2));
    end
    
    % Bootstrap the mean difference, correlation and percent correct
    experiment_condition = 1;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse);

    if excludeIncorrectTrial == 2            
        meanDiffEst_Cond1(1, ss) = sampleStat_all(1);
    else
        meanDiffEst_Cond1(1, ss) = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
    end
    meanDiffEst_Cond1(3:4, ss) = confInterval_all(:, 1);
    corr_Cond1(1, ss) = sampleStat_all(2);
    corr_Cond1(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Cond1(1, ss) = sampleStat_all(3);
    percentCorrect_Cond1(3:4, ss) = confInterval_all(:, 3);
    std_Cond1(ss) = sqrt(nanstd(estimateStim1_collapse)^2 +  nanstd(estimateStim2_collapse)^2);
    
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
    estimateStim1_collapse = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
    estimateStim2_collapse = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
    estimateStim1_collapse = estimateStim1_collapse(~isnan(estimateStim1_collapse));
    estimateStim2_collapse = estimateStim2_collapse(~isnan(estimateStim2_collapse));
    estimateStim1_collapse(estimateStim1_collapse<0) = estimateStim1_collapse(estimateStim1_collapse<0)+180;
    estimateStim2_collapse(estimateStim2_collapse<0) = estimateStim2_collapse(estimateStim2_collapse<0)+180;

    
    % Std of first and second estimate
    std_estStim1Est1 = nanstd(estimateLine1(indStim1Est1));
    std_estStim2Est1 = nanstd(estimateLine1(indStim2Est1));
    std_est1_Cond2 = sqrt((std_estStim1Est1^2 + std_estStim2Est1^2) / 2);
    std_estStim1Est2 = nanstd(estimateLine2(indStim1Est2));
    std_estStim2Est2 = nanstd(estimateLine2(indStim2Est2));
    std_est2_Cond2 = sqrt((std_estStim1Est2^2 + std_estStim2Est2^2) / 2);
    
    % Bootstrap the mean difference, correlation and percent correct
    experiment_condition = 2;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse);

    meanDiffEst_Cond2(1, ss) = sampleStat_all(1);    
    meanDiffEst_Cond2(3:4, ss) = confInterval_all(:, 1);
    corr_Cond2(1, ss) = sampleStat_all(2);
    corr_Cond2(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Cond2(1, ss) = sampleStat_all(3);
    percentCorrect_Cond2(3:4, ss) = confInterval_all(:, 3);
    std_Cond2(ss) = nanstd(estimateStim2_collapse - estimateStim1_collapse);
    
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
    estimateStim1_collapse = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
    estimateStim2_collapse = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
    estimateStim1_collapse = estimateStim1_collapse(~isnan(estimateStim1_collapse));
    estimateStim2_collapse = estimateStim2_collapse(~isnan(estimateStim2_collapse));
    
    % Std of first and second estimate
    std_estStim1Est1 = nanstd(estimateLine1(indStim1Est1));
    std_estStim2Est1 = nanstd(estimateLine1(indStim2Est1));
    std_est1_Cond3 = sqrt((std_estStim1Est1^2 + std_estStim2Est1^2) / 2);
    std_estStim1Est2 = nanstd(estimateLine2(indStim1Est2));
    std_estStim2Est2 = nanstd(estimateLine2(indStim2Est2));
    std_est2_Cond3 = sqrt((std_estStim1Est2^2 + std_estStim2Est2^2) / 2);

    % Bootstrap the mean difference, correlation and percent correct
    experiment_condition = 3;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse);

    meanDiffEst_Cond3(1, ss) = sampleStat_all(1);    
    meanDiffEst_Cond3(3:4, ss) = confInterval_all(:, 1);
    corr_Cond3(1, ss) = sampleStat_all(2);
    corr_Cond3(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Cond3(1, ss) = sampleStat_all(3);
    percentCorrect_Cond3(3:4, ss) = confInterval_all(:, 3);
    std_Cond3(ss) = nanstd(estimateStim2_collapse - estimateStim1_collapse);

    %% Save std of estimates first and second line     
    std_Cond2_split(1, ss) = std_est1_Cond2;
    std_Cond2_split(2, ss) = std_est2_Cond2;
    std_Cond3_split(1, ss) = std_est1_Cond3;
    std_Cond3_split(2, ss) = std_est2_Cond3; 
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
subplot(2, 3, 1)
hold on
minPlot_diffEst = min([meanDiffEst_Cond1(1, :) meanDiffEst_Cond2(1, :) meanDiffEst_Cond3(1, :)]) - 2;
maxPlot_diffEst = max([meanDiffEst_Cond1(1, :) meanDiffEst_Cond2(1, :) meanDiffEst_Cond3(1, :)]) + 2;
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
legend_handle = NaN(1, length(subjectAll));
for ii = 1 : length(subjectAll)
    plot([meanDiffEst_Cond1(3, ii) meanDiffEst_Cond1(4, ii)], [meanDiffEst_Cond2(1, ii) meanDiffEst_Cond2(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([meanDiffEst_Cond1(1, ii) meanDiffEst_Cond1(1, ii)], [meanDiffEst_Cond2(3, ii) meanDiffEst_Cond2(4, ii)], '-k', 'LineWidth', lineWidth)    
    legend_handle(ii) = plot(meanDiffEst_Cond1(1, ii), meanDiffEst_Cond2(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
end
if excludeIncorrectTrial == 1
    title('Correct trials only')
elseif excludeIncorrectTrial == 2
    title('All trials - Flip incorrect trials')
else
    title('All trials')
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
% legend(legend_handle, legend_meanDiff1)
xlabel('Mean difference - Cond1 (deg)')
ylabel('Mean difference - Cond2 (deg)')

subplot(2, 3, 2)
hold on
minPlot_corr = min([corr_Cond1(1, :) corr_Cond2(1, :) corr_Cond3(1, :)]) - 0.4;
maxPlot_corr = max([corr_Cond1(1, :) corr_Cond2(1, :) corr_Cond3(1, :)]) + 0.1;
plot([minPlot_corr maxPlot_corr], [minPlot_corr maxPlot_corr], 'k')
for ii = 1 : length(subjectAll)
    plot([corr_Cond1(3, ii) corr_Cond1(4, ii)], [corr_Cond2(1, ii) corr_Cond2(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([corr_Cond1(1, ii) corr_Cond1(1, ii)], [corr_Cond2(3, ii) corr_Cond2(4, ii)], '-k', 'LineWidth', lineWidth)    
    plot(corr_Cond1(1, ii), corr_Cond2(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
if excludeIncorrectTrial == 1
    title('Correct trials only')
elseif excludeIncorrectTrial == 2
    title('All trials - Flip incorrect trials')
else
    title('All trials')
end
axis([minPlot_corr, maxPlot_corr, minPlot_corr, maxPlot_corr])
axis square
% legend(legend_corr1)
xlabel('Correlation - Cond1 (deg)')
ylabel('Correlation - Cond2 (deg)')

subplot(2, 3, 3)
hold on
minPlot_pc = min([percentCorrect_Cond1(1, :) percentCorrect_Cond2(1, :) percentCorrect_Cond3(1, :)]) - 15;
maxPlot_pc = 110;
plot([minPlot_pc maxPlot_pc], [minPlot_pc maxPlot_pc])
for ii = 1 : length(subjectAll)
    plot([percentCorrect_Cond1(3, ii) percentCorrect_Cond1(4, ii)], [percentCorrect_Cond2(1, ii) percentCorrect_Cond2(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([percentCorrect_Cond1(1, ii) percentCorrect_Cond1(1, ii)], [percentCorrect_Cond2(3, ii) percentCorrect_Cond2(4, ii)], '-k', 'LineWidth', lineWidth)
    plot(percentCorrect_Cond1(1, ii), percentCorrect_Cond2(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_pc, maxPlot_pc, minPlot_pc, maxPlot_pc])
axis square
% legend(legend_pc1)
xlabel('Percent correct - Cond1 (deg)')
ylabel('Percent correct - Cond2 (deg)')

subplot(2, 3, 4)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot([meanDiffEst_Cond1(3, ii) meanDiffEst_Cond1(4, ii)], [meanDiffEst_Cond3(1, ii) meanDiffEst_Cond3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([meanDiffEst_Cond1(1, ii) meanDiffEst_Cond1(1, ii)], [meanDiffEst_Cond3(3, ii) meanDiffEst_Cond3(4, ii)], '-k', 'LineWidth', lineWidth)
    plot(meanDiffEst_Cond1(1, ii), meanDiffEst_Cond3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
if excludeIncorrectTrial == 1
    title('Correct trials only')
elseif excludeIncorrectTrial == 2
    title('All trials - Flip incorrect trials')
else
    title('All trials')
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
% legend(legend_meanDiff2)
xlabel('Mean difference - Cond1 (deg)')
ylabel('Mean difference - Cond3 (deg)')

subplot(2, 3, 5)
hold on
plot([minPlot_corr maxPlot_corr], [minPlot_corr maxPlot_corr], 'k')
for ii = 1 : length(subjectAll)
    plot([corr_Cond1(3, ii) corr_Cond1(4, ii)], [corr_Cond3(1, ii) corr_Cond3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([corr_Cond1(1, ii) corr_Cond1(1, ii)], [corr_Cond3(3, ii) corr_Cond3(4, ii)], '-k', 'LineWidth', lineWidth)    
    plot(corr_Cond1(1, ii), corr_Cond3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_corr, maxPlot_corr, minPlot_corr, maxPlot_corr])
axis square
if excludeIncorrectTrial == 1
    title('Correct trials only')
elseif excludeIncorrectTrial == 2
    title('All trials - Flip incorrect trials')
else
    title('All trials')
end
% legend(legend_corr2)
xlabel('Correlation - Cond1 (deg)')
ylabel('Correlation - Cond3 (deg)')

subplot(2, 3, 6)
hold on
plot([minPlot_pc maxPlot_pc], [minPlot_pc maxPlot_pc])
for ii = 1 : length(subjectAll)
    plot([percentCorrect_Cond1(3, ii) percentCorrect_Cond1(4, ii)], [percentCorrect_Cond3(1, ii) percentCorrect_Cond3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([percentCorrect_Cond1(1, ii) percentCorrect_Cond1(1, ii)], [percentCorrect_Cond3(3, ii) percentCorrect_Cond3(4, ii)], '-k', 'LineWidth', lineWidth)    
    plot(percentCorrect_Cond1(1, ii), percentCorrect_Cond3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
% legend(legend_pc2)
axis([minPlot_pc, maxPlot_pc, minPlot_pc, maxPlot_pc])
axis square
xlabel('Percent correct - Cond1 (deg)')
ylabel('Percent correct - Cond3 (deg)')

figure;
hold on
subplot(1, 2, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot([meanDiffEst_Cond2(3, ii) meanDiffEst_Cond2(4, ii)], [meanDiffEst_Cond3(1, ii) meanDiffEst_Cond3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([meanDiffEst_Cond2(1, ii) meanDiffEst_Cond2(1, ii)], [meanDiffEst_Cond3(3, ii) meanDiffEst_Cond3(4, ii)], '-k', 'LineWidth', lineWidth)
    plot(meanDiffEst_Cond2(1, ii), meanDiffEst_Cond3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
if excludeIncorrectTrial == 1
    title('Correct trials only')
elseif excludeIncorrectTrial == 2
    title('All trials - Flip incorrect trials')
else
    title('All trials')
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
% legend(legend_meanDiff2)
xlabel('Mean difference - Cond2 (deg)')
ylabel('Mean difference - Cond3 (deg)')

subplot(1, 2, 2)
hold on
x = 1:length(subjectAll);
for ii = 1 : length(x)
    plot(x(ii), x(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 15)
end
legend(subjectAll, 'FontSize', 20, 'Box', 'off')
for ii = 1 : length(x)
    plot(x(ii), x(ii), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', 'none', 'MarkerSize', 15)
end
xlim([10.6 11])
axis off

figure;
subplot(2, 2, 1)
hold on
minPlot = min([std_Cond2_split(1, :) std_Cond2_split(2, :)]) - 0.1;
maxPlot = max([std_Cond2_split(1, :) std_Cond2_split(2, :)]) + 0.1;
for ii = 1 : length(subjectAll)
    plot(std_Cond2_split(1, ii), std_Cond3_split(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
title('Std of firt estimate')
xlabel('Std of first estimate - Cond 2 (deg)')
ylabel('Std of first estimate - Cond 3 (deg)')
axis square

subplot(2, 2, 2)
hold on
minPlot = min([std_Cond3_split(1, :) std_Cond3_split(2, :)]) - 0.1;
maxPlot = max([std_Cond3_split(1, :) std_Cond3_split(2, :)]) + 0.1;
for ii = 1 : length(subjectAll)
    plot(std_Cond2_split(2, ii), std_Cond3_split(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
axis square
title('Std of second estimate')
xlabel('Std of second estimate - Cond 2 (deg)')
ylabel('Std of second estimate - Cond 3 (deg)')

subplot(2, 2, 3)
hold on
minPlot = min([std_Cond1 std_Cond2]) - 0.1;
maxPlot = max([std_Cond1 std_Cond2]) + 0.1;
for ii = 1 : length(subjectAll)
    plot(std_Cond1(ii), std_Cond2(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
axis square
title('Cond 2 vs Cond 1')
xlabel('Std of estimate difference - Cond 1 (deg)')
ylabel('Std of estimate difference - Cond 2 (deg)')

subplot(2, 2, 4)
hold on
minPlot = min([std_Cond1 std_Cond3]) - 0.1;
maxPlot = max([std_Cond1 std_Cond3]) + 0.1;
for ii = 1 : length(subjectAll)
    plot(std_Cond1(ii), std_Cond3(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
axis square
title('Cond 3 vs Cond 1')
xlabel('Std of estimate difference - Cond 1 (deg)')
ylabel('Std of estimate difference - Cond 3 (deg)')

%% Perform statistical tests
% Mean diff
p_meanDiff_1vs2 = signrank(meanDiffEst_Cond1(1, :), meanDiffEst_Cond2(1, :),...
                                    'method', 'exact');
fprintf('p-value mean diff cond 1 vs cond 2: %8.6f \n', p_meanDiff_1vs2)                     
p_meanDiff_1vs3= signrank(meanDiffEst_Cond1(1, :), meanDiffEst_Cond3(1, :),...
                                    'method', 'exact');
fprintf('p-value mean diff cond 1 vs cond 3: %8.6f \n', p_meanDiff_1vs3)  

% Correlation
p_corr_1vs2= signrank(corr_Cond1(1, :), corr_Cond2(1, :),...
                                    'method', 'exact');
fprintf('p-value correlation cond 1 vs cond 2: %8.6f \n', p_corr_1vs2) 
p_corr_1vs3= signrank(corr_Cond1(1, :), corr_Cond3(1, :),...
                                    'method', 'exact');
fprintf('p-value correlation cond 1 vs cond 3: %8.6f \n', p_corr_1vs3)    

% Percent correct
p_pc_1vs2= signrank(percentCorrect_Cond1(1, :), percentCorrect_Cond2(1, :),...
                                    'method', 'exact');
fprintf('p-value percent correct cond 1 vs cond 2: %8.6f \n', p_pc_1vs2) 
p_pc_1vs3= signrank(percentCorrect_Cond1(1, :), percentCorrect_Cond3(1, :),...
                                    'method', 'exact');
fprintf('p-value percent correct cond 1 vs cond 3: %8.6f \n', p_pc_1vs3)    

%% Plot individual subjects
if plotIndividualSubject
    h_currentFig = gcf;
    markerSize = 4;
    for ss = 1 : length(subjectAll)
        subject = subjectAll{ss};
        h_currentFig = figure(h_currentFig.Number+1);
        h_currentFig.Position = [300 300 1200 700];
        
        hold on

        %% Experiment 1
        if useSplit1line
            %% First line
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

            % Extract the stimulus orientation
            stimOrientation = params.lineOrientation;

            %% Second line
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

        else
            experimentName = 'HighToLow_1lineShow1';
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
            estimateStim1 = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
            estimateStim2 = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
            estimateStim1_collapse = estimateStim1(~isnan(estimateStim1));
            estimateStim2_collapse = estimateStim2(~isnan(estimateStim2));
        end

        % Analysis of percent correct
        meanEst1 = nanmean(estimateStim1_collapse);
        meanEst2 = nanmean(estimateStim2_collapse);
        criterion = linspace(meanEst1-3, meanEst2+3, 50);
        percentCorrect = NaN(1, length(criterion));
        for ii = 1 : length(criterion)
            percentCorrect(ii) = (sum(estimateStim1_collapse < criterion(ii)) + sum(estimateStim2_collapse >= criterion(ii))) / ...
                                    (sum(~isnan(estimateStim1_collapse)) + sum(~isnan(estimateStim2_collapse)));
        end

        % Find the permutation of data that has lowest correlation
        eps_correlation = 0.001;
        corr_2line = 1;
        while abs(corr_2line) > eps_correlation
            est1 = datasample(estimateStim1_collapse, length(estimateStim1_collapse), 'Replace', false);
            est2 = datasample(estimateStim2_collapse, length(estimateStim2_collapse), 'Replace', false);
            corr_2line = corr(est1, est2);
        end
        estimateStim1_collapse = est1;
        estimateStim2_collapse = est2;
        diffEst = estimateStim2_collapse - estimateStim1_collapse;
        
        % Plot
        subplot(2, 3, 1)
        hold on
        minPlot = min([estimateStim1_collapse; estimateStim2_collapse])-1;
        maxPlot = max([estimateStim1_collapse; estimateStim2_collapse])+1;
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none');          
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title('1 line')

        subplot(2, 3, 4)
        histogram(diffEst, 20, 'Normalization', 'probability')
        xlabel('Orientation difference (deg)')
        ylabel('Frequency of occurence')
        xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
        axis square
        title(['Mean :' num2str(round(meanEst2 - meanEst1, 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))])    

        %% Experiment 2
        experimentName = 'HighToLow';
        dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

        % Load data
        load(dataFile)
        trialOrder = params.trialOrder;

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
        rtLine1 = dataResponse(:, 5);
        rtLine2 = dataResponse(:, 6);

        % Extract the stimulus orientation
        stimOrientation = unique(orientationLine1);

        % Collapse subjects' estimates across presentation order
        indStim1Est1 = orientationLine1 == stimOrientation(1);
        indStim1Est2 = orientationLine2 == stimOrientation(1);
        indStim2Est1 = orientationLine1 == stimOrientation(2);
        indStim2Est2 = orientationLine2 == stimOrientation(2);
        estimateStim1_collapse = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
        estimateStim2_collapse = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
        estimateStim1_collapse = estimateStim1_collapse(~isnan(estimateStim1_collapse));
        estimateStim2_collapse = estimateStim2_collapse(~isnan(estimateStim2_collapse));
        estimateStim1_collapse(estimateStim1_collapse<0) = estimateStim1_collapse(estimateStim1_collapse<0)+180;
        estimateStim2_collapse(estimateStim2_collapse<0) = estimateStim2_collapse(estimateStim2_collapse<0)+180;        
        diffEst = estimateStim2_collapse - estimateStim1_collapse;
        percentCorrect = (sum(diffEst > 0) + sum(diffEst == 0)/2) / sum(~isnan(diffEst));
        [corrCoef, pValue] = corr(estimateStim1_collapse, estimateStim2_collapse);
        meanDiffEst_2line1(ss) = mean(diffEst);

        % Plot estimates 
        subplot(2, 3, 2)
        hold on
        minPlot = min([estimateStim1_collapse; estimateStim2_collapse])-1;
        maxPlot = max([estimateStim1_collapse; estimateStim2_collapse])+1;
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')        
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none')
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title('2 line')

        subplot(2, 3, 5)
        histogram(diffEst, 20, 'Normalization', 'probability')
        xlabel('Orientation difference (deg)')
        ylabel('Frequency of occurence')
        xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
        axis square
        title(['Mean: ' num2str(round(mean(diffEst), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1)) ])

        %% Experiment 3
        experimentName = 'HighToLow_separate';
        dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

        % Load data
        load(dataFile)
        trialOrder = params.trialOrder;

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
        rtLine1 = dataResponse(:, 5);
        rtLine2 = dataResponse(:, 6);

        % Extract the stimulus orientation
        stimOrientation = unique(orientationLine1);

        % Collapse subjects' estimates across presentation order
        indStim1Est1 = orientationLine1 == stimOrientation(1);
        indStim1Est2 = orientationLine2 == stimOrientation(1);
        indStim2Est1 = orientationLine1 == stimOrientation(2);
        indStim2Est2 = orientationLine2 == stimOrientation(2);
        estimateStim1_collapse = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
        estimateStim2_collapse = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
        indInclude = (~isnan(estimateStim1_collapse)) & (~isnan(estimateStim2_collapse));
        estimateStim1_collapse = estimateStim1_collapse(indInclude);
        estimateStim2_collapse = estimateStim2_collapse(indInclude);
        diffEst = estimateStim2_collapse - estimateStim1_collapse;
        percentCorrect = (sum(diffEst > 0) + sum(diffEst == 0)/2) / sum(~isnan(diffEst));
        [corrCoef, pValue] = corr(estimateStim1_collapse, estimateStim2_collapse);

        % Plot estimates    
        subplot(2, 3, 3)
        hold on
        minPlot = min([estimateStim1_collapse; estimateStim2_collapse])-1;
        maxPlot = max([estimateStim1_collapse; estimateStim2_collapse])+1;
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none')        
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title('2 line - redraw')

        subplot(2, 3, 6)
        histogram(diffEst, 20, 'Normalization', 'probability')
        xlabel('Orientation difference (deg)')
        ylabel('Frequency of occurence')
        xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
        axis square
        title(['Mean: ' num2str(round(mean(diffEst), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))]) 
        
        saveas(h_currentFig, [subjectAll{ss} '.pdf'])
    end
end