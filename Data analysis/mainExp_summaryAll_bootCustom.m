%%%%%%%%%%% Analyze data from the experiment with custom bootstrap function (separate) %%%%%%%%%%%%%%%
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1; % use condition with separate blocks
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
meanDiff_Cond2_split = NaN(2, length(subjectAll));
meanDiff_Cond3_split = NaN(2, length(subjectAll));
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
    std_Cond1(ss) = sqrt((nanstd(estimateStim1_collapse)^2 +  nanstd(estimateStim2_collapse)^2)/2);
    
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
    std_Cond2_split(1, ss) = sqrt((std_estStim1Est1^2 + std_estStim2Est1^2) / 2);
    std_estStim1Est2 = nanstd(estimateLine2(indStim1Est2));
    std_estStim2Est2 = nanstd(estimateLine2(indStim2Est2));
    std_Cond2_split(2, ss) = sqrt((std_estStim1Est2^2 + std_estStim2Est2^2) / 2);

    % Mean diff of 1st and 2nd estimate
    meanDiff_Cond2_split(1, ss) = nanmean(estimateLine1(indStim2Est1)) - nanmean(estimateLine1(indStim1Est1));
    meanDiff_Cond2_split(2, ss) = nanmean(estimateLine2(indStim2Est2)) - nanmean(estimateLine2(indStim1Est2));
    
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
    std_Cond2(ss) = sqrt((nanstd(estimateStim1_collapse)^2 +  nanstd(estimateStim2_collapse)^2)/2);
    
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

    % Mean diff of 1st and 2nd estimate
    meanDiff_Cond3_split(1, ss) = nanmean(estimateLine1(indStim2Est1)) - nanmean(estimateLine1(indStim1Est1));
    meanDiff_Cond3_split(2, ss) = nanmean(estimateLine2(indStim2Est2)) - nanmean(estimateLine2(indStim1Est2));
    
    % Std of first and second estimate
    std_estStim1Est1 = nanstd(estimateLine1(indStim1Est1));
    std_estStim2Est1 = nanstd(estimateLine1(indStim2Est1));
    std_Cond3_split(1, ss) = sqrt((std_estStim1Est1^2 + std_estStim2Est1^2) / 2);
    std_estStim1Est2 = nanstd(estimateLine2(indStim1Est2));
    std_estStim2Est2 = nanstd(estimateLine2(indStim2Est2));
    std_Cond3_split(2, ss) = sqrt((std_estStim1Est2^2 + std_estStim2Est2^2) / 2);

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
    std_Cond3(ss) = sqrt((nanstd(estimateStim1_collapse)^2 +  nanstd(estimateStim2_collapse)^2)/2);
end

%% Plot summary statistics of individual subject (mean diff, corr, percent correct)
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

%% Plot comparing mean diff of 2-line vs. 2-line-interupt and legend of subj
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

%% Plot analysis of std
figure;
subplot(2, 2, 1)
hold on
minPlot = min([std_Cond2_split(1, :) std_Cond2_split(2, :)]) - 0.5;
maxPlot = max([std_Cond2_split(1, :) std_Cond2_split(2, :)]) + 0.5;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
for ii = 1 : length(subjectAll)
    plot(std_Cond2_split(1, ii), std_Cond2_split(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title('2-line condition')
xlabel('Std of 1st estimate (deg)')
ylabel('Std of 2nd estimate (deg)')

subplot(2, 2, 2)
hold on
minPlot = min([std_Cond3_split(1, :) std_Cond3_split(2, :)]) - 0.5;
maxPlot = max([std_Cond3_split(1, :) std_Cond3_split(2, :)]) + 0.5;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
for ii = 1 : length(subjectAll)
    plot(std_Cond3_split(1, ii), std_Cond3_split(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title('2-line-interupt condition')
xlabel('Std of 1st estimate (deg)')
ylabel('Std of 2nd estimate (deg)')

subplot(2, 2, 3)
hold on
minPlot = min([std_Cond1 std_Cond2 std_Cond3]) - 0.5;
maxPlot = max([std_Cond1 std_Cond2 std_Cond3]) + 0.5;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
for ii = 1 : length(subjectAll)
    plot(std_Cond1(ii), std_Cond2(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title('1-line vs. 2-line')
xlabel('Std of estimate, 1-line (deg)')
ylabel('Std of estimate, 2-line (deg)')

subplot(2, 2, 4)
hold on
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
for ii = 1 : length(subjectAll)
    plot(std_Cond1(ii), std_Cond3(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title('1-line vs. 2-line interupt')
xlabel('Std of estimate, 1-line (deg)')
ylabel('Std of estimate, 2-line interupt (deg)')

%% Plot mean diff 1st vs. 2nd estimates
figure;
subplot(1, 2, 1)
hold on
minPlot = min([meanDiff_Cond2_split(1, :) meanDiff_Cond2_split(2, :)]) - 0.5;
maxPlot = max([meanDiff_Cond2_split(1, :) meanDiff_Cond2_split(2, :)]) + 0.5;
plot([minPlot maxPlot], [minPlot maxPlot])
for ii = 1 : length(subjectAll)
    plot(meanDiff_Cond2_split(1, ii), meanDiff_Cond2_split(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title('1-line vs. 2-line')
xlabel('Mean diff, 1-line (deg)')
ylabel('Mean diff, 2-line (deg)')

subplot(1, 2, 2)
hold on
minPlot = min([meanDiff_Cond3_split(1, :) meanDiff_Cond3_split(2, :)]) - 0.5;
maxPlot = max([meanDiff_Cond3_split(1, :) meanDiff_Cond3_split(2, :)]) + 0.5;
plot([minPlot maxPlot], [minPlot maxPlot])
for ii = 1 : length(subjectAll)
    plot(meanDiff_Cond3_split(1, ii), meanDiff_Cond3_split(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title('2-line vs. 2-line-interupt')
xlabel('Mean diff, 2-line (deg)')
ylabel('Mean diff, 2-line-interupt (deg)')

%% Plot correlation in cond 3 split by order of cond 2 and 3
subj_cond3_first = {'cm', 'lp', 'cr', 'zt', 'mb', 'sar'};

figure
hold on
plot([minPlot_corr maxPlot_corr], [minPlot_corr maxPlot_corr], 'k')
for ii = 1 : length(subjectAll)
    plot([corr_Cond1(3, ii) corr_Cond1(4, ii)], [corr_Cond3(1, ii) corr_Cond3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([corr_Cond1(1, ii) corr_Cond1(1, ii)], [corr_Cond3(3, ii) corr_Cond3(4, ii)], '-k', 'LineWidth', lineWidth) 
    if ~isempty(find(strcmp(subj_cond3_first, subjectAll{ii}), 1))
        plot(corr_Cond1(1, ii), corr_Cond3(1, ii), 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
    else
        plot(corr_Cond1(1, ii), corr_Cond3(1, ii), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)        
    end
end
axis([minPlot_corr, maxPlot_corr, minPlot_corr, maxPlot_corr])
xlabel('Correlation - Cond1 (deg)')
ylabel('Correlation - Cond3 (deg)')

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

% Std of estimate
p_std_1vs2= signrank(std_Cond1, std_Cond2,...
                                    'method', 'exact');
fprintf('p-value std cond 1 vs cond 2: %8.6f \n', p_std_1vs2) 
p_std_1vs3= signrank(std_Cond1, std_Cond3,...
                                    'method', 'exact');
fprintf('p-value std cond 1 vs cond 3: %8.6f \n', p_std_1vs3)    

% Mean diff 1st vs. 2nd estimate
p_meanDiff_1vs2_cond2 = signrank(meanDiff_Cond2_split(1, :), meanDiff_Cond2_split(2, :),...
                                    'method', 'exact');
fprintf('p-value mean diff 1st vs. 2nd estimate, cond 2: %8.6f \n', p_meanDiff_1vs2_cond2) 
p_meanDiff_1vs2_cond3 = signrank(meanDiff_Cond3_split(1, :), meanDiff_Cond3_split(2, :),...
                                    'method', 'exact');
fprintf('p-value mean diff 1st vs. 2nd estimate, cond 3: %8.6f \n', p_meanDiff_1vs2_cond3) 
