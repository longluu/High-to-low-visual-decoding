%%%%%%%%%%% Analyze data from the experiment with Matlab buit-in bootstrap function %%%%%%%%%%%%%%%
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
nBootstrap = 10000;
markerSize = 8;
lineWidth = 0.5;
eps_correlation = 0.0001;
plotIndividualSubject = 1;
excludeIncorrectTrial = 1; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
meanDiffEst_Exp1 = NaN(4, length(subjectAll));
meanDiffEst_Exp2 = NaN(4, length(subjectAll));
meanDiffEst_Exp3 = NaN(4, length(subjectAll));
corr_Exp1 = NaN(4, length(subjectAll));
corr_Exp2 = NaN(4, length(subjectAll));
corr_Exp3 = NaN(4, length(subjectAll));
percentCorrect_Exp1 = NaN(4, length(subjectAll));
percentCorrect_Exp2 = NaN(4, length(subjectAll));
percentCorrect_Exp3 = NaN(4, length(subjectAll));
std_Exp2 = NaN(2, length(subjectAll));
std_Exp3 = NaN(2, length(subjectAll));
pValue_meanDiff = NaN(2, length(subjectAll));
pValue_corr = NaN(2, length(subjectAll));
pValue_pc = NaN(2, length(subjectAll));
legend_meanDiff1 = cell(1, length(subjectAll));
legend_meanDiff2 = cell(1, length(subjectAll));
legend_corr1 = cell(1, length(subjectAll));
legend_corr2 = cell(1, length(subjectAll));
legend_pc1 = cell(1, length(subjectAll));
legend_pc2 = cell(1, length(subjectAll));

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
    
    % Find the permutation of data that has lowest correlation
    corr_2line = 1;
    while abs(corr_2line) > eps_correlation
        est1 = datasample(estimateStim1_collapse, length(estimateStim1_collapse), 'Replace', false);
        est2 = datasample(estimateStim2_collapse, length(estimateStim2_collapse), 'Replace', false);
        corr_2line = corr(est1, est2);
    end
    estimateStim1_collapse = est1;
    estimateStim2_collapse = est2;
    
    % Bootstrap the mean difference, correlation and percent correct
    if excludeIncorrectTrial == 0
        [confInterval_all, bootStat_all] = bootci(nBootstrap, @bootstrap_Include, [estimateStim1_collapse estimateStim2_collapse]);
    elseif excludeIncorrectTrial == 1
        [confInterval_all, bootStat_all] = bootci(nBootstrap, @bootstrap_Exclude, [estimateStim1_collapse estimateStim2_collapse]);        
    end
    meanDiffEst_Exp1(1, ss) = mean(bootStat_all(:, 1));    
    meanDiffEst_Exp1(3:4, ss) = confInterval_all(:, 1);
    corr_Exp1(1, ss) = mean(bootStat_all(:, 2));
    corr_Exp1(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Exp1(1, ss) = 100 * mean(bootStat_all(:, 3));
    percentCorrect_Exp1(3:4, ss) = 100 * confInterval_all(:, 3);
    
    
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
    diffEst_original = estimateStim2_collapse - estimateStim1_collapse;    
    
    % Std of first and second estimate
    std_estStim1Est1 = nanstd(estimateLine1(indStim1Est1));
    std_estStim2Est1 = nanstd(estimateLine1(indStim2Est1));
    std_est1_exp2 = sqrt((std_estStim1Est1^2 + std_estStim2Est1^2) / 2);
    std_estStim1Est2 = nanstd(estimateLine2(indStim1Est2));
    std_estStim2Est2 = nanstd(estimateLine2(indStim2Est2));
    std_est2_exp2 = sqrt((std_estStim1Est2^2 + std_estStim2Est2^2) / 2);
      
    % Bootstrap the mean difference, correlation and percent correct
    if excludeIncorrectTrial == 0
        [confInterval_all, bootStat_all] = bootci(nBootstrap, @bootstrap_Include, [estimateStim1_collapse estimateStim2_collapse]);
    elseif excludeIncorrectTrial == 1
        [confInterval_all, bootStat_all] = bootci(nBootstrap, @bootstrap_Exclude, [estimateStim1_collapse estimateStim2_collapse]);        
    end
    if excludeIncorrectTrial == 1
        meanDiffEst_Exp2(1, ss) = mean(estimateStim2_collapse(diffEst_original>0)) - mean(estimateStim1_collapse(diffEst_original>0));
        corr_Exp2(1, ss) = corr(estimateStim1_collapse(diffEst_original>0), estimateStim2_collapse(diffEst_original>0));        
    else
        meanDiffEst_Exp2(1, ss) = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
        corr_Exp2(1, ss) = corr(estimateStim1_collapse, estimateStim2_collapse);                
    end    
    meanDiffEst_Exp2(3:4, ss) = confInterval_all(:, 1);
    corr_Exp2(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Exp2(1, ss) = 100 * sum(estimateStim2_collapse >= estimateStim1_collapse) / length(estimateStim2_collapse);
    percentCorrect_Exp2(3:4, ss) = 100 * confInterval_all(:, 3);
    
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
    diffEst_original = estimateStim2_collapse - estimateStim1_collapse;    
    
    % Std of first and second estimate
    std_estStim1Est1 = nanstd(estimateLine1(indStim1Est1));
    std_estStim2Est1 = nanstd(estimateLine1(indStim2Est1));
    std_est1_exp3 = sqrt((std_estStim1Est1^2 + std_estStim2Est1^2) / 2);
    std_estStim1Est2 = nanstd(estimateLine2(indStim1Est2));
    std_estStim2Est2 = nanstd(estimateLine2(indStim2Est2));
    std_est2_exp3 = sqrt((std_estStim1Est2^2 + std_estStim2Est2^2) / 2);

    % Bootstrap the mean difference, correlation and percent correct
    if excludeIncorrectTrial == 0
        [confInterval_all, bootStat_all] = bootci(nBootstrap, @bootstrap_Include, [estimateStim1_collapse estimateStim2_collapse]);
    elseif excludeIncorrectTrial == 1
        [confInterval_all, bootStat_all] = bootci(nBootstrap, @bootstrap_Exclude, [estimateStim1_collapse estimateStim2_collapse]);        
    end
    if excludeIncorrectTrial == 1
        meanDiffEst_Exp3(1, ss) = mean(estimateStim2_collapse(diffEst_original>0)) - mean(estimateStim1_collapse(diffEst_original>0));
        corr_Exp3(1, ss) = corr(estimateStim1_collapse(diffEst_original>0), estimateStim2_collapse(diffEst_original>0));        
    else
        meanDiffEst_Exp3(1, ss) = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
        corr_Exp3(1, ss) = corr(estimateStim1_collapse, estimateStim2_collapse);                
    end    
    meanDiffEst_Exp3(3:4, ss) = confInterval_all(:, 1);
    corr_Exp3(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Exp3(1, ss) = 100 * sum(estimateStim2_collapse >= estimateStim1_collapse) / length(estimateStim2_collapse);
    percentCorrect_Exp3(3:4, ss) = 100 * confInterval_all(:, 3);

    %% Compute the std for 1st vs 2nd report
    std_Exp2(1, ss) = std_est1_exp2;
    std_Exp2(2, ss) = std_est2_exp2;
    std_Exp3(1, ss) = std_est1_exp3;
    std_Exp3(2, ss) = std_est2_exp3; 
    
%     %% Compute the p-value
%     pValue_meanDiff(1, ss) = sum(meanDiff_boot2 - meanDiff_boot1 > meanDiffEst_Exp2(1, ss) - meanDiffEst_Exp1(1, ss)) / nBootstrap;
%     pValue_meanDiff(2, ss) = sum(abs(meanDiff_boot3 - meanDiff_boot1) > abs(meanDiffEst_Exp3(1, ss) - meanDiffEst_Exp1(1, ss))) / nBootstrap;
%     pValue_corr(1, ss) = sum(abs(corr_boot2 - corr_boot1) > abs(corr_Exp2(1, ss) - corr_Exp1(1, ss))) / nBootstrap;
%     pValue_corr(2, ss) = sum(abs(corr_boot3 - corr_boot1) > abs(corr_Exp3(1, ss) - corr_Exp1(1, ss))) / nBootstrap;
%     pValue_pc(1, ss) = sum(abs(percentCorrect_boot2 - percentCorrect_boot1) > abs(percentCorrect_Exp2(1, ss) - percentCorrect_Exp1(1, ss))) / nBootstrap;
%     pValue_pc(2, ss) = sum(abs(percentCorrect_boot3 - percentCorrect_boot1) > abs(percentCorrect_Exp3(1, ss) - percentCorrect_Exp1(1, ss))) / nBootstrap;    
%     
%     legend_meanDiff1{ss} = num2str(round(pValue_meanDiff(1, ss), 3));
%     legend_meanDiff2{ss} = num2str(round(pValue_meanDiff(2, ss), 3));
%     legend_corr1{ss} = num2str(round(pValue_corr(1, ss), 3));
%     legend_corr2{ss} = num2str(round(pValue_corr(2, ss), 3));
%     legend_pc1{ss} = num2str(round(pValue_pc(1, ss), 3));
%     legend_pc2{ss} = num2str(round(pValue_pc(2, ss), 3));    
end

%% Plot the results
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'Purple', 'Magenta', 'Gray'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
hold on
subplot(2, 3, 1)
hold on
minPlot = min([meanDiffEst_Exp1(1, :) meanDiffEst_Exp2(1, :)]) - 2;
maxPlot = max([meanDiffEst_Exp1(1, :) meanDiffEst_Exp2(1, :)]) + 2;
legend_handle = NaN(1, length(subjectAll));
for ii = 1 : length(subjectAll)
    plot([meanDiffEst_Exp1(3, ii) meanDiffEst_Exp1(4, ii)], [meanDiffEst_Exp2(1, ii) meanDiffEst_Exp2(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([meanDiffEst_Exp1(1, ii) meanDiffEst_Exp1(1, ii)], [meanDiffEst_Exp2(3, ii) meanDiffEst_Exp2(4, ii)], '-k', 'LineWidth', lineWidth)    
    legend_handle(ii) = plot(meanDiffEst_Exp1(1, ii), meanDiffEst_Exp2(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
end
plot([minPlot maxPlot], [minPlot maxPlot])
if excludeIncorrectTrial == 1
    title('Correct trials only')
else
    title('All trials')
end
% legend(legend_handle, legend_meanDiff1)
xlabel('Mean difference - Exp1 (deg)')
ylabel('Mean difference - Exp2 (deg)')
axis equal

subplot(2, 3, 2)
hold on
minPlot = min([corr_Exp1(1, :) corr_Exp2(1, :)]) - 0.1;
maxPlot = max([corr_Exp1(1, :) corr_Exp2(1, :)]) + 0.1;
for ii = 1 : length(subjectAll)
    plot([corr_Exp1(3, ii) corr_Exp1(4, ii)], [corr_Exp2(1, ii) corr_Exp2(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([corr_Exp1(1, ii) corr_Exp1(1, ii)], [corr_Exp2(3, ii) corr_Exp2(4, ii)], '-k', 'LineWidth', lineWidth)    
    plot(corr_Exp1(1, ii), corr_Exp2(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
if excludeIncorrectTrial == 1
    title('Correct trials only')
else
    title('All trials')
end
% legend(legend_corr1)
xlabel('Correlation - Exp1 (deg)')
ylabel('Correlation - Exp2 (deg)')
axis equal

subplot(2, 3, 3)
hold on
minPlot = min([percentCorrect_Exp1(1, :) percentCorrect_Exp2(1, :)]) - 10;
maxPlot = 100;
for ii = 1 : length(subjectAll)
    plot([percentCorrect_Exp1(3, ii) percentCorrect_Exp1(4, ii)], [percentCorrect_Exp2(1, ii) percentCorrect_Exp2(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([percentCorrect_Exp1(1, ii) percentCorrect_Exp1(1, ii)], [percentCorrect_Exp2(3, ii) percentCorrect_Exp2(4, ii)], '-k', 'LineWidth', lineWidth)
    plot(percentCorrect_Exp1(1, ii), percentCorrect_Exp2(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
% legend(legend_pc1)
xlabel('Percent correct - Exp1 (deg)')
ylabel('Percent correct - Exp2 (deg)')
axis equal

subplot(2, 3, 4)
hold on
minPlot = min([meanDiffEst_Exp1(1, :) meanDiffEst_Exp3(1, :)]) - 2;
maxPlot = max([meanDiffEst_Exp1(1, :) meanDiffEst_Exp3(1, :)]) + 2;
for ii = 1 : length(subjectAll)
    plot([meanDiffEst_Exp1(3, ii) meanDiffEst_Exp1(4, ii)], [meanDiffEst_Exp3(1, ii) meanDiffEst_Exp3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([meanDiffEst_Exp1(1, ii) meanDiffEst_Exp1(1, ii)], [meanDiffEst_Exp3(3, ii) meanDiffEst_Exp3(4, ii)], '-k', 'LineWidth', lineWidth)
    plot(meanDiffEst_Exp1(1, ii), meanDiffEst_Exp3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
if excludeIncorrectTrial == 1
    title('Correct trials only')
else
    title('All trials')
end
% legend(legend_meanDiff2)
xlabel('Mean difference - Exp1 (deg)')
ylabel('Mean difference - Exp3 (deg)')
axis equal

subplot(2, 3, 5)
hold on
minPlot = min([corr_Exp1(1, :) corr_Exp3(1, :)]) - 0.1;
maxPlot = max([corr_Exp1(1, :) corr_Exp3(1, :)]) + 0.1;
for ii = 1 : length(subjectAll)
    plot([corr_Exp1(3, ii) corr_Exp1(4, ii)], [corr_Exp3(1, ii) corr_Exp3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([corr_Exp1(1, ii) corr_Exp1(1, ii)], [corr_Exp3(3, ii) corr_Exp3(4, ii)], '-k', 'LineWidth', lineWidth)    
    plot(corr_Exp1(1, ii), corr_Exp3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
if excludeIncorrectTrial == 1
    title('Correct trials only')
else
    title('All trials')
end
% legend(legend_corr2)
xlabel('Correlation - Exp1 (deg)')
ylabel('Correlation - Exp3 (deg)')
axis equal

subplot(2, 3, 6)
hold on
minPlot = min([percentCorrect_Exp1(1, :) percentCorrect_Exp3(1, :)]) - 10;
maxPlot = 100;
for ii = 1 : length(subjectAll)
    plot([percentCorrect_Exp1(3, ii) percentCorrect_Exp1(4, ii)], [percentCorrect_Exp3(1, ii) percentCorrect_Exp3(1, ii)], '-k', 'LineWidth', lineWidth)
    plot([percentCorrect_Exp1(1, ii) percentCorrect_Exp1(1, ii)], [percentCorrect_Exp3(3, ii) percentCorrect_Exp3(4, ii)], '-k', 'LineWidth', lineWidth)    
    plot(percentCorrect_Exp1(1, ii), percentCorrect_Exp3(1, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
% legend(legend_pc2)
xlabel('Percent correct - Exp1 (deg)')
ylabel('Percent correct - Exp3 (deg)')
axis equal

figure;
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
subplot(1, 2, 1)
hold on
minPlot = min([std_Exp2(1, :) std_Exp2(2, :)]) - 0.1;
maxPlot = max([std_Exp2(1, :) std_Exp2(2, :)]) + 0.1;
for ii = 1 : length(subjectAll)
    plot(std_Exp2(1, ii), std_Exp2(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
title('Experiment 2')
xlabel('Std of first estimate (deg)')
ylabel('Std of second estimate (deg)')
axis equal

subplot(1, 2, 2)
hold on
minPlot = min([std_Exp3(1, :) std_Exp3(2, :)]) - 0.1;
maxPlot = max([std_Exp3(1, :) std_Exp3(2, :)]) + 0.1;
for ii = 1 : length(subjectAll)
    plot(std_Exp3(1, ii), std_Exp3(2, ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
plot([minPlot maxPlot], [minPlot maxPlot])
title('Experiment 3')
xlabel('Std of first estimate (deg)')
ylabel('Std of second estimate (deg)')
axis equal

%% Plot individual subjects
markerSize = 4;
if plotIndividualSubject
    for ss = 1 : length(subjectAll)
        subject = subjectAll{ss};
        figure
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
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none');    
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title(['Correlation: ' num2str(round(corr(estimateStim1_collapse, estimateStim2_collapse), 2))...
            ', PC: ' num2str(round(max(percentCorrect)*100)) '%'])

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
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none')
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title(['Correlation: ' num2str(round(corrCoef, 2))...
               ', p-value: ' num2str(round(pValue, 2))...
               ', PC: ' num2str(round(percentCorrect*100)) '%'])

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
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none')
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title(['Correlation: ' num2str(round(corrCoef, 2))...
               ', p-value: ' num2str(round(pValue, 2))...
               ', PC: ' num2str(round(percentCorrect*100)) '%'])

        subplot(2, 3, 6)
        histogram(diffEst, 20, 'Normalization', 'probability')
        xlabel('Orientation difference (deg)')
        ylabel('Frequency of occurence')
        xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
        axis square
        title(['Mean: ' num2str(round(mean(diffEst), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))])    
    end
end