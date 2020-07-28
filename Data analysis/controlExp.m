%%%%%%%%%%% Analyze control experiment %%%%%%%%%%%%%%%
subjectAll = {'ll', 'ad', 'rl'};
n_subj = length(subjectAll);
meanDiffEst_same = NaN(1, length(subjectAll));
meanDiffEst_diff = NaN(1, length(subjectAll));
nBootstrap = 10000;
excludeIncorrectTrial = 2; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
alpha = 0.05;                           
                           
h1 = figure;
h2 = figure;
markerSize = 8;
lineWidth = 0.5;
hold on
estimateStim1_1line_all = cell(1, n_subj);
estimateStim2_1line_all = cell(1, n_subj);
estimateStim1_2line_all = cell(1, n_subj);
estimateStim2_2line_all = cell(1, n_subj);
meanDiffEst_Cond1 = NaN(4, length(subjectAll)+1);
meanDiffEst_Cond2 = NaN(4, length(subjectAll)+1);
corr_Cond1 = NaN(4, length(subjectAll)+1);
corr_Cond2 = NaN(4, length(subjectAll)+1);
percentCorrect_Cond1 = NaN(4, length(subjectAll)+1);
percentCorrect_Cond2 = NaN(4, length(subjectAll)+1);

%% Extract summary statistics for individual subject
for ss = 1 : n_subj
    subject = subjectAll{ss};
    session = 1;

    %% First line 49 deg
    experimentNumber = 1;
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
    estimateStim1_1line = dataResponse(:, 3);
    estimateStim1_1line(isnan(estimateStim1_1line)) = [];
    estimateStim1_1line_all{ss} = estimateStim1_1line;

    % Extract the stimulus orientation
    stimOrientation = params.lineOrientation;
    
    %% Second line 54 deg
    experimentNumber = 1;
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
    estimateStim2_1line = dataResponse(:, 4);
    estimateStim2_1line(isnan(estimateStim2_1line)) = [];
    estimateStim2_1line_all{ss} = estimateStim2_1line;
    mean_diff_1side = mean(estimateStim2_1line) - mean(estimateStim1_1line);
    mean_std_1side = mean([std(estimateStim1_1line) std(estimateStim2_1line)]);

    % Bootstrap the mean difference, correlation and percent correct
    experiment_condition = 1;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_1line, estimateStim2_1line);

    if excludeIncorrectTrial == 2            
        meanDiffEst_Cond1(1, ss) = sampleStat_all(1);
    else
        meanDiffEst_Cond1(1, ss) = mean(estimateStim2_1line) - mean(estimateStim1_1line);
    end
    meanDiffEst_Cond1(3:4, ss) = confInterval_all(:, 1);
    corr_Cond1(1, ss) = sampleStat_all(2);
    corr_Cond1(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Cond1(1, ss) = sampleStat_all(3);
    percentCorrect_Cond1(3:4, ss) = confInterval_all(:, 3);
    
    %% Show 2 lines, report 1
    experimentNumber = 1;
    experimentName = 'HighToLow_1lineShow2_1side';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    if exist(dataFile, 'file') == 2
        % Load data
        load(dataFile)

        stimulus_orientation = dataResponse(:, params.reportWhichLine);
        estimateStim1_2side = dataResponse(stimulus_orientation==49, params.reportWhichLine+2);
        estimateStim2_2side = dataResponse(stimulus_orientation==54, params.reportWhichLine+2);
        mean_diff_2side = mean(estimateStim2_2side) - mean(estimateStim1_2side);
        mean_std_2side = mean([std(estimateStim1_2side) std(estimateStim2_2side)]);
    else
        estimateStim1_2side = NaN;
        estimateStim2_2side = NaN;
        mean_diff_2side = NaN;
        mean_std_2side = NaN;
    end
    
    %% Show 2 lines, report 2
    experimentNumber = 1;
    experimentName = 'HighToLow_separate_fixed';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    % Load data
    load(dataFile)

    stimulus_orientation = [dataResponse(1, 1), dataResponse(1, 2)];
    if params.leftLine49 == 1
        estimateStim1_2line = dataResponse(:, 3);
        estimateStim2_2line = dataResponse(:, 4);
    else
        estimateStim1_2line = dataResponse(:, 4);
        estimateStim2_2line = dataResponse(:, 3);        
    end
    estimateStim1_2line_all{ss} = estimateStim1_2line;
    estimateStim2_2line_all{ss} = estimateStim2_2line;
    mean_diff_2line = mean(estimateStim2_2line) - mean(estimateStim1_2line);
    mean_std_2line = mean([std(estimateStim2_2line) std(estimateStim1_2line)]); 

    % Bootstrap the mean difference, correlation and percent correct
    experiment_condition = 2;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_2line, estimateStim2_2line);

    meanDiffEst_Cond2(1, ss) = sampleStat_all(1);    
    meanDiffEst_Cond2(3:4, ss) = confInterval_all(:, 1);
    corr_Cond2(1, ss) = sampleStat_all(2);
    corr_Cond2(3:4, ss) = confInterval_all(:, 2);
    percentCorrect_Cond2(1, ss) = sampleStat_all(3);
    percentCorrect_Cond2(3:4, ss) = confInterval_all(:, 3);
    
    %% Plot estimates
    max_x = nanmax([estimateStim1_1line; estimateStim2_1line; estimateStim1_2side; estimateStim2_2side; estimateStim1_2line; estimateStim2_2line])+1;
    min_x = nanmin([estimateStim1_1line; estimateStim2_1line; estimateStim1_2side; estimateStim2_2side; estimateStim1_2line; estimateStim2_2line])-1;
    bin_edge = min_x:max_x;

    figure(h1)
    subplot(n_subj, 3, 3*(ss-1)+1)
    hold on
    histogram(estimateStim1_1line, bin_edge)
    histogram(estimateStim2_1line, bin_edge)
    xlim([min_x, max_x])
    xlabel('Reported orientation (deg)')
    ylabel('Count')
    title(['Mean diff: ' num2str(mean_diff_1side, 3) ', Std: ' num2str(mean_std_1side, 3)])
    if ss == 1
        legend('49 deg line', '54 deg line')
    end

    subplot(n_subj, 3, 3*(ss-1)+2)
    hold on
    histogram(estimateStim1_2side, bin_edge)
    histogram(estimateStim2_2side, bin_edge)
    xlim([min_x, max_x])
    xlabel('Reported orientation (deg)')
    ylabel('Count')
    title(['Mean diff: ' num2str(mean_diff_2side, 3) ', Std: ' num2str(mean_std_2side, 3)])
    
    subplot(n_subj, 3, 3*(ss-1)+3)
    hold on
    histogram(estimateStim1_2line, bin_edge)
    histogram(estimateStim2_2line, bin_edge)
    xlim([min_x, max_x])
    xlabel('Reported orientation (deg)')
    ylabel('Count')
    title(['Mean diff: ' num2str(mean_diff_2line, 3) ', Std: ' num2str(mean_std_2line, 3)])   
    
    figure(h2)
    subplot(1, n_subj, ss)
    hold on
    plot([min_x, max_x], [min_x, max_x])
    plot(estimateStim1_2line, estimateStim2_2line, 'o',  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize-3)
    plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
    plot(mean(estimateStim1_2line), mean(estimateStim2_2line), 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')           
    xlim([min_x, max_x])
    ylim([min_x, max_x])
    xlabel('Reported orientation 49 deg line')
    ylabel('Reported orientation 54 deg line')
end

%% Extract summary statistics pooled across all subjects
% Bootstrap for 1 line condition
estimateStim1_1line_all = cell2mat(estimateStim1_1line_all);
estimateStim1_1line_all = estimateStim1_1line_all(:);
estimateStim2_1line_all = cell2mat(estimateStim2_1line_all);
estimateStim2_1line_all = estimateStim2_1line_all(:);

experiment_condition = 1;
[confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                experiment_condition, alpha, estimateStim1_1line_all, estimateStim2_1line_all);

if excludeIncorrectTrial == 2            
    meanDiffEst_Cond1(1, ss+1) = sampleStat_all(1);
else
    meanDiffEst_Cond1(1, ss+1) = mean(estimateStim2_1line_all) - mean(estimateStim1_1line_all);
end
meanDiffEst_Cond1(3:4, ss+1) = confInterval_all(:, 1);
corr_Cond1(1, ss+1) = sampleStat_all(2);
corr_Cond1(3:4, ss+1) = confInterval_all(:, 2);
percentCorrect_Cond1(1, ss+1) = sampleStat_all(3);
percentCorrect_Cond1(3:4, ss+1) = confInterval_all(:, 3);

% Bootstrap for 2 line condition
estimateStim1_2line_all = cell2mat(estimateStim1_2line_all);
estimateStim1_2line_all = estimateStim1_2line_all(:);
estimateStim2_2line_all = cell2mat(estimateStim2_2line_all);
estimateStim2_2line_all = estimateStim2_2line_all(:);

experiment_condition = 2;
[confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                experiment_condition, alpha, estimateStim1_2line_all, estimateStim2_2line_all);

if excludeIncorrectTrial == 2            
    meanDiffEst_Cond2(1, ss+1) = sampleStat_all(1);
else
    meanDiffEst_Cond2(1, ss+1) = mean(estimateStim2_2line_all) - mean(estimateStim1_2line_all);
end
meanDiffEst_Cond2(3:4, ss+1) = confInterval_all(:, 1);
corr_Cond2(1, ss+1) = sampleStat_all(2);
corr_Cond2(3:4, ss+1) = confInterval_all(:, 2);
percentCorrect_Cond2(1, ss+1) = sampleStat_all(3);
percentCorrect_Cond2(3:4, ss+1) = confInterval_all(:, 3);

%% Plot the combined subject
mean_diff = mean(estimateStim2_1line_all) - mean(estimateStim1_1line_all);
mean_diff_2line = mean(estimateStim2_2line_all) - mean(estimateStim1_2line_all);

figure
hold on
max_x = max([estimateStim1_1line_all; estimateStim2_1line_all;  estimateStim1_2line_all; estimateStim2_2line_all])+1;
min_x = min([estimateStim1_1line_all; estimateStim2_1line_all;  estimateStim1_2line_all; estimateStim2_2line_all])-1;
bin_edge = min_x:max_x;

subplot(2, 2, 1)
hold on
histogram(estimateStim1_1line_all, bin_edge)
histogram(estimateStim2_1line_all, bin_edge)
xlim([min_x, max_x])
xlabel('Reported orientation (deg)')
ylabel('Count')
title(['Mean diff: ' num2str(mean_diff, 3)])
legend('49 deg line', '54 deg line')


subplot(2, 2, 3)
hold on
histogram(estimateStim1_2line_all, bin_edge)
histogram(estimateStim2_2line_all, bin_edge)
xlim([min_x, max_x])
xlabel('Reported orientation (deg)')
ylabel('Count')
title(['Mean diff: ' num2str(mean_diff_2line, 3)])


subplot(2, 2, 2)
hold on
plot([min_x, max_x], [min_x, max_x])
plot(estimateStim1_2line_all, estimateStim2_2line_all, 'o',  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize-3)
plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_2line_all), mean(estimateStim2_2line_all), 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')           
xlim([min_x, max_x])
ylim([min_x, max_x])
xlabel('Reported orientation 49 deg line')
ylabel('Reported orientation 54 deg line')

subplot(2, 2, 4)
hold on
diff_est_2line = estimateStim2_2line_all - estimateStim1_2line_all;
min_x = min(diff_est_2line)-1;
max_x = max(diff_est_2line)+1;
bin_edge = min_x:max_x;
histogram(diff_est_2line, bin_edge)
xlim([min_x, max_x])
xlabel('Reported orientation (deg)')
ylabel('Count')

%% Plot individual subject
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ss = 1 : length(colorName)
    colorIndex(ss, :) = rgb(colorName{ss});
end

figure
hold on
subplot(2, 3, 1)
hold on
minPlot_diffEst = min([meanDiffEst_Cond1(1, :) meanDiffEst_Cond2(1, :)]) - 2;
maxPlot_diffEst = max([meanDiffEst_Cond1(1, :) meanDiffEst_Cond2(1, :)]) + 2;
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
subjectAll{4} = 'pooled';
legend_handle = NaN(1, length(subjectAll));
for ss = 1 : length(subjectAll)
    plot([meanDiffEst_Cond1(3, ss) meanDiffEst_Cond1(4, ss)], [meanDiffEst_Cond2(1, ss) meanDiffEst_Cond2(1, ss)], '-k', 'LineWidth', lineWidth)
    plot([meanDiffEst_Cond1(1, ss) meanDiffEst_Cond1(1, ss)], [meanDiffEst_Cond2(3, ss) meanDiffEst_Cond2(4, ss)], '-k', 'LineWidth', lineWidth) 
    if ss == length(subjectAll)
        legend_handle(ss) = plot(meanDiffEst_Cond1(1, ss), meanDiffEst_Cond2(1, ss), '^', 'MarkerFaceColor', colorIndex(ss, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);        
    else
        legend_handle(ss) = plot(meanDiffEst_Cond1(1, ss), meanDiffEst_Cond2(1, ss), 'o', 'MarkerFaceColor', colorIndex(ss, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    end
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
minPlot_corr = min([corr_Cond1(1, :) corr_Cond2(1, :)]) - 0.4;
maxPlot_corr = max([corr_Cond1(1, :) corr_Cond2(1, :)]) + 0.1;
plot([minPlot_corr maxPlot_corr], [minPlot_corr maxPlot_corr], 'k')
for ss = 1 : length(subjectAll)
    plot([corr_Cond1(3, ss) corr_Cond1(4, ss)], [corr_Cond2(1, ss) corr_Cond2(1, ss)], '-k', 'LineWidth', lineWidth)
    plot([corr_Cond1(1, ss) corr_Cond1(1, ss)], [corr_Cond2(3, ss) corr_Cond2(4, ss)], '-k', 'LineWidth', lineWidth)    
    if ss == length(subjectAll)
        legend_handle(ss) = plot(corr_Cond1(1, ss), corr_Cond2(1, ss), '^', 'MarkerFaceColor', colorIndex(ss, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);        
    else
        legend_handle(ss) = plot(corr_Cond1(1, ss), corr_Cond2(1, ss), 'o', 'MarkerFaceColor', colorIndex(ss, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    end
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
minPlot_pc = min([percentCorrect_Cond1(1, :) percentCorrect_Cond2(1, :)]) - 15;
maxPlot_pc = 110;
plot([minPlot_pc maxPlot_pc], [minPlot_pc maxPlot_pc])
for ss = 1 : length(subjectAll)
    plot([percentCorrect_Cond1(3, ss) percentCorrect_Cond1(4, ss)], [percentCorrect_Cond2(1, ss) percentCorrect_Cond2(1, ss)], '-k', 'LineWidth', lineWidth)
    plot([percentCorrect_Cond1(1, ss) percentCorrect_Cond1(1, ss)], [percentCorrect_Cond2(3, ss) percentCorrect_Cond2(4, ss)], '-k', 'LineWidth', lineWidth)
    if ss == length(subjectAll)
        legend_handle(ss) = plot(percentCorrect_Cond1(1, ss), percentCorrect_Cond2(1, ss), '^', 'MarkerFaceColor', colorIndex(ss, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);        
    else
        legend_handle(ss) = plot(percentCorrect_Cond1(1, ss), percentCorrect_Cond2(1, ss), 'o', 'MarkerFaceColor', colorIndex(ss, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    end
end
axis([minPlot_pc, maxPlot_pc, minPlot_pc, maxPlot_pc])
axis square
% legend(legend_pc1)
xlabel('Percent correct - Cond1 (deg)')
ylabel('Percent correct - Cond2 (deg)')

subplot(2, 3, 4)
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
