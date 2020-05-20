%%%%%%%%%%% Analyze control experiment %%%%%%%%%%%%%%%
subjectAll = {'ll', 'ad', 'rl'};
n_subj = length(subjectAll);
meanDiffEst_same = NaN(1, length(subjectAll));
meanDiffEst_diff = NaN(1, length(subjectAll));
markerSize = 6;

h1 = figure;
h2 = figure;
hold on
estimateStim1_all = cell(1, n_subj);
estimateStim2_all = cell(1, n_subj);
estimateStim1_2line_all = cell(1, n_subj);
estimateStim2_2line_all = cell(1, n_subj);

for ii = 1 : n_subj
    subject = subjectAll{ii};
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
    estimateStim1 = dataResponse(:, 3);
    estimateStim1(isnan(estimateStim1)) = [];
    estimateStim1_all{ii} = estimateStim1;

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
    estimateStim2 = dataResponse(:, 4);
    estimateStim2(isnan(estimateStim2)) = [];
    estimateStim2_all{ii} = estimateStim2;
    mean_diff_1side = mean(estimateStim2) - mean(estimateStim1);
    mean_std_1side = mean([std(estimateStim1) std(estimateStim2)]);

    %% Show 2 lines, report 1
    experimentNumber = 1;
    experimentName = 'HighToLow_1lineShow2_1side';
    dataFile = ['Data\' subject '\MainExperiment\' experimentName num2str(session) '\' experimentName '-'  num2str(experimentNumber) '.mat'];

    % Load data
    load(dataFile)

    stimulus_orientation = dataResponse(:, params.reportWhichLine);
    estimateStim1_2side = dataResponse(stimulus_orientation==49, params.reportWhichLine+2);
    estimateStim2_2side = dataResponse(stimulus_orientation==54, params.reportWhichLine+2);
    mean_diff_2side = mean(estimateStim2_2side) - mean(estimateStim1_2side);
    mean_std_2side = mean([std(estimateStim1_2side) std(estimateStim2_2side)]);

    % Split subjects' estimates based on previous trial (same or different) 
    trialOrder = params.trialOrder;
    trial_current = trialOrder(2:end);
    trial_last = trialOrder(1:end-1);
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);

    estimateStimAll = dataResponse(:, params.reportWhichLine+2);
    estimateStim1_same = estimateStimAll;
    estimateStim1_same(stimulus_orientation==54) = NaN;
    estimateStim1_same(trial_diff) = NaN;
    estimateStim1_diff = estimateStimAll;
    estimateStim1_diff(stimulus_orientation==54) = NaN;
    estimateStim1_diff(trial_same) = NaN;

    estimateStim2_same = estimateStimAll;
    estimateStim2_same(stimulus_orientation==49) = NaN;
    estimateStim2_same(trial_diff) = NaN;
    estimateStim2_diff = estimateStimAll;
    estimateStim2_diff(stimulus_orientation==49) = NaN;
    estimateStim2_diff(trial_same) = NaN;
    
    meanDiffEst_same(ii) = nanmean(estimateStim2_same) - nanmean(estimateStim1_same);
    meanDiffEst_diff(ii) = nanmean(estimateStim2_diff) - nanmean(estimateStim1_diff);
    
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
    estimateStim1_2line_all{ii} = estimateStim1_2line;
    estimateStim2_2line_all{ii} = estimateStim2_2line;
    mean_diff_2line = mean(estimateStim2_2line) - mean(estimateStim1_2line);
    mean_std_2line = mean([std(estimateStim2_2line) std(estimateStim1_2line)]); 
    
    %% Plot estimates
    max_x = max([estimateStim1; estimateStim2; estimateStim1_2side; estimateStim2_2side; estimateStim1_2line; estimateStim2_2line])+1;
    min_x = min([estimateStim1; estimateStim2; estimateStim1_2side; estimateStim2_2side; estimateStim1_2line; estimateStim2_2line])-1;
    bin_edge = min_x:max_x;

    figure(h1)
    subplot(n_subj, 3, 3*(ii-1)+1)
    hold on
    histogram(estimateStim1, bin_edge)
    histogram(estimateStim2, bin_edge)
    xlim([min_x, max_x])
    xlabel('Reported orientation (deg)')
    ylabel('Count')
    title(['Mean diff: ' num2str(mean_diff_1side, 3) ', Std: ' num2str(mean_std_1side, 3)])
    if ii == 1
        legend('49 deg line', '54 deg line')
    end

    subplot(n_subj, 3, 3*(ii-1)+2)
    hold on
    histogram(estimateStim1_2side, bin_edge)
    histogram(estimateStim2_2side, bin_edge)
    xlim([min_x, max_x])
    xlabel('Reported orientation (deg)')
    ylabel('Count')
    title(['Mean diff: ' num2str(mean_diff_2side, 3) ', Std: ' num2str(mean_std_2side, 3)])
    
    subplot(n_subj, 3, 3*(ii-1)+3)
    hold on
    histogram(estimateStim1_2line, bin_edge)
    histogram(estimateStim2_2line, bin_edge)
    xlim([min_x, max_x])
    xlabel('Reported orientation (deg)')
    ylabel('Count')
    title(['Mean diff: ' num2str(mean_diff_2line, 3) ', Std: ' num2str(mean_std_2line, 3)])   
    
    figure(h2)
    subplot(1, n_subj, ii)
    hold on
    plot([min_x, max_x], [min_x, max_x])
    plot(estimateStim1_2line, estimateStim2_2line, 'o',  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
    xlim([min_x, max_x])
    ylim([min_x, max_x])
    xlabel('Reported orientation 49 deg line')
    ylabel('Reported orientation 54 deg line')
end

%% Plot the combined subject
estimateStim1_all = cell2mat(estimateStim1_all);
estimateStim1_all = estimateStim1_all(:);
estimateStim2_all = cell2mat(estimateStim2_all);
estimateStim2_all = estimateStim2_all(:);
mean_diff = mean(estimateStim2_all) - mean(estimateStim1_all);

estimateStim1_2line_all = cell2mat(estimateStim1_2line_all);
estimateStim1_2line_all = estimateStim1_2line_all(:);
estimateStim2_2line_all = cell2mat(estimateStim2_2line_all);
estimateStim2_2line_all = estimateStim2_2line_all(:);
mean_diff_2line = mean(estimateStim2_2line_all) - mean(estimateStim1_2line_all);

figure
hold on
max_x = max([estimateStim1_all; estimateStim2_all;  estimateStim1_2line_all; estimateStim2_2line_all])+1;
min_x = min([estimateStim1_all; estimateStim2_all;  estimateStim1_2line_all; estimateStim2_2line_all])-1;
bin_edge = min_x:max_x;

subplot(2, 2, 1)
hold on
histogram(estimateStim1_all, bin_edge)
histogram(estimateStim2_all, bin_edge)
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
plot(estimateStim1_2line_all, estimateStim2_2line_all, 'o',  'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
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



%% Plot the results
% colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
%             'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
% colorIndex = NaN(length(colorName), 3);
% for ii = 1 : length(colorName)
%     colorIndex(ii, :) = rgb(colorName{ii});
% end
% legend_name = cell(1, length(subjectAll)+1);
% legend_name{1} = '';
% for ii = 1 : length(subjectAll)
%     legend_name{ii+1} = subjectAll{ii};
% end
% 
% figure
% minPlot_diffEst = min([meanDiffEst_same meanDiffEst_diff]) - 2;
% maxPlot_diffEst = max([meanDiffEst_same meanDiffEst_diff]) + 2;
% 
% hold on
% plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
% for ii = 1 : length(subjectAll)
%     plot(meanDiffEst_same(ii), meanDiffEst_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
% end
% axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
% axis square
% title('Cond 2')
% xlabel('Mean angle - same (deg)')
% ylabel('Mean angle - different (deg)')
% 
% legend(legend_name)

