%%%%%%%%%%% Analyze data from the experiment with custom bootstrap function (separate) %%%%%%%%%%%%%%%
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
lineWidth = 2;
markerSize = 8;
rt_cond1_all = cell(length(subjectAll), 2);
rt_cond2_all = cell(length(subjectAll), 2);
rt_cond3_all = cell(length(subjectAll), 2);


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
        rt_cond1_all{ss, 1} = dataResponse(:, 5);

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
        rt_cond1_all{ss, 2} = dataResponse(:, 6);
         
    else
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
    rtLine1 = dataResponse(:, 5);
    rtLine2 = dataResponse(:, 6);

    % Save the result for each subject
    rt_cond2_all{ss, 1} = rtLine1;
    rt_cond2_all{ss, 2} = rtLine2;
    
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
    rtLine1 = dataResponse(:, 5);
    rtLine2 = dataResponse(:, 6);

    % Save the result for each subject
    rt_cond3_all{ss, 1} = rtLine1;
    rt_cond3_all{ss, 2} = rtLine2;
    
end

%% Plot the average RT of one report for each condition
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
n_bins = 50;
rt1_collapse = cell2mat(rt_cond1_all);
rt2_collapse = cell2mat(rt_cond2_all);
rt3_collapse = cell2mat(rt_cond3_all);
min_x = 0;
max_x = max([rt1_collapse(:); rt2_collapse(:); rt3_collapse(:)]);

subplot(1, 3, 1)
hist(rt1_collapse(:), n_bins)
xlim([min_x max_x])
xlabel('Reaction time (s)')
ylabel('Frequency')
title(['Cond 1, Median: ' num2str(median(rt1_collapse(:)), 2) ' ,Std: ' num2str(std(rt1_collapse(:)), 2)])

subplot(1, 3, 2)
hist(rt2_collapse(:), n_bins)
xlim([min_x max_x])
xlabel('Reaction time (s)')
ylabel('Frequency')
title(['Cond 2, Median: ' num2str(median(rt2_collapse(:)), 2) ' ,Std: ' num2str(std(rt2_collapse(:)), 2)])

subplot(1, 3, 3)
hist(rt3_collapse(:), n_bins)
xlim([min_x max_x])
xlabel('Reaction time (s)')
ylabel('Frequency')
title(['Cond 3, Median: ' num2str(median(rt3_collapse(:)), 2) ' ,Std: ' num2str(std(rt3_collapse(:)), 2)])

%% Compare RT of first and second reports
% Collapse RT for each subject and report
rt_cond2_ave = NaN(length(subjectAll), 2);
rt_cond3_ave = NaN(length(subjectAll), 2);
for ss = 1 : length(subjectAll)
    rt_cond2_ave(ss, 1) = median(rt_cond2_all{ss, 1});
    rt_cond2_ave(ss, 2) = median(rt_cond2_all{ss, 2});
    
    rt_cond3_ave(ss, 1) = median(rt_cond3_all{ss, 1});   
    rt_cond3_ave(ss, 2) = median(rt_cond3_all{ss, 2});    
end

% Plot
figure
hold on
minX = floor(min([rt1_collapse(:); rt2_collapse(:); rt3_collapse(:)]))-1;
maxX = ceil(max([rt1_collapse(:); rt2_collapse(:); rt3_collapse(:)]))+1;
binEdge = minX:0.5:maxX;

subplot(2, 2, 1)
hold on
median_report1 = median(rt2_collapse(:, 1));
median_report2 = median(rt2_collapse(:, 2));
histogram(rt2_collapse(:, 1), binEdge)
histogram(rt2_collapse(:, 2), binEdge)
bin_count1 = histcounts(rt2_collapse(:, 1), binEdge);
bin_count2 = histcounts(rt2_collapse(:, 2), binEdge);
maxY = max([bin_count1 bin_count2]);
plot([median_report1 median_report1], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
plot([median_report2 median_report2], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
xlim([minX maxX])
title(['Cond 2, ' ' diff (1-2): ' num2str(round(median_report1-median_report2, 1)) ' s'])
xlabel('Reaction time (s)')
ylabel('Count')
legend('Report 1', 'Report 2') 

subplot(2, 2, 2)
hold on
median_report1 = median(rt3_collapse(:, 1));
median_report2 = median(rt3_collapse(:, 2));
histogram(rt3_collapse(:, 1), binEdge)
histogram(rt3_collapse(:, 2), binEdge)
bin_count1 = histcounts(rt3_collapse(:, 1), binEdge);
bin_count2 = histcounts(rt3_collapse(:, 2), binEdge);
maxY = max([bin_count1 bin_count2]);
plot([median_report1 median_report1], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
plot([median_report2 median_report2], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
xlim([minX maxX])
title(['Cond 3, ' ' diff (1-2): ' num2str(round(median_report1-median_report2, 1)) ' s'])
xlabel('Reaction time (s)')
ylabel('Count')
legend('Report 1', 'Report 2') 

subplot(2, 2, 3)
hold on
minRT = floor(min([rt_cond2_ave(:); rt_cond3_ave(:)]))-1;
maxRT = ceil(max([rt_cond2_ave(:); rt_cond3_ave(:)]))+1;
plot(rt_cond2_ave(:, 1), rt_cond2_ave(:, 2), 'o', ...
    'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
plot([minRT maxRT], [minRT maxRT], '--k')
axis([minRT, maxRT, minRT, maxRT])
axis square
xlabel('RT - report 1 (s)')
ylabel('RT - report 2 (s)')
title('Cond 2')
p_RT_1vs2 = signrank(rt_cond2_ave(:, 1),rt_cond2_ave(:, 2), 'method', 'exact');
fprintf('p-value RT cond 2: %8.6f \n', p_RT_1vs2) 

subplot(2, 2, 4)
hold on
plot(rt_cond3_ave(:, 1), rt_cond3_ave(:, 2), 'o', ...
    'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
plot([minRT maxRT], [minRT maxRT], '--k')
axis([minRT, maxRT, minRT, maxRT])
axis square
xlabel('RT - report 1 (s)')
ylabel('RT - report 2 (s)')
title('Cond 3')
p_RT_1vs2 = signrank(rt_cond3_ave(:, 1),rt_cond3_ave(:, 2), 'method', 'exact');
fprintf('p-value RT cond 3: %8.6f \n', p_RT_1vs2) 