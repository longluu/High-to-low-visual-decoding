%%%%%%%%%%% Compare repulsion between same and different previous trial %%%%%%%%%%%%%%%
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
session = 1;
experimentNumber = 1;
markerSize = 8;
lineWidth = 0.5;
nBootstrap = 10000;
useSplit1line = 1;
alpha = 0.05;
excludeIncorrectTrial = 2; % 0: all trials
                           % 1: correct trials
                           % 2: all trials but flip incorrect trials across diagonal line
                           
medianDiffEst_Cond1_same = NaN(1, length(subjectAll));
medianDiffEst_Cond1_diff = NaN(1, length(subjectAll));
medianDiffEst_Cond1 = NaN(1, length(subjectAll));

medianDiffEst_Cond2_same_1back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_diff_1back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_same_2back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_diff_2back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_same_3back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_diff_3back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_same_4back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_diff_4back = NaN(1, length(subjectAll));
medianDiffEst_Cond2_same_2back_cum = NaN(1, length(subjectAll));
medianDiffEst_Cond2_diff_2back_cum = NaN(1, length(subjectAll));
medianDiffEst_Cond2 = NaN(1, length(subjectAll));

medianDiffEst_Cond3_same_1back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_diff_1back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_same_2back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_diff_2back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_same_3back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_diff_3back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_same_4back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_diff_4back = NaN(1, length(subjectAll));
medianDiffEst_Cond3_same_2back_cum = NaN(1, length(subjectAll));
medianDiffEst_Cond3_diff_2back_cum = NaN(1, length(subjectAll));
medianDiffEst_Cond3 = NaN(1, length(subjectAll));

%% Compute the median angle difference for same and difference groups
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

    % Bootstrap the median difference, correlation and percent correct
    experiment_condition = 1;
    [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                    experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse);

    if excludeIncorrectTrial == 2            
        medianDiffEst_Cond1(1, ss) = sampleStat_all(1);
    else
        medianDiffEst_Cond1(1, ss) = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
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
    [medianDiffEst_Cond2(ss), medianDiffEst_Cond2_same_1back(ss), medianDiffEst_Cond2_diff_1back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
    % Split subjects' estimates based on previous 2-back trial (same or different) 
    n_back = 2;
    [~, medianDiffEst_Cond2_same_2back(ss), medianDiffEst_Cond2_diff_2back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2); 
    
    % Split subjects' estimates based on previous 3-back trial (same or different) 
    n_back = 3;
    [~, medianDiffEst_Cond2_same_3back(ss), medianDiffEst_Cond2_diff_3back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   

    % Split subjects' estimates based on previous 4-back trial (same or different) 
    n_back = 4;
    [~, medianDiffEst_Cond2_same_4back(ss), medianDiffEst_Cond2_diff_4back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   

    % Split subjects' estimates based on previous 2-back trial cumulative (same or different) 
    n_back = 12;
    [~, medianDiffEst_Cond2_same_2back_cum(ss), medianDiffEst_Cond2_diff_2back_cum(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
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
    [medianDiffEst_Cond3(ss), medianDiffEst_Cond3_same_1back(ss), medianDiffEst_Cond3_diff_1back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);   
    
    % Split subjects' estimates based on previous 2-back trial (same or different) 
    n_back = 2;
    [~, medianDiffEst_Cond3_same_2back(ss), medianDiffEst_Cond3_diff_2back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2); 
    
    % Split subjects' estimates based on previous 3-back trial (same or different) 
    n_back = 3;
    [~, medianDiffEst_Cond3_same_3back(ss), medianDiffEst_Cond3_diff_3back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
    % Split subjects' estimates based on previous 4-back trial (same or different) 
    n_back = 4;
    [~, medianDiffEst_Cond3_same_4back(ss), medianDiffEst_Cond3_diff_4back(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);  
    
    % Split subjects' estimates based on previous 2-back trial cumulative (same or different) 
    n_back = 12;
    [~, medianDiffEst_Cond3_same_2back_cum(ss), medianDiffEst_Cond3_diff_2back_cum(ss)] = split_diffEst(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2);
    
end

%% Plot the results for 2-back cum
medianDiffEst_Cond2_same_offset = medianDiffEst_Cond2_same_2back_cum - medianDiffEst_Cond1;
medianDiffEst_Cond2_diff_offset = medianDiffEst_Cond2_diff_2back_cum - medianDiffEst_Cond1;
medianDiffEst_Cond2_offset = medianDiffEst_Cond2 - medianDiffEst_Cond1;
p_crossTrial_Cond2_2back_cum = (medianDiffEst_Cond2_diff_offset - medianDiffEst_Cond2_same_offset) ./ medianDiffEst_Cond2_offset;
p_crossTrial_Cond2_2back_cum(medianDiffEst_Cond2_offset<0) = [];

medianDiffEst_Cond3_same_offset = medianDiffEst_Cond3_same_2back_cum - medianDiffEst_Cond1;
medianDiffEst_Cond3_diff_offset = medianDiffEst_Cond3_diff_2back_cum - medianDiffEst_Cond1;
medianDiffEst_Cond3_offset = medianDiffEst_Cond3 - medianDiffEst_Cond1;
p_crossTrial_Cond3_2back_cum = (medianDiffEst_Cond3_diff_offset - medianDiffEst_Cond3_same_offset) ./ medianDiffEst_Cond3_offset;
p_crossTrial_Cond3_2back_cum(medianDiffEst_Cond3_offset<0) = [];

colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([medianDiffEst_Cond2_same_2back_cum medianDiffEst_Cond2_diff_2back_cum medianDiffEst_Cond3_same_2back_cum medianDiffEst_Cond3_diff_2back_cum]) - 2;
maxPlot_diffEst = max([medianDiffEst_Cond2_same_2back_cum medianDiffEst_Cond2_diff_2back_cum medianDiffEst_Cond3_same_2back_cum medianDiffEst_Cond3_diff_2back_cum]) + 2;
hold on

subplot(1, 3, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond1_same(ii), medianDiffEst_Cond1_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 1')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond2_same_2back_cum(ii), medianDiffEst_Cond2_diff_2back_cum(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 3)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond3_same_2back_cum(ii), medianDiffEst_Cond3_diff_2back_cum(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = median(p_crossTrial_Cond2_2back_cum);
sem_p_2 = std(p_crossTrial_Cond2_2back_cum) / length(p_crossTrial_Cond2_2back_cum);
plot([0 length(p_crossTrial_Cond2_2back_cum)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2_2back_cum), p_crossTrial_Cond2_2back_cum)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = median(p_crossTrial_Cond3_2back_cum);
sem_p_3 = std(p_crossTrial_Cond3_2back_cum) / length(p_crossTrial_Cond3_2back_cum);
plot([0 length(p_crossTrial_Cond3_2back_cum)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3_2back_cum), p_crossTrial_Cond3_2back_cum)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 3')

figure
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hold on
plot([0 3], [1 1], '--k')
errorBarGraph([median_p_2; median_p_3], [median_p_2-sem_p_2; median_p_3-sem_p_3], [median_p_2+sem_p_2; median_p_3+sem_p_3], colorIndex)
box off
title('2-back cum')
ylabel('Fraction of repulsion explained by cross-trial adaptation')

%% Plot the results for 1-back
medianDiffEst_Cond2_same_offset = medianDiffEst_Cond2_same_1back - medianDiffEst_Cond1;
medianDiffEst_Cond2_diff_offset = medianDiffEst_Cond2_diff_1back - medianDiffEst_Cond1;
medianDiffEst_Cond2_offset = medianDiffEst_Cond2 - medianDiffEst_Cond1;
p_crossTrial_Cond2_1back = (medianDiffEst_Cond2_diff_offset - medianDiffEst_Cond2_same_offset) ./ medianDiffEst_Cond2_offset;
p_crossTrial_Cond2_1back(medianDiffEst_Cond2_offset<0) = [];

medianDiffEst_Cond3_same_offset = medianDiffEst_Cond3_same_1back - medianDiffEst_Cond1;
medianDiffEst_Cond3_diff_offset = medianDiffEst_Cond3_diff_1back - medianDiffEst_Cond1;
medianDiffEst_Cond3_offset = medianDiffEst_Cond3 - medianDiffEst_Cond1;
p_crossTrial_Cond3_1back = (medianDiffEst_Cond3_diff_offset - medianDiffEst_Cond3_same_offset) ./ medianDiffEst_Cond3_offset;
p_crossTrial_Cond3_1back(medianDiffEst_Cond3_offset<0) = [];

colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([medianDiffEst_Cond2_same_1back medianDiffEst_Cond2_diff_1back medianDiffEst_Cond3_same_1back medianDiffEst_Cond3_diff_1back]) - 2;
maxPlot_diffEst = max([medianDiffEst_Cond2_same_1back medianDiffEst_Cond2_diff_1back medianDiffEst_Cond3_same_1back medianDiffEst_Cond3_diff_1back]) + 2;
hold on

subplot(1, 3, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond1_same(ii), medianDiffEst_Cond1_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 1')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond2_same_1back(ii), medianDiffEst_Cond2_diff_1back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 3)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond3_same_1back(ii), medianDiffEst_Cond3_diff_1back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = median(p_crossTrial_Cond2_1back);
sem_p_2 = std(p_crossTrial_Cond2_1back) / length(p_crossTrial_Cond2_1back);
plot([0 length(p_crossTrial_Cond2_1back)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2_1back), p_crossTrial_Cond2_1back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = median(p_crossTrial_Cond3_1back);
sem_p_3 = std(p_crossTrial_Cond3_1back) / length(p_crossTrial_Cond3_1back);
plot([0 length(p_crossTrial_Cond3_1back)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3_1back), p_crossTrial_Cond3_1back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 3')

figure
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hold on
plot([0 3], [1 1], '--k')
errorBarGraph([median_p_2; median_p_3], [median_p_2-sem_p_2; median_p_3-sem_p_3], [median_p_2+sem_p_2; median_p_3+sem_p_3], colorIndex)
box off
title('1-back')
ylabel('Fraction of repulsion explained by cross-trial adaptation')

%% Plot the results for 2-back
medianDiffEst_Cond2_same_offset = medianDiffEst_Cond2_same_2back - medianDiffEst_Cond1;
medianDiffEst_Cond2_diff_offset = medianDiffEst_Cond2_diff_2back - medianDiffEst_Cond1;
medianDiffEst_Cond2_offset = medianDiffEst_Cond2 - medianDiffEst_Cond1;
p_crossTrial_Cond2_2back = (medianDiffEst_Cond2_diff_offset - medianDiffEst_Cond2_same_offset) ./ medianDiffEst_Cond2_offset;
p_crossTrial_Cond2_2back(medianDiffEst_Cond2_offset<0) = [];

medianDiffEst_Cond3_same_offset = medianDiffEst_Cond3_same_2back - medianDiffEst_Cond1;
medianDiffEst_Cond3_diff_offset = medianDiffEst_Cond3_diff_2back - medianDiffEst_Cond1;
medianDiffEst_Cond3_offset = medianDiffEst_Cond3 - medianDiffEst_Cond1;
p_crossTrial_Cond3_2back = (medianDiffEst_Cond3_diff_offset - medianDiffEst_Cond3_same_offset) ./ medianDiffEst_Cond3_offset;
p_crossTrial_Cond3_2back(medianDiffEst_Cond3_offset<0) = [];

colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([medianDiffEst_Cond2_same_2back medianDiffEst_Cond2_diff_2back medianDiffEst_Cond3_same_2back medianDiffEst_Cond3_diff_2back]) - 2;
maxPlot_diffEst = max([medianDiffEst_Cond2_same_2back medianDiffEst_Cond2_diff_2back medianDiffEst_Cond3_same_2back medianDiffEst_Cond3_diff_2back]) + 2;
hold on

subplot(1, 3, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond1_same(ii), medianDiffEst_Cond1_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 1')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond2_same_2back(ii), medianDiffEst_Cond2_diff_2back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 3)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond3_same_2back(ii), medianDiffEst_Cond3_diff_2back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = median(p_crossTrial_Cond2_2back);
sem_p_2 = std(p_crossTrial_Cond2_2back) / length(p_crossTrial_Cond2_2back);
plot([0 length(p_crossTrial_Cond2_2back)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2_2back), p_crossTrial_Cond2_2back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = median(p_crossTrial_Cond3_2back);
sem_p_3 = std(p_crossTrial_Cond3_2back) / length(p_crossTrial_Cond3_2back);
plot([0 length(p_crossTrial_Cond3_2back)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3_2back), p_crossTrial_Cond3_2back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 3')

figure
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hold on
plot([0 3], [1 1], '--k')
errorBarGraph([median_p_2; median_p_3], [median_p_2-sem_p_2; median_p_3-sem_p_3], [median_p_2+sem_p_2; median_p_3+sem_p_3], colorIndex)
box off
title('2-back')
ylabel('Fraction of repulsion explained by cross-trial adaptation')

%% Plot the results for 3-back
medianDiffEst_Cond2_same_offset = medianDiffEst_Cond2_same_3back - medianDiffEst_Cond1;
medianDiffEst_Cond2_diff_offset = medianDiffEst_Cond2_diff_3back - medianDiffEst_Cond1;
medianDiffEst_Cond2_offset = medianDiffEst_Cond2 - medianDiffEst_Cond1;
p_crossTrial_Cond2_3back = (medianDiffEst_Cond2_diff_offset - medianDiffEst_Cond2_same_offset) ./ medianDiffEst_Cond2_offset;
p_crossTrial_Cond2_3back(medianDiffEst_Cond2_offset<0) = [];

medianDiffEst_Cond3_same_offset = medianDiffEst_Cond3_same_3back - medianDiffEst_Cond1;
medianDiffEst_Cond3_diff_offset = medianDiffEst_Cond3_diff_3back - medianDiffEst_Cond1;
medianDiffEst_Cond3_offset = medianDiffEst_Cond3 - medianDiffEst_Cond1;
p_crossTrial_Cond3_3back = (medianDiffEst_Cond3_diff_offset - medianDiffEst_Cond3_same_offset) ./ medianDiffEst_Cond3_offset;
p_crossTrial_Cond3_3back(medianDiffEst_Cond3_offset<0) = [];

colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([medianDiffEst_Cond2_same_3back medianDiffEst_Cond2_diff_3back medianDiffEst_Cond3_same_3back medianDiffEst_Cond3_diff_3back]) - 2;
maxPlot_diffEst = max([medianDiffEst_Cond2_same_3back medianDiffEst_Cond2_diff_3back medianDiffEst_Cond3_same_3back medianDiffEst_Cond3_diff_3back]) + 2;
hold on

subplot(1, 3, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond1_same(ii), medianDiffEst_Cond1_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 1')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond2_same_3back(ii), medianDiffEst_Cond2_diff_3back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 3)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond3_same_3back(ii), medianDiffEst_Cond3_diff_3back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = median(p_crossTrial_Cond2_3back);
sem_p_2 = std(p_crossTrial_Cond2_3back) / length(p_crossTrial_Cond2_3back);
plot([0 length(p_crossTrial_Cond2_3back)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2_3back), p_crossTrial_Cond2_3back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = median(p_crossTrial_Cond3_3back);
sem_p_3 = std(p_crossTrial_Cond3_3back) / length(p_crossTrial_Cond3_3back);
plot([0 length(p_crossTrial_Cond3_3back)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3_3back), p_crossTrial_Cond3_3back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 3')

figure
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hold on
plot([0 3], [1 1], '--k')
errorBarGraph([median_p_2; median_p_3], [median_p_2-sem_p_2; median_p_3-sem_p_3], [median_p_2+sem_p_2; median_p_3+sem_p_3], colorIndex)
box off
title('3-back')
ylabel('Fraction of repulsion explained by cross-trial adaptation')

%% Plot the results for 4-back
medianDiffEst_Cond2_same_offset = medianDiffEst_Cond2_same_4back - medianDiffEst_Cond1;
medianDiffEst_Cond2_diff_offset = medianDiffEst_Cond2_diff_4back - medianDiffEst_Cond1;
medianDiffEst_Cond2_offset = medianDiffEst_Cond2 - medianDiffEst_Cond1;
p_crossTrial_Cond2_4back = (medianDiffEst_Cond2_diff_offset - medianDiffEst_Cond2_same_offset) ./ medianDiffEst_Cond2_offset;
p_crossTrial_Cond2_4back(medianDiffEst_Cond2_offset<0) = [];

medianDiffEst_Cond3_same_offset = medianDiffEst_Cond3_same_4back - medianDiffEst_Cond1;
medianDiffEst_Cond3_diff_offset = medianDiffEst_Cond3_diff_4back - medianDiffEst_Cond1;
medianDiffEst_Cond3_offset = medianDiffEst_Cond3 - medianDiffEst_Cond1;
p_crossTrial_Cond3_4back = (medianDiffEst_Cond3_diff_offset - medianDiffEst_Cond3_same_offset) ./ medianDiffEst_Cond3_offset;
p_crossTrial_Cond3_4back(medianDiffEst_Cond3_offset<0) = [];

colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([medianDiffEst_Cond2_same_4back medianDiffEst_Cond2_diff_4back medianDiffEst_Cond3_same_4back medianDiffEst_Cond3_diff_4back]) - 2;
maxPlot_diffEst = max([medianDiffEst_Cond2_same_4back medianDiffEst_Cond2_diff_4back medianDiffEst_Cond3_same_4back medianDiffEst_Cond3_diff_4back]) + 2;
hold on

subplot(1, 3, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond1_same(ii), medianDiffEst_Cond1_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 1')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond2_same_4back(ii), medianDiffEst_Cond2_diff_4back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

subplot(1, 3, 3)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(medianDiffEst_Cond3_same_4back(ii), medianDiffEst_Cond3_diff_4back(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('median angle - same (deg)')
ylabel('median angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = median(p_crossTrial_Cond2_4back);
sem_p_2 = std(p_crossTrial_Cond2_4back) / length(p_crossTrial_Cond2_4back);
plot([0 length(p_crossTrial_Cond2_4back)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2_4back), p_crossTrial_Cond2_4back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = median(p_crossTrial_Cond3_4back);
sem_p_3 = std(p_crossTrial_Cond3_4back) / length(p_crossTrial_Cond3_4back);
plot([0 length(p_crossTrial_Cond3_4back)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3_4back), p_crossTrial_Cond3_4back)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 3')

figure
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hold on
plot([0 3], [1 1], '--k')
errorBarGraph([median_p_2; median_p_3], [median_p_2-sem_p_2; median_p_3-sem_p_3], [median_p_2+sem_p_2; median_p_3+sem_p_3], colorIndex)
box off
title('4-back')
ylabel('Fraction of repulsion explained by cross-trial adaptation')

%% Plot the results for all back trials
p_crossTrial_1back = [p_crossTrial_Cond2_1back p_crossTrial_Cond3_1back];
p_crossTrial_2back = [p_crossTrial_Cond2_2back p_crossTrial_Cond3_2back];
p_crossTrial_3back = [p_crossTrial_Cond2_3back p_crossTrial_Cond3_3back];
p_crossTrial_4back = [p_crossTrial_Cond2_4back p_crossTrial_Cond3_4back];
p_crossTrial_2back_cum = [p_crossTrial_Cond2_2back_cum p_crossTrial_Cond3_2back_cum];

median_p_1back = median(p_crossTrial_1back);
sem_p_1back = std(p_crossTrial_1back) / length(p_crossTrial_1back);
median_p_2back = median(p_crossTrial_2back);
sem_p_2back = std(p_crossTrial_2back) / length(p_crossTrial_2back);
median_p_3back = median(p_crossTrial_3back);
sem_p_3back = std(p_crossTrial_3back) / length(p_crossTrial_3back);
median_p_4back = median(p_crossTrial_4back);
sem_p_4back = std(p_crossTrial_4back) / length(p_crossTrial_4back);
median_p_2back_cum = median(p_crossTrial_2back_cum);
sem_p_2back_cum = std(p_crossTrial_2back_cum) / length(p_crossTrial_2back_cum);

figure
hold on
plot([0 5], [1 1], '--k')
errorBarGraph([median_p_1back; median_p_2back; median_p_3back; median_p_4back; median_p_2back_cum], ...
    [median_p_1back-sem_p_1back; median_p_2back-sem_p_3back; median_p_3back-sem_p_2back; median_p_4back-sem_p_4back; median_p_2back_cum-sem_p_2back_cum], ...
    [median_p_1back+sem_p_1back; median_p_2back+sem_p_3back; median_p_3back+sem_p_2back; median_p_4back+sem_p_4back; median_p_2back_cum+sem_p_2back_cum], colorIndex)
box off
title('Collapse condition 2 and 3')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
xlabel('n-back')