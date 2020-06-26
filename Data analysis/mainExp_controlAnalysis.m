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
meanDiffEst_Cond1_same = NaN(1, length(subjectAll));
meanDiffEst_Cond1_diff = NaN(1, length(subjectAll));
meanDiffEst_Cond1 = NaN(1, length(subjectAll));
meanDiffEst_Cond2_same = NaN(1, length(subjectAll));
meanDiffEst_Cond2_diff = NaN(1, length(subjectAll));
meanDiffEst_Cond2 = NaN(1, length(subjectAll));
meanDiffEst_Cond3_same = NaN(1, length(subjectAll));
meanDiffEst_Cond3_diff = NaN(1, length(subjectAll));
meanDiffEst_Cond3 = NaN(1, length(subjectAll));

%% Compute the mean angle difference for same and difference groups
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
        
        % Bootstrap the mean difference, correlation and percent correct
        experiment_condition = 1;
        [confInterval_all, sampleStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial,...
                        experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse);

        if excludeIncorrectTrial == 2            
            meanDiffEst_Cond1(1, ss) = sampleStat_all(1);
        else
            meanDiffEst_Cond1(1, ss) = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
        end
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

        % Split subjects' estimates based on previous trial (same or different) 
        trialOrder = params.trialOrder;
        trial_current = trialOrder(2:end);
        trial_last = trialOrder(1:end-1);
        ind_same = (trial_current > 50 & trial_last > 50) | (trial_current <= 50 & trial_last <= 50);
        ind_diff = (trial_current > 50 & trial_last <= 50) | (trial_current <= 50 & trial_last > 50);
        trial_same = trial_current(ind_same);
        trial_diff = trial_current(ind_diff);

        estimateLine1_same = estimateLine1;
        estimateLine1_same(trial_diff) = NaN;
        estimateLine2_same = estimateLine2;
        estimateLine2_same(trial_diff) = NaN;
        estimateStim1_same = [estimateLine1_same(indStim1Est1); estimateLine2_same(indStim1Est2)];
        estimateStim2_same = [estimateLine2_same(indStim2Est2);estimateLine1_same(indStim2Est1)];
        estimateStim1_same = estimateStim1_same(~isnan(estimateStim1_same));
        estimateStim2_same = estimateStim2_same(~isnan(estimateStim2_same));
        estimateStim1_same(estimateStim1_same<0) = estimateStim1_same(estimateStim1_same<0)+180;
        estimateStim2_same(estimateStim2_same<0) = estimateStim2_same(estimateStim2_same<0)+180;  

        estimateLine1_diff = estimateLine1;
        estimateLine1_diff(trial_same) = NaN;
        estimateLine2_diff = estimateLine2;
        estimateLine2_diff(trial_same) = NaN;
        estimateStim1_diff = [estimateLine1_diff(indStim1Est1); estimateLine2_diff(indStim1Est2)];
        estimateStim2_diff = [estimateLine2_diff(indStim2Est2);estimateLine1_diff(indStim2Est1)];
        estimateStim1_diff = estimateStim1_diff(~isnan(estimateStim1_diff));
        estimateStim2_diff = estimateStim2_diff(~isnan(estimateStim2_diff));
        estimateStim1_diff(estimateStim1_diff<0) = estimateStim1_diff(estimateStim1_diff<0)+180;
        estimateStim2_diff(estimateStim2_diff<0) = estimateStim2_diff(estimateStim2_diff<0)+180;    

        meanDiffEst_Cond1_same(ss) = nanmean(estimateStim2_same) - nanmean(estimateStim1_same);
        meanDiffEst_Cond1_diff(ss) = nanmean(estimateStim2_diff) - nanmean(estimateStim1_diff);
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

    % Split subjects' estimates based on previous trial (same or different) 
    trialOrder = params.trialOrder;
    trial_current = trialOrder(2:end);
    trial_last = trialOrder(1:end-1);
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);
    
    estimateLine1_same = estimateLine1;
    estimateLine1_same(trial_diff) = NaN;
    estimateLine2_same = estimateLine2;
    estimateLine2_same(trial_diff) = NaN;
    estimateStim1_same = [estimateLine1_same(indStim1Est1); estimateLine2_same(indStim1Est2)];
    estimateStim2_same = [estimateLine2_same(indStim2Est2);estimateLine1_same(indStim2Est1)];
    estimateStim1_same = estimateStim1_same(~isnan(estimateStim1_same));
    estimateStim2_same = estimateStim2_same(~isnan(estimateStim2_same));
    estimateStim1_same(estimateStim1_same<0) = estimateStim1_same(estimateStim1_same<0)+180;
    estimateStim2_same(estimateStim2_same<0) = estimateStim2_same(estimateStim2_same<0)+180;  
    
    estimateLine1_diff = estimateLine1;
    estimateLine1_diff(trial_same) = NaN;
    estimateLine2_diff = estimateLine2;
    estimateLine2_diff(trial_same) = NaN;
    estimateStim1_diff = [estimateLine1_diff(indStim1Est1); estimateLine2_diff(indStim1Est2)];
    estimateStim2_diff = [estimateLine2_diff(indStim2Est2);estimateLine1_diff(indStim2Est1)];
    estimateStim1_diff = estimateStim1_diff(~isnan(estimateStim1_diff));
    estimateStim2_diff = estimateStim2_diff(~isnan(estimateStim2_diff));
    estimateStim1_diff(estimateStim1_diff<0) = estimateStim1_diff(estimateStim1_diff<0)+180;
    estimateStim2_diff(estimateStim2_diff<0) = estimateStim2_diff(estimateStim2_diff<0)+180;  

    if excludeIncorrectTrial == 0
        meanDiffEst_Cond2_same(ss) = nanmean(estimateStim2_same) - nanmean(estimateStim1_same);
        meanDiffEst_Cond2_diff(ss) = nanmean(estimateStim2_diff) - nanmean(estimateStim1_diff);
        meanDiffEst_Cond2(ss) = nanmean([estimateStim2_same; estimateStim2_diff]) ...
                                        - nanmean([estimateStim1_same; estimateStim1_diff]);
    elseif excludeIncorrectTrial == 1
        diffEst_same = estimateStim2_same - estimateStim1_same;
        diffEst_diff = estimateStim2_diff - estimateStim1_diff;
        meanDiffEst_Cond2_same(ss) = nanmean(estimateStim2_same(diffEst_same>0)) - nanmean(estimateStim1_same(diffEst_same>0));
        meanDiffEst_Cond2_diff(ss) = nanmean(estimateStim2_diff(diffEst_diff>0)) - nanmean(estimateStim1_diff(diffEst_diff>0)); 
        meanDiffEst_Cond2(ss) = nanmean([estimateStim2_same(diffEst_same>0); estimateStim2_diff(diffEst_diff>0)]) ...
                                        - nanmean([estimateStim1_same(diffEst_same>0); estimateStim1_diff(diffEst_diff>0)]);        
    else
        indFlip = (estimateStim2_same - estimateStim1_same) < 0; % flip the incorrect trials
        estStim1_swap = estimateStim1_same(indFlip);
        estStim2_swap = estimateStim2_same(indFlip);
        estimateStim1_same(indFlip) = estStim2_swap;
        estimateStim2_same(indFlip) = estStim1_swap;  
        meanDiffEst_Cond2_same(ss) = nanmean(estimateStim2_same) - nanmean(estimateStim1_same);

        indFlip = (estimateStim2_diff - estimateStim1_diff) < 0; % flip the incorrect trials
        estStim1_swap = estimateStim1_diff(indFlip);
        estStim2_swap = estimateStim2_diff(indFlip);
        estimateStim1_diff(indFlip) = estStim2_swap;
        estimateStim2_diff(indFlip) = estStim1_swap;  
        meanDiffEst_Cond2_diff(ss) = nanmean(estimateStim2_diff) - nanmean(estimateStim1_diff);    
        meanDiffEst_Cond2(ss) = nanmean([estimateStim2_same; estimateStim2_diff]) ...
                                        - nanmean([estimateStim1_same; estimateStim1_diff]);        
    end
    
    
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
    
    % Split subjects' estimates based on previous trial (same or different) 
    trialOrder = params.trialOrder;
    trial_current = trialOrder(2:end);
    trial_last = trialOrder(1:end-1);
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);
    
    estimateLine1_same = estimateLine1;
    estimateLine1_same(trial_diff) = NaN;
    estimateLine2_same = estimateLine2;
    estimateLine2_same(trial_diff) = NaN;
    estimateStim1_same = [estimateLine1_same(indStim1Est1); estimateLine2_same(indStim1Est2)];
    estimateStim2_same = [estimateLine2_same(indStim2Est2);estimateLine1_same(indStim2Est1)];
    estimateStim1_same = estimateStim1_same(~isnan(estimateStim1_same));
    estimateStim2_same = estimateStim2_same(~isnan(estimateStim2_same));
    estimateStim1_same(estimateStim1_same<0) = estimateStim1_same(estimateStim1_same<0)+180;
    estimateStim2_same(estimateStim2_same<0) = estimateStim2_same(estimateStim2_same<0)+180;  
    
    estimateLine1_diff = estimateLine1;
    estimateLine1_diff(trial_same) = NaN;
    estimateLine2_diff = estimateLine2;
    estimateLine2_diff(trial_same) = NaN;
    estimateStim1_diff = [estimateLine1_diff(indStim1Est1); estimateLine2_diff(indStim1Est2)];
    estimateStim2_diff = [estimateLine2_diff(indStim2Est2);estimateLine1_diff(indStim2Est1)];
    estimateStim1_diff = estimateStim1_diff(~isnan(estimateStim1_diff));
    estimateStim2_diff = estimateStim2_diff(~isnan(estimateStim2_diff));
    estimateStim1_diff(estimateStim1_diff<0) = estimateStim1_diff(estimateStim1_diff<0)+180;
    estimateStim2_diff(estimateStim2_diff<0) = estimateStim2_diff(estimateStim2_diff<0)+180;  
    
    if excludeIncorrectTrial == 0
        meanDiffEst_Cond3_same(ss) = nanmean(estimateStim2_same) - nanmean(estimateStim1_same);
        meanDiffEst_Cond3_diff(ss) = nanmean(estimateStim2_diff) - nanmean(estimateStim1_diff);
        meanDiffEst_Cond3(ss) = nanmean([estimateStim2_same; estimateStim2_diff]) ...
                                        - nanmean([estimateStim1_same; estimateStim1_diff]);        
    elseif excludeIncorrectTrial == 1
        diffEst_same = estimateStim2_same - estimateStim1_same;
        diffEst_diff = estimateStim2_diff - estimateStim1_diff;
        meanDiffEst_Cond3_same(ss) = nanmean(estimateStim2_same(diffEst_same>0)) - nanmean(estimateStim1_same(diffEst_same>0));
        meanDiffEst_Cond3_diff(ss) = nanmean(estimateStim2_diff(diffEst_diff>0)) - nanmean(estimateStim1_diff(diffEst_diff>0)); 
        meanDiffEst_Cond3(ss) = nanmean([estimateStim2_same(diffEst_same>0); estimateStim2_diff(diffEst_diff>0)]) ...
                                        - nanmean([estimateStim1_same(diffEst_same>0); estimateStim1_diff(diffEst_diff>0)]);                
    else
        indFlip = (estimateStim2_same - estimateStim1_same) < 0; % flip the incorrect trials
        estStim1_swap = estimateStim1_same(indFlip);
        estStim2_swap = estimateStim2_same(indFlip);
        estimateStim1_same(indFlip) = estStim2_swap;
        estimateStim2_same(indFlip) = estStim1_swap;  
        meanDiffEst_Cond3_same(ss) = nanmean(estimateStim2_same) - nanmean(estimateStim1_same);

        indFlip = (estimateStim2_diff - estimateStim1_diff) < 0; % flip the incorrect trials
        estStim1_swap = estimateStim1_diff(indFlip);
        estStim2_swap = estimateStim2_diff(indFlip);
        estimateStim1_diff(indFlip) = estStim2_swap;
        estimateStim2_diff(indFlip) = estStim1_swap;  
        meanDiffEst_Cond3_diff(ss) = nanmean(estimateStim2_diff) - nanmean(estimateStim1_diff);       
        meanDiffEst_Cond3(ss) = nanmean([estimateStim2_same; estimateStim2_diff]) ...
                                        - nanmean([estimateStim1_same; estimateStim1_diff]);        
    end
    
end

%% Compute the fraction repulsion explained by cross-trial adaptation
meanDiffEst_Cond2_same_offset = meanDiffEst_Cond2_same - meanDiffEst_Cond1;
meanDiffEst_Cond2_diff_offset = meanDiffEst_Cond2_diff - meanDiffEst_Cond1;
meanDiffEst_Cond2_offset = meanDiffEst_Cond2 - meanDiffEst_Cond1;
p_crossTrial_Cond2 = (meanDiffEst_Cond2_diff_offset - meanDiffEst_Cond2_same_offset) ./ meanDiffEst_Cond2_offset;
p_crossTrial_Cond2(meanDiffEst_Cond2_offset<0) = [];

meanDiffEst_Cond3_same_offset = meanDiffEst_Cond3_same - meanDiffEst_Cond1;
meanDiffEst_Cond3_diff_offset = meanDiffEst_Cond3_diff - meanDiffEst_Cond1;
meanDiffEst_Cond3_offset = meanDiffEst_Cond3 - meanDiffEst_Cond1;
p_crossTrial_Cond3 = (meanDiffEst_Cond3_diff_offset - meanDiffEst_Cond3_same_offset) ./ meanDiffEst_Cond3_offset;
p_crossTrial_Cond3(meanDiffEst_Cond3_offset<0) = [];

%% Plot the results
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([meanDiffEst_Cond2_same meanDiffEst_Cond2_diff meanDiffEst_Cond3_same meanDiffEst_Cond3_diff]) - 2;
maxPlot_diffEst = max([meanDiffEst_Cond2_same meanDiffEst_Cond2_diff meanDiffEst_Cond3_same meanDiffEst_Cond3_diff]) + 2;
hold on

subplot(1, 3, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(meanDiffEst_Cond1_same(ii), meanDiffEst_Cond1_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 1')
xlabel('Mean angle - same (deg)')
ylabel('Mean angle - different (deg)')

subplot(1, 3, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(meanDiffEst_Cond2_same(ii), meanDiffEst_Cond2_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('Mean angle - same (deg)')
ylabel('Mean angle - different (deg)')

subplot(1, 3, 3)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(meanDiffEst_Cond3_same(ii), meanDiffEst_Cond3_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('Mean angle - same (deg)')
ylabel('Mean angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = median(p_crossTrial_Cond2);
sem_p_2 = std(p_crossTrial_Cond2) / length(p_crossTrial_Cond2);
plot([0 length(p_crossTrial_Cond2)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2), p_crossTrial_Cond2)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = median(p_crossTrial_Cond3);
sem_p_3 = std(p_crossTrial_Cond3) / length(p_crossTrial_Cond3);
plot([0 length(p_crossTrial_Cond3)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3), p_crossTrial_Cond3)
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
ylabel('Fraction of repulsion explained by cross-trial adaptation')
