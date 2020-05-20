%%%%%%%%%%% Analyze data from the experiment with custom bootstrap function (separate) %%%%%%%%%%%%%%%
% Old: 'an', 'bl', 'ccp', 'km', 'ln', 'rl', 'sr', 'sar', 'sy', 'cr', 'cm', 'mb', 'lw', 'bg'
% New: 'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
markerSize = 8;
lineWidth = 0.5;
alpha = 0.05;
reaction_time_cond1 = cell(length(subjectAll), 2);
reaction_time_cond2 = cell(length(subjectAll), 2);
reaction_time_cond3 = cell(length(subjectAll), 2);


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
        reaction_time_cond1{ss, 1} = dataResponse(:, 5);

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
        reaction_time_cond1{ss, 2} = dataResponse(:, 6);
         
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
    orientationLine1 = dataResponse(:, 1);
    orientationLine2 = dataResponse(:, 2);
    rtLine1 = dataResponse(:, 5);
    rtLine2 = dataResponse(:, 6);

    % Extract the stimulus orientation
    stimOrientation = unique(orientationLine1);

    % Collapse subjects' estimates across presentation order
    indStim1Est1 = orientationLine1 == stimOrientation(1);
    indStim1Est2 = orientationLine2 == stimOrientation(1);
    indStim2Est1 = orientationLine1 == stimOrientation(2);
    indStim2Est2 = orientationLine2 == stimOrientation(2);
    reaction_time_cond2{ss, 1} = [rtLine1(indStim1Est1); rtLine2(indStim1Est2)];
    reaction_time_cond2{ss, 2} = [rtLine2(indStim2Est2); rtLine1(indStim2Est1)];
    
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
    rtLine1 = dataResponse(:, 5);
    rtLine2 = dataResponse(:, 6);

    % Extract the stimulus orientation
    stimOrientation = unique(orientationLine1);

    % Collapse subjects' estimates across presentation order
    indStim1Est1 = orientationLine1 == stimOrientation(1);
    indStim1Est2 = orientationLine2 == stimOrientation(1);
    indStim2Est1 = orientationLine1 == stimOrientation(2);
    indStim2Est2 = orientationLine2 == stimOrientation(2);
    reaction_time_cond3{ss, 1} = [rtLine1(indStim1Est1); rtLine2(indStim1Est2)];
    reaction_time_cond3{ss, 2} = [rtLine2(indStim2Est2); rtLine1(indStim2Est1)];
    
end

%% Plot the results
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
n_bins = 50;
rt1_collapse = cell2mat(reaction_time_cond1);
rt2_collapse = cell2mat(reaction_time_cond2);
rt3_collapse = cell2mat(reaction_time_cond3);
min_x = 0;
max_x = max([rt1_collapse(:); rt2_collapse(:); rt3_collapse(:)]);

subplot(1, 3, 1)
hist(rt1_collapse(:), n_bins)
xlim([min_x max_x])
xlabel('Reaction time (s)')
ylabel('Frequency')
title(['Mean: ' num2str(mean(rt1_collapse(:)), 2) ' ,Std: ' num2str(std(rt1_collapse(:)), 2)])

subplot(1, 3, 2)
hist(rt2_collapse(:), n_bins)
xlim([min_x max_x])
xlabel('Reaction time (s)')
ylabel('Frequency')
title(['Mean: ' num2str(mean(rt2_collapse(:)), 2) ' ,Std: ' num2str(std(rt2_collapse(:)), 2)])

subplot(1, 3, 3)
hist(rt3_collapse(:), n_bins)
xlim([min_x max_x])
xlabel('Reaction time (s)')
ylabel('Frequency')
title(['Mean: ' num2str(mean(rt3_collapse(:)), 2) ' ,Std: ' num2str(std(rt3_collapse(:)), 2)])

