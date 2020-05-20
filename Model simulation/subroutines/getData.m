function [estimate1, estimate2] = getData(subjectID, session, experimentNumber, useSplit1line)
if strcmp(subjectID, 'combined')
    subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
    estimate1 = NaN(3, length(subjectAll) * 50);
    estimate2 = NaN(3, length(subjectAll) * 50);

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

        % Accumulate the combined subject
        estimate1(1, 50*(ss-1)+1 : 50*ss) = estimateStim1_collapse;
        estimate2(1, 50*(ss-1)+1 : 50*ss) = estimateStim2_collapse;

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

        % Accumulate the combined subject
        estimate1(2, 50*(ss-1)+1 : 50*ss) = estimateStim1_collapse;
        estimate2(2, 50*(ss-1)+1 : 50*ss) = estimateStim2_collapse;

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

        % Accumulate the combined subject
        estimate1(3, 50*(ss-1)+1 : 50*ss) = estimateStim1_collapse;
        estimate2(3, 50*(ss-1)+1 : 50*ss) = estimateStim2_collapse;    
    end
else
    estimate1 = NaN(3, 50);
    estimate2 = NaN(3, 50);

    subject = subjectID;

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

    % Accumulate the combined subject
    estimate1(1, :) = estimateStim1_collapse;
    estimate2(1, :) = estimateStim2_collapse;

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

    % Accumulate the combined subject
    estimate1(2, :) = estimateStim1_collapse;
    estimate2(2, :) = estimateStim2_collapse;

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

    % Accumulate the combined subject
    estimate1(3, :) = estimateStim1_collapse;
    estimate2(3, :) = estimateStim2_collapse;    
end