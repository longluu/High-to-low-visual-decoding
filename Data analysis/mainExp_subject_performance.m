%%%%%%%%%%% Analyze eye data %%%%%%%%%%%%%%%
subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
session = 1;
experimentNumber = 1;
nBootstrapBias1line = 10000;
nPlotBootstrapSample = 200;
markerSize = 4;
rmseAll = NaN(4, length(subjectAll));
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
    estimateStim1 = dataResponse(:, 3);
    orientationLine1 = dataResponse(:, 1);

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
    estimateStim2 = dataResponse(:, 4);
    orientationLine2 = dataResponse(:, 1);
    stimOrientation = [unique(orientationLine1) unique(orientationLine2)];
    
    % Compute error
    rmse1 = sqrt((sum((estimateStim1 - stimOrientation(1)).^2) + sum((estimateStim2 - stimOrientation(2)).^2)) ...
                ./ (length(estimateStim1) + length(estimateStim2)));


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
    estimateStim1 = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
    estimateStim2 = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
    estimateStim1(isnan(estimateStim1)) = [];
    estimateStim2(isnan(estimateStim2)) = [];
    
    % Compute error
    rmse2 = sqrt((sum((estimateStim1 - stimOrientation(1)).^2) + sum((estimateStim2 - stimOrientation(2)).^2)) ...
                ./ (length(estimateStim1) + length(estimateStim2)));

    
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
    estimateStim1 = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
    estimateStim2 = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
    estimateStim1(isnan(estimateStim1)) = [];
    estimateStim2(isnan(estimateStim2)) = [];
    
    % Compute error
    rmse3 = sqrt((sum((estimateStim1 - stimOrientation(1)).^2) + sum((estimateStim2 - stimOrientation(2)).^2)) ...
                ./ (length(estimateStim1) + length(estimateStim2)));
   
   % Total error
   rmseAll(:, ss) = [rmse1; rmse2; rmse3; (rmse1+rmse2+rmse3)/3];
end

%% Plot the MSE
figure
hold on
bar(1:length(subjectAll), rmseAll')
set(gca, 'XTick', 1:length(subjectAll), 'XTickLabel', subjectAll)
xlabel('Subject')
ylabel('Root mean square error')
legend('Exp 1', 'Exp 2', 'Exp 3', 'Average')