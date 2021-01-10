subjectAll = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;

%% Plot individual subjects
if plotIndividualSubject
    h_currentFig = gcf;
    markerSize = 4;
    for ss = 1 : length(subjectAll)
        subject = subjectAll{ss};
        h_currentFig = figure(h_currentFig.Number+1);
        h_currentFig.Position = [300 300 1200 700];
        
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
        eps_correlation = 0.001;
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
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none');          
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title('1 line')

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
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')        
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none')
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title('2 line')

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
        plot([minPlot maxPlot], [minPlot maxPlot])
        plot(stimOrientation(1), stimOrientation(2), 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none')
        plot(estimateStim1_collapse, estimateStim2_collapse, 'o', 'MarkerSize', markerSize, 'MarkerFaceColor', 0.1*[1 1 1], 'MarkerEdgeColor', 'none')        
        xlabel('Reported orientation for stimulus 1')
        ylabel('Reported orientation for stimulus 2')
        axis equal
        xlim([minPlot maxPlot])
        ylim([minPlot maxPlot])
        title('2 line - redraw')

        subplot(2, 3, 6)
        histogram(diffEst, 20, 'Normalization', 'probability')
        xlabel('Orientation difference (deg)')
        ylabel('Frequency of occurence')
        xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
        axis square
        title(['Mean: ' num2str(round(mean(diffEst), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))]) 
        
        saveas(h_currentFig, [subjectAll{ss} '.pdf'])
    end
end