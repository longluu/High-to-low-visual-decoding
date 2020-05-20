%%%%%%%%%%%%%%%%%%%%% Compare model and data %%%%%%%%%%%%%%%%%%%%
% Params data
subjectID = 'combined';
useSplit1line = 1;
session = 1;
experimentNumber = 1;
markerSize = 3;
dataColor = 0.5 * [1 1 1];
stimColor = [1 0 1];
meanColor = [0 1 0];
dotMeanSize = 6;
lineWidth = 2;
alphaLevel = 1;
fontSize = 10;

%% Extract the data
[estimate1_allSubject, estimate2_allSubject] = getData(subjectID, session, experimentNumber, useSplit1line);
mean_est1 = mean(estimate1_allSubject, 2);
mean_est2 = mean(estimate2_allSubject, 2);

%% Plot the absolute judgements
figure
hold on
estimate_all = [estimate1_allSubject(:); estimate2_allSubject(:)];
minX = floor(min(estimate_all))-1;
maxX = ceil(max(estimate_all))+1;
binEdge = minX:1:maxX;

subplot(1, 3, 1)
hold on
histogram(estimate1_allSubject(1, :), binEdge)
histogram(estimate2_allSubject(1, :), binEdge)
bin_count1 = histcounts(estimate1_allSubject(1, :), binEdge);
bin_count2 = histcounts(estimate2_allSubject(1, :), binEdge);
maxY = max([bin_count1 bin_count2]);
plot([mean_est1(1) mean_est1(1)], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
plot([mean_est2(1) mean_est2(1)], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
xlim([minX maxX])
title('1 line condition')
xlabel('Reported orientation (deg)')
ylabel('Count')
legend('49^o line', '54^o line') 

subplot(1, 3, 2)
hold on
histogram(estimate1_allSubject(2, :), binEdge)
histogram(estimate2_allSubject(2, :), binEdge)
bin_count1 = histcounts(estimate1_allSubject(2, :), binEdge);
bin_count2 = histcounts(estimate2_allSubject(2, :), binEdge);
maxY = max([bin_count1 bin_count2]);
plot([mean_est1(2) mean_est1(2)], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
plot([mean_est2(2) mean_est2(2)], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
xlim([minX maxX])
title('2 line condition')
xlabel('Reported orientation (deg)')
ylabel('Count')
legend('49^o line', '54^o line') 

subplot(1, 3, 3)
hold on
histogram(estimate1_allSubject(3, :), binEdge)
histogram(estimate2_allSubject(3, :), binEdge)
bin_count1 = histcounts(estimate1_allSubject(3, :), binEdge);
bin_count2 = histcounts(estimate2_allSubject(3, :), binEdge);
maxY = max([bin_count1 bin_count2]);
plot([mean_est1(3) mean_est1(3)], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
plot([mean_est2(3) mean_est2(3)], [0 maxY], '--', 'Color', [0 0 0], 'LineWidth', lineWidth)
xlim([minX maxX])
title('2 line - redraw condition')
xlabel('Reported orientation (deg)')
ylabel('Count')
legend('49^o line', '54^o line') 