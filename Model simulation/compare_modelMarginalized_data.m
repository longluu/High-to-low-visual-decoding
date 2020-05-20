%%%%%%%%%%%%%%%%%%%%% Compare model (marginalized implementation) and data %%%%%%%%%%%%%%%%%%%%
%%%% Note that for cond 3, we don't have a marginalized implementation
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
lineWidth = 3;
alphaLevel = 1;
fontSize = 10;

% Param model
theta_stim = [49 54];
std_sensory = 2.2 * [1 1];
std_memory_1line = 4.5;
std_memory_2line = 6.2;
std_memory_redraw_1 = 6;
std_motor = 0;
n_discreteStep = 100;
n_contour = 15;
n_bins = 50;

%% Extract the data
[estimate1_allSubject, estimate2_allSubject] = getData(subjectID, session, experimentNumber, useSplit1line);

%% Model
tic
% Get the predictive distribution p(thHat|th)
[pthetahGtheta_1, pthetahGtheta_2, pthetahGtheta_3, theta] = retroBayes_marginalized(theta_stim, std_sensory, std_memory_1line, ...
    std_memory_2line, std_memory_redraw_1, n_discreteStep);
if std_motor ~= 0
    theta_centered = theta - (max(theta) + min(theta))/2;
    pMotor = normpdf(theta_centered, 0, std_motor)' * normpdf(theta_centered, 0, std_motor);
    pthetahGtheta_1 = conv2(pthetahGtheta_1, pMotor, 'same');
    pthetahGtheta_2 = conv2(pthetahGtheta_2, pMotor, 'same');
end

% Project the distribution along the diagonal axis
n_step_project = 2*n_discreteStep - 1;
axis_project_middle = n_discreteStep;
pProject_1 = NaN(1, n_step_project);
pProject_2 = NaN(1, n_step_project);
pProject_3 = NaN(1, n_step_project);
for ii = 1 : n_discreteStep-1
    % Condition 1
    pProject_1(axis_project_middle + ii) = sum(diag(pthetahGtheta_1, ii));
    pProject_1(axis_project_middle - ii) = sum(diag(pthetahGtheta_1, -ii));

    % Condition 2
    pProject_2(axis_project_middle + ii) = sum(diag(pthetahGtheta_2, ii));
    pProject_2(axis_project_middle - ii) = sum(diag(pthetahGtheta_2, -ii));

    % Condition 3
    pProject_3(axis_project_middle + ii) = sum(diag(pthetahGtheta_3, ii));
    pProject_3(axis_project_middle - ii) = sum(diag(pthetahGtheta_3, -ii));    
end
pProject_1(axis_project_middle) = sum(diag(pthetahGtheta_1));
pProject_2(axis_project_middle) = sum(diag(pthetahGtheta_2));
pProject_3(axis_project_middle) = sum(diag(pthetahGtheta_3));
pProject_1 = fliplr(pProject_1) / sum(pProject_1);
pProject_2 = fliplr(pProject_2) / sum(pProject_2);
pProject_3 = fliplr(pProject_3) / sum(pProject_3);

% Create the projection axis
[theta1, theta2] = meshgrid(theta, theta);
theta_diff = diag(rot90(theta2 - theta1));
step_axis_project = (theta_diff(2) - theta_diff(1)) / 2;
theta_augment = theta_diff(1:end-1) + step_axis_project;
theta_augment = [theta_augment; 0];
axis_project = [theta_diff theta_augment]';
axis_project = axis_project(:);
axis_project(end) = [];
toc

%% Plot the joint distribution
figure
hold on
set(gcf, 'Position', [100 50 1200 1100])

% Exp 1
% Randomize estimates 1 and 2 in condition 1 according to a predetermined order
load ind_pair_cond1_combined
estimateStim1_collapse = estimate1_allSubject(1, ind1);
estimateStim2_collapse = estimate2_allSubject(1, ind2);
diffEst = estimateStim2_collapse - estimateStim1_collapse;
percentCorrect = 100 * (sum(diffEst > 0) + sum(diffEst == 0)/2) / sum(~isnan(diffEst));
subplot(3, 3, 1)
hold on
set(gca, 'FontSize', fontSize)
minPlot = min([estimate1_allSubject(:); estimate2_allSubject(:)])-1;
maxPlot = max([estimate1_allSubject(:); estimate2_allSubject(:)])+1;
% minPlot = 30;
% maxPlot = 70;
plot([minPlot maxPlot], [minPlot maxPlot])
contour(theta, theta, pthetahGtheta_1, n_contour);
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', stimColor, 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', meanColor, 'MarkerEdgeColor', 'none')
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])

subplot(3, 3, 4)
hold on
set(gca, 'FontSize', fontSize)
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimateStim1_collapse, estimateStim2_collapse, markerSize, 'filled',...
    'MarkerFaceColor', dataColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', stimColor, 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', meanColor, 'MarkerEdgeColor', 'none')
alpha(h1, alphaLevel);
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['1 line: ' num2str(round(percentCorrect, 1)) '%'])

subplot(3, 3, 7)
hold on
set(gca, 'FontSize', fontSize)
% minHist = min(estimate2_allSubject(:) - estimate1_allSubject(:)) - 1;
% maxHist = max(estimate2_allSubject(:) - estimate1_allSubject(:)) + 1;
minHist = -40;
maxHist = 50;
histogram(diffEst, n_bins, 'Normalization', 'pdf')
plot(axis_project, pProject_1 / step_axis_project, 'r', 'LineWidth', lineWidth); 
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([minHist maxHist])
axis square
title(['Mean :' num2str(round(mean(estimateStim2_collapse) - mean(estimateStim1_collapse), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))])    

% Exp 2
estimateStim1_collapse = estimate1_allSubject(2, :);
estimateStim2_collapse = estimate2_allSubject(2, :);
diffEst = estimateStim2_collapse - estimateStim1_collapse;
percentCorrect = 100 * (sum(diffEst > 0) + sum(diffEst == 0)/2) / sum(~isnan(diffEst));
subplot(3, 3, 2)
hold on
set(gca, 'FontSize', fontSize)
plot([minPlot maxPlot], [minPlot maxPlot])
contour(theta, theta, pthetahGtheta_2, n_contour);
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', stimColor, 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', meanColor, 'MarkerEdgeColor', 'none')
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])

subplot(3, 3, 5)
hold on
set(gca, 'FontSize', fontSize)
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimateStim1_collapse, estimateStim2_collapse, markerSize, 'filled',...
    'MarkerFaceColor', dataColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', stimColor, 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', meanColor, 'MarkerEdgeColor', 'none')
alpha(h1, alphaLevel);
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['2 line: ' num2str(round(percentCorrect, 1)) '%'])

subplot(3, 3, 8)
hold on
set(gca, 'FontSize', fontSize)
histogram(diffEst, n_bins, 'Normalization', 'pdf')
plot(axis_project, pProject_2 / step_axis_project, 'r', 'LineWidth', lineWidth); 
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([minHist maxHist])
axis square
title(['Mean :' num2str(round(mean(estimateStim2_collapse) - mean(estimateStim1_collapse), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))])    

% Exp 3
estimateStim1_collapse = estimate1_allSubject(3, :);
estimateStim2_collapse = estimate2_allSubject(3, :);
diffEst = estimateStim2_collapse - estimateStim1_collapse;
percentCorrect = 100 * (sum(diffEst > 0) + sum(diffEst == 0)/2) / sum(~isnan(diffEst));
subplot(3, 3, 3)
hold on
set(gca, 'FontSize', fontSize)
plot([minPlot maxPlot], [minPlot maxPlot])
contour(theta, theta, pthetahGtheta_3, n_contour);
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', stimColor, 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', meanColor, 'MarkerEdgeColor', 'none')
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])

subplot(3, 3, 6)
hold on
set(gca, 'FontSize', fontSize)
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimateStim1_collapse, estimateStim2_collapse, markerSize, 'filled',...
    'MarkerFaceColor', dataColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', stimColor, 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', meanColor, 'MarkerEdgeColor', 'none')
alpha(h1, alphaLevel);
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['2 line (redraw): ' num2str(round(percentCorrect, 1)) '%'])

subplot(3, 3, 9)
hold on
set(gca, 'FontSize', fontSize)
histogram(diffEst, n_bins, 'Normalization', 'pdf')
plot(axis_project, pProject_3 / step_axis_project, 'r', 'LineWidth', lineWidth); 
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([minHist maxHist])
axis square
title(['Mean :' num2str(round(mean(estimateStim2_collapse) - mean(estimateStim1_collapse), 1)) ', Std: ' num2str(round(nanstd(diffEst), 1))])   

tightfig

%% Plot the marginal distribution
delta_theta = theta(2) - theta(1);
p_est1 = sum(pthetahGtheta_1, 1);
p_est1 = p_est1 / (sum(p_est1)*delta_theta);
p_est2 = sum(pthetahGtheta_1, 2);
p_est2 = p_est2 / (sum(p_est2)*delta_theta);

estimateStim1_collapse = estimate1_allSubject(1, ind1);
estimateStim2_collapse = estimate2_allSubject(1, ind2);
min_x = min([estimateStim1_collapse estimateStim2_collapse]) - 1;
max_x = max([estimateStim1_collapse estimateStim2_collapse]) + 1;
binEdge = min_x:1:max_x;

figure
hold on
histogram(estimateStim1_collapse, binEdge, 'Normalization', 'pdf')
histogram(estimateStim2_collapse, binEdge, 'Normalization', 'pdf')
plot(theta, p_est1, 'LineWidth', lineWidth)
plot(theta, p_est2, 'LineWidth', lineWidth)
bin_count1 = histcounts(estimate1_allSubject(1, :), binEdge);
bin_count2 = histcounts(estimate2_allSubject(1, :), binEdge);
maxY = max([bin_count1 bin_count2]);
xlim([min_x max_x])
title('1 line condition')
xlabel('Reported orientation (deg)')
ylabel('Count')
legend('49^o line', '54^o line') 
