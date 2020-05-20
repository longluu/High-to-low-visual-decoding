function plot_modelSampling_data(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, estimate1_allSubject, estimate2_allSubject, ...
                        theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, ...
                        markerSize, markerColor, dotMeanSize, lineWidth, n_point_kdensity, smooth_sampling, edge_bin, h_filter, n_contour, n_bins)
                    
n_trial_model = size(estimate_cond1, 2);    
n_trial_data = size(estimate1_allSubject, 2);    

% Get 2d histogram of model for contour plot
[est_hist_cond1, est_hist_cond2, est_hist_cond3_1, est_hist_cond3_2, est_hist_cond3_3] = ...
        make_2dhistogram(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, smooth_sampling, edge_bin, h_filter);    
    
figure
hold on
minPlot = 10;
maxPlot = 83;


%% Condition 1 - Model 1
subplot(3, 5, 1)
hold on
percentCorrect = 100 * sum(estimate_cond1(2, :) > estimate_cond1(1, :)) / n_trial_model;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond1, n_contour);
axis square
grid on
title(['Std sens: ' num2str(roundn(std_sensory(1), 1)) ', ' ...
    'Std mem: ' num2str(roundn(std_memory_1line, 1)) ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 6)
hold on
% Randomize estimates 1 and 2 in condition 1 according to a predetermined order
load ind_pair_cond1
estimateStim1_collapse = estimate1_allSubject(1, ind1);
estimateStim2_collapse = estimate2_allSubject(1, ind2);
percentCorrect = 100 * sum(estimateStim2_collapse > estimateStim1_collapse) / n_trial_data;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimateStim1_collapse, estimateStim2_collapse, markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimateStim1_collapse), mean(estimateStim2_collapse), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')    
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 11)
hold on
diffEst_data = estimateStim2_collapse - estimateStim1_collapse;
diffEst_mean_data = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
histogram(diffEst_data, n_bins, 'Normalization', 'pdf')
diffEst_model = estimate_cond1(2, :) - estimate_cond1(1, :);
[pDiff, x_pDiff] = ksdensity(diffEst_model, 'NumPoints', n_point_kdensity);
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst_data)-1) max(diffEst_data)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean_data, 1)) ', Std: ' num2str(roundn(nanstd(diffEst_data), 1))]) 

%% Condition 2 - Model 2a
subplot(3, 5, 2)
hold on
percentCorrect = 100 * sum(estimate_cond2(2, :) > estimate_cond2(1, :)) / n_trial_model;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond2, n_contour);
axis square
grid on
title(['Std mem: '  num2str(roundn(std_memory_2line, 1))...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 7)
hold on
percentCorrect = 100 * sum(estimate2_allSubject(2, :) > estimate1_allSubject(2, :)) / n_trial_data;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate1_allSubject(2, :), estimate2_allSubject(2, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate1_allSubject(2, :)), mean(estimate2_allSubject(2, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 12)
hold on
diffEst_data = estimate2_allSubject(2, :) - estimate1_allSubject(2, :);
diffEst_mean_data = mean(estimate2_allSubject(2, :)) - mean(estimate1_allSubject(2, :));
histogram(diffEst_data, n_bins, 'Normalization', 'pdf')
diffEst_model = estimate_cond2(2, :) - estimate_cond2(1, :);
[pDiff, x_pDiff] = ksdensity(diffEst_model, 'NumPoints', n_point_kdensity);
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst_data)-1) max(diffEst_data)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean_data, 1)) ', Std: ' num2str(roundn(nanstd(diffEst_data), 1))]) 

%% Condition 3 - Model 2b
subplot(3, 5, 3)
hold on
percentCorrect = 100 * sum(estimate_cond3_1(2, :) > estimate_cond3_1(1, :)) / n_trial_model;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond3_1, n_contour);
axis square
grid on
title(['Std mem: ' num2str(roundn(std_memory_redraw_1, 1)) ...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 8)
hold on
percentCorrect = 100 * sum(estimate2_allSubject(3, :) > estimate1_allSubject(3, :)) / n_trial_data;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate1_allSubject(3, :), estimate2_allSubject(3, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate1_allSubject(3, :)), mean(estimate2_allSubject(3, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 13)
hold on
diffEst_data = estimate2_allSubject(3, :) - estimate1_allSubject(3, :);
diffEst_mean_data = mean(estimate2_allSubject(1, :)) - mean(estimate1_allSubject(1, :));
histogram(diffEst_data, n_bins, 'Normalization', 'pdf')
diffEst_model = estimate_cond3_1(2, :) - estimate_cond3_1(1, :);
[pDiff, x_pDiff] = ksdensity(diffEst_model, 'NumPoints', n_point_kdensity);
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst_data)-1) max(diffEst_data)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean_data, 1)) ', Std: ' num2str(roundn(nanstd(diffEst_data), 1))])   

%% Condition 3 - Model 2c
subplot(3, 5, 4)
hold on
percentCorrect = 100 * sum(estimate_cond3_2(2, :) > estimate_cond3_2(1, :)) / n_trial_model;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond3_2, n_contour);
axis square
grid on
title(['Std mem: ' num2str(roundn(std_memory_redraw_2, 1)) ...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 9)
hold on
percentCorrect = 100 * sum(estimate2_allSubject(3, :) > estimate1_allSubject(3, :)) / n_trial_data;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate1_allSubject(3, :), estimate2_allSubject(3, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate1_allSubject(3, :)), mean(estimate2_allSubject(3, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 14)
hold on
diffEst_data = estimate2_allSubject(3, :) - estimate1_allSubject(3, :);
diffEst_mean_data = mean(estimate2_allSubject(1, :)) - mean(estimate1_allSubject(1, :));
histogram(diffEst_data, n_bins, 'Normalization', 'pdf')
diffEst_model = estimate_cond3_2(2, :) - estimate_cond3_2(1, :);
[pDiff, x_pDiff] = ksdensity(diffEst_model, 'NumPoints', n_point_kdensity);
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst_data)-1) max(diffEst_data)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean_data, 1)) ', Std: ' num2str(roundn(nanstd(diffEst_data), 1))])   

%% Condition 3 - Model 2d
subplot(3, 5, 5)
hold on
percentCorrect = 100 * sum(estimate_cond3_3(2, :) > estimate_cond3_3(1, :)) / n_trial_model;
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond3_3, n_contour);
axis square
grid on
title(['Std mem: ' num2str(roundn(std_memory_redraw_3, 1)) ...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 10)
hold on
percentCorrect = 100 * sum(estimate2_allSubject(3, :) > estimate1_allSubject(3, :)) / n_trial_data;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate1_allSubject(3, :), estimate2_allSubject(3, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate1_allSubject(3, :)), mean(estimate2_allSubject(3, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(3, 5, 15)
hold on
diffEst_data = estimate2_allSubject(3, :) - estimate1_allSubject(3, :);
diffEst_mean_data = mean(estimate2_allSubject(1, :)) - mean(estimate1_allSubject(1, :));
histogram(diffEst_data, n_bins, 'Normalization', 'pdf')
diffEst_model = estimate_cond3_3(2, :) - estimate_cond3_3(1, :);
[pDiff, x_pDiff] = ksdensity(diffEst_model, 'NumPoints', n_point_kdensity);
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst_data)-1) max(diffEst_data)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean_data, 1)) ', Std: ' num2str(roundn(nanstd(diffEst_data), 1))])   