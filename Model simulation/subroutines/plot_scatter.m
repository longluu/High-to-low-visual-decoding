function plot_scatter(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, ...
                        theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, markerSize, markerColor, dotMeanSize, lineWidth, n_point_kdensity)

n_trials = size(estimate_cond1, 2);                  
figure
subplot(2, 5, 1)
hold on
%     minPlot = min([memory_sample_1 memory_sample_2 estimate_1' estimate_2'])-1;
%     maxPlot = max([memory_sample_1 memory_sample_2 estimate_1' estimate_2'])+1;
minPlot = 10;
maxPlot = 83;
percentCorrect = 100 * sum(estimate_cond1(2, :) > estimate_cond1(1, :)) / n_trials;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate_cond1(1, :), estimate_cond1(2, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate_cond1(1, :)), mean(estimate_cond1(2, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')    
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['Std sensory: ' num2str(roundn(std_sensory(1), 1)) ', ' num2str(roundn(std_sensory(2), 1)) ', '...
    'Std mem: ' num2str(roundn(std_memory_1line, 1)) ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(2, 5, 6)
hold on
diffEst = estimate_cond1(2, :) - estimate_cond1(1, :);
diffEst_mean = mean(estimate_cond1(2, :)) - mean(estimate_cond1(1, :));
[pDiff, x_pDiff] = ksdensity(diffEst, 'NumPoints', n_point_kdensity);
histogram(diffEst, 30, 'Normalization', 'pdf')
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean, 1)) ', Std: ' num2str(roundn(nanstd(diffEst), 1))]) 

subplot(2, 5, 2)
hold on
percentCorrect = 100 * sum(estimate_cond2(2, :) > estimate_cond2(1, :)) / n_trials;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate_cond2(1, :), estimate_cond2(2, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate_cond2(1, :)), mean(estimate_cond2(2, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['Std mem: '  num2str(roundn(std_memory_2line, 1))...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(2, 5, 7)
hold on
diffEst = estimate_cond2(2, :) - estimate_cond2(1, :);
diffEst_mean = mean(estimate_cond2(2, :)) - mean(estimate_cond2(1, :));
[pDiff, x_pDiff] = ksdensity(diffEst, 'NumPoints', n_point_kdensity);
histogram(diffEst, 30, 'Normalization', 'pdf')
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean, 1)) ', Std: ' num2str(roundn(nanstd(diffEst), 1))]) 

subplot(2, 5, 3)
hold on
percentCorrect = 100 * sum(estimate_cond3_1(2, :) > estimate_cond3_1(1, :)) / n_trials;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate_cond3_1(1, :), estimate_cond3_1(2, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate_cond3_1(1, :)), mean(estimate_cond3_1(2, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['Std mem: ' num2str(roundn(std_memory_redraw_1, 1)) ...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(2, 5, 8)
hold on
diffEst = estimate_cond3_1(2, :) - estimate_cond3_1(1, :);
diffEst_mean = mean(estimate_cond3_1(2, :)) - mean(estimate_cond3_1(1, :));
[pDiff, x_pDiff] = ksdensity(diffEst, 'NumPoints', n_point_kdensity);
histogram(diffEst, 30, 'Normalization', 'pdf')
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean, 1)) ', Std: ' num2str(roundn(nanstd(diffEst), 1))])   

subplot(2, 5, 4)
hold on
percentCorrect = 100 * sum(estimate_cond3_2(2, :) > estimate_cond3_2(1, :)) / n_trials;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate_cond3_2(1, :), estimate_cond3_2(2, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate_cond3_2(1, :)), mean(estimate_cond3_2(2, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['Std mem: ' num2str(roundn(std_memory_redraw_2, 1)) ...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(2, 5, 9)
hold on
diffEst = estimate_cond3_2(2, :) - estimate_cond3_2(1, :);
diffEst_mean = mean(estimate_cond3_2(2, :)) - mean(estimate_cond3_2(1, :));
[pDiff, x_pDiff] = ksdensity(diffEst, 'NumPoints', n_point_kdensity);
histogram(diffEst, 30, 'Normalization', 'pdf')
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean, 1)) ', Std: ' num2str(roundn(nanstd(diffEst), 1))])   

subplot(2, 5, 5)
hold on
percentCorrect = 100 * sum(estimate_cond3_3(2, :) > estimate_cond3_3(1, :)) / n_trials;
plot([minPlot maxPlot], [minPlot maxPlot])
h1 = scatter(estimate_cond3_3(1, :), estimate_cond3_3(2, :), markerSize, 'filled',...
    'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'none'); 
plot(theta_stim(1), theta_stim(2), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [0 1 0], 'MarkerEdgeColor', 'none')
plot(mean(estimate_cond3_3(1, :)), mean(estimate_cond3_3(2, :)), 'o', 'MarkerSize', dotMeanSize, 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', 'none')        
%     alpha(h1, alphaLevel)
xlabel('Reported orientation for stimulus 1')
ylabel('Reported orientation for stimulus 2')
axis equal
xlim([minPlot maxPlot])
ylim([minPlot maxPlot])
title(['Std mem: ' num2str(roundn(std_memory_redraw_3, 1)) ...
    ', PC: ' num2str(roundn(max(percentCorrect))) '%'])

subplot(2, 5, 10)
hold on
diffEst = estimate_cond3_3(2, :) - estimate_cond3_3(1, :);
diffEst_mean = mean(estimate_cond3_3(2, :)) - mean(estimate_cond3_3(1, :));
[pDiff, x_pDiff] = ksdensity(diffEst, 'NumPoints', n_point_kdensity);
histogram(diffEst, 30, 'Normalization', 'pdf')
plot(x_pDiff, pDiff, 'r', 'LineWidth', lineWidth)
xlabel('Orientation difference (deg)')
ylabel('Frequency of occurence')
xlim([min(-2, min(diffEst)-1) max(diffEst)+1])
axis square
title(['Mean :' num2str(roundn(diffEst_mean, 1)) ', Std: ' num2str(roundn(nanstd(diffEst), 1))])   