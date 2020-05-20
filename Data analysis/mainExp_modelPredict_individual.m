%% Plot the mean estimate model vs. data
colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

% Find the best model for each subject in condition 3
mad1 = abs(data_diff_est_cond13 - model_diff_est_cond13_1);
mad2 = abs(data_diff_est_cond13 - model_diff_est_cond13_2);
mad3 = abs(data_diff_est_cond13 - model_diff_est_cond13_3);
[~, ind_best] = min([mad1; mad2; mad3], [], 1);
model_diff_est_cond13_best = NaN(size(model_diff_est_cond13_1));
model_diff_est_cond13_best(ind_best == 1) = model_diff_est_cond13_1(ind_best == 1);
model_diff_est_cond13_best(ind_best == 2) = model_diff_est_cond13_2(ind_best == 2);
model_diff_est_cond13_best(ind_best == 3) = model_diff_est_cond13_3(ind_best == 3);


figure
hold on
minPlot = min([data_diff_est_cond12(:); data_diff_est_cond13(:); model_diff_est_cond12(:); model_diff_est_cond13_1(:); model_diff_est_cond13_2(:); model_diff_est_cond13_3(:)]) - 2;
maxPlot = max([data_diff_est_cond12(:); data_diff_est_cond13(:); model_diff_est_cond12(:); model_diff_est_cond13_1(:); model_diff_est_cond13_2(:); model_diff_est_cond13_3(:)]) + 2;
data_all = [data_diff_est_cond12' repmat(data_diff_est_cond13', 1, 4)];
model_all = [model_diff_est_cond12'  model_diff_est_cond13_1'  model_diff_est_cond13_2'  model_diff_est_cond13_3'  model_diff_est_cond13_best'];
fig_title = {'2-1';'3-1';'3-2'; '3-3';'3-best'};
for kk = 1 : 5
    subplot(1, 5, kk)
    [corr_index, pvalue] = corr(data_all(:, kk), model_all(:, kk), 'Type', 'Spearman');
    mad_index = mean(abs(data_all(:, kk) - model_all(:, kk)));
    hold on
    for ii = 1 : length(subjectID)
         plot(data_all(ii, kk), model_all(ii, kk), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize);
    end
    plot([minPlot maxPlot], [minPlot maxPlot])
    axis([minPlot maxPlot minPlot maxPlot])
    title([fig_title{kk} ': Corr=' num2str(roundn(corr_index, -2)) ', p=' num2str(roundn(pvalue, -2)) ', MAD=' num2str(roundn(mad_index, -2)) ])
    xlabel('Mean difference - data (deg)')
    ylabel('Mean difference - model (deg)')
    axis square
end