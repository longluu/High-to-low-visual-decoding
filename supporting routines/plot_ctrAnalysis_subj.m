function [p_crossTrial_Cond2, p_crossTrial_Cond3] = plot_ctrAnalysis_subj(meanDiffEst_Cond1, meanDiffEst_Cond2_same, ...
                                    meanDiffEst_Cond2_diff, meanDiffEst_Cond3_same, meanDiffEst_Cond3_diff, fig_title, markerSize, nBootstrap)

meanDiffEst_Cond2_same_offset = meanDiffEst_Cond2_same - meanDiffEst_Cond1;
meanDiffEst_Cond2_diff_offset = meanDiffEst_Cond2_diff - meanDiffEst_Cond1;
p_crossTrial_Cond2 = (meanDiffEst_Cond2_diff_offset - meanDiffEst_Cond2_same_offset) ./ (meanDiffEst_Cond2_diff_offset + meanDiffEst_Cond2_same_offset);
p_crossTrial_Cond2(isnan(p_crossTrial_Cond2)) = NaN;

meanDiffEst_Cond3_same_offset = meanDiffEst_Cond3_same - meanDiffEst_Cond1;
meanDiffEst_Cond3_diff_offset = meanDiffEst_Cond3_diff - meanDiffEst_Cond1;
p_crossTrial_Cond3 = (meanDiffEst_Cond3_diff_offset - meanDiffEst_Cond3_same_offset) ./ (meanDiffEst_Cond3_diff_offset + meanDiffEst_Cond3_same_offset);
p_crossTrial_Cond3(isnan(p_crossTrial_Cond3)) = NaN;

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

subplot(1, 2, 1)
hold on
n_subject = length(meanDiffEst_Cond2_same);
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : n_subject
    plot(meanDiffEst_Cond2_same(ii), meanDiffEst_Cond2_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('mean angle - same (deg)')
ylabel('mean angle - different (deg)')

subplot(1, 2, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : n_subject
    plot(meanDiffEst_Cond3_same(ii), meanDiffEst_Cond3_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('mean angle - same (deg)')
ylabel('mean angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
nanmedian_p_2 = nanmedian(p_crossTrial_Cond2);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_Cond2);
sem_p_2 = nanstd(bootstat);
plot([0 length(p_crossTrial_Cond2)+1], [nanmedian_p_2 nanmedian_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2), p_crossTrial_Cond2)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
nanmedian_p_3 = nanmedian(p_crossTrial_Cond3);
bootstat = bootstrp(nBootstrap,@nanmedian, p_crossTrial_Cond3);
sem_p_3 = nanstd(bootstat);
plot([0 length(p_crossTrial_Cond3)+1], [nanmedian_p_3 nanmedian_p_3], '--k') 
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
errorBarGraph([nanmedian_p_2; nanmedian_p_3], [nanmedian_p_2-sem_p_2; nanmedian_p_3-sem_p_3], [nanmedian_p_2+sem_p_2; nanmedian_p_3+sem_p_3], colorIndex)
box off
title(fig_title)
ylabel('Fraction of repulsion explained by cross-trial adaptation')