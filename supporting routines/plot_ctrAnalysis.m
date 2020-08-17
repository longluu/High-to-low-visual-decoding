function [median_p_2, sem_p_2, median_p_3, sem_p_3] = ...
    plot_ctrAnalysis(diffEst_Cond1, diffEst_Cond2_same, diffEst_Cond2_diff, diffEst_Cond3_same, diffEst_Cond3_diff)

diffEst_Cond2_same_offset = diffEst_Cond2_same - diffEst_Cond1;
diffEst_Cond2_diff_offset = diffEst_Cond2_diff - diffEst_Cond1;
p_crossTrial_Cond2 = (diffEst_Cond2_diff_offset - diffEst_Cond2_same_offset) ./ (diffEst_Cond2_diff_offset + diffEst_Cond2_same_offset);
p_crossTrial_Cond2(isnan(p_crossTrial_Cond2)) = NaN;

diffEst_Cond3_same_offset = diffEst_Cond3_same - diffEst_Cond1;
diffEst_Cond3_diff_offset = diffEst_Cond3_diff - diffEst_Cond1;
p_crossTrial_Cond3 = (diffEst_Cond3_diff_offset - diffEst_Cond3_same_offset) ./ (diffEst_Cond3_diff_offset + diffEst_Cond3_same_offset);
p_crossTrial_Cond3(isnan(p_crossTrial_Cond3)) = NaN;

colorName = {'Pink', 'Brown', 'Olive', 'Teal', 'Blue', 'Black', 'Red', 'Orange', 'Yellow',...
            'Lime', 'Cyan', 'DarkViolet', 'Magenta', 'Gray', 'RosyBrown', 'PaleGreen' };
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end

figure
minPlot_diffEst = min([diffEst_Cond2_same diffEst_Cond2_diff diffEst_Cond3_same diffEst_Cond3_diff]) - 2;
maxPlot_diffEst = max([diffEst_Cond2_same diffEst_Cond2_diff diffEst_Cond3_same diffEst_Cond3_diff]) + 2;
hold on

subplot(1, 2, 1)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(diffEst_Cond2_same(ii), diffEst_Cond2_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 2')
xlabel('mean angle - same (deg)')
ylabel('mean angle - different (deg)')

subplot(1, 2, 2)
hold on
plot([minPlot_diffEst maxPlot_diffEst], [minPlot_diffEst maxPlot_diffEst], 'k')
for ii = 1 : length(subjectAll)
    plot(diffEst_Cond3_same(ii), diffEst_Cond3_diff(ii), 'o', 'MarkerFaceColor', colorIndex(ii, :), 'MarkerEdgeColor', 'none', 'MarkerSize', markerSize)
end
axis([minPlot_diffEst, maxPlot_diffEst, minPlot_diffEst, maxPlot_diffEst])
axis square
title('Cond 3')
xlabel('mean angle - same (deg)')
ylabel('mean angle - different (deg)')

figure
subplot(1, 2, 1)
hold on
median_p_2 = nanmedian(p_crossTrial_Cond2);
sem_p_2 = nanstd(p_crossTrial_Cond2) / length(p_crossTrial_Cond2);
plot([0 length(p_crossTrial_Cond2)+1], [median_p_2 median_p_2], '--k') 
bar(1:length(p_crossTrial_Cond2), p_crossTrial_Cond2)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 2')

subplot(1, 2, 2)
hold on
median_p_3 = nanmedian(p_crossTrial_Cond3);
sem_p_3 = nanstd(p_crossTrial_Cond3) / length(p_crossTrial_Cond3);
plot([0 length(p_crossTrial_Cond3)+1], [median_p_3 median_p_3], '--k') 
bar(1:length(p_crossTrial_Cond3), p_crossTrial_Cond3)
xlabel('Subject')
ylabel('Fraction of repulsion explained by cross-trial adaptation')
title('Condition 3')