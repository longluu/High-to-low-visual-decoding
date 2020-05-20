function plot_density
    figure;
    hold on
    axis_display = 0:10:70;
    

    subplot(1, 4, 1)
    hold on
    imagesc(est_hist_cond1)
    % xlim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % ylim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % set(gca, 'XTick', axis_display, 'XTickLabel', num2cell(axis_display+min(edge_bin)), ...
    %         'YTick', axis_display, 'YTickLabel', num2cell(axis_display+min(edge_bin)))
    colormap gray
    axis xy
    axis square
    axis off

    subplot(1, 4, 2)
    hold on
    imagesc(est_hist_cond2)
    % xlim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % ylim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % set(gca, 'XTick', axis_display, 'XTickLabel', num2cell(axis_display+min(edge_bin)), ...
    %         'YTick', axis_display, 'YTickLabel', num2cell(axis_display+min(edge_bin)))
    colormap gray
    axis xy
    axis square
    axis off

    subplot(1, 4, 3)
    hold on
    imagesc(est_hist_cond3_1)
    % xlim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % ylim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % set(gca, 'XTick', axis_display, 'XTickLabel', num2cell(axis_display+min(edge_bin)), ...
    %         'YTick', axis_display, 'YTickLabel', num2cell(axis_display+min(edge_bin)))
    colormap gray
    axis xy
    axis square
    axis off

    subplot(1, 4, 4)
    hold on
    imagesc(est_hist_cond3_2)
    % xlim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % ylim([min(edge_bin) max(edge_bin)]-min(edge_bin))
    % set(gca, 'XTick', axis_display, 'XTickLabel', num2cell(axis_display+min(edge_bin)), ...
    %         'YTick', axis_display, 'YTickLabel', num2cell(axis_display+min(edge_bin)))
    colormap gray
    axis xy
    axis square
    axis off