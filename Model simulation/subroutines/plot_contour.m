function plot_contour(est_hist_cond1, est_hist_cond2, est_hist_cond3_1, est_hist_cond3_2, est_hist_cond3_3, n_contour, edge_bin)
figure;
hold on
minPlot = 10;
maxPlot = 81;

subplot(1, 5, 1)
hold on
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond1, n_contour);
axis square
grid on

subplot(1, 5, 2)
hold on
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond2, n_contour);
axis square
grid on

subplot(1, 5, 3)
hold on
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond3_1, n_contour);
axis square
grid on

subplot(1, 5, 4)
hold on
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond3_2, n_contour);
axis square
grid on   

subplot(1, 5, 5)
hold on
plot([minPlot maxPlot], [minPlot maxPlot], 'k')
contour(edge_bin, edge_bin, est_hist_cond3_3, n_contour);
axis square
grid on   