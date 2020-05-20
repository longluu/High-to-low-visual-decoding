%% Parameters
% General 
flag_sampling = 1;
flag_marginal = 0;

% Model
theta_stim = [49 54];
std_sensory = [2.7 2.7];
std_memory_1line = 4;
std_memory_2line = 5.8;
std_memory_redraw_1 = 7.5;
std_memory_redraw_2 = 6;
std_memory_redraw_3 = 6;
fraction_redecode = 1;
n_trials = 200000;
n_discreteStep = 100;

% Plot
markerSize = 6;
markerColor = 0.5 * [1 1 1];
dotMeanSize = 6;
lineWidth = 3;
n_point_kdensity = 100;
smooth_sampling = 1;
edge_bin = 10:0.5:81;
h_filter = fspecial('gaussian', [3 3], 1.5);
n_contour = 15;

%% Marginalized
if flag_marginal 
    [pthetahGtheta_1, pthetahGtheta_2, pthetahGtheta_3, theta] = retroBayes_marginalized(theta_stim, std_sensory, std_memory_1line, ...
        std_memory_2line, std_memory_redraw_1, fraction_redecode, n_discreteStep);    

    %% Plot
    figure;
    minPlot = 10;
    maxPlot = 81;
    hold on

    subplot(1, 3, 1)
    hold on
    imagesc(pthetahGtheta_1)
%     xlim([minPlot maxPlot])
%     ylim([minPlot maxPlot])
    colormap gray
    axis xy
    axis square
    axis off

    subplot(1, 3, 2)
    hold on
    imagesc(pthetahGtheta_2)
%     xlim([minPlot maxPlot])
%     ylim([minPlot maxPlot])
    colormap gray
    axis xy
    axis square
    axis off

    subplot(1, 3, 3)
    hold on
%     xlim([minPlot maxPlot])
%     ylim([minPlot maxPlot])
    imagesc(pthetahGtheta_3)
    colormap gray
    axis xy
    axis square
    axis off
    
    % Plot contour
    figure;
    hold on
    n_contour = 15;

    subplot(1, 3, 1)
    hold on
    plot([minPlot maxPlot], [minPlot maxPlot], 'k')
    contour(theta, theta, pthetahGtheta_1, n_contour);
    xlim([minPlot maxPlot])
    ylim([minPlot maxPlot])
    axis square
    grid on

    subplot(1, 3, 2)
    hold on
    plot([minPlot maxPlot], [minPlot maxPlot], 'k')
    contour(theta, theta, pthetahGtheta_2, n_contour);
    xlim([minPlot maxPlot])
    ylim([minPlot maxPlot])
    axis square
    grid on

    subplot(1, 3, 3)
    hold on
    plot([minPlot maxPlot], [minPlot maxPlot], 'k')
    contour(theta, theta, pthetahGtheta_3, n_contour);
    xlim([minPlot maxPlot])
    ylim([minPlot maxPlot])
    axis square
    grid on
end

%% Sampling
if flag_sampling
    % 2-d retrospective Bayes
    [estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3] = retroBayes_sampling_2d(theta_stim, std_sensory, std_memory_1line, ...
        std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, fraction_redecode, n_trials, n_discreteStep);

    % Plot scatter
    tic
    plot_scatter(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, ...
                        theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, markerSize, markerColor, dotMeanSize, lineWidth, n_point_kdensity)
    toc
    
    % Plot contour
    [est_hist_cond1, est_hist_cond2, est_hist_cond3_1, est_hist_cond3_2, est_hist_cond3_3] = ...
            make_2dhistogram(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, smooth_sampling, edge_bin, h_filter);    
    plot_contour(est_hist_cond1, est_hist_cond2, est_hist_cond3_1, est_hist_cond3_2, est_hist_cond3_3, n_contour, edge_bin)
end