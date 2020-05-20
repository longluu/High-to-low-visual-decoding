%%%%%%%%%%%%%%%%%%%%% Compare model and data for all conditions %%%%%%%%%%%%%%%%%%%%
%%%% For condition 1 and 2, we use marginalized model fit to data.
%%%% For condition 3, we use the parameters fit to cond 1 and 2 to predict data with sampling model

%% Parameters
rerun_model = 0;

% Experiment
subjectID = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;

% Model
theta_stim = [49 54];
param_fit_all = extract_fit_param(subjectID);
n_trials = 4000000;
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

%% Run the model
if rerun_model == 1
    model_mean_est_cond1 = NaN(2, length(subjectID));
    model_mean_est_cond2 = NaN(2, length(subjectID));
    model_mean_est_cond3_1 = NaN(2, length(subjectID));
    model_mean_est_cond3_2 = NaN(2, length(subjectID));
    model_mean_est_cond3_3 = NaN(2, length(subjectID));

    for ss = 1 : length(subjectID)
        %% Extract the data
        [estimate1_allSubject, estimate2_allSubject] = getData(subjectID{ss}, session, experimentNumber, useSplit1line);

        %% Extract the model's parameters
        param_fit = param_fit_all{ss};
        std_sensory = param_fit(1) * [1 1];
        std_memory_1line = param_fit(2);
        std_memory_2line = param_fit(3);
        std_memory_redraw_1 = param_fit(3);
        std_memory_redraw_2 = param_fit(3);
        std_memory_redraw_3 = param_fit(3);
        if strcmp(subjectID{ss}, 'combined')
            const_bias = 0;
            n_bins = 60;
        else
            const_bias = param_fit(5);
            n_bins = 20;
        end

        %% Model - marginalized for Condition 1 and 2
        % Get the predictive distribution p(thHat|th)
        [pthetahGtheta_1, pthetahGtheta_2, ~, theta] = retroBayes_marginalized(theta_stim, std_sensory, std_memory_1line, ...
            std_memory_2line, std_memory_redraw_1, n_discreteStep);    

        % Shift the distribution by a constant amount
        step_theta = theta(2) - theta(1);
        const_shift = const_bias / sqrt(2);
        n_shift = round(const_shift / step_theta);
        pthetahGtheta_1 = circshift(pthetahGtheta_1, [n_shift -n_shift]);
        pthetahGtheta_2 = circshift(pthetahGtheta_2, [n_shift -n_shift]);

        % Project the distribution along the diagonal axis
        n_step_project = 2*n_discreteStep - 1;
        axis_project_middle = n_discreteStep;
        pProject_1 = NaN(1, n_step_project);
        pProject_2 = NaN(1, n_step_project);
        for ii = 1 : n_discreteStep-1
            % Condition 1
            pProject_1(axis_project_middle + ii) = sum(diag(pthetahGtheta_1, ii));
            pProject_1(axis_project_middle - ii) = sum(diag(pthetahGtheta_1, -ii));

            % Condition 2
            pProject_2(axis_project_middle + ii) = sum(diag(pthetahGtheta_2, ii));
            pProject_2(axis_project_middle - ii) = sum(diag(pthetahGtheta_2, -ii));

        end
        pProject_1(axis_project_middle) = sum(diag(pthetahGtheta_1));
        pProject_2(axis_project_middle) = sum(diag(pthetahGtheta_2));
        pProject_1 = fliplr(pProject_1) / sum(pProject_1);
        pProject_2 = fliplr(pProject_2) / sum(pProject_2);

        % Create the projection axis
        [theta1, theta2] = meshgrid(theta, theta);
        theta_diff = diag(rot90(theta2 - theta1));
        step_axis_project = (theta_diff(2) - theta_diff(1)) / 2;
        theta_augment = theta_diff(1:end-1) + step_axis_project;
        theta_augment = [theta_augment; 0];
        axis_project = [theta_diff theta_augment]';
        axis_project = axis_project(:);
        axis_project(end) = [];

        % Compute the mean
        pthetahGtheta_1_theta1 = sum(pthetahGtheta_1, 1);
        pthetahGtheta_1_theta1 = pthetahGtheta_1_theta1 / sum(pthetahGtheta_1_theta1);
        model_mean_est_cond1(1, ss) = pthetahGtheta_1_theta1 * theta';
        pthetahGtheta_1_theta2 = sum(pthetahGtheta_1, 2);
        pthetahGtheta_1_theta2 = pthetahGtheta_1_theta2 / sum(pthetahGtheta_1_theta2);
        model_mean_est_cond1(2, ss) = pthetahGtheta_1_theta2' * theta';

        pthetahGtheta_2_theta1 = sum(pthetahGtheta_2, 1);
        pthetahGtheta_2_theta1 = pthetahGtheta_2_theta1 / sum(pthetahGtheta_2_theta1);
        model_mean_est_cond2(1, ss) = pthetahGtheta_2_theta1 * theta';
        pthetahGtheta_2_theta2 = sum(pthetahGtheta_2, 2);
        pthetahGtheta_2_theta2 = pthetahGtheta_2_theta2 / sum(pthetahGtheta_2_theta2);
        model_mean_est_cond2(2, ss) = pthetahGtheta_2_theta2' * theta';

        %% Model - sampling for Condition 3
        % Get the predictive distribution p(thHat|th)
        [~, ~, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3] = retroBayes_sampling_2d(theta_stim, std_sensory, std_memory_1line, ...
            std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, n_trials, n_discreteStep);

        % Shift the estimates by a constant amount
        estimate_cond3_1(1, :) = estimate_cond3_1(1, :) - const_shift;
        estimate_cond3_1(2, :) = estimate_cond3_1(2, :) + const_shift;
        estimate_cond3_2(1, :) = estimate_cond3_2(1, :) - const_shift;
        estimate_cond3_2(2, :) = estimate_cond3_2(2, :) + const_shift;
        estimate_cond3_3(1, :) = estimate_cond3_3(1, :) - const_shift;
        estimate_cond3_3(2, :) = estimate_cond3_3(2, :) + const_shift;

        % Compute the mean
        model_mean_est_cond3_1(:, ss) = mean(estimate_cond3_1, 2);
        model_mean_est_cond3_2(:, ss) = mean(estimate_cond3_2, 2);
        model_mean_est_cond3_3(:, ss) = mean(estimate_cond3_3, 2);

        % Plot the results
        plot_model_data(pthetahGtheta_1, pthetahGtheta_2, pProject_1, pProject_2, theta, axis_project, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, estimate1_allSubject, estimate2_allSubject, ...
                                theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, ...
                                markerSize, markerColor, dotMeanSize, lineWidth, n_point_kdensity, smooth_sampling, edge_bin, h_filter, n_contour, n_bins, subjectID{ss})
    end
end
model_diff_est_cond12 = diff(model_mean_est_cond2, 1)  - diff(model_mean_est_cond1, 1);
model_diff_est_cond13_1 = diff(model_mean_est_cond3_1, 1)  - diff(model_mean_est_cond1, 1);
model_diff_est_cond13_2 = diff(model_mean_est_cond3_2, 1)  - diff(model_mean_est_cond1, 1);
model_diff_est_cond13_3 = diff(model_mean_est_cond3_3, 1)  - diff(model_mean_est_cond1, 1);

%% Analyze data
data_mean_est_cond1 = NaN(2, length(subjectID));
data_mean_est_cond2 = NaN(2, length(subjectID));
data_mean_est_cond3 = NaN(2, length(subjectID));

for ss = 1 : length(subjectID)
    % Extract the data
    [estimate1, estimate2] = getData(subjectID{ss}, session, experimentNumber, useSplit1line);
    
    % Compute the mean
    data_mean_est_cond1(1, ss) = mean(estimate1(1, :));
    data_mean_est_cond1(2, ss) = mean(estimate2(1, :));
    data_mean_est_cond2(1, ss) = mean(estimate1(2, :));
    data_mean_est_cond2(2, ss) = mean(estimate2(2, :));
    data_mean_est_cond3(1, ss) = mean(estimate1(3, :));
    data_mean_est_cond3(2, ss) = mean(estimate2(3, :));    
end
data_diff_est_cond12 = diff(data_mean_est_cond2, 1)  - diff(data_mean_est_cond1, 1);
data_diff_est_cond13 = diff(data_mean_est_cond3, 1)  - diff(data_mean_est_cond1, 1);

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
model_diff_est_cond13_best(ind_best == 2) = model_diff_est_cond13_2(ind_best == 1);
model_diff_est_cond13_best(ind_best == 3) = model_diff_est_cond13_3(ind_best == 1);


figure
hold on
minPlot = min([data_diff_est_cond12(:); data_diff_est_cond13(:); model_diff_est_cond12(:); model_diff_est_cond13_1(:); model_diff_est_cond13_2(:); model_diff_est_cond13_3(:)]) - 2;
maxPlot = max([data_diff_est_cond12(:); data_diff_est_cond13(:); model_diff_est_cond12(:); model_diff_est_cond13_1(:); model_diff_est_cond13_2(:); model_diff_est_cond13_3(:)]) + 2;
data_all = [data_diff_est_cond12' repmat(data_diff_est_cond13', 1, 4)];
model_all = [model_diff_est_cond12'  model_diff_est_cond13_1'  model_diff_est_cond13_2'  model_diff_est_cond13_3'  model_diff_est_cond13_best'];
fig_title = {'2-1';'3-1';'3-2'; '3-3';'3-best'};
for kk = 1 : 5
    subplot(1, 5, kk)
    [corr_index, pvalue] = corr(data_all(:, kk), model_all(:, kk));
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


%% Save the figures
% for ii = 2 : 16
%     file_name = ['model_fit_' subjectID{ii-1}];
%     saveas(ii, file_name, 'pdf')
% end

