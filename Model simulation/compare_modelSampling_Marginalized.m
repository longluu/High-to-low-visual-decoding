%%%%%%%%%%%%%%%%%%%%% Compare sampling and marginalized versions of the model %%%%%%%%%%%%%%%%%%%%


%% Parameters
% Model
theta_stim = [49 54];
std_sensory = [2.7 2.7];
std_memory_1line = 4;
std_memory_2line = 5;
std_memory_redraw_1 = 6;
std_memory_redraw_2 = 6;
std_memory_redraw_3 = 6;
n_trials = 1000000;
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
n_bins = 50;

%% Extract the data
[estimate1_allSubject, estimate2_allSubject] = getData(subjectID, session, experimentNumber, useSplit1line);

%% Model - sampling
% Get the predictive distribution p(thHat|th)
[estimate_cond1, estimate_cond2] = retroBayes_sampling_2d(theta_stim, std_sensory, std_memory_1line, ...
    std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, n_trials, n_discreteStep);

%% Model - marginalized
% Get the predictive distribution p(thHat|th)
[pthetahGtheta_1, pthetahGtheta_2, pthetahGtheta_3, theta] = retroBayes_marginalized(theta_stim, std_sensory, std_memory_1line, ...
    std_memory_2line, std_memory_redraw_1, n_discreteStep);    

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

%% Plot the results
tic
plot_model_sampling_marginalized(estimate_cond1, estimate_cond2, theta, pthetahGtheta_1, pthetahGtheta_2, axis_project, pProject_1, pProject_2, ...
                        theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, ...
                        markerSize, markerColor, dotMeanSize, lineWidth, n_point_kdensity, smooth_sampling, edge_bin, h_filter, n_contour, n_bins)
toc