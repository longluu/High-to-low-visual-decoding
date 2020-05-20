%%%%%%%%%%%%%%%%%%%%% Compare model (sampling implementation) and data %%%%%%%%%%%%%%%%%%%%

%% Parameters
% Experiment
subjectID = 'combined';
useSplit1line = 1;
session = 1;
experimentNumber = 1;

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

%% Model
% Get the predictive distribution p(thHat|th)
[estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3] = retroBayes_sampling_2d(theta_stim, std_sensory, std_memory_1line, ...
    std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, n_trials, n_discreteStep);


%% Plot the results
tic
plot_modelSampling_data(estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3, estimate1_allSubject, estimate2_allSubject, ...
                        theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3, ...
                        markerSize, markerColor, dotMeanSize, lineWidth, n_point_kdensity, smooth_sampling, edge_bin, h_filter, n_contour, n_bins)
toc