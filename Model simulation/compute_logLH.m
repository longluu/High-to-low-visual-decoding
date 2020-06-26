%%%%%%%%%%%%%%%%%%%%% Compute the log likelihood of model for 2 line with blank condition %%%%%%%%%%%%%%%%%%%%
%% Experiment
% Get data
subjectID = {'combined'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
[estimate1, estimate2] = getData(subjectID{1}, session, experimentNumber, useSplit1line);

% Get the empirical distribution
edge_bin = 10:1:81;
n_bin_1d = length(edge_bin);
bin_count_data = histcounts2(estimate2(3, :), estimate1(3, :), edge_bin, edge_bin);

% Plot to check
figure
hold on
imagesc(bin_count_data);
axis xy; 
colormap gray
plot(1:n_bin_1d, 1:n_bin_1d, '--white')
axis off
axis square

%% Model
% Load the model
load('modelFit-3m.mat');

% Get the model distribution
bin_count_model_sensory = histcounts2(estimate_cond3_1(2, :), estimate_cond3_1(1, :), edge_bin, edge_bin);
bin_count_model_memory = histcounts2(estimate_cond3_2(2, :), estimate_cond3_2(1, :), edge_bin, edge_bin);
bin_count_model_estimate = histcounts2(estimate_cond3_3(2, :), estimate_cond3_3(1, :), edge_bin, edge_bin);

% Plot to check
figure
hold on
for ii = 1 : 3
    subplot(1, 3, ii)
    hold on
    if ii == 1
        imagesc(bin_count_model_sensory);
    elseif ii == 2
        imagesc(bin_count_model_memory);
    elseif ii == 3
        imagesc(bin_count_model_estimate);
    end
    axis xy; 
    colormap gray
    plot(1:n_bin_1d, 1:n_bin_1d, '--white')
    axis off
    axis square
end

% Normalize the model bin count to get probability distribution
n_sample_model = size(estimate_cond3_1, 2);
p_model_sensory = bin_count_model_sensory / n_sample_model;
p_model_memory = bin_count_model_memory / n_sample_model;
p_model_estimate = bin_count_model_estimate / n_sample_model;

%% Compute the log likelihood
logLH_sensory = bin_count_data .* log(p_model_sensory);
logLH_sensory = logLH_sensory(:);
logLH_sensory(isinf(logLH_sensory) | isnan(logLH_sensory)) = [];
logLH_sensory = sum(logLH_sensory);

logLH_memory = bin_count_data .* log(p_model_memory);
logLH_memory = logLH_memory(:);
logLH_memory(isinf(logLH_memory) | isnan(logLH_memory)) = [];
logLH_memory = sum(logLH_memory);

logLH_estimate = bin_count_data .* log(p_model_estimate);
logLH_estimate = logLH_estimate(:);
logLH_estimate(isinf(logLH_estimate) | isnan(logLH_estimate)) = [];
logLH_estimate = sum(logLH_estimate);

figure
bar(1:2, [logLH_sensory-logLH_estimate logLH_memory-logLH_estimate])
set(gca, 'XTick', 1:2, 'XTickLabel', {'Sensory', 'Memory'})
xlabel('Model')
ylabel('Log likelihood relative to Estimate model')
xlim([0 3])
box off