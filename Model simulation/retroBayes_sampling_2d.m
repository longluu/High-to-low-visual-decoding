%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate the retrospective Bayesian model (1-dimension) %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [estimate_cond1, estimate_cond2, estimate_cond3_1, estimate_cond3_2, estimate_cond3_3] = retroBayes_sampling_2d(theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1, std_memory_redraw_2, std_memory_redraw_3,...
            n_trials, n_discreteStep)
tic
% Construct a grid of 2d theta
min_theta = min(theta_stim) - 3  * sqrt(max(std_sensory.^2) + std_memory_1line^2 + std_memory_2line^2 + std_memory_redraw_1^2);
max_theta = max(theta_stim) + 3  * sqrt(max(std_sensory.^2) + std_memory_1line^2 + std_memory_2line^2 + std_memory_redraw_1^2);
theta = linspace(min_theta, max_theta, n_discreteStep);
[X1,X2] = meshgrid(theta, theta);
ind_correct = X1 < X2;
ind_incorrect = X1 > X2;

% Sensory measurement of 2 lines
sensory_sample_1 = normrnd(theta_stim(1), std_sensory(1), [1 n_trials]);
sensory_sample_2 = normrnd(theta_stim(2), std_sensory(2), [1 n_trials]);

% Ordinal decision based on sensory evidence
ordinal_decision = sensory_sample_2 > sensory_sample_1;

%% Model 1: low-to-high (Condition 1)
memory_sample_1_1line = normrnd(sensory_sample_1, std_memory_1line);
memory_sample_2_1line = normrnd(sensory_sample_2, std_memory_1line);

%% Model 2a - 2 line basic (for Condition 2)
memory_sample_1_2line = normrnd(sensory_sample_1, std_memory_2line);
memory_sample_2_2line = normrnd(sensory_sample_2, std_memory_2line);
posterior_1_2line = NaN(n_trials, n_discreteStep);
posterior_2_2line = NaN(n_trials, n_discreteStep);
std_1_2line = sqrt(std_sensory(1)^2 + std_memory_2line^2);
std_2_2line = sqrt(std_sensory(2)^2 + std_memory_2line^2);

for ii = 1 : n_trials
    % Normal likelihood aroundn the sensory measurement    
    likelihood = normpdf(theta, memory_sample_2_2line(ii),std_2_2line )' *...
                  normpdf(theta, memory_sample_1_2line(ii), std_1_2line)   ;
    
    % Chop off the inconsistent part to form posterior
    posterior = likelihood;
    if ordinal_decision(ii) == 0
        posterior(ind_correct) = 0;
    else
        posterior(ind_incorrect) = 0;
    end

    % Marginal distribution
    posterior_1_2line(ii, :) = sum(posterior, 1);
    posterior_2_2line(ii, :) = sum(posterior, 2);
        
end

% Estimate
posterior_1_2line = posterior_1_2line ./ repmat(sum(posterior_1_2line, 2), 1, n_discreteStep);
posterior_2_2line = posterior_2_2line ./ repmat(sum(posterior_2_2line, 2), 1, n_discreteStep);
estimate_1_2line = posterior_1_2line * theta';
estimate_2_2line = posterior_2_2line * theta';
clear likelihood posterior posterior_1_2line posterior_2_2line

if nargout > 3
    %% Model 2b - 2 line redecode, resample from sensory sample (for Condition 3)
    memory_sample_1_2line_redraw_1 = normrnd(sensory_sample_1, std_memory_redraw_1);
    memory_sample_2_2line_redraw_1 = normrnd(sensory_sample_2, std_memory_redraw_1);
    posterior_1_2line_redraw_1 = NaN(n_trials, n_discreteStep);
    posterior_2_2line_redraw_1 = NaN(n_trials, n_discreteStep);
    std_1_2line_redraw_1 = sqrt(std_sensory(1)^2 + std_memory_redraw_1^2);
    std_2_2line_redraw_1 = sqrt(std_sensory(2)^2 + std_memory_redraw_1^2);

    for ii = 1 : n_trials
        % Normal likelihood aroundn the memory sample
        likelihood = normpdf(theta, memory_sample_2_2line_redraw_1(ii),std_2_2line_redraw_1 )' *...
                      normpdf(theta, memory_sample_1_2line_redraw_1(ii), std_1_2line_redraw_1)   ;

        % Chop off the inconsistent part to form posterior
        posterior = likelihood;
        if ordinal_decision(ii) == 0
            posterior(ind_correct) = 0;
        else
            posterior(ind_incorrect) = 0;
        end

        % Marginal distribution
        posterior_1_2line_redraw_1(ii, :) = sum(posterior, 1);
        posterior_2_2line_redraw_1(ii, :) = sum(posterior, 2);    
    end

    % Estimate
    posterior_1_2line_redraw_1 = posterior_1_2line_redraw_1 ./ repmat(sum(posterior_1_2line_redraw_1, 2), 1, n_discreteStep);
    posterior_2_2line_redraw_1 = posterior_2_2line_redraw_1 ./ repmat(sum(posterior_2_2line_redraw_1, 2), 1, n_discreteStep);
    estimate_1_redraw_1 = posterior_1_2line_redraw_1 * theta';
    estimate_2_redraw_1 = posterior_2_2line_redraw_1 * theta';
    clear likelihood posterior posterior_1_2line_redraw_1 posterior_2_2line_redraw_1 ...
            memory_sample_1_2line_redraw_1 memory_sample_2_2line_redraw_1

    %% Model 2c: 2 line redecode, resample from first memory sample (for Condition 3)
    memory_sample_1_2line_redraw_2 = normrnd(memory_sample_1_2line, std_memory_redraw_3);
    memory_sample_2_2line_redraw_2 = normrnd(memory_sample_2_2line, std_memory_redraw_3);
    posterior_1_2line_redraw_2 = NaN(n_trials, n_discreteStep);
    posterior_2_2line_redraw_2 = NaN(n_trials, n_discreteStep);
    std_1_2line_redraw_2 = sqrt(std_sensory(1)^2 + std_memory_2line^2 + std_memory_redraw_2^2);
    std_2_2line_redraw_2 = sqrt(std_sensory(2)^2 + std_memory_2line^2 + std_memory_redraw_2^2);
    for ii = 1 : n_trials
        % Model 2a - 2 line basic (for Condition 2)
        % Normal likelihood around the memory sample   
        likelihood = normpdf(theta, memory_sample_2_2line_redraw_2(ii), std_2_2line_redraw_2 )' *...
                      normpdf(theta, memory_sample_1_2line_redraw_2(ii), std_1_2line_redraw_2)   ;

        % Chop off the inconsistent part to form posterior
        posterior = likelihood;
        if ordinal_decision(ii) == 0
            posterior(ind_correct) = 0;
        else
            posterior(ind_incorrect) = 0;
        end

        % Marginal distribution
        posterior_1_2line_redraw_2(ii, :) = sum(posterior, 1);
        posterior_2_2line_redraw_2(ii, :) = sum(posterior, 2);
    end
    posterior_1_2line_redraw_2 = posterior_1_2line_redraw_2 ./ repmat(sum(posterior_1_2line_redraw_2, 2), 1, n_discreteStep);
    posterior_2_2line_redraw_2 = posterior_2_2line_redraw_2 ./ repmat(sum(posterior_2_2line_redraw_2, 2), 1, n_discreteStep);
    estimate_1_redraw_2 = posterior_1_2line_redraw_2 * theta';
    estimate_2_redraw_2 = posterior_2_2line_redraw_2 * theta';
    clear likelihood posterior posterior_1_2line_redraw_2 posterior_2_2line_redraw_2 ...
            memory_sample_1_2line_redraw_2 memory_sample_2_2line_redraw_2

    %% Model 2d: 2 line redecode, resample from first estimate (for Condition 3)
    memory_sample_1_2line_redraw_3 = normrnd(estimate_1_2line, std_memory_redraw_3);
    memory_sample_2_2line_redraw_3 = normrnd(estimate_2_2line, std_memory_redraw_3);
    posterior_1_2line_redraw_3 = NaN(n_trials, n_discreteStep);
    posterior_2_2line_redraw_3 = NaN(n_trials, n_discreteStep);
    std_1_2line_redraw_3 = sqrt(std_sensory(1)^2 + std_memory_2line^2 + std_memory_redraw_3^2);
    std_2_2line_redraw_3 = sqrt(std_sensory(2)^2 + std_memory_2line^2 + std_memory_redraw_3^2);
    for ii = 1 : n_trials
        % Model 2a - 2 line basic (for Condition 2)
        % Normal likelihood around the memory sample   
        likelihood = normpdf(theta, memory_sample_2_2line_redraw_3(ii), std_2_2line_redraw_3 )' *...
                      normpdf(theta, memory_sample_1_2line_redraw_3(ii), std_1_2line_redraw_3)   ;

        % Chop off the inconsistent part to form posterior
        posterior = likelihood;
        if ordinal_decision(ii) == 0
            posterior(ind_correct) = 0;
        else
            posterior(ind_incorrect) = 0;
        end

        % Marginal distribution
        posterior_1_2line_redraw_3(ii, :) = sum(posterior, 1);
        posterior_2_2line_redraw_3(ii, :) = sum(posterior, 2);
    end
    posterior_1_2line_redraw_3 = posterior_1_2line_redraw_3 ./ repmat(sum(posterior_1_2line_redraw_3, 2), 1, n_discreteStep);
    posterior_2_2line_redraw_3 = posterior_2_2line_redraw_3 ./ repmat(sum(posterior_2_2line_redraw_3, 2), 1, n_discreteStep);
    estimate_1_redraw_3 = posterior_1_2line_redraw_3 * theta';
    estimate_2_redraw_3 = posterior_2_2line_redraw_3 * theta';
    clear likelihood posterior posterior_1_2line_redraw_3 posterior_2_2line_redraw_3 ...
            memory_sample_1_2line_redraw_3 memory_sample_2_2line_redraw_3 memory_sample_1_2line memory_sample_2_2line
    
end
toc

%% Save the results
estimate_cond1 = [memory_sample_1_1line; memory_sample_2_1line];
estimate_cond2 = [estimate_1_2line'; estimate_2_2line'];
if nargout > 3
    half_trial = n_trials / 2;
    estimate_cond3_1 = [[estimate_1_2line(1:half_trial); estimate_1_redraw_1(half_trial+1:end)]';...
                        [estimate_2_redraw_1(1:half_trial); estimate_2_2line(half_trial+1:end)]'];
    estimate_cond3_2 = [[estimate_1_2line(1:half_trial); estimate_1_redraw_2(half_trial+1:end)]';...
                        [estimate_2_redraw_2(1:half_trial); estimate_2_2line(half_trial+1:end)]'];                
    estimate_cond3_3 = [[estimate_1_2line(1:half_trial); estimate_1_redraw_3(half_trial+1:end)]';...
                        [estimate_2_redraw_3(1:half_trial); estimate_2_2line(half_trial+1:end)]'];
end