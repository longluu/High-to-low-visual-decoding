%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate the retrospective Bayesian model (1-dimension) %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pthetahGtheta_1, pthetahGtheta_2, pthetahGtheta_3, theta] = retroBayes_marginalized(theta_stim, std_sensory, std_memory_1line, std_memory_2line, std_memory_redraw_1,...
             n_discreteStep)
        
% Define the standard deviation
if length(std_sensory) == 2
    std_sensory(2) = [];
end
std_1line = sqrt(std_sensory^2 + std_memory_1line^2); 
std_2line = sqrt(std_sensory^2 + std_memory_2line^2);

%% Condition 1
min_theta = min(theta_stim) - 3  * std_2line;
max_theta = max(theta_stim) + 3  * std_2line;
theta = linspace(min_theta, max_theta, n_discreteStep);
[X1,X2] = meshgrid(theta, theta);
ind_correct = X1 < X2;
ind_incorrect = X1 > X2;

min_m1 = theta_stim(1) - 3 * std_sensory;
max_m1 = theta_stim(1) + 3 * std_sensory;
m1 = linspace(min_m1, max_m1, n_discreteStep);

min_m2 = theta_stim(2) - 3 * std_sensory;
max_m2 = theta_stim(2) + 3 * std_sensory;
m2 = linspace(min_m2, max_m2, n_discreteStep);

pthetahGtheta_1 = normpdf(theta, theta_stim(2), std_1line )' *...
                  normpdf(theta, theta_stim(1), std_1line)   ;
              
%% Condition 2
min_mm1 = min_m1 - 3 * std_memory_2line;
max_mm1 = max_m1 + 3 * std_memory_2line;
mm1 = linspace(min_mm1, max_mm1, n_discreteStep);

min_mm2 = min_m2 - 3 * std_memory_2line;
max_mm2 = max_m2 + 3 * std_memory_2line;
mm2 = linspace(min_mm2, max_mm2, n_discreteStep);

% Inference: thetaHat_1, thetaHat_2
[MM1, MM2] = meshgrid(mm1, mm2);
MM1 = MM1(:);
MM2 = MM2(:);
posterior_1_2line_correct = NaN(length(MM1), n_discreteStep);
posterior_2_2line_correct = NaN(length(MM1), n_discreteStep);
posterior_1_2line_incorrect = NaN(length(MM1), n_discreteStep);
posterior_2_2line_incorrect = NaN(length(MM1), n_discreteStep);

for ii = 1 : length(MM1)
    % Normal likelihood around the memory sample p(mm1, mm2 | theta1, theta2)  
    likelihood = normpdf(theta, MM2(ii),std_2line )' *...
                  normpdf(theta, MM1(ii), std_2line)   ;

    % Chop off the inconsistent part to form posterior p(theta1, theta2 | mm1, mm2, CHat)  
    posterior_correct = likelihood;
    posterior_incorrect = likelihood;
    posterior_incorrect(ind_correct) = 0;
    posterior_correct(ind_incorrect) = 0;

    % Marginal distribution
    posterior_1_2line_correct(ii, :) = sum(posterior_correct, 1);
    posterior_2_2line_correct(ii, :) = sum(posterior_correct, 2);
    posterior_1_2line_incorrect(ii, :) = sum(posterior_incorrect, 1);
    posterior_2_2line_incorrect(ii, :) = sum(posterior_incorrect, 2);
end


% Compute the estimates thetaHat(mm1, mm2, Chat)
posterior_1_2line_correct = posterior_1_2line_correct ./ repmat(sum(posterior_1_2line_correct, 2), 1, n_discreteStep);
posterior_2_2line_correct = posterior_2_2line_correct ./ repmat(sum(posterior_2_2line_correct, 2), 1, n_discreteStep);
estimate_1_2line_correct = posterior_1_2line_correct * theta';
estimate_2_2line_correct = posterior_2_2line_correct * theta';

posterior_1_2line_incorrect = posterior_1_2line_incorrect ./ repmat(sum(posterior_1_2line_incorrect, 2), 1, n_discreteStep);
posterior_2_2line_incorrect = posterior_2_2line_incorrect ./ repmat(sum(posterior_2_2line_incorrect, 2), 1, n_discreteStep);
estimate_1_2line_incorrect = posterior_1_2line_incorrect * theta';
estimate_2_2line_incorrect = posterior_2_2line_incorrect * theta';

% Compute the conditional memory distribution p(mm1, mm2 | theta1, theta2, Chat)
[M1, M2] = meshgrid(m1, m2);
pmGth = normpdf(m2, theta_stim(2),std_sensory)' * normpdf(m1, theta_stim(1),std_sensory);
pmGthCcorrect = pmGth;
pmGthCcorrect(M1 > M2) = 0;
pmGthCcorrect = pmGthCcorrect(:);
pmGthCincorrect = pmGth;
pmGthCincorrect(M1 < M2) = 0;
pmGthCincorrect = pmGthCincorrect(:);
M1 = M1(:);
M2 = M2(:);
memory_distribution_correct = zeros(n_discreteStep, n_discreteStep);
memory_distribution_incorrect = zeros(n_discreteStep, n_discreteStep);

for ii = 1 : length(M1)
    % Normal distribution around the sensory sample p(mm1, mm2 | thetaHat1, thetaHat2)  
    memory_distribution_correct = memory_distribution_correct + ...
        normpdf(mm2, M2(ii),std_memory_2line )' * normpdf(mm1, M1(ii), std_memory_2line)  * pmGthCcorrect(ii) ;
    memory_distribution_incorrect = memory_distribution_incorrect + ...
        normpdf(mm2, M2(ii),std_memory_2line )' * normpdf(mm1, M1(ii), std_memory_2line)  * pmGthCincorrect(ii) ;    
end

% Marginalization by changing variable
% Correct
[f1x, f1y] = gradient(reshape(estimate_1_2line_correct, n_discreteStep, n_discreteStep), 1); 
[f2x, f2y] = gradient(reshape(estimate_2_2line_correct, n_discreteStep, n_discreteStep), 1); 
a = 1 ./ abs(f1x .* f2y - f1y .* f2x);
memory_distribution_correct = a .* memory_distribution_correct;
memory_distribution_correct = memory_distribution_correct(:);
F = scatteredInterpolant(estimate_1_2line_correct, estimate_2_2line_correct, memory_distribution_correct, 'linear', 'linear');
pthhGthCcorrect = F(X1, X2);
pthhGthCcorrect(pthhGthCcorrect<0) = 0;
pthhGthCcorrect = pthhGthCcorrect ./ sum(pthhGthCcorrect(:));

% Incorrect
[f1x, f1y] = gradient(reshape(estimate_1_2line_incorrect, n_discreteStep, n_discreteStep), 1); 
[f2x, f2y] = gradient(reshape(estimate_2_2line_incorrect, n_discreteStep, n_discreteStep), 1); 
a = 1 ./ abs(f1x .* f2y - f1y .* f2x);
memory_distribution_incorrect = a .* memory_distribution_incorrect;
memory_distribution_incorrect = memory_distribution_incorrect(:);
F = scatteredInterpolant(estimate_1_2line_incorrect, estimate_2_2line_incorrect, memory_distribution_incorrect, 'linear', 'linear');
pthhGthCincorrect = F(X1, X2);
pthhGthCincorrect(pthhGthCincorrect<0) = 0;
pthhGthCincorrect = pthhGthCincorrect ./ sum(pthhGthCincorrect(:)); 

% Compute the weight for correct and incorrect part
pNorm = normpdf(theta, theta_stim(2), std_sensory )' *...
                  normpdf(theta, theta_stim(1), std_sensory)   ;
pNorm = pNorm ./ sum(pNorm(:));
pCorrect = sum(sum(pNorm(ind_correct)));

% Final predictive distribution
if pCorrect > (1 - 1e-3)
    pthetahGtheta_2 = pthhGthCcorrect;
else
    pthetahGtheta_2 = pCorrect * pthhGthCcorrect + (1- pCorrect) * pthhGthCincorrect;
end

%% Condition 3 (pending)
if nargout == 4
    pthetahGtheta_3 = pthetahGtheta_2;
end
