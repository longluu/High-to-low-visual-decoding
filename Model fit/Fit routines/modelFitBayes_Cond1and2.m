function [fitParameter, negLogLH] = modelFitBayes_Cond1and2(SetStartPoint, initialValue, share_memory, include_motor_noise, include_const_bias, params, fileID)
%%%%%%%%%%%% Fitting the condition 1 and 2 using marginalized version

%% Initialize some variable
maxFunEval = 350;
tolFun = 1e-1;
tolX = 1e-2;
minStdSensory = 0.1;
maxStdSensory = 10;
minStdMemory1 = 0.1;
maxStdMemory1 = 10;

if share_memory
    minStdMemory2 = 0;
    maxStdMemory2 = 0;
else
    minStdMemory2 = 0.1;
    maxStdMemory2 = 10;
end
if include_motor_noise
    minStdMotor = 0.1;
    maxStdMotor = 10;
else
    minStdMotor = 0;
    maxStdMotor = 0;
end
if include_const_bias
    minConstBias = -5;
    maxConstBias = 5;
else
    minConstBias = 0;
    maxConstBias = 0;    
end
vlb = [minStdSensory minStdMemory1 minStdMemory2 minStdMotor minConstBias];
vub = [maxStdSensory maxStdMemory1 maxStdMemory2 maxStdMotor maxConstBias];    
if SetStartPoint
    searchParameter0 = initialValue;
else
    stdSensory_initial = initialValue(1) + 2*(rand-0.5);
    stdMemory1_initial =  initialValue(2) + 4*(rand-0.5);
    stdMemory2_initial =  initialValue(3) + 4*(rand-0.5);
    stdMotor_initial =  initialValue(4) + 2*(rand-0.5);
    constBias_initial =  initialValue(5) + 2*(rand-0.5);
    searchParameter0 = [stdSensory_initial stdMemory1_initial stdMemory2_initial stdMotor_initial constBias_initial];
end

%% Run the optimization
% Define the objective function
optimFun = @(searchParameter, varargin) computeNegLLH(searchParameter(1),searchParameter(2), searchParameter(3), searchParameter(4), searchParameter(5),...
    share_memory, include_motor_noise, params.theta_stim, params.n_discreteStep, params.estimate1, params.estimate2, fileID);

% Using Nealder-Mead simplex
searchOptions = struct(...
    'Display','Iter',...
    'TolFun',tolFun,...
    'TolX',tolX,...
    'MaxFunEvals',maxFunEval,...
    'FunValCheck','off');
[fitParameter, negLogLH] = fminsearchbnd(optimFun, searchParameter0, vlb, vub, searchOptions,[]);
end

function errorTotal = computeNegLLH(std_sensory, std_memory_1line, std_memory_2line, std_motor, const_bias, share_memory, include_motor_noise, theta_stim, n_discreteStep, estimate1, estimate2, fileID)

try
    %% Compute the full model prediction
    % Get the predictive distribution p(thHat|th)
    if share_memory
        std_memory_2line = std_memory_1line;
    end
    std_memory_redraw_1 = std_memory_1line;
    [pthetahGtheta_1, pthetahGtheta_2, ~, theta] = retroBayes_marginalized(theta_stim, std_sensory, std_memory_1line, ...
        std_memory_2line, std_memory_redraw_1, n_discreteStep); 
    if include_motor_noise
        theta_centered = theta - (max(theta) + min(theta))/2;
        pMotor = normpdf(theta_centered, 0, std_motor)' * normpdf(theta_centered, 0, std_motor);
        pthetahGtheta_1 = conv2(pthetahGtheta_1, pMotor, 'same');
        pthetahGtheta_2 = conv2(pthetahGtheta_2, pMotor, 'same');
    end
    
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

    % Shift the distributions by a constant amount
    n_shift = round(const_bias / step_axis_project);
    pProject_1 = circshift(pProject_1, n_shift);
    pProject_2 = circshift(pProject_2, n_shift);
    
    %% Error term for condition 1  
    epsZero = 10^(-5);
    tempEstimateModelX = axis_project;
    tempEstimateModelY = pProject_1; 
    tempEstimateModelY = tempEstimateModelY ./ diff(axis_project(1:2));
    tempEstimateDataX = estimate2(1, :) - estimate1(1, :);
    tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempEstimateDataX, 'pchip');
    tempLikelihood(tempLikelihood == 0) = epsZero; 
    logLikelihood_Cond1 = nansum(log(tempLikelihood));  
    
    %% Error term for condition 2  
    epsZero = 10^(-5);
    tempEstimateModelX = axis_project;
    tempEstimateModelY = pProject_2; 
    tempEstimateModelY = tempEstimateModelY ./ diff(axis_project(1:2));
    tempEstimateDataX = estimate2(2, :) - estimate1(2, :);
    tempLikelihood = interp1(tempEstimateModelX, tempEstimateModelY, tempEstimateDataX, 'pchip');
    tempLikelihood(tempLikelihood == 0) = epsZero; 
    logLikelihood_Cond2 = nansum(log(tempLikelihood));     
    
    errorTotal = -logLikelihood_Cond1 - logLikelihood_Cond2;
    disp(['-logLH:' num2str(round(errorTotal)) ' ' num2str(round(-logLikelihood_Cond1)) ...
            ' ' num2str(round(-logLikelihood_Cond2)) ', Params: ' num2str([roundn(std_sensory, 2) roundn(std_memory_1line, 2) roundn(std_memory_2line, 2) roundn(std_motor, 2) roundn(const_bias, 2)])])
    fprintf(fileID,'%2s %9.1f  %9.1f %9.1f %9.4f %9.4f %9.4f %9.4f %9.4f \r\n', '//', errorTotal, -logLikelihood_Cond1, -logLikelihood_Cond2,...
        std_sensory, std_memory_1line, std_memory_2line, std_motor, const_bias);
    
catch e
    keyboard
    rethrow(e)
end

end