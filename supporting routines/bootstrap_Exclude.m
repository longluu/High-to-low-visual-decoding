function [output] = bootstrap_Exclude(estMatrix)
    % Get the sample
    estimateStim1_bootstrap = estMatrix(:, 1);
    estimateStim2_bootstrap = estMatrix(:, 2);

    % Compute the statistics
    diffEst = estimateStim2_bootstrap - estimateStim1_bootstrap;
    indInclude = diffEst > 0; % include only correct trials
    meanDiff_boot = mean(estimateStim2_bootstrap(indInclude)) - mean(estimateStim1_bootstrap(indInclude));
    corr_boot = corr(estimateStim1_bootstrap(indInclude), estimateStim2_bootstrap(indInclude));
    percentCorrect_boot = sum(estimateStim2_bootstrap > estimateStim1_bootstrap) / length(estimateStim2_bootstrap);
    output = [meanDiff_boot, corr_boot, percentCorrect_boot];
end