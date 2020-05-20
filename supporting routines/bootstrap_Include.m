function [output] = bootstrap_Include(estMatrix)
    % Get the sample
    estimateStim1_bootstrap = estMatrix(:, 1);
    estimateStim2_bootstrap = estMatrix(:, 2);

    % Compute the statistics
    meanDiff_boot = mean(estimateStim2_bootstrap) - mean(estimateStim1_bootstrap);
    corr_boot = corr(estimateStim1_bootstrap, estimateStim2_bootstrap);
    percentCorrect_boot = sum(estimateStim2_bootstrap > estimateStim1_bootstrap) / length(estimateStim2_bootstrap);
    output = [meanDiff_boot, corr_boot, percentCorrect_boot];
end