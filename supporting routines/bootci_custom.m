function [confInterval_all, sampleStat_all, bootStat_all] = bootci_custom(nBootstrap, excludeIncorrectTrial, experiment_condition, alpha, estimateStim1_collapse, estimateStim2_collapse)
%% Compute the statistics on bootstrap samples
nEstStim1 = length(estimateStim1_collapse);
nEstStim2 = length(estimateStim2_collapse);
meanDiff_boot = NaN(nBootstrap, 1);
corr_boot = NaN(nBootstrap, 1);
percentCorrect_boot = NaN(nBootstrap, 1);
for ii = 1 : nBootstrap
    % Get the sample
    [estimateStim1_bootstrap, index_sample] = datasample(estimateStim1_collapse, nEstStim1);
    if experiment_condition == 1
        estimateStim2_bootstrap = datasample(estimateStim2_collapse, nEstStim2);
    else
        estimateStim2_bootstrap = estimateStim2_collapse(index_sample);        
    end
    
    % Compute percent correct first before altering the estimates
    percentCorrect_boot(ii) = 100*sum(estimateStim2_bootstrap >= estimateStim1_bootstrap) / nEstStim1;                         

    % Alter the estimates (remove or flip incorrect trials)
    diffEst = estimateStim2_bootstrap - estimateStim1_bootstrap;
    if excludeIncorrectTrial == 1
        indInclude = diffEst > 0; % include only correct trials
    else
        indInclude = 1:length(diffEst);
    end
    if excludeIncorrectTrial == 2
        indFlip = diffEst < 0; % flip the incorrect trials
        estStim1_swap = estimateStim1_bootstrap(indFlip);
        estStim2_swap = estimateStim2_bootstrap(indFlip);
        estimateStim1_bootstrap(indFlip) = estStim2_swap;
        estimateStim2_bootstrap(indFlip) = estStim1_swap;                    
    end
    
    % Compute mean and correlation
    meanDiff_boot(ii) = mean(estimateStim2_bootstrap(indInclude)) - mean(estimateStim1_bootstrap(indInclude));
    corr_boot(ii) = corr(estimateStim1_bootstrap(indInclude),...
                                estimateStim2_bootstrap(indInclude));
end   
bootStat_all = [meanDiff_boot corr_boot percentCorrect_boot];

%% Compute the statistics on actual samples
% Mean and correlation
if experiment_condition == 1
    meanDiff_sample = mean(meanDiff_boot);
    corr_sample = mean(corr_boot);
else
    if excludeIncorrectTrial == 0
        meanDiff_sample = mean(estimateStim2_collapse) - mean(estimateStim1_collapse);
        corr_sample = corr(estimateStim2_collapse, estimateStim1_collapse);
    elseif excludeIncorrectTrial == 1
        diffEst_original = estimateStim2_collapse - estimateStim1_collapse;    
        meanDiff_sample = mean(estimateStim2_collapse(diffEst_original>0)) - mean(estimateStim1_collapse(diffEst_original>0));
        corr_sample = corr(estimateStim2_collapse(diffEst_original>0), estimateStim1_collapse(diffEst_original>0));
    elseif excludeIncorrectTrial == 2
        diffEst_original = estimateStim2_collapse - estimateStim1_collapse;  
        indFlip = diffEst_original < 0; % flip the incorrect trials
        estStim1 = estimateStim1_collapse;
        estStim2 = estimateStim2_collapse;
        estStim1(indFlip) = estimateStim2_collapse(indFlip);
        estStim2(indFlip) = estimateStim1_collapse(indFlip);
        meanDiff_sample = mean(estStim2) - mean(estStim1);
        corr_sample = corr(estStim2, estStim1);        
    end
end

% Percent correct
if experiment_condition == 1
    percentCorrect_sample = mean(percentCorrect_boot);
else
    percentCorrect_sample = 100 * sum(estimateStim2_collapse >= estimateStim1_collapse) / length(estimateStim2_collapse);
end

sampleStat_all = [meanDiff_sample corr_sample percentCorrect_sample];

%% Compute confidence interval
confInterval_all(1, 1) = prctile(meanDiff_boot, alpha*100/2);
confInterval_all(2, 1) = prctile(meanDiff_boot, 100 - alpha*100/2);
confInterval_all(1, 2) = prctile(corr_boot, alpha*100/2);
confInterval_all(2, 2) = prctile(corr_boot, 100 - alpha*100/2);
confInterval_all(1, 3) = prctile(percentCorrect_boot, alpha*100/2);
confInterval_all(2, 3) = prctile(percentCorrect_boot, 100 - alpha*100/2);

% %% Compute the confidence interval
% confInterval_all = NaN(2, 3);
% for ii = 1 : 3
%     % Calculate the bias-corrected and accelerated parameters z0 and a
%     z0 = inverseNormalCDF(sum(bootStat_all(:, ii) < sampleStat_all(ii))/nReps);
% 
%     thetai = zeros(1,length(x));
%     %calculate the statistic holding one member of x out each time
%     for i=1:length(x)
%         id = [1:(i-1),(i+1):length(x)];
%         thetai(i) = myStatistic(x(id));
%     end
%     %do something related to skewness.
%     a = sum( (mean(thetai)-thetai).^3)/(6*(sum( (mean(thetai)-thetai).^2).^(3/2)));
% 
% 
%     % Calculate the 'bias-corrected and accelerated' percentiles using z0 and a
%     CIrange = 100 - alpha_level;
%     zLo = inverseNormalCDF((1-CIrange/100)/2);
%     zHi = inverseNormalCDF((1+CIrange/100)/2);
% 
%     zClo = z0 + (z0+zLo)/(1-a*(z0+zLo));
%     bcaLo = NormalCumulative(zClo,0,1);
% 
%     zChi = z0 + (z0+zHi)/(1-a*(z0+zHi));
%     bcaHi = NormalCumulative(zChi,0,1);
% 
%     CI(1) = prctile(bootstrapStat,100*bcaLo);
%     CI(2) = prctile(bootstrapStat,100*bcaHi);
% end