function [fr_allCond, se_allCond] = ...
    plot_ctrAnalysis_trial(diffEst_Cond1, diffEst_Cond2_same, diffEst_Cond2_diff, diffEst_Cond3_same, diffEst_Cond3_diff, trial_type)

%% Extract the estimate
diffEst_Cond2_same = cell2mat(diffEst_Cond2_same');
diffEst_Cond2_diff = cell2mat(diffEst_Cond2_diff');
diffEst_Cond2_same_offset = diffEst_Cond2_same - diffEst_Cond1;
diffEst_Cond2_diff_offset = diffEst_Cond2_diff - diffEst_Cond1;

diffEst_Cond3_same = cell2mat(diffEst_Cond3_same');
diffEst_Cond3_diff = cell2mat(diffEst_Cond3_diff');
diffEst_Cond3_same_offset = diffEst_Cond3_same - diffEst_Cond1;
diffEst_Cond3_diff_offset = diffEst_Cond3_diff - diffEst_Cond1;

diffEst_allCond_same_offset = [diffEst_Cond2_same_offset; diffEst_Cond3_same_offset];
diffEst_allCond_diff_offset = [diffEst_Cond2_diff_offset; diffEst_Cond3_diff_offset];

%% Bootstrap the estimates
nBootstrap = 10000;

% Compute the statistics on bootstrap samples
fr_boot_Cond2 = NaN(nBootstrap, 1); 
fr_boot_Cond3 = NaN(nBootstrap, 1); 
fr_boot_allCond = NaN(nBootstrap, 1); 

for ii = 1 : nBootstrap
    % Get the resample
    diffEst_Cond2_same_boot = datasample(diffEst_Cond2_same_offset, length(diffEst_Cond2_same_offset));
    diffEst_Cond2_diff_boot = datasample(diffEst_Cond2_diff_offset, length(diffEst_Cond2_diff_offset));
      
    diffEst_Cond3_same_boot = datasample(diffEst_Cond3_same_offset, length(diffEst_Cond3_same_offset));
    diffEst_Cond3_diff_boot = datasample(diffEst_Cond3_diff_offset, length(diffEst_Cond3_diff_offset));

    diffEst_allCond_same_boot = datasample(diffEst_allCond_same_offset, length(diffEst_allCond_same_offset));
    diffEst_allCond_diff_boot = datasample(diffEst_allCond_diff_offset, length(diffEst_allCond_diff_offset));
    
    % Compute the mean difference on resample
    fr_boot_Cond2(ii) = (nanmean(diffEst_Cond2_diff_boot) - nanmean(diffEst_Cond2_same_boot)) / (nanmean(diffEst_Cond2_diff_boot) + nanmean(diffEst_Cond2_same_boot));
    fr_boot_Cond3(ii) = (nanmean(diffEst_Cond3_diff_boot) - nanmean(diffEst_Cond3_same_boot)) / (nanmean(diffEst_Cond3_diff_boot) + nanmean(diffEst_Cond3_same_boot));
    fr_boot_allCond(ii) = (nanmean(diffEst_allCond_diff_boot) - nanmean(diffEst_allCond_same_boot)) / (nanmean(diffEst_allCond_diff_boot) + nanmean(diffEst_allCond_same_boot));    
end   

% Compute the statistics on actual samples
fr_Cond2 = (nanmean(diffEst_Cond2_diff_offset) - nanmean(diffEst_Cond2_same_offset)) / (nanmean(diffEst_Cond2_diff_offset) + nanmean(diffEst_Cond2_same_offset));
fr_Cond3 = (nanmean(diffEst_Cond3_diff_offset) - nanmean(diffEst_Cond3_same_offset)) / (nanmean(diffEst_Cond3_diff_offset) + nanmean(diffEst_Cond3_same_offset));
fr_allCond = (nanmean(diffEst_allCond_diff_offset) - nanmean(diffEst_allCond_same_offset)) / (nanmean(diffEst_allCond_diff_offset) + nanmean(diffEst_allCond_same_offset));

% Compute the standard error
se_Cond2 = std(fr_boot_Cond2);
se_Cond3 = std(fr_boot_Cond3);
se_allCond = std(fr_boot_allCond);

%% Plot the result
figure
colorName = {'Crimson', 'DarkOrange', 'Teal', 'DodgerBlue'};
colorIndex = NaN(length(colorName), 3);
for ii = 1 : length(colorName)
    colorIndex(ii, :) = rgb(colorName{ii});
end
hold on
plot([0 3], [1 1], '--k')
errorBarGraph([fr_Cond2; fr_Cond3; fr_allCond], [fr_Cond2 - se_Cond2; fr_Cond3 - se_Cond3; fr_allCond - se_allCond], ...
            [fr_Cond2 + se_Cond2; fr_Cond3 + se_Cond3; fr_allCond + se_allCond], colorIndex)
box off
title(trial_type)
ylabel('Fraction of repulsion explained by cross-trial adaptation')