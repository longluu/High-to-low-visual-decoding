function [diffEst, diffEst_same, diffEst_diff] = split_diffEst_trial(n_back, excludeIncorrectTrial, ...
        trialOrder, estimateLine1, estimateLine2, indStim1Est1, indStim1Est2, indStim2Est1, indStim2Est2)
%% Split the trials
if n_back == 1                        
    trial_current = trialOrder(2:end);
    trial_last = trialOrder(1:end-1); 
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);    
elseif n_back == 2
    trial_current = trialOrder(3:end);
    trial_last = trialOrder(1:end-2);   
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);   
elseif n_back == 3
    trial_current = trialOrder(4:end);
    trial_last = trialOrder(1:end-3);   
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);  
elseif n_back == 4
    trial_current = trialOrder(5:end);
    trial_last = trialOrder(1:end-4);   
    ind_same = (trial_current > 25 & trial_last > 25) | (trial_current <= 25 & trial_last <= 25);
    ind_diff = (trial_current > 25 & trial_last <= 25) | (trial_current <= 25 & trial_last > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);  
elseif n_back == 12
    trial_current = trialOrder(3:end);
    trial_last_1 = trialOrder(2:end-1); 
    trial_last_2 = trialOrder(1:end-2); 
    ind_same = (trial_current > 25 & trial_last_1 > 25 & trial_last_2 > 25) | (trial_current <= 25 & trial_last_1 <= 25 & trial_last_2 <= 25);
    ind_diff = (trial_current > 25 & trial_last_1 <= 25 & trial_last_2 <= 25) | (trial_current <= 25 & trial_last_1 > 25 & trial_last_2 > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);  
elseif n_back == 123
    trial_current = trialOrder(4:end);
    trial_last_1 = trialOrder(3:end-1); 
    trial_last_2 = trialOrder(2:end-2); 
    trial_last_3 = trialOrder(1:end-3);
    ind_same = (trial_current > 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25) |...
                (trial_current <= 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25);
    ind_diff = (trial_current > 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25) |...
                (trial_current <= 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);   
elseif n_back == 1234
    trial_current = trialOrder(5:end);
    trial_last_1 = trialOrder(4:end-1); 
    trial_last_2 = trialOrder(3:end-2); 
    trial_last_3 = trialOrder(2:end-3);
    trial_last_4 = trialOrder(1:end-4);
    ind_same = (trial_current > 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25 & trial_last_4 > 25) |...
                (trial_current <= 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25 & trial_last_4 <= 25);
    ind_diff = (trial_current > 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25  & trial_last_4 <= 25) |...
                (trial_current <= 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25 & trial_last_4 > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);  
elseif n_back == 12345
    trial_current = trialOrder(6:end);
    trial_last_1 = trialOrder(5:end-1); 
    trial_last_2 = trialOrder(4:end-2); 
    trial_last_3 = trialOrder(3:end-3);
    trial_last_4 = trialOrder(2:end-4);
    trial_last_5 = trialOrder(1:end-5);
    ind_same = (trial_current > 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25 & trial_last_4 > 25 & trial_last_5 > 25) |...
                (trial_current <= 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25 & trial_last_4 <= 25 & trial_last_5 <= 25);
    ind_diff = (trial_current > 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25  & trial_last_4 <= 25  & trial_last_5 <= 25) |...
                (trial_current <= 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25 & trial_last_4 > 25 & trial_last_5 > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff); 
elseif n_back == 123456
    trial_current = trialOrder(7:end);
    trial_last_1 = trialOrder(6:end-1); 
    trial_last_2 = trialOrder(5:end-2); 
    trial_last_3 = trialOrder(4:end-3);
    trial_last_4 = trialOrder(3:end-4);
    trial_last_5 = trialOrder(2:end-5);
    trial_last_6 = trialOrder(1:end-6);
    ind_same = (trial_current > 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25 & trial_last_4 > 25 & trial_last_5 > 25 & trial_last_6 > 25) |...
                (trial_current <= 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25 & trial_last_4 <= 25 & trial_last_5 <= 25 & trial_last_6 <= 25);
    ind_diff = (trial_current > 25 & trial_last_1 <= 25 & trial_last_2 <= 25 & trial_last_3 <= 25  & trial_last_4 <= 25  & trial_last_5 <= 25 & trial_last_6 <= 25) |...
                (trial_current <= 25 & trial_last_1 > 25 & trial_last_2 > 25 & trial_last_3 > 25 & trial_last_4 > 25 & trial_last_5 > 25 & trial_last_6 > 25);
    trial_same = trial_current(ind_same);
    trial_diff = trial_current(ind_diff);     
end

estimateLine1_same = estimateLine1;
indExclude = setdiff(trialOrder, trial_same);
estimateLine1_same(indExclude) = NaN;
estimateLine2_same = estimateLine2;
estimateLine2_same(indExclude) = NaN;
estimateStim1_same = [estimateLine1_same(indStim1Est1); estimateLine2_same(indStim1Est2)];
estimateStim2_same = [estimateLine2_same(indStim2Est2);estimateLine1_same(indStim2Est1)];
estimateStim1_same = estimateStim1_same(~isnan(estimateStim1_same));
estimateStim2_same = estimateStim2_same(~isnan(estimateStim2_same));
estimateStim1_same(estimateStim1_same<0) = estimateStim1_same(estimateStim1_same<0)+180;
estimateStim2_same(estimateStim2_same<0) = estimateStim2_same(estimateStim2_same<0)+180;  

estimateLine1_diff = estimateLine1;
indExclude = setdiff(trialOrder, trial_diff);
estimateLine1_diff(indExclude) = NaN;
estimateLine2_diff = estimateLine2;
estimateLine2_diff(indExclude) = NaN;
estimateStim1_diff = [estimateLine1_diff(indStim1Est1); estimateLine2_diff(indStim1Est2)];
estimateStim2_diff = [estimateLine2_diff(indStim2Est2);estimateLine1_diff(indStim2Est1)];
estimateStim1_diff = estimateStim1_diff(~isnan(estimateStim1_diff));
estimateStim2_diff = estimateStim2_diff(~isnan(estimateStim2_diff));
estimateStim1_diff(estimateStim1_diff<0) = estimateStim1_diff(estimateStim1_diff<0)+180;
estimateStim2_diff(estimateStim2_diff<0) = estimateStim2_diff(estimateStim2_diff<0)+180;   

%% Collapse all trials to get total mean difference
estimateStim1 = [estimateLine1(indStim1Est1); estimateLine2(indStim1Est2)];
estimateStim2 = [estimateLine2(indStim2Est2);estimateLine1(indStim2Est1)];
estimateStim1_collapse = estimateStim1(~isnan(estimateStim1));
estimateStim2_collapse = estimateStim2(~isnan(estimateStim2));

%% Compute the mean difference
if excludeIncorrectTrial == 0
    % Split trials
    diffEst_same = estimateStim2_same - estimateStim1_same;
    diffEst_diff = estimateStim2_diff - estimateStim1_diff;
    
    % All trials
    diffEst = estimateStim2_collapse - estimateStim1_collapse;
elseif excludeIncorrectTrial == 1
    % Split trials
    diffEst_same = estimateStim2_same - estimateStim1_same;
    diffEst_diff = estimateStim2_diff - estimateStim1_diff;
    diffEst_same(diffEst_same<0) = [];
    diffEst_diff(diffEst_diff<0) = [];     
    
    % All trials
    diffEst_original = estimateStim2_collapse - estimateStim1_collapse;    
    diffEst = diffEst_original((diffEst_original>0));    
else
    % Split trials
    indFlip = (estimateStim2_same - estimateStim1_same) < 0; % flip the incorrect trials
    estStim1_swap = estimateStim1_same(indFlip);
    estStim2_swap = estimateStim2_same(indFlip);
    estimateStim1_same(indFlip) = estStim2_swap;
    estimateStim2_same(indFlip) = estStim1_swap;  
    diffEst_same = estimateStim2_same - estimateStim1_same;

    indFlip = (estimateStim2_diff - estimateStim1_diff) < 0; % flip the incorrect trials
    estStim1_swap = estimateStim1_diff(indFlip);
    estStim2_swap = estimateStim2_diff(indFlip);
    estimateStim1_diff(indFlip) = estStim2_swap;
    estimateStim2_diff(indFlip) = estStim1_swap;  
    diffEst_diff = estimateStim2_diff - estimateStim1_diff;  
    
    % All trials
    diffEst_original = estimateStim2_collapse - estimateStim1_collapse;  
    indFlip = diffEst_original < 0; % flip the incorrect trials
    estStim1 = estimateStim1_collapse;
    estStim2 = estimateStim2_collapse;
    estStim1(indFlip) = estimateStim2_collapse(indFlip);
    estStim2(indFlip) = estimateStim1_collapse(indFlip);
    diffEst = estStim2 - estStim1;    
end
        