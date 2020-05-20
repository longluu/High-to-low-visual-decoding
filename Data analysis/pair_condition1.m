%%%%%%%%%%%%%% Determine the pairing order the estimate 1 and 2 in Condition 1 %%%%%%%%%%%%%%
subjectID = {'km', 'ccp', 'cm', 'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
eps_correlation = 0.0001;

for ii = 1 : length(subjectID)
    %% Extract the data
    [estimate1_allSubject, estimate2_allSubject] = getData(subjectID{ii}, session, experimentNumber, useSplit1line);

    %% Determine pairing order for condition 1
    % Find the permutation of data that has lowest correlation
    estimateStim1_collapse = estimate1_allSubject(1, :);
    estimateStim2_collapse = estimate2_allSubject(1, :);
    corr_2line = 1;
    while abs(corr_2line) > eps_correlation
        [est1, ind1] = datasample(estimateStim1_collapse, length(estimateStim1_collapse), 'Replace', false);
        [est2, ind2] = datasample(estimateStim2_collapse, length(estimateStim2_collapse), 'Replace', false);
        corr_2line = corr(est1, est2);
    end
    
    %% Save
    file_name = ['ind_pair_cond1_' subjectID{ii}];
    save(file_name, 'ind1', 'ind2')
end