%%%%%%%%%%%%% Run the full model fitting program %%%%%%%%%%%%%
tic
clear 
nLoops = 20;
negLLHMinLoop = [];
fitParameterAll = cell(nLoops, 1);
negLLH = NaN(nLoops, 1);

%% Choose the parameter for the model fit
% Set starting point for optimization 0: chooose randomly
%                                     1: set starting point
SetStartPoint = 0;
initialValueAll = [2.5           4.5                 5.5             3           0];
                %[stdSensory   stdMemory_1line   stdMemory_2line   motor_noise    const_bias]
share_memory = 0;
include_motor_noise = 0;
include_const_bias = 1;
subjectID = {'lp', 'if', 'cr', 'an', 'sr', 'zt', 'rl', 'lw', 'mb', 'ln', 'sar', 'bg'};
useSplit1line = 1;
session = 1;
experimentNumber = 1;
params.theta_stim = [49 54];
params.n_discreteStep = 100;

for ll = 1 : length(subjectID)
    initialValue = initialValueAll;
    
    %% Extract the data
    [estimate1, estimate2] = getData(subjectID{ll}, session, experimentNumber, useSplit1line);

    % Randomize estimates 1 and 2 in condition 1 according to a predetermined order
    file_name = ['ind_pair_cond1_' subjectID{ll}];
    load(file_name)
    estimate1(1, :) = estimate1(1, ind1);
    estimate2(1, :) = estimate2(1, ind2);
    params.estimate1 = estimate1;
    params.estimate2 = estimate2;
    
    %% Start the fitting
    currentFolder = pwd;
    fileNameRoot = 'FitResult';
    fileNumber = GetNextDataFileNumber(currentFolder, '.txt');
    fileID = fopen([fileNameRoot '-' num2str(fileNumber) '.txt'],'w');
    if isempty(subjectID{ll})
        subjectID{ll} = 'Average';
    end
    subjectName = ['Subject: ' subjectID{ll}];
    fprintf(fileID, '%27s \r\n', subjectName); 
    fprintf(fileID,'%11s %11s %11s %8s %8s %8s %8s %8s \r\n', '//-LLHTotal', '-LLH1', '-LLH2', ...
                    'stdSensory', 'stdMemory_1line', 'stdMemory_2line', 'std_motor', 'const_bias');
    for kk = 1 : nLoops
        [tempFitParameter, tempNegLLH] = modelFitBayes_Cond1and2(SetStartPoint, initialValue, share_memory, include_motor_noise, include_const_bias, params, fileID);
        fitParameterAll{kk} = tempFitParameter;
        negLLH(kk) = tempNegLLH;

        % Save the result
        fitParameter = fitParameterAll{kk};
        fprintf(fileID, ['//Iteration' '-' num2str(kk) ' \r\n']);
        fprintf(fileID,'%9.1f %9.5f %9.5f %9.5f %9.5f  %9.5f \r\n', negLLH(kk),...
            fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5));        
    end

    % Print the best params
    negLLH = negLLH(:);
    [negLLHSort,indSort] = sort(negLLH, 'ascend')
    fitParameterAll = fitParameterAll(:);
    fitParameter = fitParameterAll{indSort(1)};
    fprintf(fileID, ['//Best params' ' \r\n']);
    fprintf(fileID,'%2s %9.1f %9.5f %9.5f %9.5f %9.5f %9.5f \r\n','//', negLLH(indSort(1)),...
        fitParameter(1), fitParameter(2), fitParameter(3), fitParameter(4), fitParameter(5));
    fclose(fileID);
end
toc