function best_param_all = extract_fit_param(subjectID)
%%%%%%% Read the fit result file and extract the best parameters
best_param_all = cell(1, length(subjectID));

for ii = 1 : 16
    % Read the text file line by line and put in cell array myFile - each cell contains 1 line
    fileName = ['FitResult-' num2str(ii) '.txt'];
    fileID = fopen(fileName);
    myFile = textscan(fileID,'%s','delimiter','\n');
    myFile = myFile{1};
    
    % Extract subject's ID
    subject_name = myFile{1};
    subject_name([1:9 end]) = [];

    % Match the subject name to subjectID if any
    ind_match = find(ismember(subjectID, subject_name));
    
    % Extract the best parameters
    if ~isempty(ind_match)
        best_param = str2num(myFile{end});
        best_param_all{ind_match} = best_param(2:end);
    end
end