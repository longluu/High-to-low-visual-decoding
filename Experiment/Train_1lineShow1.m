function Train_1lineShow1(params, fileLoad)
tStartExp = GetSecs;
nTrialPerCondition = params.nTrialPerCondition;
lineOrientation = params.lineOrientation;
% setenv('DRI_PRIME','1');
Screen('Preference', 'SkipSyncTests', 0); 
Screen('Preference', 'VisualDebugLevel', 4);
Screen('Preference', 'Verbosity', 4);

% Create the data array.  
%  Column 1: orientation line 1
%  Column 2: orientation line 2
%  Column 3: subject's estimate line 1
%  Column 4: subject's estimate line 2
%  Column 5: subject's reaction time line 1
%  Column 6: subject's reaction time line 2
%  Column 7: report left line (1) or right line (2) 
nTrials = nTrialPerCondition * length(lineOrientation);
if isempty(fileLoad)
    dataResponse = NaN(nTrials, 7);
    index = 1;
    for ii = 1:length(lineOrientation)
        for jj = 1:nTrialPerCondition
            dataResponse(index,1) = lineOrientation(ii); 
            dataResponse(index,2) = lineOrientation(mod(ii, 2)+1); 
            dataResponse(index,7) = mod(jj, 2) + 1;
            index = index + 1;
        end
    end
else
    dataFile = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)], fileLoad);
    load(dataFile);
end


% Disable input to Matlab
ListenChar(2);

% Unify key names
KbName('UnifyKeyNames');

% Set up the screen and window
screens = Screen('Screens');
screenNumber = max(screens);
screen.BgColor = params.backgroundRGB;
[windowDisplay, windowRect] = Screen('OpenWindow', screenNumber, screen.BgColor);   
Screen('BlendFunction', windowDisplay, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenCenter(1), screenCenter(2)] = RectCenter(windowRect); 

% Get the size of the on screen window in pixels
[ScreenWidthPixel, ScreenHeightPixel] = Screen('WindowSize', windowDisplay);

% Retreive the maximum priority number
topPriorityLevel = MaxPriority(windowDisplay);
Priority(topPriorityLevel);

try
    %% Perform calibration
    % Set the parameters
    textFontSize = params.sizeTextPrompt;    
    textColor = params.textColor;                    % color for text and fixation
    ScreenWidthMm = params.ScreenWidthMm;                        % **the width of screen, unit mm**
    DisEye2Scr = params.DisEye2Scr;                         % **the distance from eye to screen mm*
    fixation.clr = params.fixation(1:3);                          % **color of fixation**
    fixation.sizePixel = params.fixation(end);                       % **fixation size in pixel**
    screenCenter(1) = screenCenter(1) + degToPixel(params.centerPointDeg(1), DisEye2Scr, ScreenWidthMm, ScreenWidthPixel);
    screenCenter(2) = screenCenter(2) - degToPixel(params.centerPointDeg(2), DisEye2Scr, ScreenWidthMm, ScreenWidthPixel);    
    halfWidthCaliWindowDeg = params.saccadeWindowDeg(1)/2;
    halfHeightCaliWindowDeg = params.saccadeWindowDeg(2)/2;
    calibration.PointDeg = [-halfWidthCaliWindowDeg -halfHeightCaliWindowDeg;
                        0                       -halfHeightCaliWindowDeg;
                        halfWidthCaliWindowDeg  -halfHeightCaliWindowDeg; 
                        -halfWidthCaliWindowDeg 0;
                        0                       0;
                        halfWidthCaliWindowDeg  0;
                        -halfWidthCaliWindowDeg halfHeightCaliWindowDeg;
                        0                       halfHeightCaliWindowDeg;
                        halfWidthCaliWindowDeg  halfHeightCaliWindowDeg];                        
    xCaliWindowPixel = degToPixel(calibration.PointDeg(:, 1), DisEye2Scr, ScreenWidthMm, ScreenWidthPixel) + screenCenter(1);
    yCaliWindowPixel = degToPixel(calibration.PointDeg(:, 2), DisEye2Scr, ScreenWidthMm, ScreenWidthPixel) + screenCenter(2); 
    calibration.PointPixel = [xCaliWindowPixel yCaliWindowPixel];
    calibration.timeCaliPoint = params.timeCaliPoint;
    calibration.timeFixPoint = params.timeFixPoint;

    % Prompt text
    HideCursor
    Screen('TextSize', windowDisplay, textFontSize);
    DrawFormattedText(windowDisplay, 'Left click the mouse to start the calibration', 'center', screenCenter(2), textColor);
    Screen(windowDisplay,'Flip');
    GetClicks([], 0);

    % Initialize the eye tracker
    WiseData('initialize');           
    WaitSecs(params.timeInitialize);    

    restartCalibration = 1;
    while restartCalibration
        % Display calibration target and record eye data
        % Raw eye data include:
        % rawEyeDataEyePos.L_eyeXs and Data.EyePos.L_eyeYs: x and y position of pupil(left eye)
        % rawEyeDataEyePos.R_eyeXs and Data.EyePos.R_eyeYs: x and y position of pupil(right eye)
        % rawEyeDataEyePos.L_CRXs and Data.EyePos.L_CRYs: x and y position of pupil(left eye)
        % rawEyeDataEyePos.R_CRXs and Data.EyePos.R_CRYs: x and y position of pupil(right eye)
        % rawEyeDataEyePos.R_eyeSizes: the pupil size of right eye
        % rawEyeDataEyePos.L_eyeSizes: the pupil size of left eye
        % rawEyeDataEyePos.eye_SNs: the sample number of eye data
        positionCalibrationPoint = randperm(9); % numbering scheme [1 2 3; 4 5 6; 7 8 9] 
        rawEyeData = DrawCaliSample(screenCenter, fixation, calibration, windowDisplay, textColor, positionCalibrationPoint, params.sizeTextPrompt);

        % Extract the average eye position
%         leftEye_X = NaN(9, 1);
%         leftEye_Y = NaN(9, 1);
        rightEye_X = NaN(9, 1);
        rightEye_Y = NaN(9, 1);
        rangeExtract = 600:950;
        for ii = 1 : 9
%             % Left eye
%             leftPupil_X = rawEyeData(ii).EyePos.L_eyeXs;
%             leftPupil_Y = rawEyeData(ii).EyePos.L_eyeYs;
%             leftCR_X = rawEyeData(ii).EyePos.L_CRXs;
%             leftCR_Y = rawEyeData(ii).EyePos.L_CRYs;
%             leftEye_X(positionCalibrationPoint(ii)) = mean(leftPupil_X(rangeExtract)) - mean(leftCR_X(rangeExtract));
%             leftEye_Y(positionCalibrationPoint(ii)) = mean(leftPupil_Y(rangeExtract)) - mean(leftCR_Y(rangeExtract));

            % Right eye
            rightPupil_X = rawEyeData(ii).EyePos.R_eyeXs;
            rightPupil_Y = rawEyeData(ii).EyePos.R_eyeYs;
            rightCR_X = rawEyeData(ii).EyePos.R_CRXs;
            rightCR_Y = rawEyeData(ii).EyePos.R_CRYs;
            rightEye_X(positionCalibrationPoint(ii)) = mean(rightPupil_X(rangeExtract)) - mean(rightCR_X(rangeExtract));
            rightEye_Y(positionCalibrationPoint(ii)) = mean(rightPupil_Y(rangeExtract)) - mean(rightCR_Y(rangeExtract));    
        end

        % Compute the coefficients of the mapping function between eye position and the eye gaze on monitor
        % x = a1 + a2 * X + a3 * X^2 + a4 * Y + a5 * Y^2 + a6 * X * Y
        % y = b1 + b2 * X + b3 * X^2 + b4 * Y + b5 * Y^2 + b6 * X * Y
        % x, y: the gaze position on screen
        % X, Y: the eye position (pupil-CR)
        % We first find a and b by solving E*a = m in which E = [1 X1 X1^2 Y1 Y1^2 X1*Y1; ...; 1 X9 X9^2 Y9 Y9^2 X9*Y9],
        % a = [a1 ... a5 a6], m = [x1 ... x5 x6] or [y1 ... y5 y6]
        caliPointUse = 1:9;
%         E_left = [ones(length(caliPointUse), 1) leftEye_X(caliPointUse) leftEye_X(caliPointUse).^2 ...
%                     leftEye_Y(caliPointUse) leftEye_Y(caliPointUse).^2 leftEye_X(caliPointUse).*leftEye_Y(caliPointUse)];
        E_right = [ones(length(caliPointUse), 1) rightEye_X(caliPointUse) rightEye_X(caliPointUse).^2 ...
                    rightEye_Y(caliPointUse) rightEye_Y(caliPointUse).^2 rightEye_X(caliPointUse).*rightEye_Y(caliPointUse)];

%         a_left_x = E_left\calibration.PointDeg(caliPointUse, 1);
%         a_left_y = E_left\calibration.PointDeg(caliPointUse, 2);
        a_right_x = E_right\calibration.PointDeg(caliPointUse, 1);
        a_right_y = E_right\calibration.PointDeg(caliPointUse, 2);

%         mse_left_x = sum(abs(E_left * a_left_x - calibration.PointDeg(caliPointUse, 1))) / length(caliPointUse);
%         mse_left_y = sum(abs(E_left * a_left_y - calibration.PointDeg(caliPointUse, 2))) / length(caliPointUse);
        mse_right_x = sum(abs(E_right * a_right_x - calibration.PointDeg(caliPointUse, 1))) / length(caliPointUse);
        mse_right_y = sum(abs(E_right * a_right_y - calibration.PointDeg(caliPointUse, 2))) / length(caliPointUse);

        fprintf('The mean absolute error of calibration (right-x, right-y) is %6.2f %6.2f \n', mse_right_x, mse_right_y)

        % Save data for calibration
        calibrationData.positionCalibrationPoint = positionCalibrationPoint; 
        calibrationData.rawEyeData = rawEyeData; 
        calibrationData.a_right_x = a_right_x; 
        calibrationData.a_right_y = a_right_y;
        
        % Prompt text
        textDisplay = ['The mean absolute error (right-x, right-y) is ' num2str(round(mse_right_x*100)/100) '  ' num2str(round(mse_right_y*100)/100)];        
        Screen('TextSize', windowDisplay, textFontSize-20);
        DrawFormattedText(windowDisplay, textDisplay, 'center', screenCenter(2)-100, textColor);
        DrawFormattedText(windowDisplay, 'Left click for validation, Right click to repeat calibration', 'center', screenCenter(2)+100, textColor);
        Screen(windowDisplay,'Flip');
        [~,~,~,whichButton] = GetClicks([], 0);
        if whichButton == 1
            restartCalibration = 0;
        elseif whichButton == 2
            % Quit the experiment and save everything
            error('abortCalibration');               
        end
    end

    %% Perform validation
    restartValidation = 1;
    while restartValidation
        % Choose n random validation points
        validation.nPointValidation = params.nPointValidation;
        validation.xValidation = round(rand(1,validation.nPointValidation)  * ...
                    (max(xCaliWindowPixel) - min(xCaliWindowPixel)) + min(xCaliWindowPixel));
        validation.yValidation = round(rand(1, validation.nPointValidation) * ...
                    (max(yCaliWindowPixel) - min(yCaliWindowPixel)) + min(yCaliWindowPixel));

        % Duration of one validation point
        validation.timeValiPoint = params.timeValiPoint; 

        % Draw the validation samples
        valiData = DrawValidationSample(screenCenter, fixation, validation, windowDisplay, textColor, textFontSize);

        % Extract the average eye position
%         leftEye_X = NaN(validation.nPointValidation, 1);
%         leftEye_Y = NaN(validation.nPointValidation, 1);
        rightEye_X = NaN(validation.nPointValidation, 1);
        rightEye_Y = NaN(validation.nPointValidation, 1);
        rangeExtract = 400:800;
        for ii = 1 : validation.nPointValidation
%             % Left eye
%             leftPupil_X = valiData(ii).EyePos.L_eyeXs;
%             leftPupil_Y = valiData(ii).EyePos.L_eyeYs;
%             leftCR_X = valiData(ii).EyePos.L_CRXs;
%             leftCR_Y = valiData(ii).EyePos.L_CRYs;
%             leftEye_X(ii) = mean(leftPupil_X(rangeExtract)) - mean(leftCR_X(rangeExtract));
%             leftEye_Y(ii) = mean(leftPupil_Y(rangeExtract)) - mean(leftCR_Y(rangeExtract));

            % Right eye
            rightPupil_X = valiData(ii).EyePos.R_eyeXs;
            rightPupil_Y = valiData(ii).EyePos.R_eyeYs;
            rightCR_X = valiData(ii).EyePos.R_CRXs;
            rightCR_Y = valiData(ii).EyePos.R_CRYs;
            rightEye_X(ii) = mean(rightPupil_X(rangeExtract)) - mean(rightCR_X(rangeExtract));
            rightEye_Y(ii) = mean(rightPupil_Y(rangeExtract)) - mean(rightCR_Y(rangeExtract));    
        end    

        % Compute the gaze on monitor from calibration mapping
%         E_left = [ones(validation.nPointValidation, 1) leftEye_X leftEye_X.^2 ...
%                     leftEye_Y leftEye_Y.^2 leftEye_X.*leftEye_Y];
        E_right = [ones(validation.nPointValidation, 1) rightEye_X rightEye_X.^2 ...
                    rightEye_Y rightEye_Y.^2 rightEye_X.*rightEye_Y];

        xValidationDeg = pixelToDeg(validation.xValidation, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel, screenCenter, 1);
        yValidationDeg = pixelToDeg(validation.yValidation, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel, screenCenter, 0);
%         absoluteError_xvali_left = abs(E_left * a_left_x - xValidationDeg');
%         absoluteError_yvali_left = abs(E_left * a_left_y - yValidationDeg');
        absoluteError_xvali_right = abs(E_right * a_right_x - xValidationDeg');
        absoluteError_yvali_right = abs(E_right * a_right_y - yValidationDeg');
%         mae_left_x_vali = sum(absoluteError_xvali_left) / validation.nPointValidation;
%         mae_left_y_vali = sum(absoluteError_yvali_left) / validation.nPointValidation;
        mae_right_x_vali = sum(absoluteError_xvali_right) / validation.nPointValidation;
        mae_right_y_vali = sum(absoluteError_yvali_right) / validation.nPointValidation;

        fprintf('The mean absolute error of validation (right-x, right-y) is %6.2f %6.2f \n', ...
                    mae_right_x_vali, mae_right_y_vali)
                
        % Save data for validation
        validationData.positionCalibrationPoint = [validation.xValidation; validation.yValidation]; 
        validationData.rawEyeData = valiData; 
        
        % Prompt text
        textDisplay = ['The mean absolute error (right-x, right-y) is ' num2str(round(mae_right_x_vali*100)/100) '  ' num2str(round(mae_right_y_vali*100)/100)];        
        Screen('TextSize', windowDisplay, textFontSize-20);
        DrawFormattedText(windowDisplay, textDisplay, 'center', screenCenter(2)-100, textColor);
        DrawFormattedText(windowDisplay, 'Left click for main experiment, Right click to repeat validation', 'center', screenCenter(2)+100, textColor);
        Screen(windowDisplay,'Flip');
        [~,~,~,whichButton] = GetClicks([], 0);
        if whichButton == 1
            restartValidation = 0;
        elseif whichButton == 2
            % Quit the experiment and save everything
            error('abortCalibration');               
        end
    end
    
    % Dummy code to run all the necessary functions
    GetSecs;
    KbCheck;
    WaitSecs(1);

    %% Main experiment
    % Set up the beep sound
    % Initialize Sounddriver
    InitializePsychSound(1);

    % Number of channels and Frequency of the sound
    nrchannels = 2;
    freqError = 200;
    freqResponseLeft = 400;
    freqResponseRight = 400;
    
    % How many times to we wish to play the sound
    repetitions = 1;

    % Length of the beep
    errorBeepLengthSecs = 0.5;
    responseBeepLengthSecs = 0.2;
    
    % Start immediately (0 = immediately)
    startCue = 0;

    % Should we wait for the device to really start (1 = yes)
    % INFO: See help PsychPortAudio
    waitForDeviceStart = 0;

    % Open Psych-Audio port, with the follow arguements
    % (1) [] = default sound device
    % (2) 1 = sound playback only
    % (3) 1 = default level of latency
    % (4) Requested frequency in samples per second
    % (5) 2 = stereo putput
    pahandleError = PsychPortAudio('Open', [], 1, 1, 48000, nrchannels);
    pahandleResponse = PsychPortAudio('Open', [], 1, 1, 48000, nrchannels);
    
    % Set the volume to half for this demo
    PsychPortAudio('Volume', pahandleError, 0.1);
    PsychPortAudio('Volume', pahandleResponse, 0.05);
    
    % Make a beep which we will play back to the user
    errorBeep = MakeBeep(freqError, errorBeepLengthSecs);
    responseLeftBeep = MakeBeep(freqResponseLeft, responseBeepLengthSecs);
    responseRightBeep = MakeBeep(freqResponseRight, responseBeepLengthSecs);
    
    % Convert everything to pixel
    lineDeviation = degToPixel(params.lineLocationDeg, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel);
    lineLength = degToPixel(params.lineLengthDeg, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel);
    coordinateFirstOrientation = [-lineLength*cosd(params.lineOrientation(1))/2   lineLength*cosd(params.lineOrientation(1))/2;...
                                   lineLength*sind(params.lineOrientation(1))/2  -lineLength*sind(params.lineOrientation(1))/2];
    coordinateSecondOrientation = [-lineLength*cosd(params.lineOrientation(2))/2   lineLength*cosd(params.lineOrientation(2))/2;...
                                   lineLength*sind(params.lineOrientation(2))/2  -lineLength*sind(params.lineOrientation(2))/2];
    distanceMouseInitiation = degToPixel(params.distanceMouseInitiationDeg, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel);
    pointerLengthPixel = degToPixel(params.pointerLengthDeg, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel);
    rBigCircle = min(screenCenter);
    if isempty(fileLoad)
        % Set eye data parameter
        eyeDataMain.eye_SNs = cell(1, nTrials);                     
        eyeDataMain.R_eyeXs = cell(1, nTrials);   
        eyeDataMain.R_eyeYs = cell(1, nTrials);  
        eyeDataMain.R_eyeSizes = cell(1, nTrials);  
        eyeDataMain.R_CRXs = cell(1, nTrials);  
        eyeDataMain.R_CRYs = cell(1, nTrials);  
        eyeDataMain.indicatorStimulus = cell(1, nTrials); 
        
        % New session
        trialOrder = randperm(nTrials);  
        params.trialOrder = trialOrder;
        trialIndex = 1;        
        
        % Prompt screen
        Screen('FillRect', windowDisplay, params.backgroundRGB);
        Screen('TextSize', windowDisplay, params.sizeTextPrompt);
        DrawFormattedText(windowDisplay, 'Left click to start the main experiment', 'center', screenCenter(2), params.textColor);
        Screen('Flip', windowDisplay);  
        GetClicks([], 0);
        WaitSecs(1)
            
        while trialIndex <= nTrials 
            repeatFlag = 0;
            
            % Get the orientation of the 2 lines and convert to coordinate in pixel        
            angleFirstLine = dataResponse(trialOrder(trialIndex), 1);
            if angleFirstLine == params.lineOrientation(1)
                coordinateFirst = coordinateFirstOrientation;
                coordinateSecond = coordinateSecondOrientation;
            else
                coordinateFirst = coordinateSecondOrientation;     
                coordinateSecond = coordinateFirstOrientation;         
            end

            % Get index what line to report
            whichLine = dataResponse(trialOrder(trialIndex), 7);
            
            % Fixation point
            Screen('DrawDots', windowDisplay,[screenCenter(1) screenCenter(2)], params.fixation(end), params.fixation(1:3), [], 2);            
            Screen('Flip', windowDisplay);
            WaitSecs(params.fixationDuration);
            if params.instruction
                pause
            end
            tStartFlip = GetSecs;

            % Set parameters
            tTrial = params.timeTrial;
            numDataPoint = tTrial * 1000 + 3000;  
            eye_SNs = NaN(1, numDataPoint);                     
            R_eyeXs = NaN(1, numDataPoint); 
            R_eyeYs = NaN(1, numDataPoint); 
            R_eyeSizes = NaN(1, numDataPoint); 
            R_CRXs = NaN(1, numDataPoint); 
            R_CRYs = NaN(1, numDataPoint);
            indicatorStimulus = NaN(1, numDataPoint);
            eyeVelocity = NaN(2, numDataPoint);
            eyeAcceleration = NaN(1, numDataPoint);            
            eyepos_SN = 1; 
            counterFixationInitial = 0;
            counterFixationOutsideWindow = 0;
            keepLooping = 1;
            stimNotShown = 1;
            isMouseInitiated = 0;
            isMouseSetToBigCircle = 0;
            isStartEstimate = 1;
            
            % Collect 4 data samples first (neccesary to compute the velocity)
            real_data = WiseData('get');
            last_SN = real_data(3,1); 
            while keepLooping
                real_data = WiseData('get');
                dif_SN = real_data(3,1) - last_SN;
                if dif_SN > 0
                    R_eyeXs(eyepos_SN) = real_data(1 , end);
                    R_eyeYs(eyepos_SN) = real_data( 2 , end);
                    R_eyeSizes(eyepos_SN) = real_data(3, end);
                    R_CRXs(eyepos_SN) = real_data( 5 , end);
                    R_CRYs(eyepos_SN) = real_data( 6 , end);
                    eye_SNs(eyepos_SN) = real_data(3, 1);

                    eyepos_SN = eyepos_SN + 1;
                    last_SN = real_data(3,1);
                end

                if eyepos_SN == 5 
                    break
                end
            end
            
            while keepLooping                                
                % Get eye data
                real_data = WiseData('get');
                dif_SN = real_data(3,1) - last_SN;
                if  dif_SN > 0
                    if dif_SN > 20          
                        dif_SN = 20;
                    end

                    startIndex = (41 - (2 * dif_SN - 1))/2;                    
                    for CurrentSN = startIndex : startIndex+dif_SN-1
                        % Collect the new samples, start from first to last samples
                        Col_Index = 2 * CurrentSN;    
                        R_eyeXs(eyepos_SN) = real_data( 1 , Col_Index + 1 );
                        R_eyeYs(eyepos_SN) = real_data( 2 , Col_Index + 1 );
                        R_eyeSizes(eyepos_SN) = real_data(3,Col_Index + 1);
                        R_CRXs(eyepos_SN) = real_data( 5 , Col_Index + 1 );
                        R_CRYs(eyepos_SN) = real_data( 6 , Col_Index + 1 );
                        eye_SNs(eyepos_SN) = real_data(3,1) - dif_SN + CurrentSN;

                        % Take the last 5 samples of eye position and convert to degree
                        rawEyeSignal_X = R_eyeXs(eyepos_SN-4:eyepos_SN) - R_CRXs(eyepos_SN-4:eyepos_SN);
                        rawEyeSignal_Y = R_eyeYs(eyepos_SN-4:eyepos_SN) - R_CRYs(eyepos_SN-4:eyepos_SN);
                        E_right = [ones(length(rawEyeSignal_X), 1) rawEyeSignal_X' rawEyeSignal_X'.^2 ...
                                    rawEyeSignal_Y' rawEyeSignal_Y'.^2 rawEyeSignal_X'.*rawEyeSignal_Y'];
                        eyeGaze_X_Y = E_right * [a_right_x a_right_y];

                        % Check if the eye gaze is within fixation window
                        distanceEyeToTarget = sqrt((eyeGaze_X_Y(5, 1))^2 + (eyeGaze_X_Y(5, 2))^2);                                    
                        if counterFixationInitial < params.nFixationSampleInitial
                            if distanceEyeToTarget < params.fixationWindowRadiusDeg
                               counterFixationInitial = counterFixationInitial+1;
                            else
                                counterFixationInitial = 0;
                            end
                        end
                        
%                         % Compute the velocity
%                         eyeVelocity(:, eyepos_SN) = velocityComputeOnline(eyeGaze_X_Y);
% 
%                         % Compute the acceleration
%                         eyeAcceleration(eyepos_SN) = velocityComputeOnline(eyeVelocity(2, eyepos_SN-4:eyepos_SN)');
                        
                        % If subject fixation is stable for some time,
                        % present the stimulus
                        if counterFixationInitial == params.nFixationSampleInitial
                            % Present the stimulus
                            if stimNotShown == 1
                                Screen('FillRect', windowDisplay, params.backgroundRGB);
                                Screen('DrawDots', windowDisplay,[screenCenter(1) screenCenter(2)], params.fixation(end), params.fixation(1:3), [], 2);   
                                if whichLine == 1
                                    Screen('DrawLines', windowDisplay, coordinateFirst,...
                                        params.lineWidthPixel, params.colorPointer(1, :), [screenCenter(1)-lineDeviation screenCenter(2)], 2);
                                else
                                    Screen('DrawLines', windowDisplay, coordinateSecond,...
                                        params.lineWidthPixel, params.colorPointer(2, :), [screenCenter(1)+lineDeviation screenCenter(2)], 2);
                                end
                                Screen('Flip', windowDisplay, [], [], 1);
                                
                                % Start timestamp
                                tStartStim = GetSecs;
                                indicatorStimulus(eyepos_SN) = 1;
                                stimNotShown = 0;
                            elseif GetSecs - tStartStim < params.lineDuration
                                if distanceEyeToTarget > params.fixationWindowRadiusDeg
                                   counterFixationOutsideWindow = counterFixationOutsideWindow+1;
                                else
                                    counterFixationOutsideWindow = 0;
                                end
                                if counterFixationOutsideWindow > params.nMaxGazeDeviation
                                    keepLooping = 0;
                                    repeatFlag = 1;
                                    
                                    % Warning beep
                                    PsychPortAudio('FillBuffer', pahandleError, [errorBeep; errorBeep]);
                                    PsychPortAudio('Start', pahandleError, repetitions, startCue, waitForDeviceStart);
                                    WaitSecs(params.fixationDuration);
                                end
                            elseif GetSecs - tStartStim > params.lineDuration 
                                % Time stamp for reaction time (first line)
                                if isStartEstimate == 1
                                    tStartFirstEstimate = GetSecs;
                                    indicatorStimulus(eyepos_SN) = 1;
                                    isStartEstimate = 0;
                                    Screen('FillRect', windowDisplay, params.backgroundRGB);
                                    Screen('DrawDots', windowDisplay,[screenCenter(1) screenCenter(2)], params.fixation(end), params.fixation(1:3), [], 2);   
                                    Screen('Flip', windowDisplay);
                                    
                                    % Prompt subject to estimate the first line
                                    if whichLine == 1
                                        PsychPortAudio('FillBuffer', pahandleResponse, [responseLeftBeep; responseLeftBeep]);
                                    else
                                        PsychPortAudio('FillBuffer', pahandleResponse, [responseRightBeep; responseRightBeep]);
                                    end
                                    PsychPortAudio('Start', pahandleResponse, repetitions, startCue, waitForDeviceStart);
                                    
                                    % Set the mouse to center
                                    SetMouse(screenCenter(1), screenCenter(2), windowDisplay);
                                end
                                
                                % Get the current position of the mouse
                                [xMouse, yMouse, buttonMouse] = GetMouse(windowDisplay);
                                distanceMouse = sqrt((xMouse - screenCenter(1))^2 + (yMouse - screenCenter(2))^2);
                                if distanceMouse > distanceMouseInitiation
                                    isMouseInitiated = 1;
                                end
                                
                                % Display the pointer
                                if GetSecs - tStartFlip > params.interPointerDuration
                                    if isMouseInitiated
                                        % Compute the angle of the line connecting the origin and mouse pointer
                                        thetaMouse = atan((screenCenter(2) - yMouse) / (xMouse - screenCenter(1)));

                                        if isMouseSetToBigCircle == 0
                                            % Compute the x-y coordinate of intersection of extended line and big circle
                                            if thetaMouse > 0
                                                xExtend = screenCenter(1) + rBigCircle * cos(thetaMouse);
                                            else
                                                xExtend = screenCenter(1) - rBigCircle * cos(thetaMouse);
                                            end
                                            yExtend = screenCenter(2) - rBigCircle * sin(abs(thetaMouse));

                                            % Set the mouse to the top of screen
                                            SetMouse(xExtend, yExtend, windowDisplay);

                                            % Only set the mouse to top once
                                            isMouseSetToBigCircle = 1;
                                        else                                    
                                            % Compute the pointer coordinates
                                            xDisplay = screenCenter(1) + cos(thetaMouse) * pointerLengthPixel/2;
                                            yDisplay = screenCenter(2) - sin(thetaMouse) * pointerLengthPixel/2;

                                            % Display the pointer
                                            if whichLine == 1
                                                % First line
                                                Screen('DrawLines', windowDisplay, [2*screenCenter(1)-xDisplay xDisplay; 2*screenCenter(2)-yDisplay yDisplay],...
                                                    params.lineWidthPixel, params.colorPointer(1, :), [], 2);
                                            else
                                                % Second line
                                                Screen('DrawLines', windowDisplay, [2*screenCenter(1)-xDisplay xDisplay; 2*screenCenter(2)-yDisplay yDisplay],...
                                                    params.lineWidthPixel, params.colorPointer(2, :), [], 2);
                                            end                                                
                                        end
                                    end
                                    
                                    % Draw fixation dot on top of the
                                    % pointer
                                    Screen('DrawDots', windowDisplay,[screenCenter(1) screenCenter(2)], params.fixation(end), params.fixation(1:3), [], 2);   

                                    % Get the delay
                                    Screen('Flip', windowDisplay, 1);
                                    tStartFlip = GetSecs;
                                end
                                
                                % Check if subject click mouse to stop, start the second estimate or end the trial
                                if any(buttonMouse) && exist('thetaMouse', 'var') && isMouseInitiated
                                    % Save subject's response
                                    if whichLine == 1
                                        dataResponse(trialOrder(trialIndex),3) = rad2deg(thetaMouse);
                                        dataResponse(trialOrder(trialIndex),5) = GetSecs-tStartFirstEstimate; 
                                    else
                                        dataResponse(trialOrder(trialIndex),4) = rad2deg(thetaMouse);
                                        dataResponse(trialOrder(trialIndex),6) = GetSecs-tStartFirstEstimate; 
                                    end
                                    
                                    % Wait for the click release
                                    keepLoopMouse = 1;
                                    while keepLoopMouse
                                        [~, ~, buttonMouse] = GetMouse(windowDisplay);
                                        if sum(buttonMouse) == 0
                                            keepLoopMouse = 0;
                                        end
                                    end

                                    % Stop the trial
                                    keepLooping = false;
                                end
                                

                            end
                        end                        
                    end
                    eyepos_SN = eyepos_SN + 1;  
                    last_SN = real_data(3,1);
                end

                % Check quit key
                [keyIsDown, ~, keyCode, ~] = KbCheck;
                if keyIsDown &&  sum(find(keyCode)==KbName('q'))
                    % Wait until key released
                    KbReleaseWait;

                    % Quit the experiment and save everything
                    error('abort');                                
                end                                        
            end
            
            if repeatFlag == 1
            else   
                % Record eye data
                eyeDataMain.eye_SNs{trialOrder(trialIndex)} = eye_SNs;                     
                eyeDataMain.R_eyeXs{trialOrder(trialIndex)} = R_eyeXs;   
                eyeDataMain.R_eyeYs{trialOrder(trialIndex)} = R_eyeYs;  
                eyeDataMain.R_eyeSizes{trialOrder(trialIndex)} = R_eyeSizes;  
                eyeDataMain.R_CRXs{trialOrder(trialIndex)} = R_CRXs;  
                eyeDataMain.R_CRYs{trialOrder(trialIndex)} = R_CRYs;                          
                eyeDataMain.indicatorStimulus{trialOrder(trialIndex)} = indicatorStimulus;   
                
                % Increment the trial counter
                trialIndex = trialIndex + 1;
                params.currentTrialIndex = trialIndex;
                
                % Show correct line and subject's response
                Screen('FillRect', windowDisplay, params.backgroundRGB);
                if whichLine == 1
                    Screen('DrawLines', windowDisplay, [2*screenCenter(1)-xDisplay xDisplay; 2*screenCenter(2)-yDisplay yDisplay],...
                        params.lineWidthPixel, params.colorPointer(1, :), [], 2);                
                    Screen('DrawLines', windowDisplay, coordinateFirst,...
                        params.lineWidthPixel, params.feedbackPointer, [screenCenter(1) screenCenter(2)], 2);
                else
                    Screen('DrawLines', windowDisplay, [2*screenCenter(1)-xDisplay xDisplay; 2*screenCenter(2)-yDisplay yDisplay],...
                        params.lineWidthPixel, params.colorPointer(2, :), [], 2);                                
                    Screen('DrawLines', windowDisplay, coordinateSecond,...
                        params.lineWidthPixel, params.feedbackPointer, [screenCenter(1) screenCenter(2)], 2);
                end
                Screen('DrawDots', windowDisplay,[screenCenter(1) screenCenter(2)], params.fixation(end), params.fixation(1:3), [], 2);   
                Screen('Flip', windowDisplay, 1);
                GetClicks([], 0);                 
            end     
            
            % Blank screen
            Screen('FillRect', windowDisplay, params.backgroundRGB);
            Screen('Flip', windowDisplay);
            timeWait = rand * diff(params.intertrialInterval) + params.intertrialInterval(1);
            WaitSecs(timeWait);
        end
    else
    end  
    
    % Save the data
    params.timeExp = GetSecs - tStartExp;
    if isempty(fileLoad)
        dataFolder = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)]);
        if ~exist(dataFolder,'dir')
            mkdir(dataFolder)
        end
        dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
    end
    save(dataFile,'dataResponse','params', 'eyeDataMain', 'calibrationData', 'validationData')
    
    % Prompt text
    Screen('FillRect', windowDisplay, params.backgroundRGB);
    Screen('TextSize', windowDisplay, params.sizeTextPrompt);
    DrawFormattedText(windowDisplay, 'End of block 1. Thank you for your participation!', 'center', screenCenter(2)-params.yPositionPromptText, params.textColor);
    Screen('Flip', windowDisplay);  
    KbStrokeWait   
    
    % Close the window
    sca
    ListenChar(0);
    Priority(0);
catch e  
    % Save the data
    params.timeExp = GetSecs - tStartExp;
    if isempty(fileLoad)
        dataFolder = fullfile('Data', params.subject, 'MainExperiment', [params.experimentName num2str(params.session)]);
        if ~exist(dataFolder,'dir')
            mkdir(dataFolder)
        end
        dataFile = sprintf('%s/%s-%d.mat', dataFolder, params.experimentName, GetNextDataFileNumber(dataFolder, '.mat'));    
    end
    if strcmp(e.message, 'abortCalibration')
        save(dataFile,'dataResponse','params', 'calibrationData')
    elseif strcmp(e.message, 'abortValidation')
        save(dataFile,'dataResponse','params', 'calibrationData', 'validationData')
    else
        save(dataFile,'dataResponse','params', 'eyeDataMain', 'calibrationData', 'validationData')
    end

    % Close the window
    Screen('Close')
    sca
    ListenChar(0);
    Priority(0);
    
    % Report error
    if strcmp(e.message, 'abort') || strcmp(e.message, 'abortCalibration') || strcmp(e.message, 'abortValidation')
        fprintf('- Experiment aborted, some data saved.\n');
    else      
        keyboard
    end
    
end
