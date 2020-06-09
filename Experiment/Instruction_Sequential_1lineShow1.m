function Instruction_Sequential_1lineShow1(params, fileLoad)
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
HideCursor

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
        
        % Prompt text
        textDisplay = ['The mean absolute error (right-x, right-y) is '  'xxx, xxx'];        
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
        
        % Prompt text
        textDisplay = ['The mean absolute error (right-x, right-y) is '  'xxx, xxx'];        
        Screen('TextSize', windowDisplay, textFontSize-20);
        DrawFormattedText(windowDisplay, textDisplay, 'center', screenCenter(2)-100, textColor);
        DrawFormattedText(windowDisplay, 'Left click for main experiment, Right click to repeat validation', 'center', screenCenter(2)+100, textColor);
        Screen(windowDisplay,'Flip');
        [~,~,~,whichButton] = GetClicks([], 0);
        if whichButton == 1
            restartValidation = 0;
        elseif whichButton == 2
            % Quit the experiment and save everything
            error('abortValidation');               
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
        DrawFormattedText(windowDisplay, 'Left click to start the main experiment', 'center', screenCenter(2)-params.yPositionPromptText, params.textColor);
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
            indicatorStimulus = NaN(1, numDataPoint);
            eyepos_SN = 1; 
            counterFixationInitial = 0;
            counterFixationOutsideWindow = 0;
            keepLooping = 1;
            stimNotShown = 1;
            isMouseInitiated = 0;
            isMouseSetToBigCircle = 0;
            isStartEstimate = 1;
            
            % Set the mouse to center
            SetMouse(screenCenter(1), screenCenter(2), windowDisplay);
                        
            while keepLooping                                
                % Simulate eye gaze with mouse pointer
                [xMouse, yMouse] = GetMouse(windowDisplay);

                % Check if the eye gaze is within fixation window
                distanceEyeToTarget = sqrt((xMouse - screenCenter(1))^2 + (yMouse - screenCenter(2))^2);                                    
                if counterFixationInitial < params.nFixationSampleInitial
                    if distanceEyeToTarget < 100000
                       counterFixationInitial = counterFixationInitial+1;
                    else
                        counterFixationInitial = 0;
                    end
                end


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
                            Screen('FillRect', windowDisplay, params.backgroundRGB);
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
                        if any(buttonMouse)
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
                    eyepos_SN = eyepos_SN + 1;  
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
                % Increment the trial counter
                trialIndex = trialIndex + 1;
                params.currentTrialIndex = trialIndex;
            end
            
            % Blank screen
            Screen('FillRect', windowDisplay, params.backgroundRGB);
            Screen('Flip', windowDisplay);
            timeWait = rand * diff(params.intertrialInterval) + params.intertrialInterval(1);
            WaitSecs(timeWait);
        end

    end  
        
    % Prompt text
    Screen('FillRect', windowDisplay, params.backgroundRGB);
    Screen('TextSize', windowDisplay, params.sizeTextPrompt);
    DrawFormattedText(windowDisplay, 'End of experiment. Thank you for your participation!', 'center', screenCenter(2)-params.yPositionPromptText, params.textColor);
    Screen('Flip', windowDisplay);  
    KbStrokeWait   
    params.timeExp = GetSecs - tStartExp;
    
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
