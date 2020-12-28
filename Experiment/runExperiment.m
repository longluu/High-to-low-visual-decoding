%% Define the parameters in the experiment
% Physical setup
params.ScreenWidthMm = 1210; % **the width of screen, unit mm**
params.ScreenHeightMm = 690;
params.DisEye2Scr = 560;    % **the distance from eye to screen mm*

% Stimulus
params.fixation = [255 255 255 8]; % color, size in pixel
params.fixationDuration = 0;
params.backgroundRGB = [0 0 0];	% RGB of the background.  All values are in the [0,1] range.
params.lineWidthPixel = 3; % in pixel
params.lineRGB = [255 255 255];
params.lineLocationDeg = 8; % deviation from center
params.lineLengthDeg = 6;
params.lineDuration = 1;
params.distanceMouseInitiationDeg = 6; 
params.pointerLengthDeg = 6;
params.lineFeedbackRGB = [0 255 0];
params.interPointerDuration = 0.017;
params.colorPointer = [255 0 255; 0 255 0   ];
params.feedbackPointer = [255 255 255];

% Eye gaze check
params.nFixationSampleInitial = 1000;
params.fixationWindowRadiusDeg = 3;
params.nMaxGazeDeviation = 5;

% Calibration
params.timeInitialize = 1;
params.timeCaliPoint = 1;
params.timeFixPoint = 0.8;
params.centerPointDeg = [0 -10]; 
params.saccadeWindowDeg = [20 10];

% Validation
params.nPointValidation = 10; 
params.timeValiPoint = 1; 

% Experimental condition
params.intertrialInterval = [0.3 0.6]; % time btw consecutive trials (sec)
params.lineOrientation = [49 54];
params.timeTrial = 30;
params.timeExp = [];
params.reportWhichLine = 1; % for 1 line show 2 condition only, 1: left, 2: right
params.tWaitMatch = 4.5; % for 1 line show 2 condition only, to match trial duration of 2 line condition

% Text prompt
params.textColor = [255 255 255];
params.yPositionPromptText = 150; % vertical position from center in pixel
params.sizeTextPrompt = 50;
params.sizeNumber = 80;
params.ratingScale = [1 10];
params.instruction = 0;

% Subject and session
trialCondition = 1; %##################################   
params.subject = 'xw'; %##################################  
params.session = 1; %##################################  
params.delay = [];

%% Run the driver program
fileLoad = '';
if trialCondition == 0.0
    params.nTrialPerCondition = 20;
    params.experimentName = 'Train_1lineShow1';
    Train_1lineShow1(params, fileLoad)  
elseif trialCondition == 0.01
    params.nTrialPerCondition = 20;
    params.experimentName = 'Train_1lineShow1';
    Train_1lineShow1_49(params, fileLoad) 
elseif trialCondition == 0.02
    params.nTrialPerCondition = 20;
    params.experimentName = 'Train_1lineShow1';
    Train_1lineShow1_54(params, fileLoad)     
elseif trialCondition == 0.1
    params.nTrialPerCondition = 50;
    params.experimentName = 'Instruction_1lineShow1';
    Instruction_Sequential_1lineShow1(params, fileLoad)        
elseif trialCondition == 0.2
    params.nTrialPerCondition = 50;
    params.experimentName = 'Instruction_2line';
    Instruction_Sequential_2line(params, fileLoad)  
elseif trialCondition == 0.3
    params.nTrialPerCondition = 50;
    params.experimentName = 'Instruction_2line_separate';
    Instruction_Separate_2line(params, fileLoad)      
elseif trialCondition == 1
    params.nTrialPerCondition = 50;
    params.experimentName = 'HighToLow_1lineShow2_1side';
    Main_1lineShow2_1side(params, fileLoad)   
elseif trialCondition == 1.1
    params.nTrialPerCondition = 50;
    params.experimentName = 'HighToLow_1lineShow1_49';
    Main_1lineShow1_49(params, fileLoad)
elseif trialCondition == 1.2
    params.nTrialPerCondition = 50;
    params.experimentName = 'HighToLow_1lineShow1_54';
    Main_1lineShow1_54(params, fileLoad)
elseif trialCondition == 2
    params.nTrialPerCondition = 25;
    params.experimentName = 'HighToLow';
    Main_sequentialDraw(params, fileLoad)
elseif trialCondition == 3
    params.nTrialPerCondition = 25;
    params.experimentName = 'HighToLow_separate';
    Main_separateDraw(params, fileLoad)    
end    
    
    
    
    
    
 
    