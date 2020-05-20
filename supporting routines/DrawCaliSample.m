function [Data] = DrawCaliSample(screenCenter, fixation, calibration, windowDisplay, textColor, positionCaliPoint, textFontSize)
    real_data = WiseData('get');
    isOnline = real_data(1,1);
    last_SN = real_data(3,1);
    
    numDataPoint = 9 * (calibration.timeFixPoint + calibration.timeCaliPoint) * 1000 + 100;  
    eye_SNs = NaN(1, numDataPoint);                     
    L_eyeXs = NaN(1, numDataPoint); 
    L_eyeYs = NaN(1, numDataPoint); 
    L_eyeSizes = NaN(1, numDataPoint); 
    R_eyeXs = NaN(1, numDataPoint); 
    R_eyeYs = NaN(1, numDataPoint); 
    R_eyeSizes = NaN(1, numDataPoint); 
    L_CRXs = NaN(1, numDataPoint); 
    L_CRYs = NaN(1, numDataPoint); 
    R_CRXs = NaN(1, numDataPoint); 
    R_CRYs = NaN(1, numDataPoint); 
    timeStamp.time = NaN(1, numDataPoint);
    timeStamp.sampleNumber = NaN(1, numDataPoint);

    if isOnline == 1                      % start trigger(do not change here)  
       for i_trial = 1 : 9                  
            eyepos_SN = 1;
            fixation.Loc_X = screenCenter(1);
            fixation.Loc_Y = screenCenter(2);
            
            % Present fixation dot
            Screen('DrawDots', windowDisplay, [fixation.Loc_X fixation.Loc_Y], fixation.sizePixel, fixation.clr, [], 2);
            Screen(windowDisplay, 'Flip');
            WaitSecs(calibration.timeFixPoint)
            
            % Present the calibration points
            indCaliPoint = positionCaliPoint(i_trial);
            Screen('DrawDots', windowDisplay, [calibration.PointPixel(indCaliPoint, 1) calibration.PointPixel(indCaliPoint, 2)],...
                        fixation.sizePixel, fixation.clr, [], 2);
            Screen(windowDisplay, 'Flip'); 
            
            % Start recording eye data
            tStartCaliPoint = GetSecs;
            counter = 1;
            while GetSecs - tStartCaliPoint < calibration.timeCaliPoint               
                real_data = WiseData('get');
                tempTimeStamp = GetSecs;
                dif_SN = real_data(3,1) - last_SN;
                if  dif_SN > 0
                    if dif_SN > 20          
                        dif_SN = 20;
                    end

                    for CurrentSN = 1 : dif_SN
                        Col_Index = 2 * ( dif_SN - CurrentSN ) + 2 ;    
                        L_eyeXs(eyepos_SN) = real_data( 1 , Col_Index);
                        L_eyeYs(eyepos_SN) = real_data( 2 , Col_Index);
                        L_eyeSizes(eyepos_SN) = real_data( 3 , Col_Index);
                        L_CRXs(eyepos_SN) = real_data( 5, Col_Index);
                        L_CRYs(eyepos_SN) = real_data( 6, Col_Index);

                        R_eyeXs(eyepos_SN) = real_data( 1 , Col_Index + 1 );
                        R_eyeYs(eyepos_SN) = real_data( 2 , Col_Index + 1 );
                        R_eyeSizes(eyepos_SN) = real_data(3,Col_Index + 1);
                        R_CRXs(eyepos_SN) = real_data( 5 , Col_Index + 1 );
                        R_CRYs(eyepos_SN) = real_data( 6 , Col_Index + 1 );

                        eye_SNs(eyepos_SN) = real_data(3,1) - dif_SN + CurrentSN;
                        eyepos_SN = eyepos_SN + 1;                    
                    end
                    last_SN = real_data(3,1);
                    timeStamp.sampleNumber(counter) = last_SN;
                    timeStamp.time(counter) = tempTimeStamp;
                    counter = counter + 1;
                end
            end

            L_eyeXs(eyepos_SN:end) = [];    
            L_eyeYs(eyepos_SN:end) = [];    
            L_eyeSizes(eyepos_SN:end) = [];
            R_eyeXs(eyepos_SN:end) = [];    
            R_eyeYs(eyepos_SN:end) = [];    
            R_eyeSizes(eyepos_SN:end) = [];
            eye_SNs(eyepos_SN:end) = [];    
            L_CRXs(eyepos_SN:end) =[];      
            L_CRYs(eyepos_SN:end) =[]; 
            R_CRXs(eyepos_SN:end) =[];      
            R_CRYs(eyepos_SN:end) =[];
            timeStamp.sampleNumber(counter:end) = [];
            timeStamp.time(counter:end) = [];
            
            Data(i_trial).EyePos.L_eyeXs = L_eyeXs; 
            Data(i_trial).EyePos.L_eyeYs = L_eyeYs;
            Data(i_trial).EyePos.L_eyeSizes = L_eyeSizes;
            Data(i_trial).EyePos.R_eyeXs = R_eyeXs;
            Data(i_trial).EyePos.R_eyeYs = R_eyeYs;
            Data(i_trial).EyePos.R_eyeSizes = R_eyeSizes;
            Data(i_trial).EyePos.L_CRXs = L_CRXs;
            Data(i_trial).EyePos.L_CRYs = L_CRYs;
            Data(i_trial).EyePos.R_CRXs = R_CRXs;
            Data(i_trial).EyePos.R_CRYs = R_CRYs;
            Data(i_trial).EyePos.eye_SNs = eye_SNs;
            Data(i_trial).EyePos.timeStamp.sampleNumber = timeStamp.sampleNumber;
            Data(i_trial).EyePos.timeStamp.time = timeStamp.time;            
       end
            
        Screen('TextSize', windowDisplay, textFontSize);
        DrawFormattedText(windowDisplay, 'Calibration done!', 'center', 'center', textColor);
        Screen(windowDisplay,'Flip');
        WaitSecs(1); 
    end
end