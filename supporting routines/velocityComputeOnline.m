function velocity = velocityComputeOnline(eyeGaze_X_Y)
    deltaTime = 0.001;
    velocity = (-eyeGaze_X_Y(5, :) + 8*eyeGaze_X_Y(4, :) ...
            - 8*eyeGaze_X_Y(2, :) + eyeGaze_X_Y(1, :)) / (deltaTime*12);
end