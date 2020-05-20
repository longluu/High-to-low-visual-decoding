function distanceDeg = pixelToDeg(positionPixel, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel, screenCenter, horizontalPos)
if horizontalPos == 1
    positionMm = (positionPixel - screenCenter(1)) * ScreenWidthMm / ScreenWidthPixel;
else
    positionMm = (positionPixel - screenCenter(2)) * ScreenWidthMm / ScreenWidthPixel;    
end
distanceDeg = rad2deg(atan(positionMm / DisEye2Scr));
