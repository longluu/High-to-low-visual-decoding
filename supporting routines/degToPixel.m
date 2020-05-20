function distancePixel = degToPixel(distanceDeg, DisEye2Scr, ScreenWidthMm, ScreenWidthPixel)
PixelPerMm = ScreenWidthPixel/ScreenWidthMm;     % pixels per mm
distancePixel = round(DisEye2Scr * tan(distanceDeg*pi/180) * PixelPerMm); 

