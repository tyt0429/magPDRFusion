function [g_d_az]=gyroheading(gyro,time_interval,start,g_d_az)
for i = start+1:size(gyro,1)
    if i == 1
        g_d(i,1) = -gyro(i,3) * time_interval * 180 / pi;
        g_d_az(i,1) = g_d(i,1) + 90;
        continue
    end
    
    g_d(i,1) = -gyro(i,3) * time_interval * 180/pi;
    
    if g_d(i,1) + g_d_az(i-1,1) > 360
        g_d_az(i,1) = g_d(i,1) + g_d_az(i-1,1) - 360;
        continue
    end
    
    if g_d(i,1) + g_d_az(i-1,1) < 0
        g_d_az(i,1) = g_d(i,1) + g_d_az(i-1,1) + 360;
        continue
    end
    
    g_d_az(i,1) = g_d(i,1) + g_d_az(i-1,1);
end