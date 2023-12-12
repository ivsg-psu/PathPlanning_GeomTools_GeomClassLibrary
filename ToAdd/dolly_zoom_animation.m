close all

while(1)

    for ith_frame = 1:48
        drawframe(ith_frame);
        pause(0.01);
    end
end
function drawframe(f)
    t = interp1([0 48],[0 2*pi],f);
    r = [0.528 1.84 0.223 0.664 0.689 0.0951;
        1.24 0.288 0.272 0.702 0.857 0.427;
        0.158 0.951 1.17 0.501 0.0484 0.511;
        1.37 0.911 1.22 0.54 0.665 0.656;
        0.871 0.162 2.1 0.991 1.45 0.125;
        1.57 0.489 0.39 0.989 1.38 0.53];
    
    % Make some buildings
    bar3(r,0.4,'cyan')
    ax = gca;
    
    % Ground plane
    patch([0 7 7 0],[0 0 7 7],'c')
    axis([0 7 0 7])
    light(Position=[-1.7305 -2.0113 0.2109])
    ax.DataAspectRatio=[1 1 0.75];
    set(gcf,Color='white')
    axis off
    axis vis3d
    
    % Turn on perspective
    ax.Projection='perspective';
    
    % Pick a more interesting view
    ax.CameraViewAngle = 40;
    ax.CameraPosition = [13.167 0.01 0.7388];
    distToTarget1 = ax.CameraPosition - ax.CameraTarget;
    closenessDollyZoom = 50*sin(t) + 70;
    distToTarget2 = distToTarget1 * (closenessDollyZoom/100);
    viewAngle2 = 2*atan(distToTarget1/distToTarget2*tan(ax.CameraViewAngle/2*pi/180))*180/pi;
    ax.CameraViewAngle = viewAngle2;
    ax.CameraPosition = distToTarget2 + ax.CameraTarget;
    
end