% Script to test the fcn_geometry_findIntersectionOfLineSegmentWithCircle

flag_getNewData = 0;
if flag_getNewData
    % Clear variables out of the workspace
    clearvars
    flag_getNewData = 1;
end

% Set up figure one, clearing any existing plot objects
figure(1)
clf reset
hold on
grid on
axis equal

% Define the circle
pc = [0 0];
R = 25;

% Plot the circle
fcn_geometry_plotCircle(pc,R,'b-',1)

if flag_getNewData
    % Prompt the user for two points
    title('Use mouse clicks to set the two endpoints of the line segment');
    [x,y] = ginput(2);
end
title([]);
% Plot the user line segment
plot(x,y,'k.-','markersize',10)

% Assign the user data points to point arrays
pa = [x(1) y(1)];
pb = [x(2) y(2)];

% Run the intersection solution code
[intAngles,intPoints] = fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R);

% Determine the number of intersections
Nintersections = size(intAngles,1);

if Nintersections > 0
    % Plot the results
    plot(intPoints(:,1),intPoints(:,2),'r*')
    fprintf(1,'%d intersections found between given line segment and circle\n',Nintersections);
else
    fprintf(1,'No intersections found between given line segment and circle\n');
end