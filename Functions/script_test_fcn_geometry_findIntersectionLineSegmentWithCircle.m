% Script to test the fcn_geometry_findIntersectionOfLineSegmentWithCircle

% Set flags for whether or not to use existing x,y data to run the check
% and whether to get new data (if requested) with the mouse or simply run
% the static test cases
flag_getNewData = 1;
flag_getDataViaMouse = 0;

if flag_getNewData
    % Clear variables out of the workspace
    clearvars -except flag*
% Check to make sure the get new data flag was not unset by accident (for
% example if the script was not previously run) and reset it if the x or y
% data are missing
elseif ~exist('x','var') || ~exist('y','var')
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

if flag_getNewData && flag_getDataViaMouse
    % Prompt the user for two points
    title('Use mouse clicks to set the two endpoints of the line segment');
    [x,y] = ginput(2);
elseif flag_getNewData
  x = [2; 10];
  y = [2; 5];
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
    fprintf(1,'%d intersection(s) found between given line segment and circle\n',Nintersections);
else
    fprintf(1,'No intersections found between given line segment and circle\n');
end