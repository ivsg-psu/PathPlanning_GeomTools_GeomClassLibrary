% script_test_fcn_geometry_plotArc
% Tests fcn_geometry_plotArc

% Revision history:
%      2024_02_12 - S . Brennan
%      -- Edited for new function using plotCircle as starter

close all

%% BASIC example for one arc plotting from 0 to 90, counter-clockwise
figNum = 1;
figure(figNum); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example for one arc plotting from 0 to 90, clockwise
figNum = 101;
figure(figNum); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
flag_arc_is_counterclockwise = 0;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);


%% BASIC example for one arc plotting from 0 to 90, negative
figNum = 2;
figure(figNum); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = -90 * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example for one arc plotting from -90 to 90, positive
figNum = 3;
figure(figNum); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example for one arc plotting from 90 to -90, positive
% Must add 360 degrees to end point 

figNum = 4;
figure(figNum); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians =  90 * pi/180;
end_angle_in_radians   =  (-90 + 360) * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example for multiple arcs
figNum = 5;
figure(figNum); clf;

centers = [1 3; 2 4];


radii = [2; 3];
start_angle_in_radians = [0; 180] * pi/180;
end_angle_in_radians = [90; 200] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example for multiple arcs, with mixed start
figNum = 501;
figure(figNum); clf;

centers = [1 3; 2 4];


radii = [2; 3];
start_angle_in_radians = [0; 180] * pi/180;
end_angle_in_radians = [90; -200] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);





%% BASIC example - pass in color string
figNum = 6;
figure(figNum); clf;

centers = [1 2];
radii = 3;
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = 'm';

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);


%% BASIC example - pass in point attributes
figNum = 7;
figure(figNum); clf;

centers    = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = 'r.';

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);


%% BASIC example - show that can pass in color index
figNum = 8;
figure(figNum); clf;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = [0.5 0.5 1]; % A light blue

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);


%% BASIC example - show that can pass in full complex string
figNum = 9;
figure(figNum); clf;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = sprintf(' ''.'',''Color'',[0 0.5 0],''MarkerSize'', 20');

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example - show that multiple arcs work with negative endpoints
figNum = 901;
figure(figNum); clf;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90-360; 200-360; 320-360] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = sprintf(' ''.'',''Color'',[0 0.5 0],''MarkerSize'', 20');

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);


%% BASIC example 6 - change colors on plotting
figNum = 10;
figure(figNum); clf;
hold on;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];


for i_arc=1:length(centers(:,1))
    fcn_geometry_plotArc(centers(i_arc,:),radii(i_arc), start_angle_in_radians(i_arc,:), end_angle_in_radians(i_arc,:),(flag_arc_is_counterclockwise), (degree_step), [0  0 0.3*i_arc], figNum);
end

%% BASIC example - show that can change number of points in the arc via degree_step
figNum = 11;
figure(figNum); clf;

% Show defaults in blue
centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = [];
format = sprintf(' ''.'',''Color'',[0 0 1],''MarkerSize'', 20');

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

% Plot sparse in red
centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = 10;
format = sprintf(' ''.'',''Color'',[1 0 0],''MarkerSize'', 40');

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format), figNum);

%% BASIC example - show that can can use fcn_geometry_plotArc to generate arc data, in matrix or cell arrays, even without plotting 
% set figNum to empty after clearing the figure
figNum = 12;
figure(figNum); clf;

centers = [3 4];
radii = 2; 
start_angle_in_radians = 45 * pi/180;
end_angle_in_radians = 135 * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = []; % Default is 1 degree
format = [];
figNum = [];

arc_points_matrix = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format),(figNum));

% Pull multiple arcs at the same time
centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
degree_step = 5;
format = [];
figNum = [];

arc_points_cell_array = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format),(figNum));


figNum = 12;
figure(figNum); clf;
figure(figNum);
hold on; grid on;
axis equal
plot(arc_points_matrix(:,1),arc_points_matrix(:,2),'b.-','MarkerSize',20);
plot(arc_points_cell_array{1}(:,1),arc_points_cell_array{1}(:,2),'k.-','MarkerSize',20);


%% Test of fast implementation mode 

centers = [3 4];
radii = 2; 
start_angle_in_radians = 45 * pi/180;
end_angle_in_radians = 135 * pi/180;
flag_arc_is_counterclockwise = [];
degree_step = []; % Default is 1 degree
format = [];
figNum = [];

arc_points_matrix = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format),(figNum));


% Perform the calculation in slow mode
figNum = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    arc_points_matrix = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format),(figNum));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
figNum = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    arc_points_matrix = fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (flag_arc_is_counterclockwise), (degree_step), (format),(figNum));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_plotArc:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Speed test

%% Break cases follow
% - these are ones that intentionally crash the code by passing invalid
% arguments
if 1==0
%% BREAK CASES 1 - break on centers
figNum = 999;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4];
fcn_geometry_plotArc(centers,radii,[0.1 0.1 1])

end
