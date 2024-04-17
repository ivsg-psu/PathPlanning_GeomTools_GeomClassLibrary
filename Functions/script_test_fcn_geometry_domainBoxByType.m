% script_test_fcn_geometry_domainBoxByType
% Tests fcn_geometry_domainBoxByType

% Revision history:
% 2024_02_12 - S . Brennan
% -- Edited for new function using plotCircle as starter
% 2024_04_11 - S . Brennan
% -- Added assertions

close all

rng(1);

%% check line boxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  _      _
% | |    (_)
% | |     _ _ __   ___  ___
% | |    | | '_ \ / _ \/ __|
% | |____| | | | |  __/\__ \
% |______|_|_| |_|\___||___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASIC example to produce domain for one line, typical situation
fig_num = 1;
figure(fig_num); clf;

unit_line_projection_vector = fcn_geometry_calcUnitVector([1 1]);
base_point_on_line = [0 0];
transverse_distance_to_lowest_point = -1;
transverse_distance_to_highest_point = 5;
distance_from_line_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'line',...
    unit_line_projection_vector, base_point_on_line, ...
    [transverse_distance_to_lowest_point, transverse_distance_to_highest_point], ...
    distance_from_line_to_boundary,...
    (fig_num));


assert(issimplified(domain_box));

%% BASIC example to show it works even if transverse distances in wrong order
fig_num = 2;
figure(fig_num); clf;

unit_line_projection_vector = fcn_geometry_calcUnitVector([1 1]);
base_point_on_line = [0 0];
transverse_distance_to_lowest_point  = 1;
transverse_distance_to_highest_point = -5;
distance_from_line_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'line',...
    unit_line_projection_vector, base_point_on_line, ...
    [transverse_distance_to_lowest_point, transverse_distance_to_highest_point], ...
    distance_from_line_to_boundary,...
    (fig_num));


assert(issimplified(domain_box));

%% BASIC example to show it works with negative distances
fig_num = 3;
figure(fig_num); clf;

unit_line_projection_vector = fcn_geometry_calcUnitVector([1 1]);
base_point_on_line = [0 0];
transverse_distance_to_lowest_point  = 1;
transverse_distance_to_highest_point = -5;
distance_from_line_to_boundary = -0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'line',...
    unit_line_projection_vector, base_point_on_line, ...
    [transverse_distance_to_lowest_point, transverse_distance_to_highest_point], ...
    distance_from_line_to_boundary,...
    (fig_num));


assert(issimplified(domain_box));



%% BASIC example - show that can calculate line results without a figure

unit_line_projection_vector = fcn_geometry_calcUnitVector([1 1]);
base_point_on_line = [0 0];
transverse_distance_to_lowest_point  = 1;
transverse_distance_to_highest_point = -5;
distance_from_line_to_boundary = -0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'line',...
    unit_line_projection_vector, base_point_on_line, ...
    [transverse_distance_to_lowest_point, transverse_distance_to_highest_point], ...
    distance_from_line_to_boundary,...
    (-1));


assert(issimplified(domain_box));
%% check arc boxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     /\
%    /  \   _ __ ___ ___
%   / /\ \ | '__/ __/ __|
%  / ____ \| | | (__\__ \
% /_/    \_\_|  \___|___/
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BASIC example to produce domain for one arc, typical situation
fig_num = 20;
figure(fig_num); clf;

circleCenter = [1 3];
circleRadius = 2;
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
degree_step = 1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

assert(issimplified(domain_box));

%% BASIC example to produce domain for one arc, uncertainty produces negative radii
fig_num = 21;
figure(fig_num); clf;

circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
degree_step = 1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 2;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

assert(~issimplified(domain_box));

%% BASIC example for one arc plotting from 0 to 90, negative
fig_num = 22;
figure(fig_num); clf;

circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = -90 * pi/180;
degree_step = -1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

assert(issimplified(domain_box));

%% BASIC example for one arc plotting from -90 to 90, positive
fig_num = 23;
figure(fig_num); clf;

circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

assert(issimplified(domain_box));

%% BASIC example for one arc plotting from 90 to -90, positive
% Must add 360 degrees to end point 

fig_num = 24;
figure(fig_num); clf;

circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians =  90 * pi/180;
end_angle_in_radians   =  (-90 + 360) * pi/180;
degree_step = 1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

assert(issimplified(domain_box));

%% BASIC example - show that can change number of points in the arc via degree_step
fig_num = 25;
figure(fig_num); clf;

circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 30;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

assert(issimplified(domain_box));

%% BASIC example - show that can calculate arc results without a figure

circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 30;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    ([]));

assert(issimplified(domain_box));

%% Speed test
% Set defaults
circleCenter = [1 3];
circleRadius = 2; 
start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 30;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;



% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;

    domain_box1 = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));
    
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;

    domain_box2 = fcn_geometry_domainBoxByType(...
        'arc',...
        circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
        (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

assert(isequal(domain_box1,domain_box2));
assert(averageTimeSlow>averageTimeFast);

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_domainBoxByType:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Break cases follow
% - these are ones that intentionally crash the code by passing invalid
% arguments
if 1==0
%% BREAK CASES 1 - break on centers
fig_num = 999;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4];
fcn_geometry_domainBoxByType(centers,radii,[0.1 0.1 1])

end
