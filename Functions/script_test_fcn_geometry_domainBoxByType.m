% script_test_fcn_geometry_domainBoxByType
% Tests fcn_geometry_domainBoxByType

% Revision history:
%      2024_02_12 - S . Brennan
%      -- Edited for new function using plotCircle as starter

close all

%% BASIC example to produce domain for one arc, typical situation
fig_num = 1;
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

%% BASIC example to produce domain for one arc, uncertainty produces negative radii
fig_num = 1;
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

%% BASIC example for one arc plotting from 0 to 90, negative
fig_num = 2;
figure(fig_num); clf;

start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = -90 * pi/180;
degree_step = -1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

%% BASIC example for one arc plotting from -90 to 90, positive
fig_num = 3;
figure(fig_num); clf;


start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));
%% BASIC example for one arc plotting from 90 to -90, positive
% Must add 360 degrees to end point 

fig_num = 4;
figure(fig_num); clf;

start_angle_in_radians =  90 * pi/180;
end_angle_in_radians   =  (-90 + 360) * pi/180;
degree_step = 1;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

%% BASIC example - show that can change number of points in the arc via degree_step
fig_num = 5;
figure(fig_num); clf;

start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 30;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    (fig_num));

%% BASIC example - show that can calculate results without a figure

start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = 30;
angles = (start_angle_in_radians:(degree_step*pi/180):end_angle_in_radians)';
distance_from_circle_to_boundary = 0.5;

domain_box = fcn_geometry_domainBoxByType(...
    'arc',...
    circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
    ([]));

%% Speed test
% Set defaults
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

    domain_box = fcn_geometry_domainBoxByType(...
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

    domain_box = fcn_geometry_domainBoxByType(...
        'arc',...
        circleCenter, circleRadius, angles, distance_from_circle_to_boundary, ...
        (fig_num));

    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

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
