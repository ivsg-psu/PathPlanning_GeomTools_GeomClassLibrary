% script_test_fcn_geometry_plotArc
% Tests fcn_geometry_plotArc

% Revision history:
%      2024_02_12 - S . Brennan
%      -- Edited for new function using plotCircle as starter

close all

%% BASIC example for one arc plotting from 0 to 90, positive
fig_num = 1;
figure(fig_num); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);


%% BASIC example for one arc plotting from 0 to 90, negative
fig_num = 2;
figure(fig_num); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = -90 * pi/180;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);

%% BASIC example for one arc plotting from -90 to 90, positive
fig_num = 3;
figure(fig_num); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians = -90 * pi/180;
end_angle_in_radians   =  90 * pi/180;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);

%% BASIC example for one arc plotting from 90 to -90, positive
% Must add 360 degrees to end point 

fig_num = 4;
figure(fig_num); clf;

centers = [1 3];
radii = 2; 
start_angle_in_radians =  90 * pi/180;
end_angle_in_radians   =  (-90 + 360) * pi/180;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);

%% BASIC example for multiple arcs
fig_num = 5;
figure(fig_num); clf;

centers = [1 3; 2 4];
radii = [2; 3];
start_angle_in_radians = [0; 180] * pi/180;
end_angle_in_radians = [90; 200] * pi/180;
degree_step = [];
format = [];

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);



%% BASIC example - pass in color string
fig_num = 6;
figure(fig_num); clf;

centers = [1 2];
radii = 3;
start_angle_in_radians = 0 * pi/180;
end_angle_in_radians = 90 * pi/180;
degree_step = [];
format = 'm';

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);


%% BASIC example - pass in point attributes
fig_num = 7;
figure(fig_num); clf;

centers    = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
degree_step = [];
format = 'r.';

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);


%% BASIC example - show that can pass in color index
fig_num = 8;
figure(fig_num); clf;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
degree_step = [];
format = [0.5 0.5 1]; % A light blue

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);


%% BASIC example - show that can pass in full complex string
fig_num = 9;
figure(fig_num); clf;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;
degree_step = [];
format = sprintf(' ''.'',''Color'',[0.5 1 0.5],''MarkerSize'', 20''');

fcn_geometry_plotArc(centers, radii, start_angle_in_radians, end_angle_in_radians, (degree_step), (format), fig_num);


%% BASIC example 6 - change colors on plotting
fig_num = 6;
figure(fig_num); clf;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
start_angle_in_radians = [0; 180; 270] * pi/180;
end_angle_in_radians = [90; 200; 320] * pi/180;

for i_arc=1:length(centers(:,1))
    fcn_geometry_plotArc(centers(i_arc,:),radii(i_arc), start_angle_in_radians(i_arc,:), end_angle_in_radians(i_arc,:),[0.3*i_arc 0.3*i_arc 1],fig_num);
end

%% Speed test

%% Break cases follow
% - these are ones that intentionally crash the code by passing invalid
% arguments
if 1==0
%% BREAK CASES 1 - break on centers
fig_num = 999;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4];
fcn_geometry_plotArc(centers,radii,[0.1 0.1 1])

end
