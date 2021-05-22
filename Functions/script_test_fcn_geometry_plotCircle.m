% script_test_fcn_geometry_plotCircle
% Tests fcn_geometry_plotCircle

% Revision history:
%      2021_05_22:
%      -- Edited for new function

close all

%% BASIC example for one circle
fig_num = 1;
figure(fig_num); axis square; grid minor;

center = [1 3];
radius = [2];
fcn_geometry_plotCircle(center,radius,fig_num);


%% BASIC example for multiple circles
fig_num = 2;
figure(fig_num); axis square; grid minor;

centers = [1 3; 2 4];
radii = [2; 3];
fcn_geometry_plotCircle(centers,radii,fig_num);

