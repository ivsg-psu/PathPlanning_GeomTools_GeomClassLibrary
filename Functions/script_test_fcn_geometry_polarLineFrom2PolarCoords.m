% script_test_fcn_geometry_polarLineFrom2PolarCoords
% Exercises the function: fcn_geometry_polarLineFrom2PolarCoords

% Revision history:
% 2021_05_27
% -- wrote the code

close all;
clc;



%% script_test_fcn_geometry_polarLineFrom2PolarCoords
fig_num = 1;

start_point_cart = [1 1];
end_point_cart   = [4 3];

[theta1,r1] = cart2pol(start_point_cart(:,1),start_point_cart(:,2));
[theta2,r2]   = cart2pol(end_point_cart(:,1),end_point_cart(:,2));
points = [theta1,r1;theta2,r2];

% Calculate the line
[phi,rho] = fcn_geometry_polarLineFrom2PolarCoords(points,fig_num);
axis([0 5 0 5]);
assert(isequal(round(phi,4),-0.9828));
assert(isequal(round(rho,4),-0.2774));
