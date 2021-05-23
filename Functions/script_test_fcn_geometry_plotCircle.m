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
radius = [2]; %#ok<*NBRAK>
fcn_geometry_plotCircle(center,radius,[],fig_num);


%% BASIC example for multiple circles
fig_num = 2;
figure(fig_num); axis square; grid minor;

centers = [1 3; 2 4];
radii = [2; 3];
fcn_geometry_plotCircle(centers,radii,[],fig_num);


%% BASIC example 3
fig_num = 3;

centers = [1 2];
radii = 3;
fcn_geometry_plotCircle(centers,radii,'b',fig_num)

%% BASIC example 4
fig_num = 4;

centers    = [1 2; 2 4; 3 5]; 
radii = [3; 4; 5];
fcn_geometry_plotCircle(centers,radii,'r.',fig_num)

%% BASIC example 5
fig_num = 5;

centers  = [1 2; 2 4; 3 5]; 
radii = [3; 4; 5];

% Do a light blue 
fcn_geometry_plotCircle(centers,radii,[0.5 0.5 1],fig_num)

%% BASIC example 6
fig_num = 6;

centers  = [1 2; 2 4; 3 5]; 
radii = [3; 4; 5];

for i_circle=1:length(centers(:,1))
    fcn_geometry_plotCircle(centers(i_circle,:),radii(i_circle),[0.3*i_circle 0.3*i_circle 1],fig_num);
end

%% Break cases follow 
% - these are ones that intentionally crash the code by passing invalid
% arguments
if 1==0
%% BREAK CASES 1 - break on centers
fig_num = 999;
 
centers  = [1 2; 2 4; 3 5]; 
radii = [3; 4];
fcn_geometry_plotCircle(centers,radii,[0.1 0.1 1])

end
