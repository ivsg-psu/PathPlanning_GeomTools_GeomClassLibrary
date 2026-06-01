% script_test_fcn_geometry_plotCircle
% Tests fcn_geometry_plotCircle

% Revision history:
%      2021_05_22:
%      -- Edited for new function

close all

%% BASIC example for one circle
figNum = 1;
figure(figNum); axis square; grid minor;

center = [1 3];
radius = [2]; %#ok<*NBRAK>
fcn_geometry_plotCircle(center,radius,[],figNum);

%% BASIC example for multiple circles
figNum = 2;
figure(figNum); axis square; grid minor;

centers = [1 3; 2 4];
radii = [2; 3];
fcn_geometry_plotCircle(centers,radii,[],figNum);


%% BASIC example 3
figNum = 3;

centers = [1 2];
radii = 3;
fcn_geometry_plotCircle(centers,radii,'b',figNum);

%% BASIC example 4
figNum = 4;

centers    = [1 2; 2 4; 3 5];
radii = [3; 4; 5];
fcn_geometry_plotCircle(centers,radii,'r.',figNum);


%% BASIC example 5
figNum = 5;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];

% Do a light blue
fcn_geometry_plotCircle(centers,radii,[0.5 0.5 1],figNum);

%% BASIC example 6 - many circles
figNum = 6;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4; 5];

for i_circle=1:length(centers(:,1))
    fcn_geometry_plotCircle(centers(i_circle,:),radii(i_circle),[0.3*i_circle 0.3*i_circle 1],figNum);
end

%% BASIC example 7 - complex plot string
figNum = 7;
figure(figNum); clf;

centers    = [1 2];
radii = [3];
fcn_geometry_plotCircle(centers,radii,sprintf(' ''-'',''Color'',[0.6 0 0],''LineWidth'',2 '),figNum);



%% Break cases follow
% - these are ones that intentionally crash the code by passing invalid
% arguments
if 1==0
%% BREAK CASES 1 - break on centers
figNum = 999;

centers  = [1 2; 2 4; 3 5];
radii = [3; 4];
fcn_geometry_plotCircle(centers,radii,[0.1 0.1 1])

end
