% script_test_fcn_geometry_euclideanPointsToPointsDistance
% Tests fcn_geometry_euclideanPointsToPointsDistance

% Revision History:
% 2021-05-28 - S. Brennan
% -- revised function to prep for geometry class 
% -- rewrote function to use vector sum
% -- added plotting option
% 2021-06-05
% -- fixed comments, added debugging option


%% BASIC example - single points in 2D
fig_num = 1;

pt1 = [1 1];  
pt2 = [2 3];
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

%% BASIC example - two points in 2D
fig_num = 2;

pt1 = [1 1; 0 0];  
pt2 = [2 3; 4 0] ;
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

%% BASIC example - many points in 2D
fig_num = 3;

pt1 = rand(5,2);  
pt2 = rand(5,2);  
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);


%% BASIC example - single points in 3D
fig_num = 31;

pt1 = [1 1 0];  
pt2 = [2 3 2];
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

%% BASIC example - two points in 3D
fig_num = 32;

pt1 = [1 1 0; 0 0 1];  
pt2 = [2 3 4; 4 0 2] ;
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

%% BASIC example - multiple points in 3D
fig_num = 33;

pt1 = [-1 1 0; 0 0 1; -3 -2 -4];  
pt2 = [2 3 4; 4 0 2; -5 3 -2] ;
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);

%% BASIC example - multiple points in 3D
fig_num = 34;

pt1 = rand(5,3);
pt2 = rand(5,3);
dist=fcn_geometry_euclideanPointsToPointsDistance(pt1,pt2,fig_num);
