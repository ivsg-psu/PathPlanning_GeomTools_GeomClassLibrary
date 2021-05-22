% script_test_fcn_geometry_findTangentPointFromPointToCircle.m
% Tests fcn_geometry_findTangentPointFromPointToCircle

% Revision history:
%      2021_04_22:
%      -- first write of the code copying functionality from fcn_FastestTraversal_checkInputsToFunctions

%% BASIC example for one circle and one point (positive cross product)
fig_num = 1;
centers = [0 0];
radii = 1;
points = [2 3];
cross_prod = 1; % Keep the postive cross product
points_tangent = fcn_geometry_findTangentPointFromPointToCircle(...
    centers,radii,points,cross_prod,fig_num); %#ok<*NASGU>

%% BASIC example for one circle and one point (negative cross product)
fig_num = 2;
centers = [0 0];
radii = 1;
points = [2 3];
cross_prod = -3;
points_tangent = fcn_geometry_findTangentPointFromPointToCircle(...
    centers,radii,points,cross_prod,fig_num); %#ok<*NASGU>

%% ADVANCED example that uses vectors of centers and points
fig_num = 3;
centers = [0 0; 1 4];
radii = [1; 0.5];
points = [2 3; 3 4];
cross_prod = [-3; 3];
points_tangent = fcn_geometry_findTangentPointFromPointToCircle(...
    centers,radii,points,cross_prod,fig_num);

%% ADVANCED example that lets user select a point among circles, and move it around
fig_num = 3;

centers = [0 0; 1 4; 4 2];
radii = [1; 2; 0.5];
points = [-2*ones(length(radii),1) 4*ones(length(radii),1)];
cross_prod = ones(length(radii),1);
points_tangent = fcn_geometry_findTangentPointFromPointToCircle(...
    centers,radii,points,cross_prod,fig_num);
% Loop until right button is hit
button = 1;
while sum(button) <=1   % read ginputs until a mouse right-button occurs
    % Get a new point and redo plot
    [xp,yp,button] = ginput(1);
    points = [xp*ones(length(centers(:,1)),1),...
        yp*ones(length(centers(:,1)),1)];
    points_tangent = fcn_geometry_findTangentPointFromPointToCircle(...
        centers,radii,points,cross_prod,fig_num);
end