% script_test_fcn_geometry_findTangentPointsFromPointToCircle.m
% Tests fcn_geometry_findTangentPointsFromPointToCircle

% Revision history:
%      2021_04_22:
%      -- first write of the code copying functionality from fcn_FastestTraversal_checkInputsToFunctions

%% BASIC example for one circle and one point
fig_num = 1;
centers = [0 0];
radii = 1;
points = [2 3];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num); %#ok<*NASGU>

%% ADVANCED example that uses vectors of centers and points
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [2 3; 3 4];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);

%% ADVANCED example that has one point too close to the center
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [0.5 0.5; 3 4];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);

%% ADVANCED example that has one point too close to the center
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [1 2; 1.5 4];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);



%% ADVANCED example that lets user select a point among circles, and move it around
fig_num = 999;

centers = [0 0; 1 4; 4 2];
radii = [1; 2; 0.5];
points = [-2*ones(length(radii),1) 4*ones(length(radii),1)];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);
% Loop until right button is hit
button = 1;
while sum(button) <=1   % read ginputs until a mouse right-button occurs
    % Get a new point and redo plot
    [xp,yp,button] = ginput(1);
    points = [xp*ones(length(centers(:,1)),1),...
        yp*ones(length(centers(:,1)),1)];
    points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
        centers,radii,points,fig_num);
end

%% FAIL CONDITIONS
if 1==0
    %% ADVANCED example that has one point on the center
    fig_num = 2;
    centers = [0 0; 1 4];
    radii = [1; 1];
    points = [1 1; 1 4];
    points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
        centers,radii,points,fig_num);
    
end