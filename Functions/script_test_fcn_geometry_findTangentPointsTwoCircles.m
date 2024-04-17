%% script_test_fcn_geometry_findTangentPointsTwoCircles
% Tests fcn_geometry_findTangentPointsTwoCircles

% Revision history:
% 2021_04_22 - S. Brennan
% -- first write of the code copying functionality from fcn_FastestTraversal_checkInputsToFunctions
% 2024_04_16 - S. Brennan
% -- added assertions

close all

%% BASIC example - find all the points for one circle
fig_num = 1;
figure(fig_num); clf;

centers_start = [0 0];
centers_end   = [0 4];
radii_start   = 0.3; 
radii_end     = 0.2;
flag_inside_or_outside = 0;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);
 


% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs
assert(isequal(round(points_tangent_start,4),[0.2976,0.0375;-0.2976,0.0375;0.2999,0.0075;-0.2999,0.0075]));
assert(isequal(round(points_tangent_end,4),[-0.1984,3.9750;0.1984,3.9750;0.1999,4.0050;-0.1999,4.005]));

%% BASIC example - find all the points for two circles
fig_num = 2;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = 0;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);
 

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only inside points
fig_num = 3;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = -1;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);


% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only outside points
fig_num = 4;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = 1;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only outside points with equal radii
fig_num = 5;
figure(fig_num); clf;

centers_start = [0 0; 1 1; 2 2; -2 2; 2 -2; -2 -2];
centers_end   = [0 4; 5 1; 4 4; -4 4; 4 -4; -4 -4];
radii_start   = [0.3; 0.3; 0.3;  0.3;  0.3;  0.3;];
radii_end     = [0.3; 0.3; 0.3;  0.3;  0.3;  0.3;];
flag_inside_or_outside = 1;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only outside points, with mixed equal/unequal
% radii
fig_num = 6;
figure(fig_num); clf;

centers_start = [0 0; 1 1; 2 2; -2 2; 2 -2; -2 -2];
centers_end   = [0 4; 5 1; 4 4; -4 4; 4 -4; -4 -4];
radii_start   = [0.3; 0.3; 0.3;  0.3;  0.3;  0.3;];
radii_end     = [0.6; 0.3; 0.8;  0.3;  0.9;  0.2;];
flag_inside_or_outside = 1;
voting_points_start = [];
voting_points_end   = [];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only outside points with voting to force only
% one set to be returned.
fig_num = 7;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = 1;
voting_points_start = [-0.3 0; 0.7 5];
voting_points_end   = [-0.3 4; 3.7 8];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only inside points with voting to force only
% one set to be returned.
fig_num = 8;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = -1;
voting_points_start = [0.3 0; 2 5];
voting_points_end   = [-0.3 4; 3.7 8];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% BASIC example - find only inside points with voting to force only
% one set to be returned. 
fig_num = 9;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = -1;
voting_points_start = [0.3 0; 2 5];
voting_points_end   = [-0.3 4; 3.7 8];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);

% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs


%% ADVANCED example - find only inside points with voting to force only
% one set to be returned, but give a bad vote (this throws a warning
% and gives the wrong answer!)

fig_num = 9999;
figure(fig_num); clf;

centers_start = [0 0; 1 5];
centers_end   = [0 4; 4 8];
radii_start   = [0.3; 0.5];
radii_end     = [0.2; 0.7];
flag_inside_or_outside = -1;
voting_points_start = [0.3 0; 2 5];
voting_points_end   = [-0.3 4; 5 8];

[points_tangent_start, points_tangent_end] = ...
    fcn_geometry_findTangentPointsTwoCircles(...
    centers_start,...
    centers_end,...
    radii_start,...
    radii_end,...
    flag_inside_or_outside,...
    voting_points_start,voting_points_end,...
    fig_num);


% Check variable sizes
assert(length(points_tangent_start(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_start(:,1))>=1); % Does it have 1 or more rows?
assert(length(points_tangent_end(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent_end(:,1))>=1); % Does it have 1 or more rows?

% Check outputs

