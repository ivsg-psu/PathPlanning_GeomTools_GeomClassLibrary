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

assert(isequal(round(points_tangent,4),[0.9533,-0.3022;-0.6456,0.7637]));

%% ADVANCED example that uses vectors of centers and points
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [2 3; 3 4];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);

assert(length(points_tangent(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent(:,1))>=1); % Does it have 1 or more rows?

%% ADVANCED example that has one point too close to the center
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [0.5 0.5; 3 4];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);

assert(length(points_tangent(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent(:,1))>=1); % Does it have 1 or more rows?

%% ADVANCED example that has one point too close to the center
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [1 2; 1.5 4];
points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points,fig_num);

assert(length(points_tangent(1,:))==2); % Does it have 2 columns?
assert(length(points_tangent(:,1))>=1); % Does it have 1 or more rows?

%% Test of fast implementation mode 

centers = [0 0];
radii = 1;
points = [2 3];

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points, (fig_num));
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
    points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
    centers,radii,points, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);



%% ADVANCED example that lets user select a point among circles, and move it around
% enable_advanced_example = false; % flag advanced example off for non-interactive execution
% if enable_advanced_example
%     fig_num = 999;
% 
%     centers = [0 0; 1 4; 4 2];
%     radii = [1; 2; 0.5];
%     points = [-2*ones(length(radii),1) 4*ones(length(radii),1)];
%     points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
%         centers,radii,points,fig_num);
%     % Loop until right button is hit
%     button = 1;
%     while sum(button) <=1   % read ginputs until a mouse right-button occurs
%         % Get a new point and redo plot
%         [xp,yp,button] = ginput(1);
%         points = [xp*ones(length(centers(:,1)),1),...
%             yp*ones(length(centers(:,1)),1)];
%         points_tangent = fcn_geometry_findTangentPointsFromPointToCircle(...
%             centers,radii,points,fig_num);
%     end
% end

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
