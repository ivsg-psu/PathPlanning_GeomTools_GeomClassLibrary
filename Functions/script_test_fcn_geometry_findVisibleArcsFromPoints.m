% script_test_fcn_geometry_findVisibleArcsFromPoints.m
% Tests fcn_geometry_findVisibleArcsFromPoints

% Revision history:
%      2021_05_22:
%      -- first write of the code,
%      -- copying functionality from
%      -- script_test_fcn_geometry_findTangentPointsFromPointToCircle

%% BASIC example for one circle and one point
fig_num = 1;
centers = [0 0];
radii = 1;
points = [2 3];
visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(...
    centers,radii,points,fig_num);

fcn_summarize(visible_arc_angles,...
    centers,...
    radii,...
    points);

assert(isequal(round(visible_arc_angles,4),2.5559));

%% ADVANCED example that uses vectors of centers and points
fig_num = 2;
centers = [0 0; 1 4];
radii = [1; 1];
points = [2 3; 3 4];
visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(...
    centers,radii,points,fig_num);

fcn_summarize(visible_arc_angles,...
    centers,...
    radii,...
    points);

%% ADVANCED example that has one point too close to the center
fig_num = 3;
centers = [0 0; 1 4];
radii = [1; 1];
points = [0.5 0.5; 3 4];
visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(...
    centers,radii,points,fig_num);

fcn_summarize(visible_arc_angles,...
    centers,...
    radii,...
    points);

%% ADVANCED example that has one point too close to the center
fig_num = 4;
centers = [0 0; 1 4];
radii = [1; 1];
points = [1 2; 1.5 4];
visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(...
    centers,radii,points,fig_num);

fcn_summarize(visible_arc_angles,...
    centers,...
    radii,...
    points);



%% ADVANCED example that lets user select a point among circles, and move it around
enable_advanced_example = false; % flag advanced example off for non-interactive execution
if enable_advanced_example
    fig_num = 999;

    centers = [0 0; 1 4; 4 2];
    radii = [1; 2; 0.5];
    points = [-2*ones(length(radii),1) 4*ones(length(radii),1)];
    points_tangent = fcn_geometry_findVisibleArcsFromPoints(...
        centers,radii,points,fig_num);
    % Loop until right button is hit
    button = 1;
    while sum(button) <=1   % read ginputs until a mouse right-button occurs
        % Get a new point and redo plot
        [xp,yp,button] = ginput(1);
        points = [xp*ones(length(centers(:,1)),1),...
            yp*ones(length(centers(:,1)),1)];
        points_tangent = fcn_geometry_findVisibleArcsFromPoints(...
            centers,radii,points,fig_num);
    end
end

%% FAIL CONDITIONS
if 1==0
    %% ADVANCED example that has one point on the center
    fig_num = 2;
    centers = [0 0; 1 4];
    radii = [1; 1];
    points = [1 1; 1 4];
    visible_arc_angles = fcn_geometry_findVisibleArcsFromPoints(...
        centers,radii,points,fig_num);

    fcn_summarize(visible_arc_angles,...
        centers,...
        radii,...
        points);

end



%% Function start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcn_summarize(...
    visible_arc_angles,...
    centers,...
    radii,...
    points)

for i=1:length(visible_arc_angles(:,1))
    fprintf(1,'\n\nCenters: %.2f %.2f\n',centers(i,1),centers(i,2));
    fprintf(1,'Raddi: %.2f ',radii(i,1));
    fprintf(1,'Points: %.2f %.2f\n',points(i,1),points(i,2));
    fprintf(1,'Calculated visibility angle: %.2f\n',visible_arc_angles(i,1)*180/pi);

end % Ends for loop
end % Ends function
