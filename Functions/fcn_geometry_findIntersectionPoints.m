function intersection_point = fcn_geometry_findIntersectionPoints(breakPointsCell, fig_num)

% This function is still under development. This function will/may be
% merged into fcn_geometry_findBreakPoints function. 

% Need to update the instructions

% Feb 21, 2024 - Aneesh Batchu
% Wrote the code originally

% Fit parameters of Line Segment
%
% breakPointsCell{}.fitParameters(1) = unit_projection_vector_x 
% breakPointsCell{}.fitParameters(2) = unit_projection_vector_y 
% breakPointsCell{}.fitParameters(3) = base_point_x 
% breakPointsCell{}.fitParameters(4) = base_point_y
% breakPointsCell{}.fitParameters(5) = station_distance_min 
% breakPointsCell{}.fitParameters(6) = station_distance_max

% Creating a matrix to store the best fit parameters
fitParametersMatrix = zeros([length(breakPointsCell),6]);

for i = 1: length(breakPointsCell)
    fitParametersMatrix(i,:) = breakPointsCell{i}.fitParameters;
    
end

% Finding the slopes of the lines
slopesLines = fitParametersMatrix(:,2)./fitParametersMatrix(:,1);
% Finding the y-intercepts of the lines 
yInterceptsLines = fitParametersMatrix(:,4) - slopesLines.*fitParametersMatrix(:,3);

% Line: ax+by+c = 0 ---> -mx+y-yIntercept = 0
% 
% a = -slopesLines; b = 1; c = -yInterceptsLines
%
% -mx + b = c
% AA = [-m1  1; -m2 1]
% bb = [c1; c2]

% Solving for intersecting point 
AA = [-1.*(slopesLines), ones(size(slopesLines))]; 
bb = (yInterceptsLines);

% Finding the intersection point
intersection_point = AA\bb;

% Printing the point
fprintf('Intersection Point: (%.4f, %.4f)\n', intersection_point(1), intersection_point(2));

% Plotting the intersection point using the best fit parameters
figure(fig_num)
hold on 

% Define x range for plotting
x_range = linspace(-10, 10, 100);

% Calculate y values for each line
y1 = slopesLines(1,1) * x_range + yInterceptsLines(1,1);
y2 = slopesLines(2,1) * x_range + yInterceptsLines(2,1);

% Plot the lines
% figure(888)
plot(x_range, y1, 'k--', 'LineWidth', 2);
plot(x_range, y2, 'b--', 'LineWidth', 2); % plot second line in red



end