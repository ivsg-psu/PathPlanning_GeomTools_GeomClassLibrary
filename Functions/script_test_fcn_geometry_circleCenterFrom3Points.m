% script_test_fcn_geometry_circleCenterFrom3Points - this is is a script written
% to test the function: fcn_geometry_circleCenterFrom3Points.m.
%
% Revision history:
% 2020_03_20 - started writing the function and script
% 2020_05_22 - added more comments
% 2021_05_23 
% -- merged previous function into geometry class
% -- automated input argument checking
% -- changed from x,y separate inputs into points inputs
close all;

%% BASIC example 1
geometry_circleCenterFrom3Points_case1 = matlab.unittest.TestCase.forInteractiveUse;
fig_num = 1;
points = [0 0; 1 4; 0.5 -1];
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);

fcn_summarize(...     
    centers,...
    radii,...
    points);
assertEqual(geometry_circleCenterFrom3Points_case1, round(centers,4), [3.6667,1.2083])
assertEqual(geometry_circleCenterFrom3Points_case1, round(radii,4), 3.8606)

%% ADVANCED example that uses vectors of x and y
fig_num = 100;
points = [0 0; 1 4; 0.5 -1; -1 4];
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);


fcn_summarize(...     
    centers,...
    radii,...
    points);
 
%% ADVANCED example that lets user select N points
fig_num = 1000;
figure(fig_num); clf; grid on; axis equal;

points = ginput; % Get arbitrary N points until user hits return
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);

fcn_summarize(...     
    centers,...
    radii,...
    points);
       

%% ADVANCED example that uses vectors of x and y
fig_num = 101;
points = [0 0; 0.5 4; 1 -1; 4 -3; 6 2; 7 -2; 9 3; 11 3; 15 -0.5];
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);
plot(points(:,1),points(:,2),'r-');

fcn_summarize(...     
    centers,...
    radii,...
    points);

%% BASIC example that gives same result as above points
fig_num = 102;
points = [0 0; 0.5 4; 1 -1; 4 -3; 6 2; 7 -2; 9 3; 11 3; 15 -0.5];
hold on
figure(1); clf;
for i=1:length(points(:,1))-2
    fcn_geometry_circleCenterFrom3Points(points(i:i+2,:),fig_num);
    plot(points(:,1),points(:,2),'r-');
    pause;
end


%% ADVANCED example that lets user select N points
fig_num = 1001;
x = [1; 2];
y = [1; 2];
button = 1;
figure(fig_num); clf; grid on; axis equal;
while sum(button) <=1   % read ginputs until a mouse right-button occurs   
    % Get a new point and redo plot
    [x(end+1),y(end+1),button] = ginput(1); %#ok<SAGROW>
    points = [x, y];
    fcn_geometry_circleCenterFrom3Points(points,fig_num);     
end


%% ADVANCED example that lets user select N points
fig_num = 1002;
x = [1; 2; 3];
y = [1; 2; 3];
figure(fig_num); clf; grid on; axis equal;

% The following tests the 3 input form:
button = 1;
while sum(button) <=1   % read ginputs until a mouse right-button occurs
    % Shift points up to prep for next input
    x(1:end-1) = x(2:end);
    y(1:end-1) = y(2:end);
    
    % Get a new point and redo plot
    [x(end),y(end),button] = ginput(1);
    points = [x,y];
    fcn_geometry_circleCenterFrom3Points(points,fig_num);
end


%% Function start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fcn_summarize(...
    centers,...
    radii,...
    points)

for i=1:length(centers(:,1))
    fprintf(1,'\n\n');
    fprintf(1,'Circle %.0d\n',i);
    for j=i:i+2
        fprintf(1,'Points %.0d: %.2f %.2f\n',j, points(j,1),points(j,2));
    end
    fprintf(1,'Centers: %.2f %.2f\n',centers(i,1),centers(i,2));
    fprintf(1,'Radii: %.2f \n',radii(i,1));
end % Ends for loop
end % Ends function
