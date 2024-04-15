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
fig_num = 1;
points = [0 0; 1 4; 0.5 -1];
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);

% fcn_summarize(...
%     centers,...
%     radii,...
%     points);
assert(isequal(round(centers,4), [3.6667,1.2083]));
assert(isequal(round(radii,4), 3.8606));

%% ADVANCED example that uses vectors of x and y
fig_num = 100;
points = [0 0; 1 4; 0.5 -1; -1 4];
[centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);


% fcn_summarize(...
%     centers,...
%     radii,...
%     points);

assert(isequal(round(centers(1,:),4), [3.6667    1.2083]));
assert(isequal(round(centers(2,:),4), [0.0000    1.5750]));
assert(isequal(round(radii(1,1),4), 3.8606));
assert(isequal(round(radii(2,1),4), 2.6231));

%% ADVANCED example that lets user select N points
% enable_advanced_example = false; % flag advanced example off for non-interactive execution
% if enable_advanced_example
%     fig_num = 1000;
%     figure(fig_num); clf; grid on; axis equal;
% 
%     points = ginput; % Get arbitrary N points until user hits return
%     [centers,radii] = fcn_geometry_circleCenterFrom3Points(points,fig_num);
% 
%     fcn_summarize(...
%         centers,...
%         radii,...
%         points);
% end


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
end


%% ADVANCED example that lets user select N points
% enable_advanced_example = false; % flag advanced example off for non-interactive execution
% if enable_advanced_example
%     fig_num = 1001;
%     x = [1; 2];
%     y = [1; 2];
%     button = 1;
%     figure(fig_num); clf; grid on; axis equal;
%     while sum(button) <=1   % read ginputs until a mouse right-button occurs
%         % Get a new point and redo plot
%         [x(end+1),y(end+1),button] = ginput(1); %#ok<SAGROW>
%         points = [x, y];
%         fcn_geometry_circleCenterFrom3Points(points,fig_num);
%     end
% end


%% ADVANCED example that lets user select N points
% enable_advanced_example = false; % flag advanced example off for non-interactive execution
% if enable_advanced_example
%     fig_num = 1002;
%     x = [1; 2; 3];
%     y = [1; 2; 3];
%     figure(fig_num); clf; grid on; axis equal;
% 
%     % The following tests the 3 input form:
%     button = 1;
%     while sum(button) <=1   % read ginputs until a mouse right-button occurs
%         % Shift points up to prep for next input
%         x(1:end-1) = x(2:end);
%         y(1:end-1) = y(2:end);
% 
%         % Get a new point and redo plot
%         [x(end),y(end),button] = ginput(1);
%         points = [x,y];
%         fcn_geometry_circleCenterFrom3Points(points,fig_num);
%     end
% end

%% BASIC example that tests the 3 separate point input styles
fig_num = 103;
points1 = [0 0; 5 0];
points2 = [-2 2; 7 2];
points3 = [0 4; 5 4];
fcn_geometry_circleCenterFrom3Points(points1, points2, points3,fig_num);


%% BASIC example that tests flag_max_speed

points1 = [0 0; 5 0];
points2 = [-2 2; 7 2];
points3 = [0 4; 5 4];


% Perform the calculation in slow mode
REPS = 1000; 
minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    fcn_geometry_circleCenterFrom3Points(points1, points2, points3,[]);   
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    fcn_geometry_circleCenterFrom3Points(points1, points2, points3,-1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_circleCenterFrom3Points:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);



%% FAIL CASES
if 1==0
    %% fails because flag_max_speed is set to value other than -1
    fig_num = -2;
    points1 = [0 0; 5 0];
    points2 = [-2 2; 7 2];
    points3 = [0 4; 5 4];
    flag_max_speed = 0;

    fcn_geometry_circleCenterFrom3Points(points1, points2, points3,fig_num);
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
