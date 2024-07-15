%% script_test_fcn_geometry_fillColorFromNumberOrName
% Exercises the function: fcn_geometry_fillColorFromNumberOrName

% Revision history:
% 2024_04_11 - S. Brennan
% -- wrote the code

close all;


%% Basic call - integer color number

color_vector = fcn_geometry_fillColorFromNumberOrName(2);

% Check the output type and size
assert(isequal(size(color_vector),[1 3]));


%% Basic call - fractional color number with current color map

Npoints = 100;
x_data = linspace(-2,2,Npoints)';
y_data = x_data*3 + 7;
max_y = max(y_data);
min_y = min(y_data);
points = [x_data y_data];

fractional_y_data = (y_data - min_y)/(max_y - min_y);
color_vector = fcn_geometry_fillColorFromNumberOrName(fractional_y_data);

% Check the output type and size
assert(isequal(size(color_vector),[Npoints 3]));

% Make the plot
fig_num = 1;
figure(fig_num);
clf;
hold on;
grid on;
axis equal;

for ith_point = 1:Npoints
    plot(points(ith_point,1),points(ith_point,2),'.','Color',color_vector(ith_point,:));
end

%% Basic call - fractional color number with red to green color map

Npoints = 100;
x_data = linspace(-2,2,Npoints)';
y_data = x_data*3 + 7;
max_y = max(y_data);
min_y = min(y_data);
points = [x_data y_data];

fractional_y_data = (y_data - min_y)/(max_y - min_y);
color_vector = fcn_geometry_fillColorFromNumberOrName(fractional_y_data,[],'redtogreen');

% Check the output type and size
assert(isequal(size(color_vector),[Npoints 3]));

% Make the plot
fig_num = 2;
figure(fig_num);
clf;
hold on;
grid on;
axis equal;

for ith_point = 1:Npoints
    plot(points(ith_point,1),points(ith_point,2),'.','Color',color_vector(ith_point,:));
end

%% Basic call - fractional color number with custom colormap

Npoints = 100;
x_data = linspace(-2,2,Npoints)';
y_data = x_data*3 + 7;
max_y = max(y_data);
min_y = min(y_data);
points = [x_data y_data];

fractional_y_data = (y_data - min_y)/(max_y - min_y);
color_vector = fcn_geometry_fillColorFromNumberOrName(fractional_y_data,[],'jet');

% Check the output type and size
assert(isequal(size(color_vector),[Npoints 3]));

% Make the plot
fig_num = 3;
figure(fig_num);
clf;
hold on;
grid on;
axis equal;

for ith_point = 1:Npoints
    plot(points(ith_point,1),points(ith_point,2),'.','Color',color_vector(ith_point,:));
end
%% Basic call - string color number

color_vector = fcn_geometry_fillColorFromNumberOrName(2,'regression arc');

% Check the output type and size
assert(isequal(size(color_vector),[1 3]));
assert(isequal(color_vector,[1 0 0]));

%% Basic call - string color number

color_vector = fcn_geometry_fillColorFromNumberOrName(2,'vector regression segment fit');

% Check the output type and size
assert(isequal(size(color_vector),[1 3]));
assert(isequal(color_vector,[0 0 1]));

%% Test of fast mode

% Perform the calculation in slow mode
fig_num = [];
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    color_vector = fcn_geometry_fillColorFromNumberOrName(2,'vector regression segment fit',fig_num);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
minTimeFast = Inf;
tic;
for i=1:REPS
    tstart = tic;
    color_vector = fcn_geometry_fillColorFromNumberOrName(2,'vector regression segment fit',[],fig_num);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fillColorFromNumberOrName:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

