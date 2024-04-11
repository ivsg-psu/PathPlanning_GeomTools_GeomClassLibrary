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
    color_vector = fcn_geometry_fillColorFromNumberOrName(2,'vector regression segment fit',fig_num);
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

