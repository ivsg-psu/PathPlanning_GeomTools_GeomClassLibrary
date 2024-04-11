% script_test_fcn_geometry_fitSlopeInterceptNPoints
% Exercises the function: fcn_geometry_fitSlopeInterceptNPoints
% Revision history:
% 2020_06_25
% -- wrote the code
% 2021_05_24
% -- revised name to fcn_geometry_fitSlopeInterceptNPoints

close all;

fig_num = 1;

%% Test 1: a basic test
points = [2 3; 4 5];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
assert(isequal(slope,1));
assert(isequal(intercept,1));

%% Test 2: a vertical line
fig_num = fig_num + 1;
points = [2 0; 2 2];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 3: many points randomly determined
fig_num = fig_num + 1;
slope = -3;
intercept = 4;
Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.2*randn(Npoints,1);
points = [x_data,y_data];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 4: many vertical points
fig_num = fig_num + 1;
Npoints = 1000;
x_data = 2*ones(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 5: a singular situation (determinant is zero - which gives b = 0)
fig_num = fig_num + 1;
points = [6 4; 3 2];
[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test of fast mode
% Fill in data to test
slope = -3;
intercept = 4;
Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.2*randn(Npoints,1);
points = [x_data,y_data];


% Perform the calculation in slow mode
REPS = 10; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points, []);
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
REPS = 1000; minTimeFast = Inf; nsum = 10;
tic;
for i=1:REPS
    tstart = tic;
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points, -1);
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitSlopeInterceptNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

%% Debugging - This is a fail case
fig_num = 999;
points = [
   9.991259411578580   0.200000000000000
  10.001803670408686   0.300000000000000
   9.982375541376497   0.400000000000000
   9.977758789496919   0.500000000000000
  10.008682631868705   0.700000000000000
   9.988763772969360   0.800000000000000
  10.003402293646682   1.000000000000000
  10.011493908019808   1.200000000000000
   9.984955058249714   1.300000000000000
   9.995849995905278   1.400000000000000
  10.006337144959792   1.500000000000000
   9.982672092520335   1.600000000000000
  10.007135029591263   1.700000000000000
   9.998450373432329   1.800000000000000
  10.005582924433416   1.900000000000000
   9.979543989497849   2.000000000000000
   9.974383738360030   2.100000000000000
  10.001778511906023   2.500000000000000
   9.988736840713367   2.600000000000000
   9.989763064150692   2.800000000000000
  10.009874363357989   3.000000000000000
   9.984802416100781   3.100000000000000
   9.998229855546459   3.300000000000000
   9.984344000425939   3.400000000000000
   9.975609845675763   3.500000000000000
   9.991928970032030   3.600000000000000
   9.990306079204805   3.900000000000000
   9.983588133844744   4.000000000000000
  10.005622755909043   4.100000000000000
   9.981054848371585   4.200000000000000
   9.993381147363705   4.300000000000000
  10.007368013343468   4.400000000000000
  10.007946121327633   4.500000000000000
   9.990592139085731   4.600000000000000
  10.009923513002883   4.700000000000000
   9.971033609842408   4.900000000000000
    ];

[slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,999);

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_fitSlopeInterceptNPoints(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end