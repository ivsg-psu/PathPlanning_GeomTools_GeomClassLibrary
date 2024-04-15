% script_test_fcn_geometry_fitVectorToNPoints
% Exercises the function: fcn_geometry_fitVectorToNPoints

% Revision history:
% 2021_05_24 - S. Brennan
% -- wrote the code
% -- revised from fcn_geometry_fitSlopeInterceptNPoints
% 2024_04_14 - S. Brennan
% -- added assertion testing

close all;


%% Test 1: a basic test
% -x + y + 1 = 0;
fig_num = 1;
points = [2 1; 3 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[0.5000,-0.5000]));
assert(isequal(round(unit_vector,4),[0.7071,0.7071]));

%% Test 2: a basic test where root point is at origin
fig_num = 2;
% -x + y = 0;
points = [1 1; 2 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[0 0]));
assert(isequal(round(unit_vector,4),[0.7071,0.7071]));

%% Test 3: a basic test where not pointed in [1 1] direction
fig_num = 11;
points = [2 2; 0 4];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[2 2]));
assert(isequal(round(unit_vector,4),[-0.7071,0.7071]));

%% Test 4: pointed opposite direction to show direction does not matter
% Same results as test 2
fig_num = 12;
points = [2 2; 1 1];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[0 0]));
assert(isequal(round(unit_vector,4),[0.7071,0.7071]));

%% Test 3.1: horizontal line
% WARNING: an give either direction!
fig_num = 31;
% y = 2;

points = [2 2; 5 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[0 2]));
assert(isequal(round(unit_vector,4),[-1 0]));

%% Test 3.2: horizontal line passing through origin
fig_num = 32;
% y = 2;

points = [2 0; 5 0];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));


%% Test 3.3: vertical line
fig_num = 33;
% x = 2;

points = [2 0; 2 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[2 0]));
assert(isequal(round(unit_vector,4),[0 1]));


%% Test 3.4: vertical line passing through origin
fig_num = 34;
% x = 2;

points = [0 0; 0 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(round(root_point,4),[0 0]));
assert(isequal(round(unit_vector,4),[0 1]));



%% Test 4: many points randomly determined
close all;

fig_num = 40;
A = -3;
B = 2;
C = 4;

Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*(-A/B) + (-C/B) + 0.2*randn(Npoints,1);
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(size(root_point),[1 2]));
assert(isequal(size(unit_vector),[1 2]));

%% Test 4.1: many horizontal points randomly determined
close all;

fig_num = 41;
A = 0;
B = 2;
C = 4;

Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*(-A/B) + (-C/B) + 0.2*randn(Npoints,1);
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(size(root_point),[1 2]));
assert(isequal(size(unit_vector),[1 2]));

%% Test 4.2: many vertical points randomly determined
fig_num = 42;

Npoints = 1000;
x_data = 2*ones(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(size(root_point),[1 2]));
assert(isequal(size(unit_vector),[1 2]));


%% Test 4.3: many vertical points
fig_num = 43;

Npoints = 1000;
x_data = 2*ones(Npoints,1)+ 0.2*randn(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(size(root_point),[1 2]));
assert(isequal(size(unit_vector),[1 2]));


%% Test 5: a singular situation (determinant is zero - which gives errors in normal regression - can only be solved in this vector method)
fig_num = 50;
points = [6 4; 3 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(size(root_point),[1 2]));
assert(isequal(size(unit_vector),[1 2]));


%% Noisy vertical line 
% this was giving wrong answers with linear regression and was why this
% function was written
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
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'Figure: %.0d - Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',fig_num, root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

assert(isequal(size(root_point),[1 2]));
assert(isequal(size(unit_vector),[1 2]));

%% Test of fast implementation mode 

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

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points, (fig_num));
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
    [root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_fitVectorToNPoints:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);

assert(averageTimeFast<averageTimeSlow);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
    fprintf(1,'Root point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end
