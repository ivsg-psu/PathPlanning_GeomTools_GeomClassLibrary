% script_test_fcn_geometry_calcOrthogonalVector
% Exercises the function: fcn_geometry_calcOrthogonalVector

% Revision history:
% 2024_01_24 - S. Brennan
% -- wrote the code


close all;
clc;


%% Test 1: a basic test in 2D
fig_num = 1;
figure(fig_num);
clf;

input_vectors = [3 3]; 
seed_points = [];
unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num)); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors.*unit_vectors,2);
assert(all(abs(dot_product_sums)<(eps*100)));

%% Test 2: many vectors in 2D
fig_num = 2;
figure(fig_num);
clf;

input_vectors = randn(5,2); 
seed_points = [];
unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num)); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors.*unit_vectors,2);
assert(all(abs(dot_product_sums)<(eps*100)));

%% Test 3: a basic test in 3D
fig_num = 3;
figure(fig_num);
clf;

input_vectors = [1/2^0.5 1/2^0.5 0]; 
seed_points = [];
unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num)); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors*unit_vectors',2);
assert(all(abs(dot_product_sums)<(eps*100)));

%% Test 4: a basic test in 3D, many points
fig_num = 4;
figure(fig_num);
clf;

step = 0.05;

input_vectors = (step:step:1)'.*[3 2 4]; 
seed_points = [];
unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num)); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors*unit_vectors',2);
assert(all(abs(dot_product_sums)<(eps*100)));


%% Test 5: testing seed points
fig_num = 5;
figure(fig_num);
clf;


input_vectors = [3 3 3]; 
seed_points = [1 0 0];
unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num)); 

% Check that they are all unit length
length_errors = ones(length(unit_vectors(:,1)),1) - sum(unit_vectors.^2,2).^0.5;
assert(all(abs(length_errors)<(eps*100)));

% Check that dot products are zero
dot_product_sums = sum(input_vectors*unit_vectors',2);
assert(all(abs(dot_product_sums)<(eps*100)));


%% Testing fast mode
% Perform the calculation in slow mode
fig_num = [];

step = 0.05;

input_vectors = (step:step:1)'.*[3 2 4]; 
seed_points = [];

REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num));
    telapsed = toc(tstart);
    minTimeSlow = min(telapsed,minTimeSlow);
end
averageTimeSlow = toc/REPS;

% Perform the operation in fast mode
fig_num = -1;
REPS = 1000; minTimeFast = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    unit_vectors = fcn_geometry_calcOrthogonalVector(input_vectors, seed_points, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_calcOrthogonalVector:\n');
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
    [slope,intercept] = fcn_geometry_calcOrthogonalVector(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end