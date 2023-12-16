% script_test_fcn_geometry_calcUnitVector
% Exercises the function: fcn_geometry_calcUnitVector
% Revision history:
% 2023_12_14
% -- wrote the code

close all;
clc;


%% Test 1: a basic test
fig_num = 1;
input_vectors = [3 3]; 

unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num); %#ok<NASGU>


%% Test 2: many vectors
fig_num = 2;
input_vectors = randn(10,2); 
unit_vectors = fcn_geometry_calcUnitVector(input_vectors, fig_num);

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [slope,intercept] = fcn_geometry_calcUnitVector(points,fig_num);
    fprintf(1,'\n\nSlope is: %.2f, Intercept is: %.2f\n',slope,intercept);
end