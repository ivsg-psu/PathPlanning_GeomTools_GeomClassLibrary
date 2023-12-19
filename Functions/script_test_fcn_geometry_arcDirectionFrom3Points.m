% script_test_fcn_geometry_arcDirectionFrom3Points - this is is a script written
% to test the function: fcn_geometry_arcDirectionFrom3Points.m.
%

% Revision history:
% 2023_12_19 - sbrennan@psu.edu
% -- original write of the code
close all;

%% BASIC example 1 - positive
fig_num = 1;

points1 = [0 0];
points2 = [1 4];
points3 = [0 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,1))

%% BASIC example 2 - negative
fig_num = 2;

points1 = [0 0];
points2 = [-1 4];
points3 = [0 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,-1))
%% BASIC example 2 - stacked
fig_num = 3;

points1 = [0 0; 3 3;];
points2 = [-1 4; 4 4];
points3 = [0 5; 3 5];
is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3, fig_num);

assert(isequal(is_counterClockwise,-1))



%% Fail cases follow
if 1==0
    %% FAIL 1: points2 not same length
    points1 = [2 3; 3 4];
    points2 = [4 5];
    points3 = [5 6; 7 8];
    is_counterClockwise = fcn_geometry_arcDirectionFrom3Points(points1, points2, points3,fig_num);

end



