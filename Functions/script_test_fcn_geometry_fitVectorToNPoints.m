% script_test_fcn_geometry_fitVectorToNPoints
% Exercises the function: fcn_geometry_fitVectorToNPoints
% Revision history:
% 2021_05_24
% -- wrote the code
% -- revised from fcn_geometry_fitSlopeInterceptNPoints
%
close all;
clc;


%% Test 1: a basic test
% -x + y + 1 = 0;
fig_num = 1;
points = [2 1; 3 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 2: a basic test
fig_num = 2;
% -x + y = 0;
points = [1 1; 2 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 1: a basic test
fig_num = 11;
points = [2 7; 5 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 1: a basic test
fig_num = 12;
points = [5 2; 8 7];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));


%% Test 3.1: horizontal line
fig_num = 31;
% y = 2;

points = [2 2; 5 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 3.1: horizontal line
fig_num = 32;
% y = 2;

points = [2 0; 5 0];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));


%% Test 3: vertical line
fig_num = 3;
% x = 2;

points = [2 0; 2 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));


%% Test 3: vertical line
fig_num = 3;
% x = 2;

points = [0 0; 0 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));


%% Test 4: many points randomly determined
close all;

fig_num = 4;
A = -3;
B = 2;
C = 4;

Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*(-A/B) + (-C/B) + 0.2*randn(Npoints,1);
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 5: many points randomly determined
close all;

fig_num = 4;
A = 0;
B = 2;
C = 4;

Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*(-A/B) + (-C/B) + 0.2*randn(Npoints,1);
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 4: many vertical points
fig_num = fig_num + 1;

Npoints = 1000;
x_data = 2*ones(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Test 4.2: many vertical points
fig_num = 42;

Npoints = 1000;
x_data = 2*ones(Npoints,1)+ 0.2*randn(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];

[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));


%% Test 5: a singular situation (determinant is zero - which gives b = 0)
fig_num = fig_num + 1;
points = [6 4; 3 2];
[root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));

%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_fitVectorToNPoints(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end