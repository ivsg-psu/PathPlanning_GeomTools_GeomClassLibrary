% script_test_fcn_geometry_fitVectorToPlane
% Exercises the function: fcn_geometry_fitVectorToPlane
% Revision history:
% 2021_05_24
% -- wrote the code
% -- revised from fcn_geometry_fitSlopeInterceptNPoints
%
close all;

%% Test 1: a basic test of generic plane

fig_num = 1;
figure(fig_num);
clf;
rng(1823);

true_ABC = fcn_geometry_calcUnitVector([ 5 3 1]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
x = points(:,1);
y = points(:,2);
true_z = points(:,3);

z = true_z;

[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);


plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,4),round(unit_vector,4)));
assert(isequal(round(true_root,4),round(vector_root,4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),round(0.000,4)));
assert(isequal(round(plane_distances,4),round(plane_distance_from_origin*ones(length(x),1),4)));


%% Test 2: a basic test of perfect XY plane

fig_num = 2;
figure(fig_num);
clf;
rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([0 0 1]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
x = points(:,1);
y = points(:,2);
true_z = points(:,3);

z = true_z;

[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);

plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,4),round(unit_vector,4)));
assert(isequal(round(true_root,4),round(vector_root,4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),round(0.000,4)));
assert(isequal(round(plane_distances,4),round(plane_distance_from_origin*ones(length(x),1),4)));

%% Test 3: a basic test of perfect XZ plane

fig_num = 3;
figure(fig_num);
clf;
rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([0 1 0]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
x = points(:,1);
y = points(:,2);
true_z = points(:,3);

z = true_z;

[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);

plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,4),round(unit_vector,4)));
assert(isequal(round(true_root,4),round(vector_root,4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),round(0.000,4)));
assert(isequal(round(plane_distances,4),round(plane_distance_from_origin*ones(length(x),1),4)));

%% Test 4: a basic test of perfect YZ plane

fig_num = 4;
figure(fig_num);
clf;
rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([1 0 0]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
x = points(:,1);
y = points(:,2);
true_z = points(:,3);

z = true_z;

[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);

plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,4),round(unit_vector,4)));
assert(isequal(round(true_root,4),round(vector_root,4)));
assert(isequal(round(standard_deviation_in_plane_orthogonals,4),round(0.000,4)));
assert(isequal(round(plane_distances,4),round(plane_distance_from_origin*ones(length(x),1),4)));

%% Test 21: a basic test of noisy XY plane

fig_num = 21;
figure(fig_num);
clf;
rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([0 0 1]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
true_x = points(:,1);
true_y = points(:,2);
true_z = points(:,3);

x = true_x + 0.01*randn(N_points,1);
y = true_y + 0.01*randn(N_points,1);
z = true_z + 0.01*randn(N_points,1);


[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);

plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,2),round(unit_vector,2)));
assert(isequal(round(true_root,1),round(vector_root,1)));
assert(standard_deviation_in_plane_orthogonals<3*0.02);
assert(length(plane_distances(:,1))==length(x));

%% Test 31: a basic test of noisy XZ plane

fig_num = 3;
figure(fig_num);
clf;
rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([0 1 0]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

true_x = points(:,1);
true_y = points(:,2);
true_z = points(:,3);

x = true_x + 0.01*randn(N_points,1);
y = true_y + 0.01*randn(N_points,1);
z = true_z + 0.01*randn(N_points,1);

[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);

plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,2),round(unit_vector,2)));
assert(isequal(round(true_root,1),round(vector_root,1)));
assert(standard_deviation_in_plane_orthogonals<3*0.02);
assert(length(plane_distances(:,1))==length(x));

%% Test 41: a basic test of noisy YZ plane

fig_num = 4;
figure(fig_num);
clf;
rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([1 0 0]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
true_x = points(:,1);
true_y = points(:,2);
true_z = points(:,3);

x = true_x + 0.01*randn(N_points,1);
y = true_y + 0.01*randn(N_points,1);
z = true_z + 0.01*randn(N_points,1);

[vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],fig_num);

plane_distance_from_origin = sum(true_ABC.*points(1,:),2);
true_root = true_ABC*plane_distance_from_origin;

% Print results
fprintf(1,'\n\nResults of plane fitting in figure %.0d: \n',fig_num);
fprintf(1,'True values: \n');
table_data = [true_ABC true_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Fitted values: \n');
table_data = [unit_vector vector_root];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

fprintf(1,'Error in fit: \n');
table_data = [abs(unit_vector-true_ABC) abs(vector_root-true_root)];
header_strings = [{'A'}, {'B'},{'C'},{'Root X'},{'Root Y'},{'Root Z'}]; % Headers for each column
formatter_strings = [{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'},{'%.6f'}]; % How should each column be printed?
N_chars = [12, 12, 12, 12, 12, 12]; % Specify spaces for each column
fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings, N_chars);

assert(isequal(round(true_ABC,2),round(unit_vector,2)));
assert(isequal(round(true_root,1),round(vector_root,1)));
assert(standard_deviation_in_plane_orthogonals<3*0.02);
assert(length(plane_distances(:,1))==length(x));

%% Test of fast implementation mode 

rng(1823);
warning('on','backtrace');

true_ABC = fcn_geometry_calcUnitVector([1 0 0]);
true_base_point = [1 2 3];

N_points = 100;

direction_ortho = fcn_geometry_calcOrthogonalVector(ones(N_points,1)*true_ABC, [], -1);
points = true_base_point + randn(N_points,1).*direction_ortho;

% For debugging
% figure(38383);
% clf;
% hold on;
% axis equal;
% grid on;
% view(3);
% plot3(points(:,1),points(:,2),points(:,3),'k.');
    

% Fill in x, y, and z
x = points(:,1);
y = points(:,2);
true_z = points(:,3);

z = true_z;

% Perform the calculation in slow mode
fig_num = [];
REPS = 100; minTimeSlow = Inf; 
tic;
for i=1:REPS
    tstart = tic;
    [vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],(fig_num));
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
    [vector_root, unit_vector, standard_deviation_in_plane_orthogonals, plane_distances] = fcn_geometry_fitVectorToPlane([x y z],(fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'\n\nComparison of fast and slow modes of fcn_geometry_fitVectorToPlane:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.5f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.5f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.5f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.5f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);
%% Fail conditions
if 1==0
    %% FAIL 1: points not long enough
    points = [2 3];
    [root_point, unit_vector] = fcn_geometry_fitVectorToPlane(points,fig_num);
    fprintf(1,'\n\nRoot point is: %.2f %.2f, Unit vector is: %.2f %.2f\n',root_point(1,1),root_point(1,2),unit_vector(1,1),unit_vector(1,2));
end
