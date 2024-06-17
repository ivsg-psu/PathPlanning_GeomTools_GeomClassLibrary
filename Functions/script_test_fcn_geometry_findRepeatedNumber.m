% script_test_fcn_geometry_findRepeatedNumber
% Exercises the function: fcn_geometry_findRepeatedNumber
% Revision history:
% 2024_6_17
% -- wrote the code

%% Test 1

points = [1, 1, 1, 2, 4, 4];
output = fcn_geometry_findRepeatedNumber(points);
display(output)
assert(isequal(length(output),4))

%%  Test 2 
num_points = 50;
range = [1,10];
points = randi(range,num_points,1);
output = fcn_geometry_findRepeatedNumber(points);
disp(output)
assert(isequal(length(output),10));