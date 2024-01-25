% script_test_fcn_geometry_findVectorsInAlignment
% Tests fcn_geometry_findVectorsInAlignment

% Revision History
%
% 2024_01_24 - Aneeh Batchu
% -- wrote the script
%

% To-DO
% -- write more fail conditions


%% Clear Workspace
clc
close all

%% Basic Example to test if the input vectors are in the vicinity of reference vector

% NOTE
%
% Unit Reference vector is shown in GREEN color
% Unit Input Vectors are shown in RED
% Vectors that are aligned with Unit Reference vector are shown in BLUE

inputVectors = [0.5 0.5 0.5];

% Reference Vector
refVector = [0.5 0.5 0.5];

% This is the tolerance limit to determine if a vector is in alignment with
% the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than or equal to tolerance limit, the
% input vector is considered to be in the vicinity of the reference vector
tolerance = 0.1;

[dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance, 1);

% The unit vector of input vector and reference vector are same. Therefore,
% the euclidean distance between both the vectors is zero.
assert(isequal(dist_btw_unit_refVectAndInputVec, 0));
assert(isequal(vectorsCloseToRef, inputVectors))


%% Basic Example to test if the input vectors are in alignment of reference vector 

inputVectors = [1 2 3];

% Reference Vector
refVector = [0.5 0.5 0.5];

% This is the tolerance limit to determine if a vector is in alignment with
% the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than or equal to tolerance limit, the
% input vector is considered to be in the vicinity of the reference vector
tolerance = 0.1;

[dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance,2);

unit_inputVectors = fcn_geometry_calcUnitVector(inputVectors, -1);
unit_refVector = fcn_geometry_calcUnitVector(refVector, -1); 

dist_calculated = sum((unit_inputVectors-unit_refVector).^2,2).^0.5;
% The distance is calculated manually to check the assertion
assert(isequal(dist_btw_unit_refVectAndInputVec, dist_calculated))

% The answer should be a zero*3 double empty matrix. Since, the input
% vector is not in alignment with the reference vector
assert(isequal(vectorsCloseToRef,  0 * zeros(0, 3)))

%% Basic Example to test if the input vectors are in alignment with reference vector

inputVectors = [1 1 1; 2 2 2; 1 3 1; 2 2 1; 8 -3 2; -3 -3 -3; 5 -3 2];

% Reference Vector
refVector = [0.5 0.5 0.5];

% This is the tolerance limit to determine if a vector is in alignment with
% the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than or equal to tolerance limit, the
% input vector is considered to be in the vicinity of the reference vector
tolerance = 0.5;

[dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance, 3);

assert(~isequal(dist_btw_unit_refVectAndInputVec, 0));

% Based on dist_btw_unit_refVectAndInputVec, this assertion is verfied
vectors_aligned = [1     1     1;
                   2     2     2;
                   2     2     1];

assert(isequal(vectorsCloseToRef, vectors_aligned));

% disp(dist_btw_unit_refVectAndInputVec);
% disp(vectorsCloseToRef);

%% Basic Example to test if the function works for 2D

inputVectors = [1 1];

% Reference Vector
refVector = [0.5 1];

% This is the tolerance limit to determine if a vector is in alignment with
% the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than or equal to tolerance limit, the
% input vector is considered to be in the vicinity of the reference vector
tolerance = 0.5;

[dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance);

unit_inputVectors = fcn_geometry_calcUnitVector(inputVectors, -1);
unit_refVector = fcn_geometry_calcUnitVector(refVector, -1); 

dist_calculated = sum((unit_inputVectors-unit_refVector).^2,2).^0.5;
% The distance is calculated manually to check the assertion
assert(isequal(dist_btw_unit_refVectAndInputVec, dist_calculated));

assert(isequal(vectorsCloseToRef, inputVectors));


%% Basic Example to test if the function works for 4D

inputVectors = [1 1 1 0.5];

% Reference Vector
refVector = [0.5 1 1 1];

% This is the tolerance limit to determine if a vector is in alignment with
% the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than or equal to tolerance limit, the
% input vector is considered to be in the vicinity of the reference vector
tolerance = 0.5;

[dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance);

unit_inputVectors = fcn_geometry_calcUnitVector(inputVectors, -1);
unit_refVector = fcn_geometry_calcUnitVector(refVector, -1); 

dist_calculated = sum((unit_inputVectors-unit_refVector).^2,2).^0.5;
% The distance is calculated manually to check the assertion
assert(isequal(dist_btw_unit_refVectAndInputVec, dist_calculated));

assert(isequal(vectorsCloseToRef, inputVectors));

%% Testing fast mode
% Perform the calculation in slow mode
fig_num = [];
REPS = 1000; minTimeSlow = Inf;
tic;
for i=1:REPS
    tstart = tic;
    [dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance, (fig_num));
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
    [dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance, (fig_num));
    telapsed = toc(tstart);
    minTimeFast = min(telapsed,minTimeFast);
end
averageTimeFast = toc/REPS;

fprintf(1,'Comparison of fast and slow modes of fcn_geometry_calcUnitVector:\n');
fprintf(1,'N repetitions: %.0d\n',REPS);
fprintf(1,'Slow mode average speed per call (seconds): %.8f\n',averageTimeSlow);
fprintf(1,'Slow mode fastest speed over all calls (seconds): %.8f\n',minTimeSlow);
fprintf(1,'Fast mode average speed per call (seconds): %.8f\n',averageTimeFast);
fprintf(1,'Fast mode fastest speed over all calls (seconds): %.8f\n',minTimeFast);
fprintf(1,'Average ratio of fast mode to slow mode (unitless): %.3f\n',averageTimeSlow/averageTimeFast);
fprintf(1,'Fastest ratio of fast mode to slow mode (unitless): %.3f\n',minTimeSlow/minTimeFast);


%% Fail Conditions

if 1 == 0

    inputVectors = [1 1];
    % Reference Vector
    refVector = [0.5 1 1];
    tolerance = 0.5;
    [dist_btw_unit_refVectAndInputVec, vectorsCloseToRef] = fcn_geometry_findVectorsInAlignment(refVector, inputVectors, tolerance);


end

