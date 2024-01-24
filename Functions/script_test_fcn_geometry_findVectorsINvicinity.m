% script_fcn_fcn_geometry_findVectorsINvicinity
% Tests fcn_geometry_findVectorsINvicinity

% Revision History
%
% 2024_01_24 - Aneeh Batchu
% -- wrote the script

%% Basic Example to test if the input vectors are in the vicinity of reference vector

inputVectors = [1 1 1; 2 2 2; 1 3 1; 2 2 1; 8 -3 2; -3 -3 -3; 5 -3 2];

% Reference Vector
refVector = [0.5 0.5 0.5];

% This is the tolerance limit to determine if a vector is in the vicinity
% of the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than tolerance limit, the input vector
% is considered to be in the vicinity of the reference vector.
tolerance = 0.5;

[vectorsCloseTOref, SameDirection_vectorsCloseTOref] = fcn_geometry_findVectorsINvicinity(inputVectors, refVector, tolerance);

disp(vectorsCloseTOref);
disp(SameDirection_vectorsCloseTOref);

%% Basic Example to test if the input vectors are in the vicinity of reference vector when tolerance is more than the previous case

inputVectors = [1 1 1; 2 2 2; 1 3 1; 2 2 1; 8 -3 2; -3 -3 -3; 5 -3 2];

% Reference Vector
refVector = [0.5 0.5 0.5];

% This is the tolerance limit to determine if a vector is in the vicinity
% of the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than tolerance limit, the input vector
% is considered to be in the vicinity of the reference vector.
tolerance = 2;

[vectorsCloseTOref, SameDirection_vectorsCloseTOref] = fcn_geometry_findVectorsINvicinity(inputVectors, refVector, tolerance);

disp(vectorsCloseTOref);
disp(SameDirection_vectorsCloseTOref);

%% Example: RefVector is [-1 -1 -1]

inputVectors = [1 1 1; 2 2 2; 1 3 1; 2 2 1; 8 -3 2; -3 -3 -3; 5 -3 2];

% Reference Vector
refVector = [-1 -1 -1];

% This is the tolerance limit to determine if a vector is in the vicinity
% of the reference vector or not. Here, the euclidean distance of the unit
% vectors of input vectors and the unit vector of reference vector is
% computed. If the distance is less than tolerance limit, the input vector
% is considered to be in the vicinity of the reference vector.
tolerance = 2;

[vectorsCloseTOref, SameDirection_vectorsCloseTOref] = fcn_geometry_findVectorsINvicinity(inputVectors, refVector, tolerance);

disp(vectorsCloseTOref);
disp(SameDirection_vectorsCloseTOref);



