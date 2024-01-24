% scripts_test_fcn_geometry_separatePointsIntoGrids.m

% This is a script to exercise the function: fcn_geometry_separatePointsIntoGrids.m
% This function was written on  2024_01_22 by V. Wagh (vbw5054@psu.edu)

% Revision history:
% 2024_01_22
% -- first write of the code

close all;

%% Basic Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ____            _        ______                           _
%  |  _ \          (_)      |  ____|                         | |
%  | |_) | __ _ ___ _  ___  | |__  __  ____ _ _ __ ___  _ __ | | ___
%  |  _ < / _` / __| |/ __| |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \
%  | |_) | (_| \__ \ | (__  | |____ >  < (_| | | | | | | |_) | |  __/
%  |____/ \__,_|___/_|\___| |______/_/\_\__,_|_| |_| |_| .__/|_|\___|
%                                                      | |
%                                                      |_|
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
%% BASIC example 1: only X
inputPoints = [3.2 7; 4.1 3; 1.2 2.6];
gridSize = 2;
gridBoundaries = [0 6];
[gridIndices,gridDomains,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries);

assert(isequal(length(gridIndices),3))
assert(isequal(length(gridDomains),3))
assert(isequal(length(gridCenters),3))

assert(isequal(gridIndices{1},3))
assert(isequal(gridIndices{2},1))
assert(isequal(gridIndices{3},2))
%% BASIC example 2: XY
inputPoints = [2.1 4.3; 1.2 3.2; 3.3 1.1];
gridSize = 2;
gridBoundaries = [0 6; 0 6];
[gridIndices,gridDomains,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries);

assert(isequal(length(gridIndices),9))
assert(isequal(length(gridDomains),9))
assert(isequal(length(gridCenters),9))

assert(isequal(gridIndices{2},2))
assert(isequal(gridIndices{4},3))
assert(isequal(gridIndices{6},1))
%% BASIC example 3: XYZ
inputPoints = [2.1 4.3; 1.2 3.2; 3.3 1.1];
gridSize = 2;
gridBoundaries = [0 6; 0 6; 0 6];
[gridIndices,gridDomains,gridCenters] = fcn_geometry_separatePointsIntoGrids(inputPoints, gridSize, gridBoundaries);

assert(isequal(length(gridIndices),27))
assert(isequal(length(gridDomains),27))
assert(isequal(length(gridCenters),27))
%%